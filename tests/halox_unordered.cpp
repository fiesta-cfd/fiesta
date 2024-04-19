#include <cassert>
#include "fiesta.hpp"
#include "input.hpp"
#include <vector>
#include "mpi.hpp"
#include "rkfunction.hpp"
#include "cart3d.hpp"
#include "bc.hpp"
#include <iostream>
#include "log2.hpp"
#include <cstdlib>

int main(int argc, char* argv[]) {
  {
    struct inputConfig cf;

    cf.ndim=3;
    cf.glbl_nci=6;
    cf.glbl_ncj=6;
    cf.glbl_nck=6;
    cf.ng=3;
    cf.xPer=1;
    cf.yPer=1;
    cf.zPer=1;
    cf.nvt=5;
    cf.nv=5;

    cf.ceq=false;
    cf.noise=false;
    cf.visc=false;
    cf.buoyancy=false;

    cf.xProcs=2;
    cf.yProcs=2;
    cf.zProcs=2;

    cf.dx=1;
    cf.dy=1;
    cf.dz=1;

    cf.R = 1;
    cf.ns=1;
    cf.speciesName= {"TestAir"};
    cf.gamma = {1.0};
    cf.M = {1.0};
    cf.mu = {1.0};

    cf.mpiScheme = atoi(argv[1]);

    Kokkos::InitArguments kokkosArgs;
    kokkosArgs.ndevices = 1;
    MPI_Init(NULL,NULL);
    int temp_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &temp_rank);
    Log::Logger(5,0,temp_rank);
    Kokkos::initialize(kokkosArgs);

    Log::debug("AT A");
    mpi_init(cf);

    Log::debug("AT B");
    rk_func *f;
    f = new cart3d_func(cf);

    //cf.mpiScheme=4;

    Log::debug("AT C");
    if (cf.mpiScheme == 1)
      cf.m = std::make_shared<copyHaloExchange>(cf,f->var);
    else if (cf.mpiScheme == 2)
      cf.m = std::make_shared<packedHaloExchange>(cf,f->var);
    else if (cf.mpiScheme == 3)
      cf.m = std::make_shared<directHaloExchange>(cf,f->var);

    Log::debug("MPI Scheme: {}",cf.mpiScheme);

    FS4DH varH = Kokkos::create_mirror_view(f->var);

    // x-y
    for (int i=0;i<9;++i)
      for (int j=0;j<9;++j)
        for (int k=0;k<9;++k)
          varH(i,j,k,0) = 10+cf.rank;
    for (int i=3;i<6;++i)
      for (int j=3;j<6;++j)
        for (int k=3;k<6;++k)
          varH(i,j,k,0) = cf.rank;

    Log::debug("varH(0,0,0,0)={}",varH(0,0,0,0));
    Log::debug("varH(1,1,1,0)={}",varH(1,1,1,0));
    Log::debug("varH(2,2,2,0)={}",varH(2,2,2,0));

    Log::debug("varH(6,6,6,0)={}",varH(6,6,6,0));
    Log::debug("varH(7,7,7,0)={}",varH(7,7,7,0));
    Log::debug("varH(8,8,8,0)={}",varH(8,8,8,0));

    Kokkos::deep_copy(f->var,varH);
    applyBCs(cf,f);
    Kokkos::deep_copy(varH,f->var);

    //Log::debugAll("{} {} {} {} {} {}",cf.xMinus,cf.xPlus,cf.yMinus,cf.yPlus,cf.zMinus,cf.zPlus);
    //Log::debugAll("cf.ngi={}, cf.ngj={}, cf.ngk={}",cf.ngi,cf.ngj,cf.ngk);

    if (cf.rank==0){
      // neighbors
      assert(cf.xMinus==4);
      assert(cf.xPlus==4);
      assert(cf.yMinus==2);
      assert(cf.yPlus==2);
      assert(cf.zMinus==1);
      assert(cf.zPlus==1);

      // x
      assert(varH(0,4,4,0)==4);
      assert(varH(1,4,4,0)==4);
      assert(varH(2,4,4,0)==4);
      assert(varH(6,4,4,0)==4);
      assert(varH(7,4,4,0)==4);
      assert(varH(8,4,4,0)==4);

      // y
      assert(varH(4,0,4,0)==2);
      assert(varH(4,1,4,0)==2);
      assert(varH(4,2,4,0)==2);
      assert(varH(4,6,4,0)==2);
      assert(varH(4,7,4,0)==2);
      assert(varH(4,8,4,0)==2);

      // z
      assert(varH(4,4,0,0)==1);
      assert(varH(4,4,1,0)==1);
      assert(varH(4,4,2,0)==1);
      assert(varH(4,4,6,0)==1);
      assert(varH(4,4,7,0)==1);
      assert(varH(4,4,8,0)==1);
    }
  
  }
  Kokkos::finalize();
  MPI_Finalize();
  return 0;
}
