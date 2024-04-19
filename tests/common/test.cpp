#include "fiesta.hpp"
#include "input.hpp"
#include <vector>
#ifdef HAVE_MPI
#include "mpi.hpp"
#endif
#include "rkfunction.hpp"
#include "cart3d.hpp"
#include "bc.hpp"
#include <iostream>
#include "debug.hpp"

rk_func * initTests(struct inputConfig& cf){
  cf.ndim=3;
  cf.glbl_nci=3;
  cf.glbl_ncj=3;
  cf.glbl_nck=3;
  cf.nci=3;
  cf.ncj=3;
  cf.nck=3;
  cf.ngi=9;
  cf.ngj=9;
  cf.ngk=9;
  cf.ni=4;
  cf.nj=4;
  cf.nk=4;
  cf.ng=3;
  cf.xPer=0;
  cf.yPer=0;
  cf.zPer=0;
  cf.nvt=5;
  cf.nv=5;

  cf.xMinus=-1;
  cf.xPlus=-1;
  cf.yMinus=-1;
  cf.yPlus=-1;
  cf.zMinus=-1;
  cf.zPlus=-1;

  cf.ceq=false;
  cf.noise=false;
  cf.visc=false;
  cf.buoyancy=false;

  cf.xProcs=1;
  cf.yProcs=1;
  cf.zProcs=1;

  cf.dx=1;
  cf.dy=1;
  cf.dz=1;

  cf.R = 1;
  cf.ns=1;
  cf.speciesName= {"TestAir"};
  cf.gamma = {1.0};
  cf.M = {1.0};
  cf.mu = {1.0};

  Kokkos::InitArguments kokkosArgs;
  kokkosArgs.ndevices = 1;
#ifdef HAVE_MPI
  MPI_Init(NULL,NULL);
  int temp_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &temp_rank);
  Log::Logger(5,0,temp_rank);
  Log::debug("TEST HAVE MPI");
#else
  Log::Logger(5,0,0);
#endif
  Kokkos::initialize(kokkosArgs);

  Log::debug("AT A");
#ifdef HAVE_MPI
  Log::debug("TEST MPI INIT");
  mpi_init(cf);
#endif

  Log::debug("AT B");
  rk_func *f;
  f = new cart3d_func(cf);

  Log::debug("AT C");

#ifdef HAVE_MPI
  cf.mpiScheme=1;
  //cf.m = std::make_shared<mpiBuffers>(cf);
  if (cf.mpiScheme == 1)
    //cf.m = new copyHaloExchange(cf, f->var);
    cf.m = std::make_shared<copyHaloExchange>(cf,f->var);
  else if (cf.mpiScheme == 2)
    //cf.m = new packedHaloExchange(cf, f->var);
    cf.m = std::make_shared<packedHaloExchange>(cf,f->var);
  else if (cf.mpiScheme == 3)
    //cf.m = new directHaloExchange(cf, f->var);
    cf.m = std::make_shared<directHaloExchange>(cf,f->var);
#endif

  Log::debug("AT D");
  return f;
}
