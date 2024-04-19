#include "catch2/catch_all.hpp"
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

TEST_CASE( "Outflow BC Test", "[outflow_bc]" ) {
  std::cout << "BEGIN TEST INIT\n";
  struct inputConfig cf;
  cf.xPer=0;
  cf.yPer=0;
  cf.zPer=0;
  rk_func *f = initTests(cf);

  cf.bcR=BCType::outflow;
  cf.bcL=BCType::outflow;
  cf.bcB=BCType::outflow;
  cf.bcT=BCType::outflow;
  cf.bcH=BCType::outflow;
  cf.bcF=BCType::outflow;

  FS4DH varH = Kokkos::create_mirror_view(f->var);

  // x-y
  varH(3,3,3,0) = 1;
  varH(4,3,3,0) = 2;
  varH(5,3,3,0) = 3;
  varH(3,4,3,0) = 4;
  varH(4,4,3,0) = 5;
  varH(5,4,3,0) = 6;
  varH(3,5,3,0) = 7;
  varH(4,5,3,0) = 8;
  varH(5,5,3,0) = 9;

  Kokkos::deep_copy(f->var,varH);
  applyBCs(cf,f);
  Kokkos::deep_copy(varH,f->var);

  REQUIRE(varH(4,0,3,0)==8);
  REQUIRE(varH(0,4,3,0)==6);
  REQUIRE(varH(4,8,3,0)==2);
  REQUIRE(varH(8,4,3,0)==4);
  REQUIRE(varH(0,0,3,0)==9);
  REQUIRE(varH(8,0,3,0)==7);
  REQUIRE(varH(0,8,3,0)==3);
  REQUIRE(varH(8,8,3,0)==1);

  // y-z
  varH(3,3,3,0) = 1;
  varH(3,4,3,0) = 2;
  varH(3,5,3,0) = 3;
  varH(3,3,4,0) = 4;
  varH(3,4,4,0) = 5;
  varH(3,5,4,0) = 6;
  varH(3,3,5,0) = 7;
  varH(3,4,5,0) = 8;
  varH(3,5,5,0) = 9;

  Kokkos::deep_copy(f->var,varH);
  applyBCs(cf,f);
  Kokkos::deep_copy(varH,f->var);

  REQUIRE(varH(3,4,0,0)==8);
  REQUIRE(varH(3,0,4,0)==6);
  REQUIRE(varH(3,4,8,0)==2);
  REQUIRE(varH(3,8,4,0)==4);
  REQUIRE(varH(3,0,0,0)==9);
  REQUIRE(varH(3,8,0,0)==7);
  REQUIRE(varH(3,0,8,0)==3);
  REQUIRE(varH(3,8,8,0)==1);
}

TEST_CASE( "Refelctive BC Test", "[reflective_bc]" ) {
  std::cout << "BEGIN TEST INIT\n";
  struct inputConfig cf;
  cf.xPer=0;
  cf.yPer=0;
  cf.zPer=0;
  rk_func *f = initTests(cf);
  cf.bcR=BCType::reflective;
  cf.bcL=BCType::reflective;
  cf.bcB=BCType::reflective;
  cf.bcT=BCType::reflective;
  cf.bcH=BCType::reflective;
  cf.bcF=BCType::reflective;

  FS4DH varH = Kokkos::create_mirror_view(f->var);

  // x-y
  varH(3,3,3,0) = 1; varH(3,3,3,1) = 1;
  varH(4,3,3,0) = 2; varH(4,3,3,1) = 2;
  varH(5,3,3,0) = 3; varH(5,3,3,1) = 3;
  varH(3,4,3,0) = 4; varH(3,4,3,1) = 4;
  varH(4,4,3,0) = 5; varH(4,4,3,1) = 5;
  varH(5,4,3,0) = 6; varH(5,4,3,1) = 6;
  varH(3,5,3,0) = 7; varH(3,5,3,1) = 7;
  varH(4,5,3,0) = 8; varH(4,5,3,1) = 8;
  varH(5,5,3,0) = 9; varH(5,5,3,1) = 9;

  Kokkos::deep_copy(f->var,varH);
  applyBCs(cf,f);
  Kokkos::deep_copy(varH,f->var);

  REQUIRE(varH(4,0,3,0)== 8); REQUIRE(varH(4,0,3,1)==-8);
  REQUIRE(varH(0,4,3,0)==-6); REQUIRE(varH(0,4,3,1)== 6);
  REQUIRE(varH(4,8,3,0)== 2); REQUIRE(varH(4,8,3,1)==-2);
  REQUIRE(varH(8,4,3,0)==-4); REQUIRE(varH(8,4,3,1)== 4);
  REQUIRE(varH(0,0,3,0)==-9); REQUIRE(varH(0,0,3,1)==-9);
  REQUIRE(varH(8,0,3,0)==-7); REQUIRE(varH(8,0,3,1)==-7);
  REQUIRE(varH(0,8,3,0)==-3); REQUIRE(varH(0,8,3,1)==-3);
  REQUIRE(varH(8,8,3,0)==-1); REQUIRE(varH(8,8,3,1)==-1);

  // y-z
  varH(3,3,3,1) = 1; varH(3,3,3,2) = 1;
  varH(3,4,3,1) = 2; varH(3,4,3,2) = 2;
  varH(3,5,3,1) = 3; varH(3,5,3,2) = 3;
  varH(3,3,4,1) = 4; varH(3,3,4,2) = 4;
  varH(3,4,4,1) = 5; varH(3,4,4,2) = 5;
  varH(3,5,4,1) = 6; varH(3,5,4,2) = 6;
  varH(3,3,5,1) = 7; varH(3,3,5,2) = 7;
  varH(3,4,5,1) = 8; varH(3,4,5,2) = 8;
  varH(3,5,5,1) = 9; varH(3,5,5,2) = 9;

  Kokkos::deep_copy(f->var,varH);
  applyBCs(cf,f);
  Kokkos::deep_copy(varH,f->var);

  REQUIRE(varH(3,4,0,1)== 8); REQUIRE(varH(3,4,0,2)==-8);
  REQUIRE(varH(3,0,4,1)==-6); REQUIRE(varH(3,0,4,2)== 6);
  REQUIRE(varH(3,4,8,1)== 2); REQUIRE(varH(3,4,8,2)==-2);
  REQUIRE(varH(3,8,4,1)==-4); REQUIRE(varH(3,8,4,2)== 4);
  REQUIRE(varH(3,0,0,1)==-9); REQUIRE(varH(3,0,0,2)==-9);
  REQUIRE(varH(3,8,0,1)==-7); REQUIRE(varH(3,8,0,2)==-7);
  REQUIRE(varH(3,0,8,1)==-3); REQUIRE(varH(3,0,8,2)==-3);
  REQUIRE(varH(3,8,8,1)==-1); REQUIRE(varH(3,8,8,2)==-1);
}

TEST_CASE( "Hydrostatic BC Test", "[hydrostatic_bc]" ) {
  std::cout << "BEGIN TEST INIT\n";
  {
  struct inputConfig cf;
  cf.xPer=0;
  cf.yPer=0;
  cf.zPer=0;
  rk_func *f = initTests(cf);

  FS4DH varH = Kokkos::create_mirror_view(f->var);
  varH(3,3,4,0) = 0;
  varH(4,3,4,0) = 0;
  varH(5,3,4,0) = 0;
  varH(3,4,4,0) = 0;
  varH(4,4,4,0) = 0;
  varH(5,4,4,0) = 0;
  varH(3,5,4,0) = 0;
  varH(4,5,4,0) = 0;
  varH(5,5,4,0) = 0;

  varH(3,3,4,1) = 0;
  varH(4,3,4,1) = 0;
  varH(5,3,4,1) = 0;
  varH(3,4,4,1) = 0;
  varH(4,4,4,1) = 0;
  varH(5,4,4,1) = 0;
  varH(3,5,4,1) = 0;
  varH(4,5,4,1) = 0;
  varH(5,5,4,1) = 0;

  varH(3,3,4,2) = 0;
  varH(4,3,4,2) = 0;
  varH(5,3,4,2) = 0;
  varH(3,4,4,2) = 0;
  varH(4,4,4,2) = 0;
  varH(5,4,4,2) = 0;
  varH(3,5,4,2) = 0;
  varH(4,5,4,2) = 0;
  varH(5,5,4,2) = 0;

  varH(3,3,4,3) = 1;
  varH(4,3,4,3) = 2;
  varH(5,3,4,3) = 3;
  varH(3,4,4,3) = 4;
  varH(4,4,4,3) = 5;
  varH(5,4,4,3) = 6;
  varH(3,5,4,3) = 7;
  varH(4,5,4,3) = 8;
  varH(5,5,4,3) = 9;

  varH(3,3,4,4) = 1;
  varH(4,3,4,4) = 1;
  varH(5,3,4,4) = 1;
  varH(3,4,4,4) = 1;
  varH(4,4,4,4) = 1;
  varH(5,4,4,4) = 1;
  varH(3,5,4,4) = 1;
  varH(4,5,4,4) = 1;
  varH(5,5,4,4) = 1;

  Kokkos::deep_copy(f->var,varH);

  cf.bcR=BCType::outflow;
  cf.bcL=BCType::outflow;
  cf.bcB=BCType::hydrostatic;
  cf.bcT=BCType::hydrostatic;
  cf.bcH=BCType::outflow;
  cf.bcF=BCType::outflow;
  //cf.ndim=2;
  applyBCs(cf,f);

  Kokkos::deep_copy(varH,f->var);

  //e
  REQUIRE(varH(4,2,4,3)==-1);
  REQUIRE(varH(4,0,4,3)==-7);
  REQUIRE(varH(0,4,4,3)==6);
  REQUIRE(varH(4,8,4,3)==17);
  REQUIRE(varH(8,4,4,3)==4);

  REQUIRE(varH(0,0,4,3)==-6);
  REQUIRE(varH(8,0,4,3)==-8);
  REQUIRE(varH(0,8,4,3)==18);
  REQUIRE(varH(8,8,4,3)==16);
  }
  Kokkos::finalize();
#ifdef HAVE_MPI
  MPI_Finalize();
#endif
}

TEST_CASE( "Periodic BC Test", "[periodic_bc]" ) {
  std::cout << "BEGIN TEST INIT\n";
  {
  struct inputConfig cf;
  cf.ndim=3;
  cf.glbl_nci=3;
  cf.glbl_ncj=3;
  cf.glbl_nck=10;
  cf.nci=3;
  cf.ncj=3;
  cf.nck=10;
  cf.ngi=9;
  cf.ngj=9;
  cf.ngk=16;
  cf.ni=4;
  cf.nj=4;
  cf.nk=11;
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
  std::cout << "TEST HAVE MPI\n";
  MPI_Init(NULL,NULL);
#endif
  Kokkos::initialize(kokkosArgs);

#ifdef HAVE_MPI
  std::cout << "DOING MPI_INIT\n";
  mpi_init(cf);
#endif

  rk_func *f;
  f = new cart3d_func(cf);


#ifdef HAVE_MPI
  cf.mpiScheme=1;
  cf.m = std::make_shared<orderedHostHaloExchange>(cf,f->var);
#endif

  FS4DH varH = Kokkos::create_mirror_view(f->var);
  varH(4,4,3,0)  = 50;
  varH(4,4,4,0)  = 51;
  varH(4,4,5,0)  = 52;
  varH(4,4,6,0)  = 53;
  varH(4,4,7,0)  = 54;
  varH(4,4,8,0)  = 55;
  varH(4,4,9,0)  = 56;
  varH(4,4,10,0) = 57;
  varH(4,4,11,0) = 58;
  varH(4,4,12,0) = 59;

  Kokkos::deep_copy(f->var,varH);

  cf.bcR=BCType::reflective;
  cf.bcL=BCType::reflective;
  cf.bcB=BCType::reflective;
  cf.bcT=BCType::reflective;
  cf.bcH=BCType::reflective;
  cf.bcF=BCType::reflective;
  applyBCs(cf,f);

  std::cout << cf.ngi << ", " << cf.ngj << ", " << cf.ngk << "\n";
  std::cout << cf.zPlus << " : " << cf.zMinus << "\n";

  Kokkos::deep_copy(varH,f->var);

  //e
  std::cout << "Left:  "
            << varH(4,4, 0,0) << ", "
            << varH(4,4, 1,0) << ", "
            << varH(4,4, 2,0) << ", "
            << varH(4,4, 3,0) << ", "
            << varH(4,4, 4,0) << ", "
            << varH(4,4, 5,0) << ", "
            << varH(4,4, 6,0) << ", "
            << varH(4,4, 7,0) << ", "
            << varH(4,4, 8,0) << ", "
            << varH(4,4, 9,0) << ", "
            << varH(4,4,10,0) << ", "
            << varH(4,4,11,0) << ", "
            << varH(4,4,12,0) << ", "
            << varH(4,4,13,0) << ", "
            << varH(4,4,14,0) << ", "
            << varH(4,4,15,0) << "\n";

  //std::cout << "Right: "
  //          << varH(4,4,7,0) << ", "
  //          << varH(4,4,8,0) << ", "
  //          << varH(4,4,9,0) << "\n";

  REQUIRE(varH(4,4,0,0)==57);
  REQUIRE(varH(4,4,1,0)==58);
  REQUIRE(varH(4,4,2,0)==59);
  
  REQUIRE(varH(4,4,13,0)==50);
  REQUIRE(varH(4,4,14,0)==51);
  REQUIRE(varH(4,4,15,0)==52);
  }
  Kokkos::finalize();
#ifdef HAVE_MPI
  MPI_Finalize();
#endif
}
