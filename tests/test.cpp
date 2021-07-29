#include "catch2/catch_all.hpp"
#include "fiesta.hpp"
#include "input.hpp"
#include <vector>
#include "mpi.hpp"
#include "rkfunction.hpp"
#include "cart2d.hpp"
#include "bc.hpp"
#include <iostream>

rk_func * initTests(struct inputConfig& cf){
  cf.ndim=3;
  cf.glbl_nci=3;
  cf.glbl_ncj=3;
  cf.glbl_nck=3;
  cf.ng=3;
  //cf.xPer=0;
  //cf.yPer=0;
  //cf.zPer=0;
  cf.nvt=4;
  cf.nv=4;

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
  MPI_Init(NULL,NULL);
  Kokkos::initialize(kokkosArgs);


  mpi_init(cf);

  rk_func *f;
  f = new cart2d_func(cf);

  cf.m = std::make_shared<mpiBuffers>(cf);

  return f;
}

TEST_CASE( "Outflow BC Test", "[outflow_bc]" ) {
  struct inputConfig cf;
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
  struct inputConfig cf;
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
  {
  struct inputConfig cf;
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
  MPI_Finalize();
}

TEST_CASE( "Periodic BC Test", "[periodic_bc]" ) {
  {
  struct inputConfig cf;
  cf.xPer=1;
  cf.yPer=1;
  cf.zPer=1;
  rk_func *f = initTests(cf);

  FS4DH varH = Kokkos::create_mirror_view(f->var);
  varH(3,3,4,0) = 1;
  varH(4,3,4,0) = 2;
  varH(5,3,4,0) = 3;
  varH(3,4,4,0) = 4;
  varH(4,4,4,0) = 5;
  varH(5,4,4,0) = 6;
  varH(3,5,4,0) = 7;
  varH(4,5,4,0) = 8;
  varH(5,5,4,0) = 9;

  varH(4,3,3,1) = 1;
  varH(4,3,4,1) = 2;
  varH(4,3,5,1) = 3;
  varH(4,4,3,1) = 4;
  varH(4,4,4,1) = 5;
  varH(4,4,5,1) = 6;
  varH(4,5,3,1) = 7;
  varH(4,5,4,1) = 8;
  varH(4,5,5,1) = 9;

  Kokkos::deep_copy(f->var,varH);

  cf.bcR=BCType::reflective;
  cf.bcL=BCType::reflective;
  cf.bcB=BCType::outflow;
  cf.bcT=BCType::outflow;
  cf.bcH=BCType::outflow;
  cf.bcF=BCType::outflow;
  applyBCs(cf,f);

  std::cout << cf.xPlus << " : " << cf.xMinus << "\n";

  Kokkos::deep_copy(varH,f->var);

  //e
  REQUIRE(varH(0,4,4,0)==4);
  REQUIRE(varH(1,4,4,0)==5);
  REQUIRE(varH(2,4,4,0)==6);
  
  REQUIRE(varH(6,4,4,0)==4);
  REQUIRE(varH(7,4,4,0)==5);
  REQUIRE(varH(8,4,4,0)==6);

  REQUIRE(varH(4,4,0,1)==4);
  REQUIRE(varH(4,4,1,1)==5);
  REQUIRE(varH(4,4,2,1)==6);
  
  REQUIRE(varH(4,4,6,1)==4);
  REQUIRE(varH(4,4,7,1)==5);
  REQUIRE(varH(4,4,8,1)==6);
  }
  Kokkos::finalize();
  MPI_Finalize();
}
