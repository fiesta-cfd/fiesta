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
#include <cassert>
#include "debug.hpp"
#include "test.hpp"

int main(){
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
  assert(varH(4,2,4,3)==-1);
  assert(varH(4,0,4,3)==-7);
  assert(varH(0,4,4,3)==6);
  assert(varH(4,8,4,3)==17);
  assert(varH(8,4,4,3)==4);

  assert(varH(0,0,4,3)==-6);
  assert(varH(8,0,4,3)==-8);
  assert(varH(0,8,4,3)==18);
  assert(varH(8,8,4,3)==16);
  }
  Kokkos::finalize();
#ifdef HAVE_MPI
  MPI_Finalize();
#endif


  return 0;
}
