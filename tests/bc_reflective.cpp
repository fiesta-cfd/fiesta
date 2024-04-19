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
  {
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

    assert(varH(4,0,3,0)== 8); assert(varH(4,0,3,1)==-8);
    assert(varH(0,4,3,0)==-6); assert(varH(0,4,3,1)== 6);
    assert(varH(4,8,3,0)== 2); assert(varH(4,8,3,1)==-2);
    assert(varH(8,4,3,0)==-4); assert(varH(8,4,3,1)== 4);
    assert(varH(0,0,3,0)==-9); assert(varH(0,0,3,1)==-9);
    assert(varH(8,0,3,0)==-7); assert(varH(8,0,3,1)==-7);
    assert(varH(0,8,3,0)==-3); assert(varH(0,8,3,1)==-3);
    assert(varH(8,8,3,0)==-1); assert(varH(8,8,3,1)==-1);

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

    assert(varH(3,4,0,1)== 8); assert(varH(3,4,0,2)==-8);
    assert(varH(3,0,4,1)==-6); assert(varH(3,0,4,2)== 6);
    assert(varH(3,4,8,1)== 2); assert(varH(3,4,8,2)==-2);
    assert(varH(3,8,4,1)==-4); assert(varH(3,8,4,2)== 4);
    assert(varH(3,0,0,1)==-9); assert(varH(3,0,0,2)==-9);
    assert(varH(3,8,0,1)==-7); assert(varH(3,8,0,2)==-7);
    assert(varH(3,0,8,1)==-3); assert(varH(3,0,8,2)==-3);
    assert(varH(3,8,8,1)==-1); assert(varH(3,8,8,2)==-1);
  }
  Kokkos::finalize();
#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return 0;
}

