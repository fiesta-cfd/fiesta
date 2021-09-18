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

    assert(varH(4,0,3,0)==8);
    assert(varH(0,4,3,0)==6);
    assert(varH(4,8,3,0)==2);
    assert(varH(8,4,3,0)==4);
    assert(varH(0,0,3,0)==9);
    assert(varH(8,0,3,0)==7);
    assert(varH(0,8,3,0)==3);
    assert(varH(8,8,3,0)==1);

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

    assert(varH(3,4,0,0)==8);
    assert(varH(3,0,4,0)==6);
    assert(varH(3,4,8,0)==2);
    assert(varH(3,8,4,0)==4);
    assert(varH(3,0,0,0)==9);
    assert(varH(3,8,0,0)==7);
    assert(varH(3,0,8,0)==3);
    assert(varH(3,8,8,0)==1);
  }
  Kokkos::finalize();
#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return 0;
}
