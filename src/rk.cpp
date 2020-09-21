#include "rk.hpp"
#include "fiesta.hpp"
#include "Kokkos_Core.hpp"
#include "debug.hpp"
#include "input.hpp"
#include "rkfunction.hpp"
#include "bc.hpp"

void rkAdvance(struct inputConfig &cf, class rk_func *f){
  typedef Kokkos::MDRangePolicy<Kokkos::Rank<3>> policy_1;
  // apply boundary conditions
  applyBCs(cf, f->var);

  // First Stage Compute
  f->compute();

  // assign temporary variables
  FS4D mytmp = f->tmp1;
  FS4D myvar = f->var;
  FS4D mydvar = f->dvar;

  // First stage update
  f->timers["rk"].reset();
  Kokkos::parallel_for(
      "Loop1", policy_1({0, 0, 0}, {cf.ngi, cf.ngj, cf.ngk}),
      KOKKOS_LAMBDA(const int i, const int j, const int k) {
        for (int v=0; v<cf.nvt; ++v){
          mytmp(i, j, k, v) = myvar(i, j, k, v);
          myvar(i, j, k, v) =
            myvar(i, j, k, v) + 0.5 * cf.dt * mydvar(i, j, k, v);
        }
      });
  Kokkos::fence();
  f->timers["rk"].accumulate();

  // apply boundary conditions
  applyBCs(cf, f->var);

  // Second stage compute
  f->compute();

  // assign temporary variables
  // mytmp = f->tmp1;
  myvar = f->var;
  mydvar = f->dvar;

  // Second stage update
  f->timers["rk"].reset();
  Kokkos::parallel_for(
      "Loop2", policy_1({0, 0, 0}, {cf.ngi, cf.ngj, cf.ngk}),
      KOKKOS_LAMBDA(const int i, const int j, const int k) {
        for (int v=0; v<cf.nvt; ++v){
          myvar(i, j, k, v) = mytmp(i, j, k, v) + cf.dt * mydvar(i, j, k, v);
        }
      });
  Kokkos::fence();
  f->timers["rk"].accumulate();
}
