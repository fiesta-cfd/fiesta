/*
  Copyright 2019-2021 The University of New Mexico

  This file is part of FIESTA.
  
  FIESTA is free software: you can redistribute it and/or modify it under the
  terms of the GNU Lesser General Public License as published by the Free
  Software Foundation, either version 3 of the License, or (at your option) any
  later version.
  
  FIESTA is distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
  A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
  details.
  
  You should have received a copy of the GNU Lesser General Public License
  along with FIESTA.  If not, see <https://www.gnu.org/licenses/>.
*/

#include "rkfunction.hpp"
#include "rk.hpp"
#include "kokkosTypes.hpp"
#include "Kokkos_Core.hpp"
#include "debug.hpp"
#include "input.hpp"
#include "bc.hpp"

void rkAdvance(struct inputConfig &cf, class std::unique_ptr<class rk_func>&f){
  typedef Kokkos::MDRangePolicy<Kokkos::Rank<3>> policy_1;
  // First Stage Compute
  f->compute();

  // assign temporary variables
  FS4D mytmp = f->tmp1;
  FS4D myvar = f->var;
  FS4D mydvar = f->dvar;
  FSCAL dt = cf.dt;
  FSCAL nvt = cf.nvt;

  // First stage update
  f->timers["rk"].reset();
  Kokkos::parallel_for(
      "Loop1", policy_1({0, 0, 0}, {cf.ngi, cf.ngj, cf.ngk}),
      KOKKOS_LAMBDA(const int i, const int j, const int k) {
        for (int v=0; v<nvt; ++v){
          mytmp(i, j, k, v) = myvar(i, j, k, v);
          myvar(i, j, k, v) =
            myvar(i, j, k, v) + 0.5 * dt * mydvar(i, j, k, v);
        }
      });
  Kokkos::fence();
  f->timers["rk"].accumulate();

  // Second stage compute
  f->compute();

  // assign temporary variables
  myvar = f->var;
  mydvar = f->dvar;

  // Second stage update
  f->timers["rk"].reset();
  Kokkos::parallel_for(
      "Loop2", policy_1({0, 0, 0}, {cf.ngi, cf.ngj, cf.ngk}),
      KOKKOS_LAMBDA(const int i, const int j, const int k) {
        for (int v=0; v<nvt; ++v){
          myvar(i, j, k, v) = mytmp(i, j, k, v) + dt * mydvar(i, j, k, v);
        }
      });
  Kokkos::fence();
  f->timers["rk"].accumulate();
}
