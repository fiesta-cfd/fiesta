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

#include "Kokkos_Core.hpp"
#include "debug.hpp"
#include "kokkosTypes.hpp"
#include "input.hpp"
#ifndef NOMPI
#include "mpi.hpp"
#include "mpi.h"
#endif

struct bc_L {
  FS4D u;
  int n, ng, bc_type;

  bc_L(int n_, int ng_, int bc_type_, FS4D u_)
      : n(n_), ng(ng_), bc_type(bc_type_), u(u_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int j, const int k, const int v) const {
    for (int i = 0; i < ng; ++i)
      if (v == 0 && bc_type == 1)
        u(n - i - 1, j, k, v) = -u(n + i, j, k, v);
      else
        u(n - i - 1, j, k, v) = u(n + i, j, k, v);
  }
};

struct bc_R {
  FS4D u;
  int n, ng, bc_type;

  bc_R(int n_, int ng_, int bc_type_, FS4D u_)
      : n(n_), ng(ng_), bc_type(bc_type_), u(u_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int j, const int k, const int v) const {
    for (int i = 0; i < ng; ++i)
      if (v == 0 && bc_type == 1)
        u(n + i, j, k, v) = -u(n - i - 1, j, k, v);
      else
        u(n + i, j, k, v) = u(n - i - 1, j, k, v);
  }
};

struct bc_B {
  FS4D u;
  int n, ng, bc_type;

  bc_B(int n_, int ng_, int bc_type_, FS4D u_)
      : n(n_), ng(ng_), bc_type(bc_type_), u(u_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, const int k, const int v) const {
    for (int j = 0; j < ng; ++j)
      if (v == 1 && bc_type == 1)
        u(i, n - j - 1, k, v) = -u(i, n + j, k, v);
      else
        u(i, n - j - 1, k, v) = u(i, n + j, k, v);
  }
};

struct bc_T {
  FS4D u;
  int n, ng, bc_type;

  bc_T(int n_, int ng_, int bc_type_, FS4D u_)
      : n(n_), ng(ng_), bc_type(bc_type_), u(u_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, const int k, const int v) const {
    for (int j = 0; j < ng; ++j)
      if (v == 1 && bc_type == 1)
        u(i, n + j, k, v) = -u(i, n - j - 1, k, v);
      else
        u(i, n + j, k, v) = u(i, n - j - 1, k, v);
  }
};

struct bc_H {
  FS4D u;
  int n, ng, bc_type;

  bc_H(int n_, int ng_, int bc_type_, FS4D u_)
      : n(n_), ng(ng_), bc_type(bc_type_), u(u_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, const int j, const int v) const {
    for (int k = 0; k < ng; ++k)
      if (v == 2 && bc_type == 1)
        u(i, j, n - k - 1, v) = -u(i, j, n + k, v);
      else
        u(i, j, n - k - 1, v) = u(i, j, n + k, v);
  }
};

struct bc_F {
  FS4D u;
  int n, ng, bc_type;

  bc_F(int n_, int ng_, int bc_type_, FS4D u_)
      : n(n_), ng(ng_), bc_type(bc_type_), u(u_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, const int j, const int v) const {
    for (int k = 0; k < ng; ++k)
      if (v == 2 && bc_type == 1)
        u(i, j, n + k, v) = -u(i, j, n - k - 1, v);
      else
        u(i, j, n + k, v) = u(i, j, n - k - 1, v);
  }
};

struct bc_xPer {
  FS4D u;
  int ng, nci;

  bc_xPer(int ng_, int nci_, FS4D u_) : ng(ng_), nci(nci_), u(u_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, const int j, const int k, const int v) const {
    double tmp;
    for (int g = 1; g < ng + 1; ++g) {
      tmp = u(i - g, j, k, v);
      u(i - g, j, k, v) = u(nci + g, j, k, v);
      u(nci + g, j, k, v) = tmp;
    }
  }
};

struct bc_yPer {
  FS4D u;
  int ng, ncj;

  bc_yPer(int ng_, int ncj_, FS4D u_) : ng(ng_), ncj(ncj_), u(u_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, const int j, const int k, const int v) const {
    double tmp;
    for (int g = 1; g < ng + 1; ++g) {
      tmp = u(i, j - g, k, v);
      u(i, j - g, k, v) = u(i, ncj + g, k, v);
      u(i, ncj + g, k, v) = tmp;
    }
  }
};

struct bc_zPer {
  FS4D u;
  int ng, nck;

  bc_zPer(int ng_, int nck_, FS4D u_) : ng(ng_), nck(nck_), u(u_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, const int j, const int k, const int v) const {
    double tmp;
    for (int g = 1; g < ng + 1; ++g) {
      tmp = u(i, j, k - g, v);
      u(i, j, k - g, v) = u(i, j, nck + g, v);
      u(i, j, nck + g, v) = tmp;
    }
  }
};

void applyBCs(struct inputConfig cf, class rk_func *f) {

  typedef Kokkos::MDRangePolicy<Kokkos::Rank<3>> policy_bl;
  typedef Kokkos::MDRangePolicy<Kokkos::Rank<4>> policy_bl4;

#ifndef NOMPI
  f->timers["halo"].reset();
  cf.m->haloExchange();
  f->timers["halo"].accumulate();
  f->timers["bc"].reset();
#else
  f->timers["bc"].reset();
  if (cf.xPer == 1)
    Kokkos::parallel_for(
        policy_bl4({0, 0, 0, 0}, {cf.nci, cf.ngj, cf.ngk, cf.nvt}),
        bc_xPer(cf.ng, cf.nci, f->var));
  if (cf.yPer == 1)
    Kokkos::parallel_for(
        policy_bl4({0, 0, 0, 0}, {cf.ngi, cf.ncj, cf.ngk, cf.nvt}),
        bc_yPer(cf.ng, cf.ncj, f->var));
  if (cf.ndim == 3 && cf.zPer == 1)
    Kokkos::parallel_for(
        policy_bl4({0, 0, 0, 0}, {cf.ngi, cf.ngj, cf.nck, cf.nvt}),
        bc_zPer(cf.ng, cf.nck, f->var));

#endif

  if (cf.xMinus < 0) {
    Kokkos::parallel_for(policy_bl({0, 0, 0}, {cf.ngj, cf.ngk, cf.nvt}),
                         bc_L(cf.ng, cf.ng, cf.bcL, f->var));
  }

  if (cf.xPlus < 0) {
    Kokkos::parallel_for(policy_bl({0, 0, 0}, {cf.ngj, cf.ngk, cf.nvt}),
                         bc_R(cf.ng + cf.nci, cf.ng, cf.bcR, f->var));
  }

  if (cf.yMinus < 0) {
    Kokkos::parallel_for(policy_bl({0, 0, 0}, {cf.ngi, cf.ngk, cf.nvt}),
                         bc_B(cf.ng, cf.ng, cf.bcB, f->var));
  }

  if (cf.yPlus < 0) {
    Kokkos::parallel_for(policy_bl({0, 0, 0}, {cf.ngi, cf.ngk, cf.nvt}),
                         bc_T(cf.ng + cf.ncj, cf.ng, cf.bcT, f->var));
  }

  if (cf.ndim == 3) {
    if (cf.zMinus < 0) {
      Kokkos::parallel_for(policy_bl({0, 0, 0}, {cf.ngi, cf.ngj, cf.nvt}),
                           bc_H(cf.ng, cf.ng, cf.bcH, f->var));
    }

    if (cf.zPlus < 0) {
      Kokkos::parallel_for(policy_bl({0, 0, 0}, {cf.ngi, cf.ngj, cf.nvt}),
                           bc_F(cf.ng + cf.nck, cf.ng, cf.bcF, f->var));
    }
  }
  f->timers["bc"].accumulate();
}
