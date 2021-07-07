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

#include <cstdio>

struct bc_gen {
  FS4D u;
  int ihat, jhat, khat, ng, type, nv;

  bc_gen(int type_, int ng_, int i_, int j_, int k_, int nv_, FS4D u_)
      : type(type_), ng(ng_), ihat(i_), jhat(j_), khat(k_), nv(nv_), u(u_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, const int j, const int k) const {

    if (type == 0){
      for (int v=0; v<nv; ++v){
        for (int n=0; n<ng; ++n){
          u(i+ihat*(n+1),j+jhat*(n+1),k+khat*(n+1),v) = u(i-ihat*n,j-jhat*n,k-khat*n,v);
        }
      }
    }
    if (type == 1){
      double reflect=1.0;
      for (int v=0; v<nv; ++v){
        reflect=1.0;
        if(ihat != 0 && v==0) reflect=-1.0;
        if(jhat != 0 && v==1) reflect=-1.0;
        if(khat != 0 && v==2) reflect=-1.0;
        for (int n=0; n<ng; ++n){
          u(i+ihat*(n+1),j+jhat*(n+1),k+khat*(n+1),v) = reflect*u(i-ihat*n,j-jhat*n,k-khat*n,v);
        }
      }
    }
    if (type == 3){
      // Do normal extension first (to get momentums and other variables
      for (int v=0; v<nv; ++v){
        for (int n=0; n<ng; ++n){
          u(i+ihat*(n+1),j+jhat*(n+1),k+khat*(n+1),v) = u(i-ihat*n,j-jhat*n,k-khat*n,v);
        }
      }
      // correct energy variable with linear pressure extension
      int ig,jg,kg; //coordinates of ghost points
      int i1,j1,k1; //coordinates of real points
      int i2,j2,k2; //coordinates of real points
      for (int n=0; n<ng; ++n){
        ig = i+ihat*(n+1);
        jg = j+jhat*(n+1);
        kg = k+khat*(n+1);
        i1= ig-ihat;
        j1= jg-jhat;
        k1= kg-khat;
        i2= ig-2*ihat;
        j2= jg-2*jhat;
        k2= kg-2*khat;

        //compute kinetic energies
        double ke = 0.5*(u(ig,jg,kg,0)*u(ig,jg,kg,0) + u(ig,jg,kg,1)*u(ig,jg,kg,1))/u(ig,jg,kg,3);
        double ke1= 0.5*(u(i1,j1,k1,0)*u(i1,j1,k1,0) + u(i1,j1,k1,1)*u(i1,j1,k1,1))/u(i1,j1,k1,3);
        double ke2= 0.5*(u(i2,j2,k2,0)*u(i2,j2,k2,0) + u(i2,j2,k2,1)*u(i2,j2,k2,1))/u(i2,j2,k2,3);

        // extrapolate pressure
        u(ig,jg,kg,2)= 2*(u(i1,j1,k1,2)-ke1) - (u(i2,j2,k2,2)-ke2) + ke;
      }
    }
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
    Kokkos::parallel_for(policy_bl({cf.ng, 0, 0}, {cf.ng+1, cf.ngj, cf.ngk}),
                         bc_gen(cf.bcL, cf.ng, -1, 0, 0, cf.nvt, f->var));
  }

  if (cf.xPlus < 0) {
    Kokkos::parallel_for(policy_bl({cf.ng+cf.nci-1, 0, 0}, {cf.ng+cf.nci, cf.ngj, cf.ngk}),
                         bc_gen(cf.bcR, cf.ng, 1, 0, 0, cf.nvt, f->var));
  }

  if (cf.yMinus < 0) {
    Kokkos::parallel_for(policy_bl({0, cf.ng, 0}, {cf.ngi, cf.ng+1, cf.ngk}),
                         bc_gen(cf.bcB, cf.ng, 0, -1, 0, cf.nvt, f->var));
  }

  if (cf.yPlus < 0) {
    Kokkos::parallel_for(policy_bl({0, cf.ng+cf.ncj-1, 0}, {cf.ngi, cf.ng+cf.ncj, cf.ngk}),
                         bc_gen(cf.bcT, cf.ng, 0, 1, 0, cf.nvt, f->var));
  }

  if (cf.ndim == 3) {
    if (cf.zMinus < 0) {
      Kokkos::parallel_for(policy_bl({0, 0, cf.ng}, {cf.ngi, cf.ngj, cf.ng+1}),
                           bc_gen(cf.bcH, cf.ng, 0, 0, -1, cf.nvt, f->var));
    }

    if (cf.zPlus < 0) {
      Kokkos::parallel_for(policy_bl({0, 0, cf.ng+cf.nck-1}, {cf.ngi, cf.ngj, cf.ng+cf.nck}),
                           bc_gen(cf.bcF, cf.ng, 0, 0, 1, cf.nvt, f->var));
    }
  }
  f->timers["bc"].accumulate();
}
