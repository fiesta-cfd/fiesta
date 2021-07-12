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

#include "status.hpp"
#include "input.hpp"
#include <iomanip>
#include <locale>
#ifndef NOMPI
#include "mpi.h"
#endif
#include "output.hpp"
#include "rkfunction.hpp"
#include <cmath>
#include "fmt/core.h"
#include "pretty.hpp"

using namespace std;

struct maxVarFunctor2d {

  FS4D var;
  int v;

  maxVarFunctor2d(FS4D var_, int v_) : var(var_), v(v_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, const int j, double &lmax) const {

    double s = var(i, j, 0, v);

    if (abs(s) > lmax)
      lmax = s;
  }
};

struct minVarFunctor2d {

  FS4D var;
  int v;

  minVarFunctor2d(FS4D var_, int v_) : var(var_), v(v_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, const int j, double &lmin) const {

    double s = var(i, j, 0, v);

    if (abs(s) < lmin)
      lmin = s;
  }
};

struct minVarFunctor3d {
  FS4D var;
  int v;

  minVarFunctor3d(FS4D var_, int v_) : var(var_), v(v_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, const int j, const int k, double &lmin) const {
    double s = var(i, j, k, v);
    if (abs(s) < lmin)
      lmin = s;
  }
};

struct maxVarFunctor3d {

  FS4D var;
  int v;

  maxVarFunctor3d(FS4D var_, int v_) : var(var_), v(v_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, const int j, const int k, double &lmax) const {

    double s = var(i, j, k, v);

    if (abs(s) > lmax)
      lmax = s;
  }
};

bool isBad(double val){
      if ((isnormal(val) || val == 0) && (val > -1.0e+200 && val < 1.0e+200))
        return false;
      else
        return true;
}

bool isConcern(double val){
      if (val < -1.0e+16 || val > 1.0e+16)
        return true;
      else
        return false;
}

void statusCheck(int cFlag, struct inputConfig cf, rk_func *f, double time, fiestaTimer &wall, fiestaTimer &sim) {
  double max[cf.nvt];
  double min[cf.nvt];
  ansiColors c(cFlag);
  string smin,smax;

  if (cf.rank == 0) {
    cf.log->message("[{}] Reporting Status",cf.t);
    cout << fmt::format("{: >8}Timestep:  {}{}{}/{}{}{} ({}{:.0f}%{})\n","",
        c(magenta),cf.t,c(reset),c(magenta),cf.tend,c(reset),c(green),100.0*(double)cf.t/(double)cf.tend,c(reset));
    cout << fmt::format("{: >8}Sim Time:  {}{:.2g}{}s\n","",c(magenta),cf.time,c(reset));
    cout << fmt::format("{: >8}Wall Time: {}{:.0f}{}s\n","",c(magenta),wall.check(),c(reset));
    cout << fmt::format("{: >8}ETR:       {}{:.0f}{}s\n","", c(magenta),cf.nt*sim.check()/(cf.t-1-cf.tstart)-sim.check(),c(reset));

    cout << fmt::format("{: <8}{: <16}{: >11}{: >11}\n","","","Min","Max");
  }

  // var
  if (cf.ndim == 2) {
    policy_f cell_pol = policy_f({cf.ng, cf.ng}, {cf.ngi - cf.ng, cf.ngj - cf.ng});
    for (int v = 0; v < cf.nvt; ++v) {
      Kokkos::parallel_reduce(cell_pol, maxVarFunctor2d(f->var, v), Kokkos::Max<double>(max[v]));
      Kokkos::parallel_reduce(cell_pol, minVarFunctor2d(f->var, v), Kokkos::Min<double>(min[v]));
      #ifndef NOMPI
        MPI_Allreduce(&max[v], &max_recv[v], 1, MPI_DOUBLE, MPI_MAX, cf.comm);
        MPI_Allreduce(&min[v], &min_recv[v], 1, MPI_DOUBLE, MPI_MIN, cf.comm);
        max[v] = max_recv[v];
        min[v] = min_recv[v];
      #endif
    }
  } else {
    policy_f3 cell_pol = policy_f3({cf.ng, cf.ng, cf.ng}, {cf.ngi - cf.ng, cf.ngj - cf.ng, cf.ngk-cf.ng});
    for (int v = 0; v < cf.nvt; ++v) {
      Kokkos::parallel_reduce(cell_pol, maxVarFunctor3d(f->var, v), Kokkos::Max<double>(max[v]));
      Kokkos::parallel_reduce(cell_pol, minVarFunctor3d(f->var, v), Kokkos::Min<double>(min[v]));
      #ifndef NOMPI
        MPI_Allreduce(&max[v], &max_recv[v], 1, MPI_DOUBLE, MPI_MAX, cf.comm);
        MPI_Allreduce(&min[v], &min_recv[v], 1, MPI_DOUBLE, MPI_MIN, cf.comm);
        max[v] = max_recv[v];
        min[v] = min_recv[v];
      #endif
    }
  }

  if (cf.rank == 0) {
    for (int v = 0; v < cf.nvt; ++v) {
      if (isBad(min[v]))
        smin = format("{}{:>11.2e}{}",c(red),min[v],c(reset));
      else if (isConcern(min[v]))
        smin = format("{}{:>11.2e}{}",c(yellow),min[v],c(reset));
      else
        smin = format("{}{:>11.2e}{}",c(magenta),min[v],c(reset));

      if (isBad(max[v]))
        smax = format("{}{:>11.2e}{}",c(red),max[v],c(reset));
      else if (isConcern(max[v]))
        smax = format("{}{:>11.2e}{}",c(yellow),max[v],c(reset));
      else
        smax = format("{}{:>11.2e}{}",c(magenta),max[v],c(reset));

      cout << fmt::format("{: >8}{: <16}{: >11}{: >11}\n","",f->varNames[v],smin,smax);
    }
  }

  // varx
  if (cf.ndim == 2) {
    policy_f cell_pol = policy_f({cf.ng, cf.ng}, {cf.ngi - cf.ng, cf.ngj - cf.ng});
    for (size_t v = 0; v < f->varxNames.size(); ++v) {
      Kokkos::parallel_reduce(cell_pol, maxVarFunctor2d(f->varx, v), Kokkos::Max<double>(max[v]));
      Kokkos::parallel_reduce(cell_pol, minVarFunctor2d(f->varx, v), Kokkos::Min<double>(min[v]));
      #ifndef NOMPI
        MPI_Allreduce(&max[v], &max_recv[v], 1, MPI_DOUBLE, MPI_MAX, cf.comm);
        MPI_Allreduce(&min[v], &min_recv[v], 1, MPI_DOUBLE, MPI_MIN, cf.comm);
        max[v] = max_recv[v];
        min[v] = min_recv[v];
      #endif
    }
  } else {
    policy_f3 cell_pol = policy_f3({cf.ng,cf.ng,cf.ng}, {cf.ngi - cf.ng, cf.ngj - cf.ng, cf.ngk-cf.ng});
    for (size_t v = 0; v < f->varxNames.size(); ++v) {
      Kokkos::parallel_reduce(cell_pol, maxVarFunctor3d(f->varx, v), Kokkos::Max<double>(max[v]));
      Kokkos::parallel_reduce(cell_pol, minVarFunctor3d(f->varx, v), Kokkos::Min<double>(min[v]));
      #ifndef NOMPI
        MPI_Allreduce(&max[v], &max_recv[v], 1, MPI_DOUBLE, MPI_MAX, cf.comm);
        MPI_Allreduce(&min[v], &min_recv[v], 1, MPI_DOUBLE, MPI_MIN, cf.comm);
        max[v] = max_recv[v];
        min[v] = min_recv[v];
      #endif
    }
  }

  if (cf.rank == 0) {
    for (size_t v = 0; v < f->varxNames.size(); ++v) {
      if (isBad(min[v]))
        smin = format("{}{:>11.2e}{}",c(red),min[v],c(reset));
      else if (isConcern(min[v]))
        smin = format("{}{:>11.2e}{}",c(yellow),min[v],c(reset));
      else
        smin = format("{}{:>11.2e}{}",c(magenta),min[v],c(reset));

      if (isBad(max[v]))
        smax = format("{}{:>11.2e}{}",c(red),max[v],c(reset));
      else if (isConcern(max[v]))
        smax = format("{}{:>11.2e}{}",c(yellow),max[v],c(reset));
      else
        smax = format("{}{:>11.2e}{}",c(magenta),max[v],c(reset));

      cout << fmt::format("{: >8}{: <16}{: >11}{: >11}\n","",f->varxNames[v],smin,smax);
    }
  }
}
