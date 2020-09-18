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

using namespace std;

struct maxVarFunctor2d {

  FS4D var;
  int v;

  maxVarFunctor2d(FS4D var_, int v_) : var(var_), v(v_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, const int j, double &lmax) const {

    double s = var(i, j, 0, v);

    if (s > lmax)
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

    if (s < lmin)
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
    if (s < lmin)
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

    if (s > lmax)
      lmax = s;
  }
};

void statusCheck(int cFlag, struct inputConfig cf, rk_func *f, double time,
                 fiestaTimer &wall, fiestaTimer &sim) {
  //policy_f cell_pol =
  //    policy_f({cf.ng, cf.ng}, {cf.ngi - cf.ng, cf.ngj - cf.ng});

  double max[cf.nvt];
  double min[cf.nvt];
  string vname;

  if (cf.rank == 0) {
    cout << c(cFlag, YEL) << "    Status: " << c(cFlag, NON) << endl;
 
    cout << "      Time Step:       " << c(cFlag, CYA) 
                                      << setprecision(0) << scientific
                                      << (float)cf.t
                                      << c(cFlag, NON)
         << "/" 
                                      << c(cFlag, CYA)
                                      << setprecision(0) << scientific
                                      << (float)cf.tend
                                      << c(cFlag, NON)
         << endl;

    cout << "      Simulation Time: " << c(cFlag, CYA) << setprecision(2)
         << scientific << time << c(cFlag, NON) << endl;
    cout << "      Wall Time:       " << c(cFlag, CYA) << wall.checkf()
         << c(cFlag, NON) << endl;
    cout << "      ETR:             " << c(cFlag, CYA)
         << sim.formatTime((double)cf.nt*sim.check()/(cf.t-1-cf.tstart)-sim.check())
         << c(cFlag, NON) << endl;
  }

  // Print Header
  if (cf.rank == 0) {
    cout << "      " << setw(13) << " " << setw(11) << right << "Min"
         << setw(11) << right << "Max" << endl;
    cout << "      "
         << "---------------------------------" << endl;
  }

  // var
  if (cf.ndim == 2) {
    policy_f cell_pol =
        policy_f({cf.ng, cf.ng}, {cf.ngi - cf.ng, cf.ngj - cf.ng});
    for (int v = 0; v < cf.nvt; ++v) {
      Kokkos::parallel_reduce(cell_pol, maxVarFunctor2d(f->var, v),
                              Kokkos::Max<double>(max[v]));
      Kokkos::parallel_reduce(cell_pol, minVarFunctor2d(f->var, v),
                              Kokkos::Min<double>(min[v]));

      #ifndef NOMPI
        MPI_Allreduce(&max[v], &max[v], 1, MPI_DOUBLE, MPI_MAX, cf.comm);
        MPI_Allreduce(&min[v], &min[v], 1, MPI_DOUBLE, MPI_MIN, cf.comm);
      #endif
    }
  } else {
    policy_f3 cell_pol = policy_f3({cf.ng,      cf.ng,          cf.ng},
                               {cf.ngi - cf.ng, cf.ngj - cf.ng, cf.ngk-cf.ng});
    for (int v = 0; v < cf.nvt; ++v) {
      Kokkos::parallel_reduce(cell_pol, maxVarFunctor3d(f->var, v),
                              Kokkos::Max<double>(max[v]));
      Kokkos::parallel_reduce(cell_pol, minVarFunctor3d(f->var, v),
                              Kokkos::Min<double>(min[v]));
      #ifndef NOMPI
        MPI_Allreduce(&max[v], &max[v], 1, MPI_DOUBLE, MPI_MAX, cf.comm);
        MPI_Allreduce(&min[v], &min[v], 1, MPI_DOUBLE, MPI_MIN, cf.comm);
      #endif
    }
  }

  if (cf.rank == 0) {
    // Check for nans and infs and print table
    for (int v = 0; v < cf.nvt; ++v) {
      stringstream ss, smax, smin;

      vname = f->varNames[v];

      if ((isnormal(min[v]) || min[v] == 0) &&
          (min[v] > -1.0e+200 && min[v] < 1.0e+200))
        smin << c(cFlag, CYA) << setw(11) << right << setprecision(2)
             << scientific << min[v] << c(cFlag, NON);
      else
        smax << c(cFlag, RED) << setw(11) << right << setprecision(2)
             << scientific << min[v] << c(cFlag, NON);

      if ((isnormal(max[v]) || max[v] == 0) &&
          (max[v] > -1.0e+200 && max[v] < 1.0e+200))
        smax << c(cFlag, CYA) << setw(11) << right << setprecision(2)
             << scientific << max[v] << c(cFlag, NON);
      else
        smax << c(cFlag, RED) << setw(11) << right << setprecision(2)
             << scientific << max[v] << c(cFlag, NON);

      cout << "      " << setw(13) << left << vname << smin.str() << smax.str()
           << endl;
    }
  }

  // varx
  if (cf.ndim == 2) {
    policy_f cell_pol =
        policy_f({cf.ng, cf.ng}, {cf.ngi - cf.ng, cf.ngj - cf.ng});
    for (int v = 0; v < cf.nvt; ++v) {
      Kokkos::parallel_reduce(cell_pol, maxVarFunctor2d(f->varx, v),
                              Kokkos::Max<double>(max[v]));
      Kokkos::parallel_reduce(cell_pol, minVarFunctor2d(f->varx, v),
                              Kokkos::Min<double>(min[v]));

      #ifndef NOMPI
        MPI_Allreduce(&max[v], &max[v], 1, MPI_DOUBLE, MPI_MAX, cf.comm);
        MPI_Allreduce(&min[v], &min[v], 1, MPI_DOUBLE, MPI_MIN, cf.comm);
      #endif
    }
  } else {
    policy_f3 cell_pol = policy_f3({cf.ng,      cf.ng,          cf.ng},
                               {cf.ngi - cf.ng, cf.ngj - cf.ng, cf.ngk-cf.ng});
    for (int v = 0; v < cf.nvt; ++v) {
      Kokkos::parallel_reduce(cell_pol, maxVarFunctor3d(f->varx, v),
                              Kokkos::Max<double>(max[v]));
      Kokkos::parallel_reduce(cell_pol, minVarFunctor3d(f->varx, v),
                              Kokkos::Min<double>(min[v]));
      #ifndef NOMPI
        MPI_Allreduce(&max[v], &max[v], 1, MPI_DOUBLE, MPI_MAX, cf.comm);
        MPI_Allreduce(&min[v], &min[v], 1, MPI_DOUBLE, MPI_MIN, cf.comm);
      #endif
    }
  }

  if (cf.rank == 0) {

    // Check for nans and infs and print table
    for (int v = 0; v < f->varxNames.size(); ++v) {
      stringstream ss, smax, smin;

      vname = f->varxNames[v];

      if ((isnormal(min[v]) || min[v] == 0) &&
          (min[v] > -1.0e+200 && min[v] < 1.0e+200))
        smin << c(cFlag, CYA) << setw(11) << right << setprecision(2)
             << scientific << min[v] << c(cFlag, NON);
      else
        smax << c(cFlag, RED) << setw(11) << right << setprecision(2)
             << scientific << min[v] << c(cFlag, NON);

      if ((isnormal(max[v]) || max[v] == 0) &&
          (max[v] > -1.0e+200 && max[v] < 1.0e+200))
        smax << c(cFlag, CYA) << setw(11) << right << setprecision(2)
             << scientific << max[v] << c(cFlag, NON);
      else
        smax << c(cFlag, RED) << setw(11) << right << setprecision(2)
             << scientific << max[v] << c(cFlag, NON);

      cout << "      " << setw(13) << left << vname << smin.str() << smax.str()
           << endl;
    }
  }
}
