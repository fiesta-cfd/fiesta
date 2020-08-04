#include "status.hpp"
#include "input.hpp"
#include <iomanip>
#include <locale>
#ifndef NOMPI
#include "mpi.h"
#endif
#include "output.hpp"
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

void statusCheck(int cFlag, struct inputConfig cf, FS4D var, int t, double time,
                 fiestaTimer &wall) {
  policy_f cell_pol =
      policy_f({cf.ng, cf.ng}, {cf.ngi - cf.ng, cf.ngj - cf.ng});

  //    double myMax[cf.nvt];
  double max[cf.nvt];
  //    double myMin[cf.nvt];
  double min[cf.nvt];
  string vname;

  if (cf.rank == 0) {
    cout << c(cFlag, YEL) << "    Status: " << c(cFlag, NON) << endl;
    cout << "      Time Step:       " << c(cFlag, CYA) << t << c(cFlag, NON)
         << "/" << c(cFlag, CYA) << cf.tend << c(cFlag, NON) << endl;
    cout << "      Simulation Time: " << c(cFlag, CYA) << setprecision(5)
         << scientific << time << c(cFlag, NON) << endl;
    cout << "      Wall Time:       " << c(cFlag, CYA) << wall.checkf()
         << c(cFlag, NON) << endl;
  }

  if (cf.ndim == 2) {
    policy_f cell_pol =
        policy_f({cf.ng, cf.ng}, {cf.ngi - cf.ng, cf.ngj - cf.ng});
    for (int v = 0; v < cf.nvt; ++v) {
      Kokkos::parallel_reduce(cell_pol, maxVarFunctor2d(var, v),
                              Kokkos::Max<double>(max[v]));
      Kokkos::parallel_reduce(cell_pol, minVarFunctor2d(var, v),
                              Kokkos::Min<double>(min[v]));

#ifndef NOMPI
      MPI_Allreduce(&max[v], &max[v], 1, MPI_DOUBLE, MPI_MAX, cf.comm);
      MPI_Allreduce(&min[v], &min[v], 1, MPI_DOUBLE, MPI_MIN, cf.comm);
#endif
    }
  } else {
    policy_f3 cell_pol =
        policy_f3({cf.ng, cf.ng, cf.ng},
                  {cf.ngi - cf.ng, cf.ngj - cf.ng, cf.ngk - cf.ng});
    for (int v = 0; v < cf.nvt; ++v) {
      Kokkos::parallel_reduce(cell_pol, maxVarFunctor3d(var, v),
                              Kokkos::Max<double>(max[v]));
      Kokkos::parallel_reduce(cell_pol, minVarFunctor3d(var, v),
                              Kokkos::Min<double>(min[v]));
#ifndef NOMPI
      MPI_Allreduce(&max[v], &max[v], 1, MPI_DOUBLE, MPI_MAX, cf.comm);
      MPI_Allreduce(&min[v], &min[v], 1, MPI_DOUBLE, MPI_MIN, cf.comm);
#endif
    }
  }

  // Test inf, nan and zero
  // max[1] = max[1]/0.0;
  // min[2] = min[2]/-0.0;
  // min[0] = nan("");
  // max[3] = 0.0;

  // Print Header
  if (cf.rank == 0) {
    cout << "      " << setw(13) << " " << setw(11) << right << "Min"
         << setw(11) << right << "Max" << endl;
    cout << "      "
         << "---------------------------------" << endl;

    // Check for nans and infs and print table
    for (int v = 0; v < cf.nvt; ++v) {
      stringstream ss, smax, smin;

      if (cf.ndim == 2) {
        if (v == 0)
          vname = "X-Momentum: ";
        else if (v == 1)
          vname = "Y-Momentum: ";
        else if (v == 2)
          vname = "Energy: ";
        else if (v > 2 && v < cf.nv) {
          ss << "Density " << (v - 2) << ": ";
          vname = ss.str();
        } else {
          ss << "C-" << (v - cf.nv) << ": ";
          vname = ss.str();
        }
      } else {
        if (v == 0)
          vname = "X-Momentum: ";
        else if (v == 1)
          vname = "Y-Momentum: ";
        else if (v == 2)
          vname = "Z-Momentum: ";
        else if (v == 3)
          vname = "Energy: ";
        else if (v > 3 && v < cf.nv) {
          ss << "Density " << (v - 3) << ": ";
          vname = ss.str();
        } else {
          ss << "C-" << (v - cf.nv) << ": ";
          vname = ss.str();
        }
      }

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
