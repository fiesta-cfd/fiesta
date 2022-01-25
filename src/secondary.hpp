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

#ifndef SECONDARY_H
#define SECONDARY_H

#include "kokkosTypes.hpp"

struct copyExtraVars2D {
  FS4D varx;
  FS3D vel;
  FS2D rho, p, T;

  copyExtraVars2D(FS4D varx_, FS3D v_, FS2D p_, FS2D rho_, FS2D T_)
      : varx(varx_), vel(v_), p(p_), rho(rho_), T(T_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, const int j) const {
    varx(i,j,0,0) = vel(i,j,0);
    varx(i,j,0,1) = vel(i,j,1);
    varx(i,j,0,2) = p(i,j);
    varx(i,j,0,3) = T(i,j);
    varx(i,j,0,4) = rho(i,j);
  }
};
struct copyExtraVars3D {
  FS4D varx;
  FS4D vel;
  FS3D rho, p, T;

  copyExtraVars3D(FS4D varx_, FS4D v_, FS3D p_, FS3D rho_, FS3D T_)
      : varx(varx_), vel(v_), p(p_), rho(rho_), T(T_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, const int j, const int k) const {
    varx(i,j,k,0) = vel(i,j,k,0);
    varx(i,j,k,1) = vel(i,j,k,1);
    varx(i,j,k,2) = vel(i,j,k,2);
    varx(i,j,k,3) =   p(i,j,k);
    varx(i,j,k,4) =   T(i,j,k);
    varx(i,j,k,5) = rho(i,j,k);
  }
};
struct calculateRhoPT2D {
private:
  FS4D var;
  FS2D rho, p, T;
  FS1D cd;

public:
  calculateRhoPT2D(const FS4D& var_, FS2D& p_, FS2D& rho_, FS2D& T_, const FS1D& cd_)
      : var(var_), p(p_), rho(rho_), T(T_), cd(cd_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, const int j) const {

    int ns = (int)cd(0);
    int nv = (int)cd(4);
    FSCAL gamma, gammas, Rs;
    FSCAL Cp = 0;
    FSCAL Cv = 0;

    rho(i, j) = 0.0;

    // Total Density for this cell
    for (int s = 0; s < ns; ++s) {
      rho(i, j) = rho(i, j) + var(i, j, 0, 3 + s);
    }

    // Calculate mixture ratio of specific heats
    for (int s = 0; s < ns; ++s) {
      gammas = cd(6 + 3 * s);
      Rs = cd(6 + 3 * s + 1);

      // accumulate mixture heat capacity by mass fraction weights
      Cp =
          Cp + (var(i, j, 0, 3 + s) / rho(i, j)) * (gammas * Rs / (gammas - 1));
      Cv = Cv + (var(i, j, 0, 3 + s) / rho(i, j)) * (Rs / (gammas - 1));
    }
    gamma = Cp / Cv;

    // calculate pressure assuming perfect gas
    p(i, j) =
        (gamma - 1) * (var(i, j, 0, 2) -
                       (0.5 / rho(i, j)) * (var(i, j, 0, 0) * var(i, j, 0, 0) +
                                            var(i, j, 0, 1) * var(i, j, 0, 1)));

    // calcualte cell temperature
    T(i, j) = p(i, j) / ((Cp - Cv) * rho(i, j));
  }
};

struct calculateRhoPT3D {
  FS4D var;
  FS3D p, rho, T;
  FS1D cd;

  calculateRhoPT3D(FS4D var_, FS3D p_, FS3D rho_, FS3D T_, FS1D cd_)
      : var(var_), p(p_), rho(rho_), T(T_), cd(cd_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, const int j, const int k) const {

    int ns = (int)cd(0);
    int nv = (int)cd(4);
    FSCAL gamma, gammas, Rs;
    FSCAL Cp = 0;
    FSCAL Cv = 0;

    rho(i, j, k) = 0.0;

    for (int s = 0; s < ns; ++s) {
      rho(i, j, k) = rho(i, j, k) + var(i, j, k, 4 + s);
    }

    for (int s = 0; s < ns; ++s) {
      gammas = cd(6 + 3 * s);
      Rs = cd(6 + 3 * s + 1);

      Cp = Cp +
           (var(i, j, k, 4 + s) / rho(i, j, k)) * (gammas * Rs / (gammas - 1));
      Cv = Cv + (var(i, j, k, 4 + s) / rho(i, j, k)) * (Rs / (gammas - 1));
    }

    gamma = Cp / Cv;

    p(i, j, k) = (gamma - 1) *
                 (var(i, j, k, 3) -
                  (0.5 / rho(i, j, k)) * (var(i, j, k, 0) * var(i, j, k, 0) +
                                          var(i, j, k, 1) * var(i, j, k, 1) +
                                          var(i, j, k, 2) * var(i, j, k, 2)));
    // calcualte cell temperature
    T(i, j,k) = p(i, j,k) / ((Cp - Cv) * rho(i, j,k));
  }
};

#endif
