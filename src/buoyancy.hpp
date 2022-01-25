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

#ifndef bouyancy_hpp
#define bouyancy_hpp
struct computeBuoyancy3D {
  FS4D dvar;
  FS4D var,varx;
  FS3D rho;
  FSCAL g;
  FSCAL rhoRef;

  computeBuoyancy3D(FS4D dvar_, FS4D var_, FS4D varx_, FS3D rho_, FSCAL g_, FSCAL r_)
      : dvar(dvar_), var(var_), varx(varx_), rho(rho_), g(g_), rhoRef(r_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, const int j, const int k) const {

    FSCAL v = var(i, j, k, 1) / rho(i,j,k);
    FSCAL rhop = rho(i,j,k) - rhoRef;
    //FSCAL eps = 1e-6;

    //if (rhop >= eps || rhop <= -eps) {
      FSCAL f = -g * rhop;
      dvar(i, j, k, 1) += f;
      dvar(i, j, k, 3) += v * f;
      //varx(i,j,k,21) = f;
      //varx(i,j,k,22) = v*f;
    //}
  }
};

struct computeBuoyancy {
  FS4D dvar;
  FS4D var;
  FS2D rho;
  FSCAL g;
  FSCAL rhoRef;

  computeBuoyancy(FS4D dvar_, FS4D var_, FS2D rho_, FSCAL g_, FSCAL r_)
      : dvar(dvar_), var(var_), rho(rho_), g(g_), rhoRef(r_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, const int j) const {

    //FSCAL v = var(i, j, 0, 1) / rho(i, j);
    FSCAL v = var(i, j, 0, 1) / var(i,j,0,3);
    FSCAL rhop = rho(i, j) - rhoRef;
    FSCAL eps = 1e-6;

    if (rhop >= eps || rhop <= -eps) {
      FSCAL f = -g * rhop;
      dvar(i, j, 0, 1) += f;
      dvar(i, j, 0, 2) += abs(v) * f;
    }
  }
};
#endif
