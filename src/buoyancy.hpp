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
struct computeBuoyancy {
  FS4D dvar;
  FS4D var;
  FS2D rho;

  computeBuoyancy(FS4D dvar_, FS4D var_, FS2D rho_)
      : dvar(dvar_), var(var_), rho(rho_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, const int j) const {

    double v = var(i, j, 0, 1) / rho(i, j);
    double rhop = rho(i, j) - 1.0;
    double eps = 1e-6;

    if (rhop >= eps || rhop <= -eps) {
      double f = -9.81 * rhop;
      ;

      dvar(i, j, 0, 1) += f;
      dvar(i, j, 0, 2) += v * f;

      //            printf("%f\n",f);
    }
  }
};
#endif
