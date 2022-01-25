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

#ifndef CART3D_H
#define CART3D_H

#include "Kokkos_Core.hpp"
#include "kokkosTypes.hpp"
#include "input.hpp"
#include "rkfunction.hpp"

class cart3d_func : public rk_func {

public:
  cart3d_func(struct inputConfig &cf);
  void compute();
  void preStep();
  void postStep();
  void preSim();
  void postSim();

  FS3D p;       // Pressure
  FS3D T;       // Temperature
  FS3D rho;     // Total Density
  FS4D vel;     // Velocity
  FS3D qx;      // Heat Fluxes in X direction
  FS3D qy;      // Heat Fluxes in Y direction
  FS3D qz;      // Heat Fluxes in Z direction
  FS3D fluxx;   // Weno Fluxes in X direction
  FS3D fluxy;   // Weno Fluxes in Y direction
  FS3D fluxz;   // Weno Fluxes in Z direction
  FS5D stressx; // Stress tensor on X faces
  FS5D stressy; // Stress tensor on Y faces
  FS5D stressz; // Stress tensor on Z faces
  FS4D gradRho; // Density Gradient array
  FS4D cFlux;
  FS6D mFlux;
  FS3D_I noise;
  FS1D cd; // Device configuration array
  FSCAL dxmag;
};

#endif
