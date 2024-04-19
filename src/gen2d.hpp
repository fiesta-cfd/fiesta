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

#ifndef GEN2D_H
#define GEN2D_H

#include "kokkosTypes.hpp"
#include "input.hpp"
#include "noise.hpp"
#include "rkfunction.hpp"
#ifdef HAVE_MPI
#include "mpi.hpp"
#endif
#include "advect.hpp"
#include "flux.hpp"
#include "metric.hpp"
#include "presgrad.hpp"
#include "secondary.hpp"
#include "velocity.hpp"

class gen2d_func : public rk_func {

public:
  gen2d_func(struct inputConfig &cf_);

  void compute();
  void preStep();
  void postStep();
  void preSim();
  void postSim();

  FS2D p;       // Pressure
  FS2D T;       // Temperature
  FS3D tvel;    // Transformed velocity
  FS2D rho;     // Total Density
  FS2D fluxx;   // Weno Fluxes in X direction
  FS2D fluxy;   // Weno Fluxes in Y direction
  FS2D_I noise; // Noise indicator array
  FS1D cd;      // Device configuration array
  FS4D metrics; // jacobian metrics

#ifdef HAVE_MPI
  FS4D ls, lr, rs, rr, bs, br, ts, tr;
  FS4DH lsH, lrH, rsH, rrH, bsH, brH, tsH, trH;
#endif
};

#endif
