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

#ifndef FIESTA_INPUT_H
#define FIESTA_INPUT_H

#include "writer.hpp"
#include "Kokkos_Core.hpp"
#include "kokkosTypes.hpp"
#include "lua.hpp"
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#ifndef NOMPI
#include "mpi.h"
#endif
#include <string>
#include "timer.hpp"
//#include "rkfunction.hpp"
#include <map>
#include <vector>

// Lua error function
void error(lua_State *L, const char *fmt, ...);

// Lua get boolean value
int getglobbool(lua_State *L, const char *var);

// Lua get integer value
int getglobint(lua_State *L, const char *var);

// Lua get double value
double getglobdbl(lua_State *L, const char *var);

// configuration structure
struct inputConfig {
  fiestaTimer totalTimer;
  fiestaTimer initTimer;
  fiestaTimer simTimer;
  fiestaTimer loadTimer;
  fiestaTimer gridTimer;
  fiestaTimer writeTimer;

  class writer *w;
  class mpiHaloExchange *m;

  int colorFlag,timeFormat;
  std::string inputFname;
  std::string title;
  std::string terrainName;
  int xmp, ymp, zmp;
  int ndim;
  int glbl_ni, glbl_nj, glbl_nk;
  int glbl_nci, glbl_ncj, glbl_nck;
  std::vector<std::string> speciesName;
  double *gamma;
  double *M;
  double *mu;
  double R;
  int scheme;
  int visc;
  int gravity;
  double g_accel;
  double *g_vec;
  double dt, dx, dy, dz;
  int xProcs, yProcs, zProcs, numProcs;
  int xPlus, yPlus, zPlus;
  int xMinus, yMinus, zMinus;
  int rank;
  int nci, ncj, nck;
  int nt, ni, nj, nk, nv, ns, nvt;
  int iStart, jStart, kStart;
  int iEnd, jEnd, kEnd;
#ifndef NOMPI
  int mpiScheme;
  MPI_Comm comm;
#endif
  int cF, cB, cZ;
  int ng, ngi, ngj, ngk;
  int xPer, yPer, zPer;
  int restart;
  std::string restartName;
  std::string pathName;
  int tstart, tend;
  double time;
  int ceq,st;
  double kap, eps, alpha, beta, betae;
  int bcL, bcR, bcB, bcT, bcH, bcF;
  double n_dh, n_coff, n_eta;
  int noise, n_nt;
  int t;
  int grid;
  double h, tdx, tdy;

  int globalGridDims[3];
  int globalCellDims[3];
  int localGridDims[3];
  int localCellDims[3];
  int subdomainOffset[3];

  int out_freq, stat_freq, write_freq, restart_freq;
};

struct commandArgs {
  int versionFlag;
  int colorFlag;
  int timeFormat;
  int numThreads;
  int numDevices;
  std::string fileName;
};

struct commandArgs getCommandlineOptions(int argc, char **argv);
struct inputConfig executeConfiguration(struct commandArgs cargs);

int loadInitialConditions(struct inputConfig cf,  FS4D &v, FS4D &g);
int loadGrid(struct inputConfig cf, FS4D &v);

#endif // FIESTA_INPUT_HPP
