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

#include "input.hpp"
#include "Kokkos_Core.hpp"
#include "debug.hpp"
#include "kokkosTypes.hpp"
#include "string.h"
#include <getopt.h>
#include <iostream>
#include <string>
#include <regex>
#include "luaReader.hpp"
#include "hdf.hpp"
#include <sys/types.h>
#include <sys/stat.h>

struct commandArgs getCommandlineOptions(int argc, char **argv){

  // create command argumet structure
  struct commandArgs cArgs;

  // Defaults
  cArgs.fileName = "fiesta.lua";
  cArgs.timeFormat = 2;
  cArgs.versionFlag = 0;
  cArgs.colorFlag = 0;
  cArgs.numDevices = 1;
  cArgs.numThreads = 1;

  // create options
  static struct option long_options[] = {
      {"version", no_argument, NULL, 'v'},
      {"color", optional_argument, NULL, 'c'},
      {"time-format", optional_argument, NULL, 't'},
      {"kokkos-ndevices", optional_argument, NULL, 'k'},
      {"kokkos-num-devices", optional_argument, NULL, 'k'},
      {"kokkos-threads", optional_argument, NULL, 'k'},
      {"kokkos-help", optional_argument, NULL, 'k'},
      {"kokkos-numa", optional_argument, NULL, 'k'},
      {NULL, 0, NULL, 0}};

  std::string copt;
  int c = 1;
  int opt_index;
  while ((c = getopt_long(argc, argv, "vctn:", long_options, &opt_index)) != -1) {
    switch (c) {
    case 'c':
      if (optarg)
        copt = std::string(optarg);
      else
        copt = "auto";
      break;
    case 'n':
#ifdef HAVE_CUDA
      cArgs.numDevices = atoi(optarg);
#elif HAVE_OPENMP
      cArgs.numThreads = atoi(optarg);
#endif
      break;
    case 'v':
      cArgs.versionFlag = 1;
      break;
    case 't':
      if (optarg)
        cArgs.timeFormat = atoi(optarg);
      else
        cArgs.timeFormat = 0;
      break;
    }
  }

  if (copt.compare("off") == 0)
    cArgs.colorFlag = 0;
  else if (copt.compare("on") == 0)
    cArgs.colorFlag = 1;
  else if (copt.compare("auto") == 0)
    cArgs.colorFlag = 2;

  if (optind < argc)
    cArgs.fileName = std::string(argv[optind]);
  else
    cArgs.fileName = std::string("fiesta.lua");

  // If the --version command line option is present, then just print version
  // and compilation information and exit.
  if (cArgs.versionFlag) {
    cout << "Fiesta" << endl;
    cout << "Version:    '" << FIESTA_VERSION << "'" << endl;
    cout << "Build Type: '" << FIESTA_OPTIONS << "'" << endl;
    cout << "Build Time: '" << FIESTA_BTIME << "'" << endl;
    exit(EXIT_SUCCESS);
  }

  return cArgs;
}

struct inputConfig executeConfiguration(struct commandArgs cargs) {
  struct inputConfig cf;
  cf.colorFlag = cargs.colorFlag;
  cf.timeFormat = cargs.timeFormat;

  luaReader L(cargs.fileName);

  // Required Parameters
  L.get("nt",    cf.nt);
  L.get("ni",    cf.glbl_nci);
  L.get("nj",    cf.glbl_ncj);
  L.get("ns",    cf.ns);
  L.get("dt",    cf.dt);
  L.get("R",     cf.R);
  L.get("title", cf.title);
  L.getArray("species_names",cf.speciesName,cf.ns);

  // Defaultable Parameters
  L.get("tstart",   cf.tstart,  0);
  L.get("time",     cf.time,    0.0);
  L.get("ceq",      cf.ceq,     0);
  L.get("noise",    cf.noise,   0);
  L.get("buoyancy", cf.gravity, 0);
  L.get("ndim",     cf.ndim,    2);
  L.get("visc",     cf.visc,    0);
  L.get("procsx",   cf.xProcs,  1);
  L.get("procsy",   cf.yProcs,  1);
  L.get("ng",       cf.ng,      3);
  L.get("xPer",     cf.xPer,    0);
  L.get("yPer",     cf.yPer,    0);
  L.get("bcXmin",   cf.bcL,     0);
  L.get("bcXmax",   cf.bcR,     0);
  L.get("bcYmin",   cf.bcB,     0);
  L.get("bcYmax",   cf.bcT,     0);
  L.get("restart",  cf.restart, 0);
  L.get("out_freq",     cf.out_freq,    0);
  L.get("write_freq",   cf.write_freq,  0);
  L.get("restart_freq", cf.restart_freq,0);
  L.get("stat_freq",    cf.stat_freq,   0);
  L.get("restartName",  cf.restartName, "restart-0000000.h5");
  L.get("terrainName", cf.terrainName, "terrain.h5");
  L.get("pathName",     cf.pathName, ".");

  // Check if pathName is accessable
  struct stat st;
  if(stat(cf.pathName.c_str(),&st) != 0){
    printf("Cannot access: %s\n",cf.pathName.c_str());
    exit(EXIT_FAILURE);
  }else if(st.st_mode & S_IFDIR == 0){
    printf("Cannot access: %s\n",cf.pathName.c_str());
    exit(EXIT_FAILURE);
  }

  std::string scheme, grid, mpi;
  L.get("scheme", scheme,"weno5");
  L.get("grid",   grid,"cartesian");
  L.get("mpi",    mpi, "host");

  // Dependent Parameters
  if (cf.ceq == 1) {
    L.get("kappa",  cf.kap);
    L.get("epsilon",cf.eps);
    L.get("alpha",  cf.alpha);
    L.get("beta",   cf.beta);
    L.get("betae",  cf.betae);
    L.get("st",     cf.st,10);
  }
  if (cf.noise == 1) {
    L.get("n_dh",  cf.n_dh);
    L.get("n_eta", cf.n_eta);
    L.get("n_coff",cf.n_coff);
    L.get("n_nt",  cf.n_nt,1);
  }
  if (cf.ndim == 3) {
    L.get("nk",     cf.glbl_nck);
    L.get("procsz", cf.zProcs,1);
    L.get("zPer",   cf.zPer,0);
    L.get("bcZmin", cf.bcH,0);
    L.get("bcZmaz", cf.bcF,0);
  }
  if (grid.compare("cartesian") == 0) {
    cf.grid = 0;
    L.get("dx",cf.dx);
    L.get("dy",cf.dy);
    if (cf.ndim == 3)
      L.get("dz",cf.dz);
  }
  if (grid.compare("terrain") == 0) {
    if (cf.ndim == 3){
      cf.grid = 2;
      cf.dx=1.0;
      cf.dy=1.0;
      cf.dz = 1.0;
      L.get("tdx",cf.tdx);
      L.get("tdy",cf.tdy);
      L.get("h",cf.h);
    } else {
      printf("ndim must be equal to 3 for terrain");
      exit(EXIT_FAILURE);
    }
  }
  // Array Parameters
  cf.M = (double *)malloc(cf.ns * sizeof(double));
  cf.gamma = (double *)malloc(cf.ns * sizeof(double));
  cf.mu = (double *)malloc(cf.ns * sizeof(double));
  L.getArray("gamma",cf.gamma,cf.ns);
  L.getArray("M",cf.M,cf.ns);
  L.getArray("mu",cf.mu,cf.ns);

  // Close Lua File
  L.close();

  // Save input file name
  cf.inputFname = cargs.fileName;

  // Set Scheme Number from Name
  cf.scheme = 1;
  if (scheme.compare("weno5") == 0)
    cf.scheme = 1;
  if (scheme.compare("centered4") == 0)
    cf.scheme = 2;
  if (scheme.compare("quick") == 0)
    cf.scheme = 3;

#ifndef NOMPI
  // MPI halo exchange strategy from name
  if (mpi.compare("host") == 0)
    cf.mpiScheme = 1;
  else if (mpi.compare("gpu-aware") == 0)
    cf.mpiScheme = 2;
  else if (mpi.compare("gpu-type") == 0)
    cf.mpiScheme = 3;
  else {
      printf("Invalid MPI communication scheme.\n");
      exit(EXIT_FAILURE);
  }
#endif

  // Set Grid Options from Grid String
  if (grid.compare("generalized") == 0) {
    cf.grid = 1;
    cf.dx = 1.0;
    cf.dy = 1.0;
    if (cf.ndim == 3)
      cf.dz = 1.0;
  }

  // Calculate time indices
  cf.t = cf.tstart;
  cf.tend = cf.tstart + cf.nt;


  // Set defaults for 2d grids
  if (cf.ndim == 2) {
    cf.nv = 3 + cf.ns;
    cf.dz = cf.dx;
    cf.zProcs = 1;
    cf.glbl_nck = 1;
    cf.bcH = 0;
    cf.bcF = 0;
    cf.zPer = 0;
  }else{
    cf.nv = 4 + cf.ns;
  }

  // nvt = Number of Variables Total (including c variables)
  cf.nvt = cf.nv;
  if (cf.ceq == 1)
    cf.nvt += 5;

  /* calculate number of grid vertices from number of cells */
  cf.glbl_ni = cf.glbl_nci + 1;
  cf.glbl_nj = cf.glbl_ncj + 1;
  if (cf.ndim == 3)
    cf.glbl_nk = cf.glbl_nck + 1;
  else
    cf.glbl_nk = 1;

  /* Set MPI defaults for 1 processor*/
  cf.rank = 0;
  cf.numProcs = 1;
  cf.nci = cf.glbl_nci;
  cf.ncj = cf.glbl_ncj;
  cf.nck = cf.glbl_nck;

  cf.ngi = cf.nci + 2 * cf.ng;
  cf.ngj = cf.ncj + 2 * cf.ng;
  if (cf.ndim == 3)
    cf.ngk = cf.nck + 2 * cf.ng;
  else
    cf.ngk = 1;

  cf.iStart = 0;
  cf.iEnd = cf.nci;
  cf.jStart = 0;
  cf.jEnd = cf.ncj;
  cf.kStart = 0;
  cf.kEnd = cf.nck;

  cf.ni = cf.nci + 1;
  cf.nj = cf.ncj + 1;
  if (cf.ndim == 3)
    cf.nk = cf.nck + 1;
  else
    cf.nk = 1;

  cf.xMinus = -1 + cf.xPer;
  cf.xPlus = -1 + cf.xPer;
  cf.yMinus = -1 + cf.yPer;
  cf.yPlus = -1 + cf.yPer;
  cf.zMinus = -1 + cf.zPer;
  cf.zPlus = -1 + cf.zPer;

  return cf;
}

int loadInitialConditions(struct inputConfig cf, FS4D &deviceV, FS4D &deviceG) {

  //int ii,jj,kk;
  double x,y,z;
  FS4DH hostV = Kokkos::create_mirror_view(deviceV);
  FS4DH hostG = Kokkos::create_mirror_view(deviceG);

  Kokkos::deep_copy(hostG,deviceG);

  luaReader L(cf.inputFname);

  for (int v = 0; v < cf.nv; ++v) {
    if (cf.ndim == 3) {
      for (int k = cf.ng; k < cf.nck + cf.ng; ++k) {
        for (int j = cf.ng; j < cf.ncj + cf.ng; ++j) {
          for (int i = cf.ng; i < cf.nci + cf.ng; ++i) {
            x = 0;
            y = 0;
            z = 0;
            for (int ix=0; ix<2; ++ix){
              for (int iy=0; iy<2; ++iy){
                for (int iz=0; iz<2; ++iz){
                  x += hostG(i+ix-cf.ng,j+iy-cf.ng,k+iz-cf.ng,0);
                  y += hostG(i+ix-cf.ng,j+iy-cf.ng,k+iz-cf.ng,1);
                  z += hostG(i+ix-cf.ng,j+iy-cf.ng,k+iz-cf.ng,2);
                }
              }
            }
            x = x/8.0;
            y = y/8.0;
            z = z/8.0;
            hostV(i, j, k, v) = L.call("f",4,x,y,z,(double)v);
          }
        }
      }
    } else {
      for (int j = cf.ng; j < cf.ncj + cf.ng; ++j) {
        for (int i = cf.ng; i < cf.nci + cf.ng; ++i) {
            x = 0;
            y = 0;
            for (int ix=0; ix<2; ++ix){
              for (int iy=0; iy<2; ++iy){
                x += hostG(i+ix-cf.ng,j+iy-cf.ng,0,0);
                y += hostG(i+ix-cf.ng,j+iy-cf.ng,0,1);
              }
            }
            x = x/4.0;
            y = y/4.0;
            hostV(i, j, 0, v) = L.call("f",4,x,y,0.0,(double)v);
        }
      }
    }
  }

  L.close();

  Kokkos::deep_copy(deviceV, hostV);

  return 0;
}

void error(lua_State *L, const char *fmt, ...) {
  va_list argp;
  va_start(argp, fmt);
  vfprintf(stderr, fmt, argp);
  va_end(argp);
  lua_close(L);
  exit(EXIT_FAILURE);
}

int loadGrid(struct inputConfig cf, FS4D &deviceV) {

  FS4DH hostV = Kokkos::create_mirror_view(deviceV);

  if (cf.grid == 1) {

    luaReader L(cf.inputFname);

    int ii,jj,kk;
    for (int v = 0; v < cf.ndim; ++v) {
      for (int k = 0; k < cf.nk; ++k) {
        for (int j = 0; j < cf.nj; ++j) {
          for (int i = 0; i < cf.ni; ++i) {
            ii = cf.iStart + i;
            jj = cf.jStart + j;
            kk = cf.kStart + k;
            hostV(i, j, k, v) = L.call("g",4,(double)ii,(double)jj,(double)kk,(double)v);
          }
        }
      }
    }
    L.close();
    Kokkos::deep_copy(deviceV, hostV);
  } else if (cf.grid == 2) {
    cf.w->readTerrain(cf,deviceV);
  } else {
    for (int k = 0; k < cf.nk; ++k) {
      for (int j = 0; j < cf.nj; ++j) {
        for (int i = 0; i < cf.ni; ++i) {
          hostV(i, j, k, 0) = (cf.iStart + i) * cf.dx;
          hostV(i, j, k, 1) = (cf.jStart + j) * cf.dy;
          hostV(i, j, k, 2) = (cf.kStart + k) * cf.dz;
        }
      }
    }
    Kokkos::deep_copy(deviceV, hostV);
  }


  return 0;
}
