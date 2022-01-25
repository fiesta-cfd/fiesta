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
#include <sys/types.h>
#include <sys/stat.h>
#include "unistd.h"
#include "log2.hpp"
#include "bc.hpp"
#include <filesystem>

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
  cArgs.verbosity = 3;

  // create options
  static struct option long_options[] = {
      {"version", no_argument, NULL, 'V'},
      {"verbosity", optional_argument, NULL, 'v'},
      {"color", optional_argument, NULL, 'c'},
      {"kokkos-ndevices", optional_argument, NULL, 'k'},
      {"kokkos-num-devices", optional_argument, NULL, 'k'},
      {"kokkos-threads", optional_argument, NULL, 'k'},
      {"kokkos-help", optional_argument, NULL, 'k'},
      {"kokkos-numa", optional_argument, NULL, 'k'},
      {NULL, 0, NULL, 0}};

  std::string copt;
  int c = 1;
  int opt_index;
  while ((c = getopt_long(argc, argv, "Vv:ctn:", long_options, &opt_index)) != -1) {
    switch (c) {
    case 'c':
      if (optarg)
        copt = std::string(optarg);
      else
        copt = "on";
      break;
    case 'n':
#ifdef HAVE_CUDA
      cArgs.numDevices = atoi(optarg);
#elif HAVE_OPENMP
      cArgs.numThreads = atoi(optarg);
#endif
      break;
    case 'v':
      if (optarg)
        cArgs.verbosity = atoi(optarg);
      else
        cArgs.verbosity = 3;
      break;
    case 'V':
      cArgs.versionFlag = 1;
      break;
    }
  }

  cArgs.colorFlag = false;
  if (copt.compare("on") == 0)
    cArgs.colorFlag = true;
  else if (copt.compare("auto") == 0 && isatty(fileno(stdout)))
      cArgs.colorFlag = true;

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

inputConfig::~inputConfig(){}

//int Log::verbosity;
//ansiColors *Log::c;
//int Log::rank;
//Kokkos::Timer *Log::timer;

void executeConfiguration(struct inputConfig &cf, struct commandArgs cargs){
  cf.colorFlag = cargs.colorFlag;
  cf.timeFormat = cargs.timeFormat;
  cf.verbosity = cargs.verbosity;

  luaReader L(cargs.fileName,"fiesta");

  // Required Parameters
  L.get({"title"}, cf.title);

  L.get({"time","nt"},cf.nt,0);
  L.get({"time","dt"},cf.dt);
  L.get({"time","start_index"},   cf.tstart,  0);
  L.get({"time","start_time"},     cf.time,    0.0);
  // Calculate time indices
  cf.t = cf.tstart;
  if(cf.nt==0){
    cf.tinterval = false;
    L.get({"time","tend"},cf.tend,0);
    if(cf.tend==0){
      Log::error("Must specify either 'fiesta.time.nt' or 'fiesta.time.tend'");
      exit(EXIT_FAILURE);
    }
    cf.nt=cf.tend-cf.tstart;
  }else{
    cf.tinterval=true;
    cf.tend = cf.tstart + cf.nt;
  }


  // Defaultable Parameters
  L.get({"eos","R"},        cf.R,       8.314462);
  //L.get("cequations"},      cf.ceq,     0);
  //L.get("noise"},    cf.noise,   0);
  L.get({"grid","ndim"},     cf.ndim,    2);
  L.get({"viscosity","enabled"},     cf.visc,    false);
  L.get({"ng"},       cf.ng,      3);
  L.get({"bc","xperiodic"},     cf.xPer,    false);
  L.get({"bc","yperiodic"},     cf.yPer,    false);

  L.get({"progress","frequency"},     cf.out_freq,    0);
  L.get({"write_frequency"},   cf.write_freq,  0);
  //L.get({"restart_frequency"}, cf.restart_freq,0);
  L.get({"status","frequency"},    cf.stat_freq,   0);

  L.get({"terrain","name"}, cf.terrainName, std::string("terrain.h5"));

  L.get({"restart","enabled"}, cf.restart, false);
  L.get({"restart","name"},    cf.restartName, std::string("restart-0000000.h5"));
  L.get({"restart","path"},    cf.pathName, std::string("."));
  L.get({"restart","reset"}, cf.restartReset, false);
  L.get({"restart","auto"}, cf.autoRestart, false);
  L.get({"restart","auto_name"}, cf.autoRestartName, std::string("restart-auto"));
  L.get({"restart","time_remaining"}, cf.restartTimeRemaining, 0);
  if(cf.autoRestart){
    cf.restart=true;
    cf.restartName=format("{}/{}.h5",cf.pathName,cf.autoRestartName);
    Log::debug("Auto Restart Name: '{}'",cf.restartName);
    if (!std::filesystem::exists(cf.restartName)){
      Log::warning("Auto Restart '{}' not found. Initiating fresh run.",cf.autoRestartName);
      cf.restart=false;
    }
    L.get({"restart","frequency"}, cf.restart_freq,cf.tend);
  }else{
    cf.autoRestartName="restart";
    L.get({"restart","frequency"}, cf.restart_freq,0);
  }

  std::string scheme, grid, mpi;
  L.get({"mpi","type"}, mpi, std::string("host"));
  L.get({"hdf5","chunk"}, cf.chunkable, false);
  L.get({"hdf5","compress"}, cf.compressible, false);

  if (!cf.chunkable && cf.compressible){
    cf.compressible = false;
    Log::warning("HDF5 Compression requires chunking.  Disabling.");
  }

  L.get({"advection_scheme"}, scheme, std::string("weno5"));
  L.get({"grid","type"},   grid, std::string("cartesian"));

  //vector<FSCAL> dx;
  //L.getArray("dx",dx,cf.ndim);
  L.get({"grid","dx"},cf.dxvec,cf.ndim);

  vector<size_t> ni;
  //L.getArray("ni",ni,cf.ndim);
  L.get({"grid","ni"},ni,cf.ndim);
  cf.glbl_nci = ni[0];
  cf.glbl_ncj = ni[1];
  if (cf.ndim == 3) 
    cf.glbl_nck = ni[2];
  else
    cf.glbl_nck = 1.0;

  vector<size_t> procs;
  //L.getArray("procs",procs,cf.ndim);
  L.get({"mpi","procs"},procs,cf.ndim);
  cf.xProcs=procs[0];
  cf.yProcs=procs[1];
  if (cf.ndim == 3) 
    cf.zProcs=procs[2];
  else
    cf.zProcs=1;

  std::string bcname;

  if (!cf.xPer){
    L.get({"bc","xmin"},bcname);
    cf.bcL=parseBC(bcname);
    L.get({"bc","xmax"},bcname);
    cf.bcR=parseBC(bcname);
  }
  if (!cf.yPer){
    L.get({"bc","ymin"},bcname);
    cf.bcB=parseBC(bcname);
    L.get({"bc","ymax"},bcname);
    cf.bcT=parseBC(bcname);
  }

  // Dependent Parameters
  //
  L.get({"buoyancy","enabled"}, cf.buoyancy, false);
  if (cf.buoyancy){
    L.get({"buoyancy","acceleration"},cf.gAccel,9.81);
    L.get({"buoyancy","rho_reference"},cf.rhoRef,0.0);
  }
  std::string test;
  L.get({"ceq","enabled"},cf.ceq,false);
  if(cf.ceq){
    L.get({"ceq","kappa"},cf.kap,1.23456);
    L.get({"ceq","epsilon"},cf.eps);
    L.get({"ceq","alpha"},cf.alpha);
    L.get({"ceq","beta"},cf.beta);
    L.get({"ceq","betae"},cf.betae);
    L.get({"ceq","st"},cf.st,0);
  }

  L.get({"noise","enabled"},cf.noise,false);
  if(cf.noise){
    L.get({"noise","dh"},cf.n_dh);
    L.get({"noise","eta"},cf.n_eta);
    L.get({"noise","coff"},cf.n_coff);
    L.get({"noise","nt"},cf.n_nt,1);
    L.get({"noise","mode"},cf.n_mode,1);
  }
  if (cf.ndim == 3) {
    L.get({"bc","zperiodic"},   cf.zPer,false);
    if(!cf.zPer){
      L.get({"bc","zmin"}, bcname);
      cf.bcH=parseBC(bcname);
      L.get({"bc","zmax"}, bcname);
      cf.bcF=parseBC(bcname);
    }
  }
  if (grid.compare("cartesian") == 0) {
    cf.grid = 0;
    cf.dx=cf.dxvec[0];
    cf.dy=cf.dxvec[1];
    if (cf.ndim == 3)
      cf.dz=cf.dxvec[2];
  }
  if (grid.compare("terrain") == 0) {
    if (cf.ndim == 3){
      cf.grid = 2;
      cf.dx=1.0;
      cf.dy=1.0;
      cf.dz = 1.0;
      L.get({"tdx"},cf.tdx);
      L.get({"tdy"},cf.tdy);
      L.get({"h"},cf.h);
    } else {
      printf("ndim must be equal to 3 for terrain");
      exit(EXIT_FAILURE);
    }
  }
  L.getSpeciesData(cf);

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

#ifdef HAVE_MPI
  // MPI halo exchange strategy from name
  if (mpi.compare("host") == 0)
    cf.mpiScheme = 1;
  else if (mpi.compare("gpu-aware") == 0)
    cf.mpiScheme = 2;
  else if (mpi.compare("gpu-type") == 0)
    cf.mpiScheme = 3;
  else if (mpi.compare("gpu-aware-ordered") == 0)
    cf.mpiScheme = 4;
  else if (mpi.compare("gpu-aware-unordered") == 0)
    cf.mpiScheme = 5;
  else if (mpi.compare("host-ordered") == 0)
    cf.mpiScheme = 6;
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


  // Set defaults for 2d grids
  if (cf.ndim == 2) {
    cf.nv = 3 + cf.ns;
    cf.dz = cf.dx;
    cf.zProcs = 1;
    cf.glbl_nck = 1;
    cf.bcH = BCType::outflow;
    cf.bcF = BCType::outflow;
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

  // initialize signal flags to inactive
  cf.restartFlag=0;
  cf.exitFlag=0;

  cf.xmp = 0;
  cf.ymp = 0;
  cf.zmp = 0;
#ifndef HAVE_MPI
  cf.globalGridDims.push_back(cf.glbl_ni);
  cf.globalGridDims.push_back(cf.glbl_nj);
  cf.globalCellDims.push_back(cf.glbl_nci);
  cf.globalCellDims.push_back(cf.glbl_ncj);
  cf.localGridDims.push_back(cf.ni);
  cf.localGridDims.push_back(cf.nj);
  cf.localCellDims.push_back(cf.nci);
  cf.localCellDims.push_back(cf.ncj);
  cf.subdomainOffset.push_back(cf.iStart);
  cf.subdomainOffset.push_back(cf.jStart);
  if(cf.ndim==3){
    cf.globalGridDims.push_back(cf.glbl_nk);
    cf.globalCellDims.push_back(cf.glbl_nck);
    cf.localGridDims.push_back(cf.nk);
    cf.localCellDims.push_back(cf.nck);
    cf.subdomainOffset.push_back(cf.kStart);
  }
#endif

  cf.ioThisStep = false;
}

int loadInitialConditions(struct inputConfig cf, FS4D &deviceV, FS4D &deviceG) {
  FSCAL x,y,z;
  FS4DH hostV = Kokkos::create_mirror_view(deviceV);
  FS4DH hostG = Kokkos::create_mirror_view(deviceG);
  if(cf.grid!=0){
    Kokkos::deep_copy(hostG,deviceG);
  }

  luaReader L(cf.inputFname,"fiesta");

  for (int v = 0; v < cf.nv; ++v) {
    if (cf.ndim == 3) {
      for (int k = cf.ng; k < cf.nck + cf.ng; ++k) {
        for (int j = cf.ng; j < cf.ncj + cf.ng; ++j) {
          for (int i = cf.ng; i < cf.nci + cf.ng; ++i) {
            if(cf.grid==0){
              // compute cell center coordintes for uniform grids
              x=cf.dx*(i-cf.ng+cf.subdomainOffset[0]) + 0.5*cf.dx;
              y=cf.dy*(j-cf.ng+cf.subdomainOffset[1]) + 0.5*cf.dy;
              z=cf.dz*(k-cf.ng+cf.subdomainOffset[2]) + 0.5*cf.dz;
            }else{
              // average cell nodes to compute cell center coordinate for non-uniform grids
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
            }
            hostV(i, j, k, v) = L.call("initial_conditions",4,x,y,z,(FSCAL)v);
          }
        }
      }
    } else {
      for (int j = cf.ng; j < cf.ncj + cf.ng; ++j) {
        for (int i = cf.ng; i < cf.nci + cf.ng; ++i) {
            if(cf.grid==0){
              // compute cell center coordintes for uniform grids
              x=cf.dx*(i-cf.ng+cf.subdomainOffset[0]) + 0.5*cf.dx;
              y=cf.dy*(j-cf.ng+cf.subdomainOffset[1]) + 0.5*cf.dy;
            }else{
              // average cell nodes to compute cell center coordinate for non-uniform grids
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
            }
            hostV(i, j, 0, v) = L.call("initial_conditions",4,x,y,0.0,(FSCAL)v);
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

    luaReader L(cf.inputFname,"fiesta");

    int ii,jj,kk;
    for (int v = 0; v < cf.ndim; ++v) {
      for (int k = 0; k < cf.nk; ++k) {
        for (int j = 0; j < cf.nj; ++j) {
          for (int i = 0; i < cf.ni; ++i) {
            ii = cf.iStart + i;
            jj = cf.jStart + j;
            kk = cf.kStart + k;
            hostV(i, j, k, v) = L.call("initialize_grid",4,(FSCAL)ii,(FSCAL)jj,(FSCAL)kk,(FSCAL)v);
          }
        }
      }
    }
    L.close();
    Kokkos::deep_copy(deviceV, hostV);
  } else if (cf.grid == 2) {
    cf.w->readTerrain(cf,deviceV);
  }


  return 0;
}
