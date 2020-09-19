#include "input.hpp"
#include "Kokkos_Core.hpp"
#include "debug.hpp"
#include "fiesta.hpp"
#include "particle.hpp"
#include "string.h"
#include <getopt.h>
#include <iostream>
#include <string>
#include <regex>
#include "luaReader.hpp"

struct commandArgs getCommandlineOptions(int argc, char **argv){

  // create command argumet structure
  struct commandArgs cArgs;

  cArgs.fileName = "fiesta.lua";
  cArgs.timeFormat = 2;
  cArgs.versionFlag = 0;
  cArgs.colorFlag = 0;

  // create options
  static struct option long_options[] = {
      {"version", no_argument, NULL, 'v'},
      {"color", optional_argument, NULL, 'c'},
      {"time-format", optional_argument, NULL, 't'},
      {"kokkos-ndevices", optional_argument, NULL, 'n'},
      {"kokkos-num-devices", optional_argument, NULL, 'n'},
      {"kokkos-threads", optional_argument, NULL, 'n'},
      {"kokkos-help", optional_argument, NULL, 'n'},
      {"kokkos-numa", optional_argument, NULL, 'n'},
      {NULL, 0, NULL, 0}};

  std::string copt;
  int c = 1;
  int opt_index;
  while ((c = getopt_long(argc, argv, "", long_options, &opt_index)) != -1) {
    switch (c) {
    case 'c':
      if (optarg)
        copt = std::string(optarg);
      else
        copt = "auto";
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

  return cArgs;
}

struct inputConfig executeConfiguration(struct commandArgs cargs) {
  struct inputConfig cf;

  luaReader L(cargs.fileName);

  // Required Parameters
  cf.nt       = L.getDouble("nt");
  cf.glbl_nci = L.getInt("ni");
  cf.glbl_ncj = L.getInt("nj");
  cf.ns       = L.getInt("ns");
  cf.dt       = L.getDouble("dt");
  cf.R        = L.getDouble("R");
  cf.title    = L.getString("title");
  L.getStrings("species_names",cf.ns,cf.speciesName);

  // Defaultable Parameters
  cf.tstart   = L.getInt("tstart",0);
  cf.time     = L.getDouble("time",0.0);
  cf.ceq      = L.getBool("ceq",0);
  cf.particle = L.getBool("particle",0);
  cf.noise    = L.getBool("noise",0);
  cf.gravity  = L.getBool("buoyancy",0);
  cf.ndim     = L.getDouble("ndim",2);
  cf.visc     = L.getBool("visc",0);
  cf.xProcs   = L.getInt("procsx",1);
  cf.yProcs   = L.getInt("procsy",1);
  cf.ng       = L.getInt("ng",3);
  cf.xPer     = L.getInt("xPer",0);
  cf.yPer     = L.getInt("yPer",0);
  cf.bcL      = L.getInt("bcXmin",0);
  cf.bcR      = L.getInt("bcXmax",0);
  cf.bcB      = L.getInt("bcYmin",0);
  cf.bcT      = L.getInt("bcYmax",0);
  cf.restart  = L.getBool("restart",0);
  cf.out_freq        = L.getInt("out_freq",0);
  cf.write_freq      = L.getInt("write_freq",0);
  cf.restart_freq    = L.getInt("restart_freq",0);
  cf.stat_freq       = L.getInt("stat_freq",0);
  cf.restartName     = L.getString("restartName","restart-0000000.h5");
  std::string scheme = L.getString("scheme","weno5");
  std::string grid   = L.getString("grid","cartesian");


  // Dependent Parameters
  if (cf.ceq == 1) {
    cf.kap = L.getDouble("kappa");
    cf.eps = L.getDouble("epsilon");
    cf.alpha = L.getDouble("alpha");
    cf.beta = L.getDouble("beta");
    cf.betae = L.getDouble("betae");
    cf.st = L.getInt("st",10);
  }
  if (cf.noise == 1) {
    cf.n_dh = L.getDouble("n_dh");
    cf.n_eta = L.getDouble("n_eta");
    cf.n_coff = L.getDouble("n_coff");
    cf.n_nt = L.getInt("n_nt",1);
  }
  if (cf.ndim == 3) {
    cf.glbl_nck = L.getInt("nk");
    cf.zProcs = L.getInt("procsz",1);
    cf.zPer = L.getInt("zPer",0);
    cf.bcH = L.getInt("bcZmin",0);
    cf.bcF = L.getInt("bcZmax",0);
  }
  if (cf.particle == 1) {
    cf.p_np = L.getInt("p_np");
  } else {
    cf.p_np = 0;
  }
  if (grid.compare("cartesian") == 0) {
    cf.grid = 0;
    cf.dx = L.getDouble("dx");
    cf.dy = L.getDouble("dy");
    if (cf.ndim == 3)
      cf.dz = L.getDouble("dz");
  }

  // Dependent Array Parameters
  cf.M = (double *)malloc(cf.ns * sizeof(double));
  cf.gamma = (double *)malloc(cf.ns * sizeof(double));
  cf.mu = (double *)malloc(cf.ns * sizeof(double));
  L.getDoubles("gamma",cf.ns,cf.gamma);
  L.getDoubles("M",cf.ns,cf.M);
  L.getDoubles("mu",cf.ns,cf.mu);

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

  // Set Grid Options from Grid String
  if (grid.compare("generalized") == 0) {
    cf.grid = 1;
    cf.dx = 1.0;
    cf.dy = 1.0;
    if (cf.ndim == 3)
      cf.dz = 1.0;
  }

  // Calculate time indices
  cf.t        = cf.tstart;
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
  }

  // Calculate total number of flow variables
  cf.nv = 4 + cf.ns;

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

int loadInitialConditions(struct inputConfig cf, const FS4D deviceV) {

  int ii,jj,kk;
  FS4DH hostV = Kokkos::create_mirror_view(deviceV);

  luaReader L(cf.inputFname);

  for (int v = 0; v < cf.nv; ++v) {
    if (cf.ndim == 3) {
      for (int k = cf.ng; k < cf.nck + cf.ng; ++k) {
        for (int j = cf.ng; j < cf.ncj + cf.ng; ++j) {
          for (int i = cf.ng; i < cf.nci + cf.ng; ++i) {
            ii = i - cf.ng;
            jj = j - cf.ng;
            kk = k - cf.ng;
            hostV(i, j, k, v) = L.call("f",4,ii,jj,kk,v);
          }
        }
      }
    } else {
      for (int j = cf.ng; j < cf.ncj + cf.ng; ++j) {
        for (int i = cf.ng; i < cf.nci + cf.ng; ++i) {
            ii = i - cf.ng;
            jj = j - cf.ng;
            kk = 0;
            hostV(i, j, 0, v) = L.call("f",4,ii,jj,kk,v);
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

int loadGrid(struct inputConfig cf, const FS4D deviceV) {

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
            hostV(i, j, k, v) = L.call("g",4,ii,jj,kk,v);
          }
        }
      }
    }
    L.close();
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
  }
  Kokkos::deep_copy(deviceV, hostV);

  return 0;
}
int loadParticles(struct inputConfig cf, const FSP2D deviceV) {

  if (cf.particle == 1) {
    FSP2DH particlesH = Kokkos::create_mirror_view(deviceV);

    double z;
    luaReader L(cf.inputFname);

    for (int p = 0; p < cf.p_np; ++p) {
      for (int v = 0; v < cf.ndim+2; ++v) {
        z=L.call("p",2,p,v);

        particlesH(p).state = 1;
        particlesH(p).u = 0.0;
        particlesH(p).v = 0.0;
        if (v == 0)
          particlesH(p).x = z;
        if (v == 1)
          particlesH(p).y = z;
        if (v == 2)
          particlesH(p).r = z;
        if (v == 3)
          particlesH(p).m = z;
      }
    }

    L.close();

    Kokkos::deep_copy(deviceV, particlesH);
  }

  return 0;
}
