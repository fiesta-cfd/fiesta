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

// Lua error function
void error(lua_State *L, const char *fmt, ...) {
  va_list argp;
  va_start(argp, fmt);
  vfprintf(stderr, fmt, argp);
  va_end(argp);
  lua_close(L);
  exit(EXIT_FAILURE);
}

// Lua get boolean value
int getglobbool(lua_State *L, const char *var, int defaultable, int default_value) {
  int result;
  int isnil;
  lua_getglobal(L, var);

  if (lua_isnoneornil(L,-1)){
    if (defaultable)
      return default_value;
    else
      error(L,"Error Reading Input File: Could not read value for '%s', a boolean value was expected.\n",var);
  }else{
    if ( lua_isstring(L,-1) ) {
      const char *iresult;
      iresult = lua_tostring(L,-1);

      std::string str(iresult);
      std::regex rgxtrue(R"(\.?(true|on|enable|enabled|1)\.?)",std::regex_constants::icase);
      std::regex rgxfalse(R"(\.?(false|off|disable|disabled|0)\.?)",std::regex_constants::icase);

      if ( std::regex_match(str,rgxtrue) ){
        result = 1;
      } else if ( std::regex_match(str,rgxfalse) ) {
        result = 0;
      } else {
        error(L,"Error Reading Input File: Unknown value for '%s', a boolean value was expected.\n",var);
      }
    } else {
      result = lua_toboolean(L,-1);
    }
  }
  lua_pop(L, 1);
  return result;
}

// Lua get integer value
int getglobint(lua_State *L, const char *var, int defaultable, int default_value) {
  int isnum, result;
  lua_getglobal(L, var);
  if (lua_isnoneornil(L, -1)){
    if (defaultable)
      return default_value;
    else
      error(L, "Error Reading Input File: Could not read value for '%s', an integer was expected.\n", var);
  }
  result = (int)lua_tointegerx(L, -1, &isnum);
  if (!isnum)
    error(L, "Error Reading Input File: Could not read value for '%s', an integer was expected.\n", var);
  lua_pop(L, 1);
  return result;
}

// Lua get double value
double getglobdbl(lua_State *L, const char *var, int defaultable, double default_value) {
  int isnum;
  double result;
  lua_getglobal(L, var);
  if (lua_isnoneornil(L, -1)){
    if (defaultable)
      return default_value;
    else
      error(L, "Error Reading Input File: Could not read value for '%s', an integer was expected.\n", var);
  }
  result = (double)lua_tonumberx(L, -1, &isnum);
  if (!isnum)
    error(L, "Error Reading Input File: Could not read value for '%s', an integer was expected.\n", var);
  lua_pop(L, 1);
  return result;
}

// Lua get string
std::string getglobstr(lua_State *L, const char *var, int defaultable, std::string default_value) {
  // size_t len;
  const char *result;
  result = (char *)malloc(32 * sizeof(char));
  lua_getglobal(L, var);
  if (lua_isnoneornil(L, -1)){
    if (defaultable)
      return default_value;
    else
      error(L, "Error Reading Input File: Could not read value for '%s', an string was expected.\n", var);
  }
  result = lua_tostring(L, -1);
  return std::string(result);
}

struct commandArgs getCommandlineOptions(int argc, char **argv){

  // create command argumet structure
  struct commandArgs cArgs;

  cArgs.fileName = "fiesta.lua";
  cArgs.timeFormat = 0;
  cArgs.versionFlag = 0;
  cArgs.colorFlag = 0;

  // create options
  static struct option long_options[] = {
      {"version", no_argument, NULL, 'v'},
      {"color", optional_argument, NULL, 'c'},
      {"decimal-time", optional_argument, NULL, 't'},
      {"kokkos-ndevices", optional_argument, NULL, 'n'},
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
        cArgs.timeFormat = 2;
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

  cf.inputFname = cargs.fileName;

  /* Create Lua State */
  lua_State *L = luaL_newstate(); // Opens Lua
  luaL_openlibs(L);               // opens the standard libraries

  /* Open and run Lua configuration file */
  if (luaL_loadfile(L, cf.inputFname.c_str()) || lua_pcall(L, 0, 0, 0))
    error(L, "Cannot run config file: %s\n", lua_tostring(L, -1));

  /* get global variables from Lua results */
  cf.ndim = getglobint(L, "ndim",1,2);
  cf.nt = getglobint(L, "nt",0,0);
  cf.glbl_nci = getglobint(L, "ni",0,0);
  cf.glbl_ncj = getglobint(L, "nj",0,0);
  cf.ns = getglobint(L, "ns",0,0);
  cf.dt = getglobdbl(L, "dt",0,0.0);
  cf.R = getglobdbl(L, "R",0,0.0);
  cf.visc = getglobbool(L, "visc",1,0);
  cf.xProcs = getglobint(L, "procsx",1,1);
  cf.yProcs = getglobint(L, "procsy",1,1);
  cf.ng = 3;
  cf.xPer = getglobint(L, "xPer",1,0);
  cf.yPer = getglobint(L, "yPer",1,0);
  cf.bcL = getglobint(L, "bcXmin",1,0);
  cf.bcR = getglobint(L, "bcXmax",1,0);
  cf.bcB = getglobint(L, "bcYmin",1,0);
  cf.bcT = getglobint(L, "bcYmax",1,0);
  cf.restart = getglobbool(L, "restart",1,0);
  cf.restartName = getglobstr(L, "restartName",1,"restart-0000000.h5");
  std::string scheme = getglobstr(L, "scheme",1,"weno5");
  cf.title = getglobstr(L, "title",0,"none");
  std::string grid = getglobstr(L, "grid",1,"cartesian");
  cf.tstart = getglobint(L, "tstart",1,0);
  cf.t = cf.tstart;
  cf.time = getglobdbl(L, "time",1,0.0);
  cf.ceq = getglobbool(L, "ceq",1,0);
  if (cf.ceq == 1) {
    cf.kap = getglobdbl(L, "kappa",0,0.0);
    cf.eps = getglobdbl(L, "epsilon",0,0.0);
    cf.alpha = getglobdbl(L, "alpha",0,0.0);
    cf.beta = getglobdbl(L, "beta",0,0.0);
    cf.betae = getglobdbl(L, "betae",0,0.0);
  }
  cf.noise = getglobbool(L, "noise",1,0);
  if (cf.noise == 1) {
    cf.n_dh = getglobdbl(L, "n_dh",0,0.0);
    cf.n_eta = getglobdbl(L, "n_eta",0,0.0);
    cf.n_coff = getglobdbl(L, "n_coff",0,0.0);
    cf.n_nt = getglobint(L, "n_nt",1,1);
  }
  if (cf.ndim == 3) {
    cf.glbl_nck = getglobint(L, "nk",0,0);
    cf.zProcs = getglobint(L, "procsz",1,1);
    cf.zPer = getglobint(L, "zPer",1,0);
    cf.bcH = getglobint(L, "bcZmin",1,0);
    cf.bcF = getglobint(L, "bcZmax",1,0);
  }
  cf.particle = getglobbool(L, "particle",1,0);
  if (cf.particle == 1) {
    cf.p_np = getglobint(L, "p_np",0,0);
  } else {
    cf.p_np = 0;
  }

  cf.gamma = (double *)malloc(cf.ns * sizeof(double));
  cf.M = (double *)malloc(cf.ns * sizeof(double));
  cf.g_vec = (double *)malloc(3 * sizeof(double));
  cf.tend = cf.tstart + cf.nt;


  cf.scheme = 1;
  if (scheme.compare("weno5") == 0)
    cf.scheme = 1;
  if (scheme.compare("centered4") == 0)
    cf.scheme = 2;
  if (scheme.compare("quick") == 0)
    cf.scheme = 3;

  if (grid.compare("generalized") == 0) {
    cf.grid = 1;
    cf.dx = 1.0;
    cf.dy = 1.0;
    if (cf.ndim == 3)
      cf.dz = 1.0;
  }
  if (grid.compare("cartesian") == 0) {
    cf.grid = 0;
    cf.dx = getglobdbl(L, "dx",0,0.0);
    cf.dy = getglobdbl(L, "dy",0,0.0);
    if (cf.ndim == 3)
      cf.dz = getglobdbl(L, "dz",0,0.0);
  }

  int isnum;

  lua_getglobal(L, "gamma");
  for (int s = 0; s < cf.ns; ++s) {
    lua_pushnumber(L, s + 1);
    lua_gettable(L, -2);
    cf.gamma[s] = (double)lua_tonumberx(L, -1, &isnum);
    lua_pop(L, 1);
  }
  lua_getglobal(L, "M");
  for (int s = 0; s < cf.ns; ++s) {
    lua_pushnumber(L, s + 1);
    lua_gettable(L, -2);
    cf.M[s] = (double)lua_tonumberx(L, -1, &isnum);
    lua_pop(L, 1);
  }
  if (cf.visc == 1) {
    cf.mu = (double *)malloc(cf.ns * sizeof(double));
    lua_getglobal(L, "mu");
    for (int s = 0; s < cf.ns; ++s) {
      lua_pushnumber(L, s + 1);
      lua_gettable(L, -2);
      cf.mu[s] = (double)lua_tonumberx(L, -1, &isnum);
      lua_pop(L, 1);
    }
  } else {
    cf.mu = (double *)malloc(cf.ns * sizeof(double));
    for (int s = 0; s < cf.ns; ++s)
      cf.mu[s] = 0.0;
  }

  cf.gravity = getglobbool(L, "buoyancy",1,0);
  // cf.gravity = 0;
  // if (cf.gravity == 1){
  //     cf.g_accel     = getglobdbl (L, "g_accel",1,9.8);
  //     lua_getglobal(L, "g_vec");
  //     for (int d=0; d<3; ++d){
  //         lua_pushnumber(L, d+1);
  //         lua_gettable(L, -2);
  //         cf.g_vec[d] = (double)lua_tonumberx(L, -1, &isnum);
  //         lua_pop(L,1);
  //     }
  // }

  cf.out_freq = getglobint(L, "out_freq",1,0);
  cf.write_freq = getglobint(L, "write_freq",1,0);
  cf.restart_freq = getglobint(L, "restart_freq",1,0);
  cf.stat_freq = getglobint(L, "stat_freq",1,0);

  cf.nv = 4 + cf.ns;

  /* Done with Lua */
  lua_close(L);

  if (cf.ndim == 2) {
    cf.nv = 3 + cf.ns;
    cf.dz = cf.dx;
    cf.zProcs = 1;
    cf.glbl_nck = 1;
    cf.bcH = 0;
    cf.bcF = 0;
    cf.zPer = 0;
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

int loadInitialConditions(struct inputConfig cf, const FS4D deviceV) {

  double z;
  int isnum;

  lua_State *L = luaL_newstate(); // Opens Lua
  luaL_openlibs(L);               // opens the standard libraries

  /* Open and run Lua configuration file */
  if (luaL_loadfile(L, cf.inputFname.c_str()) || lua_pcall(L, 0, 0, 0))
    error(L, "Cannot run config file: %s\n", lua_tostring(L, -1));

  FS4DH hostV = Kokkos::create_mirror_view(deviceV);
  for (int v = 0; v < cf.nv; ++v) {
    if (cf.ndim == 3) {
      for (int k = cf.ng; k < cf.nck + cf.ng; ++k) {
        for (int j = cf.ng; j < cf.ncj + cf.ng; ++j) {
          for (int i = cf.ng; i < cf.nci + cf.ng; ++i) {
            int ii = i - cf.ng;
            int jj = j - cf.ng;
            int kk = k - cf.ng;
            lua_getglobal(L, "f");
            lua_pushnumber(L, cf.iStart + ii);
            lua_pushnumber(L, cf.jStart + jj);
            lua_pushnumber(L, cf.kStart + kk);
            lua_pushnumber(L, v);
            if (lua_pcall(L, 4, 1, 0) != LUA_OK)
              error(L, "error running function 'f': %s\n", lua_tostring(L, -1));
            z = lua_tonumberx(L, -1, &isnum);
            if (!isnum)
              error(L, "function 'f' should return a number");
            lua_pop(L, 1);

            hostV(i, j, k, v) = z;
          }
        }
      }
    } else {
      for (int j = cf.ng; j < cf.ncj + cf.ng; ++j) {
        for (int i = cf.ng; i < cf.nci + cf.ng; ++i) {
          int ii = i - cf.ng;
          int jj = j - cf.ng;
          lua_getglobal(L, "f");
          lua_pushnumber(L, cf.iStart + ii);
          lua_pushnumber(L, cf.jStart + jj);
          lua_pushnumber(L, cf.kStart);
          lua_pushnumber(L, v);
          if (lua_pcall(L, 4, 1, 0) != LUA_OK)
            error(L, "error running function 'f': %s\n", lua_tostring(L, -1));
          z = lua_tonumberx(L, -1, &isnum);
          if (!isnum)
            error(L, "function 'f' should return a number");
          lua_pop(L, 1);

          hostV(i, j, 0, v) = z;
        }
      }
    }
  }

  lua_close(L);

  Kokkos::deep_copy(deviceV, hostV);

  return 0;
}

int loadGrid(struct inputConfig cf, const FS4D deviceV) {

  FS4DH hostV = Kokkos::create_mirror_view(deviceV);

  if (cf.grid == 1) {

    double z;
    int isnum;

    lua_State *L = luaL_newstate(); // Opens Lua
    luaL_openlibs(L);               // opens the standard libraries

    /* Open and run Lua configuration file */
    if (luaL_loadfile(L, cf.inputFname.c_str()) || lua_pcall(L, 0, 0, 0))
      error(L, "Cannot run config file: %s\n", lua_tostring(L, -1));

    for (int v = 0; v < cf.ndim; ++v) {
      for (int k = 0; k < cf.nk; ++k) {
        for (int j = 0; j < cf.nj; ++j) {
          for (int i = 0; i < cf.ni; ++i) {
            lua_getglobal(L, "g");
            lua_pushnumber(L, cf.iStart + i);
            lua_pushnumber(L, cf.jStart + j);
            lua_pushnumber(L, cf.kStart + k);
            lua_pushnumber(L, v);
            if (lua_pcall(L, 4, 1, 0) != LUA_OK)
              error(L, "error running function 'f': %s\n", lua_tostring(L, -1));
            z = lua_tonumberx(L, -1, &isnum);
            if (!isnum)
              error(L, "function 'f' should return a number");
            lua_pop(L, 1);

            hostV(i, j, k, v) = z;
          }
        }
      }
    }

    lua_close(L);

  } else {
    for (int k = 0; k < cf.nk; ++k) {
      for (int j = 0; j < cf.nj; ++j) {
        for (int i = 0; i < cf.ni; ++i) {
          hostV(i, j, k, 0) = (cf.iStart + i) * cf.dx;
          hostV(i, j, k, 1) = (cf.jStart + j) * cf.dy;
          hostV(i, j, k, 2) = (cf.kStart + k) * cf.dz;
          // printf("(%d, %d, %d) : (%.1f, %.1f,
          // %.1f)\n",i,j,k,hostV(i,j,k,0),hostV(i,j,k,1),hostV(i,j,k,2));
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
    int isnum;

    lua_State *L = luaL_newstate(); // Opens Lua
    luaL_openlibs(L);               // opens the standard libraries

    /* Open and run Lua configuration file */
    if (luaL_loadfile(L, cf.inputFname.c_str()) || lua_pcall(L, 0, 0, 0))
      error(L, "Cannot run config file: %s\n", lua_tostring(L, -1));
    for (int p = 0; p < cf.p_np; ++p) {
      for (int v = 0; v < cf.ndim; ++v) {
        lua_getglobal(L, "p");
        lua_pushnumber(L, p);
        lua_pushnumber(L, v);
        if (lua_pcall(L, 2, 1, 0) != LUA_OK)
          error(L, "error running function 'f': %s\n", lua_tostring(L, -1));
        z = lua_tonumberx(L, -1, &isnum);
        if (!isnum)
          error(L, "function 'p' should return a number");
        lua_pop(L, 1);

        particlesH(p).state = 1;
        if (v == 0)
          particlesH(p).x = z;
        if (v == 1)
          particlesH(p).y = z;
      }
    }

    Kokkos::deep_copy(deviceV, particlesH);
  }

  return 0;
}
