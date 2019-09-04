#ifndef INPUT_H
#define INPUT_H

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
extern "C" {
#include "lua.h"
#include "lauxlib.h"
#include "lualib.h"
}
#include <mpi.h>

// Lua error function
void error(lua_State *L, const char *fmt, ...);

//Lua get boolean value
int getglobbool(lua_State *L, const char *var);

//Lua get integer value
int getglobint(lua_State *L, const char *var);

//Lua get double value
double getglobdbl(lua_State *L, const char *var);

//configuration structure
struct inputConfig {
    int glbl_ni,glbl_nj,glbl_nk;
    int glbl_nci,glbl_ncj,glbl_nck;
    double gamma;
    double dx,dy,dz;
    int xProcs,yProcs,zProcs,numProcs;
    int xPlus,yPlus,zPlus;
    int xMinus,yMinus,zMinus;
    int rank;
    int nci,ncj,nck;
    int ni,nj,nk;
    int iStart,jStart,kStart;
    int iEnd,jEnd,kEnd;
    MPI_Comm comm;
    int cF,cB,cZ;
};

struct inputConfig executeConfiguration(char * fname);

#endif