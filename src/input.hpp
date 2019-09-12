#ifndef INPUT_H
#define INPUT_H

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include "lua.hpp"
#include "Kokkos_Core.hpp"
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
    char inputFname[32];
    int glbl_ni,glbl_nj,glbl_nk;
    int glbl_nci,glbl_ncj,glbl_nck;
    double gamma;
    double dt,dx,dy,dz;
    int xProcs,yProcs,zProcs,numProcs;
    int xPlus,yPlus,zPlus;
    int xMinus,yMinus,zMinus;
    int rank;
    int nci,ncj,nck;
    int nt,ni,nj,nk,nv;
    int iStart,jStart,kStart;
    int iEnd,jEnd,kEnd;
    MPI_Comm comm;
    int cF,cB,cZ;
    int ng,ngi,ngj,ngk;

    int out_freq;
};

struct inputConfig executeConfiguration(char * fname);

int loadInitialConditions(struct inputConfig cf, const Kokkos::View<double****> v);

#endif
