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
    int ndim;
    int glbl_ni,glbl_nj,glbl_nk;
    int glbl_nci,glbl_ncj,glbl_nck;
    double *gamma;
    double *M;
    double R;
    int visc;
    double dt,dx,dy,dz;
    int xProcs,yProcs,zProcs,numProcs;
    int xPlus,yPlus,zPlus;
    int xMinus,yMinus,zMinus;
    int rank;
    int nci,ncj,nck;
    int nt,ni,nj,nk,nv,ns;
    int iStart,jStart,kStart;
    int iEnd,jEnd,kEnd;
    MPI_Comm comm;
    int cF,cB,cZ;
    int ng,ngi,ngj,ngk;
    int xPer,yPer,zPer;
    int restart;
    const char * sfName;
    int tstart;
    double time;
    double ceq,kap,eps,alpha,beta,betae;
    int bcL,bcR,bcB,bcT,bcH,bcF;

    int out_freq, write_freq, restart_freq;
};

struct inputConfig executeConfiguration(char * fname);

int loadInitialConditions(struct inputConfig cf, const Kokkos::View<double****> v);

#endif
