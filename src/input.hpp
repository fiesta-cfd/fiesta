#ifndef INPUT_H
#define INPUT_H

#include "fiesta.hpp"
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include "lua.hpp"
#include "Kokkos_Core.hpp"
#ifndef NOMPI
#include <mpi.h>
#endif
#include <string>

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
    std::string inputFname;
    std::string title;
    int xmp,ymp,zmp;
    int ndim;
    int glbl_ni,glbl_nj,glbl_nk;
    int glbl_nci,glbl_ncj,glbl_nck;
    double *gamma;
    double *M;
    double *mu;
    double R;
    int scheme;
    int visc;
    int gravity;
    double g_accel;
    double* g_vec;
    double dt,dx,dy,dz;
    int xProcs,yProcs,zProcs,numProcs;
    int xPlus,yPlus,zPlus;
    int xMinus,yMinus,zMinus;
    int rank;
    int nci,ncj,nck;
    int nt,ni,nj,nk,nv,ns,nvt;
    int iStart,jStart,kStart;
    int iEnd,jEnd,kEnd;
#ifndef NOMPI
    MPI_Comm comm;
#endif
    int cF,cB,cZ;
    int ng,ngi,ngj,ngk;
    int xPer,yPer,zPer;
    int restart;
    const char * sfName;
    int tstart,tend;
    double time;
    double ceq,kap,eps,alpha,beta,betae;
    int bcL,bcR,bcB,bcT,bcH,bcF;
    double n_dh,n_coff,n_eta;
    int noise,n_nt;
    int particle, p_np;
    int t;
    int grid;

    int out_freq, stat_freq, write_freq, restart_freq;
};

void getCommandlineOptions(int argc, char** argv, int &vFlag, int &cFlag, std::string &fName);
struct inputConfig executeConfiguration(std::string fName);

int loadInitialConditions(struct inputConfig cf, const FS4D v);
int loadGrid(struct inputConfig cf, const FS4D v);

#endif
