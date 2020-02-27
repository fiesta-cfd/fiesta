#include "input.hpp"
#include "Kokkos_Core.hpp"
#include "debug.hpp"
#include "string.h"

// Lua error function
void error(lua_State *L, const char *fmt, ...){
    va_list argp;
    va_start(argp, fmt);
    vfprintf(stderr, fmt, argp);
    va_end(argp);
    lua_close(L);
    exit(EXIT_FAILURE);
}

//Lua get boolean value
int getglobbool(lua_State *L, const char *var) {
    int result;
    lua_getglobal(L, var);
    result = (int)lua_toboolean(L, -1);
    lua_pop(L,1);
    return result;
}

//Lua get integer value
int getglobint(lua_State *L, const char *var) {
    int isnum, result;
    lua_getglobal(L, var);
    result = (int)lua_tointegerx(L, -1, &isnum);
    if (!isnum)
        error(L, "'%s' should be an integer\n", var);
    lua_pop(L,1);
    return result;
}

//Lua get double value
double getglobdbl(lua_State *L, const char *var) {
    int isnum;
    double result;
    lua_getglobal(L, var);
    result = (double)lua_tonumberx(L, -1, &isnum);
    if (!isnum)
        error(L, "'%s' should be a double\n", var);
    lua_pop(L,1);
    return result;
}

//Lua get string
char * getglobstr(lua_State *L, const char * var){
    //size_t len;
    const char *iresult;
    char *result;
    result = (char *)malloc(32*sizeof(char));
    lua_getglobal(L, var);
    iresult = lua_tostring(L, -1);
    strcpy(result,iresult);
    //result = lua_tolstring(L, -1, &len);
    return result;
}

struct inputConfig executeConfiguration(char * fname){

    struct inputConfig myConfig;
    /* Create Lua State */
    lua_State *L = luaL_newstate(); //Opens Lua
    luaL_openlibs(L);               //opens the standard libraries

    /* Open and run Lua configuration file */
    if (luaL_loadfile(L,fname) || lua_pcall(L,0,0,0))
        error(L, "Cannot run config file: %s\n", lua_tostring(L, -1));

    /* get global variables from Lua results */
    myConfig.ndim        = getglobint (L, "ndim" );
    myConfig.nt          = getglobint (L, "nt" );
    myConfig.glbl_nci    = getglobint (L, "ni" );
    myConfig.glbl_ncj    = getglobint (L, "nj" );
    myConfig.ns          = getglobint (L, "ns" );
    myConfig.dt          = getglobdbl (L, "dt" );
    myConfig.dx          = getglobdbl (L, "dx" );
    myConfig.dy          = getglobdbl (L, "dy" );
    myConfig.R           = getglobdbl (L, "R" );
    myConfig.visc        = getglobdbl (L, "visc" );
    myConfig.xProcs      = getglobint (L, "procsx");
    myConfig.yProcs      = getglobint (L, "procsy");
    //myConfig.ng          = getglobint (L, "ng" );
    myConfig.ng = 3;
    myConfig.xPer        = getglobint (L, "xPer" );
    myConfig.yPer        = getglobint (L, "yPer" );
    myConfig.bcL         = getglobint (L, "bcXmin" );
    myConfig.bcR         = getglobint (L, "bcXmax" );
    myConfig.bcB         = getglobint (L, "bcYmin" );
    myConfig.bcT         = getglobint (L, "bcYmax" );
    myConfig.restart     = getglobint (L, "restart");
    myConfig.sfName      = getglobstr (L, "restartName");
    myConfig.tstart      = getglobint (L, "tstart");
    myConfig.time        = getglobdbl (L, "time");
    myConfig.ceq         = getglobdbl (L, "ceq");
    myConfig.kap         = getglobdbl (L, "kappa");
    myConfig.eps         = getglobdbl (L, "epsilon");
    myConfig.alpha       = getglobdbl (L, "alpha");
    myConfig.beta        = getglobdbl (L, "beta");
    myConfig.betae       = getglobdbl (L, "betae");
    if (myConfig.ndim == 3){
        myConfig.glbl_nck    = getglobint (L, "nk" );
        myConfig.dz          = getglobdbl (L, "dz" );
        myConfig.zProcs      = getglobint (L, "procsz");
        myConfig.zPer        = getglobint (L, "zPer" );
        myConfig.bcH         = getglobint (L, "bcZmin" );
        myConfig.bcF         = getglobint (L, "bcZmax" );
    }

    myConfig.gamma = (double*)malloc(myConfig.ns*sizeof(double));
    myConfig.M = (double*)malloc(myConfig.ns*sizeof(double));
    myConfig.tend = myConfig.tstart+myConfig.nt;

    int isnum;

    lua_getglobal(L, "gamma");
    for (int s=0; s<myConfig.ns; ++s){
        lua_pushnumber(L, s+1);
        lua_gettable(L, -2);
        myConfig.gamma[s] = (double)lua_tonumberx(L, -1, &isnum);
        lua_pop(L,1);
    }
    lua_getglobal(L, "M");
    for (int s=0; s<myConfig.ns; ++s){
        lua_pushnumber(L, s+1);
        lua_gettable(L, -2);
        myConfig.M[s] = (double)lua_tonumberx(L, -1, &isnum);
        lua_pop(L,1);
    }

    myConfig.out_freq    = getglobint (L, "out_freq");
    myConfig.write_freq    = getglobint (L, "write_freq");
    myConfig.restart_freq    = getglobint (L, "restart_freq");

    myConfig.nv = 4 + myConfig.ns;

    snprintf(myConfig.inputFname,32,"%s",fname);

    /* Done with Lua */
    lua_close(L);

    if (myConfig.ndim == 2){
        myConfig.nv = 3 + myConfig.ns;
        myConfig.dz = myConfig.dx;
        myConfig.zProcs = 1;
        myConfig.glbl_nck = 1;
        myConfig.bcH = 0;
        myConfig.bcF = 0;
        myConfig.zPer = 0;
    }

    /* calculate number of nodes from number of cells */
    myConfig.glbl_ni = myConfig.glbl_nci + 1;
    myConfig.glbl_nj = myConfig.glbl_ncj + 1;
    myConfig.glbl_nk = myConfig.glbl_nck + 1;

    return myConfig;
}

int loadInitialConditions(struct inputConfig cf, const Kokkos::View<double****> deviceV){
    
    double z;
    int isnum;
    
    lua_State *L = luaL_newstate(); //Opens Lua
    luaL_openlibs(L);               //opens the standard libraries

    /* Open and run Lua configuration file */
    if (luaL_loadfile(L,cf.inputFname) || lua_pcall(L,0,0,0))
        error(L, "Cannot run config file: %s\n", lua_tostring(L, -1));

    typename Kokkos::View<double****>::HostMirror hostV = Kokkos::create_mirror_view(deviceV);
    for (int v=0; v<cf.nv; ++v){
        if (cf.ndim == 3){
            for (int k=cf.ng; k<cf.nck+cf.ng; ++k){
                for (int j=cf.ng; j<cf.ncj+cf.ng; ++j){
                    for (int i=cf.ng; i<cf.nci+cf.ng; ++i){
                        int ii = i - cf.ng;
                        int jj = j - cf.ng;
                        int kk = k - cf.ng;
                        lua_getglobal(L,"f");
                        lua_pushnumber(L,cf.iStart+ii);
                        lua_pushnumber(L,cf.jStart+jj);
                        lua_pushnumber(L,cf.kStart+kk);
                        lua_pushnumber(L,v);
                        if (lua_pcall(L,4,1,0) != LUA_OK)
                            error(L, "error running function 'f': %s\n",lua_tostring(L, -1));
                        z = lua_tonumberx(L,-1,&isnum);
                        if (!isnum)
                            error(L, "function 'f' should return a number");
                        lua_pop(L,1);
                        
                        hostV(i,j,k,v) = z;
                    }
                }
            }
        }else{
            for (int j=cf.ng; j<cf.ncj+cf.ng; ++j){
                for (int i=cf.ng; i<cf.nci+cf.ng; ++i){
                    int ii = i - cf.ng;
                    int jj = j - cf.ng;
                    lua_getglobal(L,"f");
                    lua_pushnumber(L,cf.iStart+ii);
                    lua_pushnumber(L,cf.jStart+jj);
                    lua_pushnumber(L,cf.kStart);
                    lua_pushnumber(L,v);
                    if (lua_pcall(L,4,1,0) != LUA_OK)
                        error(L, "error running function 'f': %s\n",lua_tostring(L, -1));
                    z = lua_tonumberx(L,-1,&isnum);
                    if (!isnum)
                        error(L, "function 'f' should return a number");
                    lua_pop(L,1);
                    
                    hostV(i,j,0,v) = z;
                }
            }
        }
    }
    
    lua_close(L);
    
    Kokkos::deep_copy(deviceV,hostV);
    //printf("cells: %f, %f, %f\n",hostV(2001,3,0,0),hostV(2002,3,0,0),hostV(2003,3,0,0));
    
    return 0;
}
