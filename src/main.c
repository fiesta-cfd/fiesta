#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include "lua.h"
#include "lauxlib.h"
#include "lualib.h"

void error(lua_State *L, const char *fmt, ...){
    va_list argp;
    va_start(argp, fmt);
    vfprintf(stderr, fmt, argp);
    va_end(argp);
    lua_close(L);
    exit(EXIT_FAILURE);
}

int getglobbool(lua_State *L, const char *var) {
    int isnum, result;
    lua_getglobal(L, var);
    result = (int)lua_toboolean(L, -1);
    lua_pop(L,1);
    return result;
}

int getglobint(lua_State *L, const char *var) {
    int isnum, result;
    lua_getglobal(L, var);
    result = (int)lua_tointegerx(L, -1, &isnum);
    if (!isnum)
        error(L, "'%s' should be an integer\n", var);
    lua_pop(L,1);
    return result;
}

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

int main(){
    char buff[256];
    int ni,nj,nk,nci,ncj,nck;
    int idx;
    double gamma;
    double dx,dy,dz;
    lua_State *L = luaL_newstate(); //Opens Lua
    //luaL_openlibs(L);               //opens the standard libraries

    if (luaL_loadfile(L,"input.lua") || lua_pcall(L,0,0,0))
        error(L, "Cannot run config file: %s\n", lua_tostring(L, -1));
    nci   = getglobint (L, "ni" );
    ncj   = getglobint (L, "nj" );
    nck   = getglobint (L, "nk" );
    dx    = getglobdbl (L, "dx" );
    dy    = getglobdbl (L, "dy" );
    dz    = getglobdbl (L, "dz" );
    gamma = getglobdbl (L, "gamma" );

    ni = nci + 1;
    nj = ncj + 1;
    nk = nck + 1;

    printf("Gamma = %4.2f\n",gamma);
    printf("ni = %d, dx = %f\n",ni,dx);
    printf("nj = %d, dy = %f\n",nj,dy);
    printf("nk = %d, dz = %f\n",nk,dz);

    double *x = malloc(ni*nj*nk*sizeof(double));
    double *y = malloc(ni*nj*nk*sizeof(double));
    double *z = malloc(ni*nj*nk*sizeof(double));
    double *v = malloc(ni*nj*nk*sizeof(double));

    for (int k=1; k<nk; ++k){
        for (int j=1; j<nj; ++j){
            for (int i=1; i<ni; ++i){
                idx = (ni*nj)*k+ni*j+i;
                //printf("(%d,%d,%d), %d/%d\n",i,j,k,idx,ni*nj*nk);
                x[idx] = i*dx;
                y[idx] = j*dy;
                z[idx] = k*dz;
                v[idx] = i*j*k/(ni*nj*nk);
            }
        }
    }

    printf("%d : %f, %f, %f\n",ni*nj*nk-1,x[ni*nj*nk-1],y[ni*nj*nk-1],z[ni*nj*nk-1]);

    lua_close(L);
    return 0;
}
