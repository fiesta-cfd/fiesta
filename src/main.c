#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include "lua.h"
#include "lauxlib.h"
#include "lualib.h"
#include "cgnslib.h"
#include <mpi.h>

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
    MPI_Init(NULL,NULL);

    char buff[256];
    int ni,nj,nk,nci,ncj,nck;
    int idx;
    double gamma;
    double dx,dy,dz;

    int procsx,procsy,procsz;

    int index_file, icelldim, iphysdim, index_base, index_zone, index_coord, index_flow, index_field;
    cgsize_t isize[3][3];

    lua_State *L = luaL_newstate(); //Opens Lua
    //luaL_openlibs(L);               //opens the standard libraries

    if (luaL_loadfile(L,"input.lua") || lua_pcall(L,0,0,0))
        error(L, "Cannot run config file: %s\n", lua_tostring(L, -1));
    nci    = getglobint (L, "ni" );
    ncj    = getglobint (L, "nj" );
    nck    = getglobint (L, "nk" );
    dx     = getglobdbl (L, "dx" );
    dy     = getglobdbl (L, "dy" );
    dz     = getglobdbl (L, "dz" );
    gamma  = getglobdbl (L, "gamma" );
    procsx = getglobint (L, "procsx");
    procsy = getglobint (L, "procsy");
    procsz = getglobint (L, "procsz");

    ni = nci + 1;
    nj = ncj + 1;
    nk = nck + 1;

    int rank,numprocs;
    int dims[3] = {procsz,procsx,procsy};
    int periods[3] = {0,0,0};
    int coords[3];
    int left,right,top,bottom,front,back;

    MPI_Comm cartcomm;
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0){
        printf("Gamma = %4.2f\n",gamma);
        printf("ni = %d, dx = %f\n",ni,dx);
        printf("nj = %d, dy = %f\n",nj,dy);
        printf("nk = %d, dz = %f\n",nk,dz);
    }

    printf("I am %d of %d\n",rank,numprocs);

    double *x = malloc(ni*nj*nk*sizeof(double));
    double *y = malloc(ni*nj*nk*sizeof(double));
    double *z = malloc(ni*nj*nk*sizeof(double));
    double *v = malloc(ni*nj*nk*sizeof(double));


    for (int k=0; k<nk; ++k){
        for (int j=0; j<nj; ++j){
            for (int i=0; i<ni; ++i){
                idx = (ni*nj)*k+ni*j+i;
                x[idx] = i*dx;
                y[idx] = j*dy;
                z[idx] = k*dz;
                v[idx] = (double)((i+1)*(j+1)*(k+1))/(double)(ni*nj*nk);
                //printf("(%d,%d,%d): %f, %f, %f, %f\n",i,j,k,x[idx],y[idx],z[idx],v[idx]);
            }
        }
    }

    lua_close(L);

    if (rank == 0){
        if (cg_open("grid.cgns",CG_MODE_WRITE,&index_file)) cg_error_exit();

        cg_base_write(index_file,"Base",3,3,&index_base);
        
        isize[0][0] = ni;
        isize[0][1] = nj;
        isize[0][2] = nk;

        isize[1][0] = isize[0][0] - 1;
        isize[1][1] = isize[0][1] - 1;
        isize[1][2] = isize[0][2] - 1;
        
        isize[2][0] = 0;
        isize[2][1] = 0;
        isize[2][2] = 0;

        cg_zone_write(index_file,index_base,"Zone 1",*isize,CG_Structured,&index_zone);

        cg_coord_write(index_file, index_base, index_zone, CG_RealDouble, "CoordinateX",x,&index_coord);
        cg_coord_write(index_file, index_base, index_zone, CG_RealDouble, "CoordinateY",y,&index_coord);
        cg_coord_write(index_file, index_base, index_zone, CG_RealDouble, "CoordinateZ",z,&index_coord);

        cg_close(index_file);

        if (cg_open("grid.cgns",CG_MODE_MODIFY,&index_file)) cg_error_exit();

        cg_sol_write(index_file,index_base,index_zone,"Sollution 1",CG_Vertex,&index_flow);
        cg_field_write(index_file,index_base,index_zone,index_flow,CG_RealDouble,"Density",v,&index_field);

        cg_close(index_file);
    }

    MPI_Finalize();
    return 0;
}
