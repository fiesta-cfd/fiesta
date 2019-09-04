#include "input.hpp"
#include "mpi_init.hpp"
#include "cgns.hpp"
#include "Kokkos_Core.hpp"
#include <mpi.h>

#define MYDBG printf("%s:%d\n",__FILE__,__LINE__);
#define MYDBG0 if(rank==0) printf("%s:%d\n",__FILE__,__LINE__);

void fnExit1(void){
    Kokkos::finalize();
}

int main(int argc, char* argv[]){
    MPI_Init(NULL,NULL);
    Kokkos::initialize(argc, argv);

    atexit(fnExit1);

    int idx;

    struct inputConfig cf;

    //cf = executeConfiguration("input.lua");
    cf = executeConfiguration(argv[1]);

    cf = mpi_init(cf);

    if (cf.rank == 0){
        printf("Gamma = %4.2f\n",cf.gamma);
        printf("glbl_ni = %d, dx = %f\n",cf.glbl_ni,cf.dx);
        printf("glbl_nj = %d, dy = %f\n",cf.glbl_nj,cf.dy);
        printf("glbl_nk = %d, dz = %f\n",cf.glbl_nk,cf.dz);
    }

    Kokkos::View<double*[2]> a("testView",5);

    

    //printf("I am %3d of %3d: (%2d,%2d,%2d,%2d,%2d,%2d), (%d,%d,%d), (%2d,%2d,%2d), (%2d,%2d,%2d), (%2d,%2d,%2d,%2d,%2d,%2d)\n"
    //        ,rank,numprocs,left,right,bottom,top,back,front
    //        ,coords[0],coords[1],coords[2]
    //        ,nci,ncj,nck
    //        ,ni,nj,nk
    //        ,starti,endi,startj,endj,startk,endk);

    /* allocate grid coordinate and flow variables */
    double *x = (double*)malloc(cf.ni*cf.nj*cf.nk*sizeof(double));
    double *y = (double*)malloc(cf.ni*cf.nj*cf.nk*sizeof(double));
    double *z = (double*)malloc(cf.ni*cf.nj*cf.nk*sizeof(double));
    double *v = (double*)malloc(cf.nci*cf.ncj*cf.nck*sizeof(double));


    /* calculate rank local grid coordinates */
    for (int k=0; k<cf.nk; ++k){
        for (int j=0; j<cf.nj; ++j){
            for (int i=0; i<cf.ni; ++i){
                idx = (cf.ni*cf.nj)*k+cf.ni*j+i;
                x[idx] = (cf.iStart + i)*cf.dx;
                y[idx] = (cf.jStart + j)*cf.dy;
                z[idx] = (cf.kStart + k)*cf.dz;
            }
        }
    }

    //cf = writeGrid(cf,x,y,z,"grid.cgns");

    /* contrive sample flow variable */
    for (int k=0; k<cf.nck; ++k){
        for (int j=0; j<cf.ncj; ++j){
            for (int i=0; i<cf.nci; ++i){
                idx = (cf.nci*cf.ncj)*k+cf.nci*j+i;
                //v[idx] = (double)((starti+i+1)*(startj+j+1)*(startk+k+1))/(double)(ni*nj*nk);
                v[idx] = cf.rank;
            }
        }
    }
    
    writeSolution(cf,x,y,z,v,1,0.65);
    writeSolution(cf,x,y,z,v,2,0.85);

    MPI_Finalize();
    //Kokkos::finalize();
    return 0;
}
