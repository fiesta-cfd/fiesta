//#include "input.h"
#include "mpi_init.h"
#include "cgnslib.h"
#include "pcgnslib.h"
#include <mpi.h>

#define MYDBG printf("%s:%d\n",__FILE__,__LINE__);
#define MYDBG0 if(rank==0) printf("%s:%d\n",__FILE__,__LINE__);

int main(){
    MPI_Init(NULL,NULL);

    int idx;

    int index_file, icelldim, iphysdim, index_base, index_zone, index_coord, index_flow, index_field;
    cgsize_t isize[3][3];

    struct inputConfig cf;

    cf = executeConfiguration("input.lua");

    cf = mpi_init(cf);

    if (cf.rank == 0){
        printf("Gamma = %4.2f\n",cf.gamma);
        printf("glbl_ni = %d, dx = %f\n",cf.glbl_ni,cf.dx);
        printf("glbl_nj = %d, dy = %f\n",cf.glbl_nj,cf.dy);
        printf("glbl_nk = %d, dz = %f\n",cf.glbl_nk,cf.dz);
    }

    //printf("I am %3d of %3d: (%2d,%2d,%2d,%2d,%2d,%2d), (%d,%d,%d), (%2d,%2d,%2d), (%2d,%2d,%2d), (%2d,%2d,%2d,%2d,%2d,%2d)\n"
    //        ,rank,numprocs,left,right,bottom,top,back,front
    //        ,coords[0],coords[1],coords[2]
    //        ,nci,ncj,nck
    //        ,ni,nj,nk
    //        ,starti,endi,startj,endj,startk,endk);

    /* allocate grid coordinate and flow variables */
    double *x = malloc(cf.ni*cf.nj*cf.nk*sizeof(double));
    double *y = malloc(cf.ni*cf.nj*cf.nk*sizeof(double));
    double *z = malloc(cf.ni*cf.nj*cf.nk*sizeof(double));
    double *v = malloc(cf.nci*cf.ncj*cf.nck*sizeof(double));


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


    int Cx, Cy, Cz;
    /* calculate rank local start and end indices for CGNS file shapes */
    cgsize_t start[3] = {cf.iStart+1,cf.jStart+1,cf.kStart+1};
    cgsize_t end[3] = {cf.iEnd+1,cf.jEnd+1,cf.kEnd+1};
    cgsize_t endc[3] = {cf.iEnd,cf.jEnd,cf.kEnd};

    /* specify overall file shape for cgns */
    isize[0][0] = cf.glbl_ni;
    isize[0][1] = cf.glbl_nj;
    isize[0][2] = cf.glbl_nk;

    isize[1][0] = isize[0][0] - 1;
    isize[1][1] = isize[0][1] - 1;
    isize[1][2] = isize[0][2] - 1;
    
    isize[2][0] = 0;
    isize[2][1] = 0;
    isize[2][2] = 0;

    cgp_mpi_comm(cf.comm);
    
    /* create cgns file and write grid coordinates in parallel*/
    if (cgp_open("pout.cgns", CG_MODE_WRITE, &index_file) ||
        cg_base_write(index_file,"Base",3,3,&index_base) ||
        cg_zone_write(index_file,index_base,"Zone 1",*isize,CG_Structured,&index_zone))
        cgp_error_exit();

    if (cgp_coord_write(index_file,index_base,index_zone,CG_RealDouble, "CoordinateX", &Cx) ||
        cgp_coord_write(index_file,index_base,index_zone,CG_RealDouble, "CoordinateY", &Cy) ||
        cgp_coord_write(index_file,index_base,index_zone,CG_RealDouble, "CoordinateZ", &Cz))
        cgp_error_exit();

    if (cgp_coord_write_data(index_file,index_base,index_zone,Cx, start, end, x) ||
        cgp_coord_write_data(index_file,index_base,index_zone,Cy, start, end, y) ||
        cgp_coord_write_data(index_file,index_base,index_zone,Cz, start, end, z))
        cgp_error_exit();

    cgp_close(index_file);

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


    int index_sol;

    /* open cgns file and write cell centered flow variable */
    if (cgp_open("pout.cgns", CG_MODE_MODIFY, &index_file))
        cgp_error_exit();

    if (cg_sol_write(index_file,index_base,index_zone,"FlowSolution",CG_CellCenter, &index_sol))
        cgp_error_exit();

    if (cgp_field_write(index_file,index_base,index_zone,index_sol,CG_RealDouble,"Density",&index_flow))
        cgp_error_exit();

    if (cgp_field_write_data(index_file,index_base,index_zone,index_sol,index_flow,start,endc,v))
        cgp_error_exit();

    cgp_close(index_file);

    MPI_Finalize();
    return 0;
}
