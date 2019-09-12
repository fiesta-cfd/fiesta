#include "mpi_init.hpp"
#include "lsdebug.hpp"

struct inputConfig mpi_init(struct inputConfig cf){

    int rem;

    int dims[3] = {cf.xProcs,cf.yProcs,cf.zProcs};
    int periods[3] = {0,0,0};
    int coords[3];

    /* Get basic MPI parameters */
    MPI_Comm_size(MPI_COMM_WORLD, &cf.numProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &cf.rank);

    /* Create cartesian topology and get rank dimensions and neighbors */
    MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods, 0, &cf.comm);
    MPI_Comm_rank(cf.comm, &cf.rank);
    MPI_Cart_coords(cf.comm, cf.rank, 3, coords);
    MPI_Cart_shift(cf.comm, 0, 1, &cf.xMinus, &cf.xPlus);
    MPI_Cart_shift(cf.comm, 1, 1, &cf.yMinus, &cf.yPlus);
    MPI_Cart_shift(cf.comm, 2, 1, &cf.zMinus, &cf.zPlus);

    /* Distribute grid cells to mpi ranks including uneven remainders */
    rem = cf.glbl_nci % cf.xProcs;
    cf.nci = floor(cf.glbl_nci/cf.xProcs);
    if (coords[0] < rem){
        cf.nci = cf.nci + 1;
        cf.iStart = cf.nci*coords[0];
    }else{
        cf.iStart = rem*(cf.nci+1) + (coords[0]-rem)*cf.nci;
    }
    cf.iEnd = cf.iStart + cf.nci;
    cf.ni = cf.nci + 1;
    cf.ngi = cf.nci + 2*cf.ng;

    rem = cf.glbl_ncj % cf.yProcs;
    cf.ncj = floor(cf.glbl_ncj/cf.yProcs);
    if (coords[1] < rem){
        cf.ncj = cf.ncj + 1;
        cf.jStart = cf.ncj*coords[1];
    }else{
        cf.jStart = rem*(cf.ncj+1) + (coords[1]-rem)*cf.ncj;
    }
    cf.jEnd = cf.jStart + cf.ncj;
    cf.nj = cf.ncj + 1;
    cf.ngj = cf.ncj + 2*cf.ng;

    rem = cf.glbl_nck % cf.zProcs;
    cf.nck = floor(cf.glbl_nck/cf.zProcs);
    if (coords[2] < rem){
        cf.nck = cf.nck + 1;
        cf.kStart = cf.nck*coords[2];
    }else{
        cf.kStart = rem*(cf.nck+1) + (coords[2]-rem)*cf.nck;
    }
    cf.kEnd = cf.kStart + cf.nck;
    cf.nk = cf.nck + 1;
    cf.ngk = cf.nck + 2*cf.ng;

    return cf;
}

void haloExchange(struct inputConfig cf, Kokkos::View<double****> deviceV){
    typename Kokkos::View<double****>::HostMirror hostV = Kokkos::create_mirror_view(deviceV);
    Kokkos::deep_copy(hostV,deviceV);

    MYDBG0
    double *leftOut   = (double*)malloc(cf.ng*cf.ngj*cf.ngk*cf.nv*sizeof(double));
    double *leftIn    = (double*)malloc(cf.ng*cf.ngj*cf.ngk*cf.nv*sizeof(double));
    double *rightOut  = (double*)malloc(cf.ng*cf.ngj*cf.ngk*cf.nv*sizeof(double));
    double *rightIn   = (double*)malloc(cf.ng*cf.ngj*cf.ngk*cf.nv*sizeof(double));
    
    double *bottomOut = (double*)malloc(cf.ngi*cf.ng*cf.ngk*cf.nv*sizeof(double));
    double *bottomIn  = (double*)malloc(cf.ngi*cf.ng*cf.ngk*cf.nv*sizeof(double));
    double *topOut    = (double*)malloc(cf.ngi*cf.ng*cf.ngk*cf.nv*sizeof(double));
    double *topIn     = (double*)malloc(cf.ngi*cf.ng*cf.ngk*cf.nv*sizeof(double));

    double *backOut   = (double*)malloc(cf.ngi*cf.ngj*cf.ng*cf.nv*sizeof(double));
    double *backIn    = (double*)malloc(cf.ngi*cf.ngj*cf.ng*cf.nv*sizeof(double));
    double *frontOut  = (double*)malloc(cf.ngi*cf.ngj*cf.ng*cf.nv*sizeof(double));
    double *frontIn   = (double*)malloc(cf.ngi*cf.ngj*cf.ng*cf.nv*sizeof(double));
    MYDBG0

    int idx;
    MYDBG0
    for (int v=0; v<cf.nv; ++v){
        for (int k=0; k<cf.ngk; ++k){
            for (int j=0; j<cf.ngj; ++j){
                for (int i=cf.ng; i<cf.ng+cf.ng; ++i){
                    idx = v*cf.ng*cf.ngj*cf.ngk + k*cf.ng*cf.ngj + j*cf.ng + (i-cf.ng);
                    leftOut[idx]  = hostV(i,j,k,v);
                    rightOut[idx] = hostV(i+cf.nci-cf.ng,j,k,v);
                }
            }
        }
    }
    MYDBG0
    for (int v=0; v<cf.nv; ++v){
        for (int k=0; k<cf.ngk; ++k){
            for (int j=cf.ng; j<cf.ng+cf.ng; ++j){
                for (int i=0; i<cf.ngi; ++i){
                    idx = v*cf.ngi*cf.ng*cf.ngk + k*cf.ngi*cf.ng + (j-cf.ng)*cf.ngi + i;
                    bottomOut[idx]  = hostV(i,j,k,v);
                    topOut[idx] = hostV(i+cf.nci-cf.ng,j,k,v);
                }
            }
        }
    }
    MYDBG0
    for (int v=0; v<cf.nv; ++v){
        for (int k=cf.ng; k<cf.ng+cf.ng; ++k){
            for (int j=0; j<cf.ngj; ++j){
                for (int i=0; i<cf.ngi; ++i){
                    idx = v*cf.ngi*cf.ngj*cf.ng + (k-cf.ng)*cf.ngi*cf.ngj + j*cf.ngi + i;
                    backOut[idx]  = hostV(i,j,k,v);
                    frontOut[idx] = hostV(i+cf.nci-cf.ng,j,k,v);
                }
            }
        }
    }
    MYDBG0


}
