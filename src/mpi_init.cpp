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
    

    int idx;
    MPI_Request reqs[12];

    
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
    
    MPI_Barrier(cf.comm);
    
    //send and recieve left
    MPI_Isend(leftOut, cf.ng*cf.ngj*cf.ngk*cf.nv, MPI_DOUBLE, cf.xMinus,           0, cf.comm, &reqs[0]);
    MPI_Irecv(leftIn,  cf.ng*cf.ngj*cf.ngk*cf.nv, MPI_DOUBLE, cf.xMinus, MPI_ANY_TAG, cf.comm, &reqs[1]);

    //send and recieve right
    MPI_Isend(rightOut, cf.ng*cf.ngj*cf.ngk*cf.nv, MPI_DOUBLE, cf.xPlus,           0, cf.comm, &reqs[2]);
    MPI_Irecv(rightIn,  cf.ng*cf.ngj*cf.ngk*cf.nv, MPI_DOUBLE, cf.xPlus, MPI_ANY_TAG, cf.comm, &reqs[3]);
    
    //send and recieve bottom
    MPI_Isend(bottomOut, cf.ngi*cf.ng*cf.ngk*cf.nv, MPI_DOUBLE, cf.yMinus,           0, cf.comm, &reqs[4]);
    MPI_Irecv(bottomIn,  cf.ngi*cf.ng*cf.ngk*cf.nv, MPI_DOUBLE, cf.yMinus, MPI_ANY_TAG, cf.comm, &reqs[5]);
    
    //send and recieve top
    MPI_Isend(topOut, cf.ngi*cf.ng*cf.ngk*cf.nv, MPI_DOUBLE, cf.yPlus,           0, cf.comm, &reqs[6]);
    MPI_Irecv(topIn,  cf.ngi*cf.ng*cf.ngk*cf.nv, MPI_DOUBLE, cf.yPlus, MPI_ANY_TAG, cf.comm, &reqs[7]);
    
    //send and recieve back
    MPI_Isend(backOut, cf.ngi*cf.ngj*cf.ng*cf.nv, MPI_DOUBLE, cf.zMinus,           0, cf.comm, &reqs[8]);
    MPI_Irecv(backIn,  cf.ngi*cf.ngj*cf.ng*cf.nv, MPI_DOUBLE, cf.zMinus, MPI_ANY_TAG, cf.comm, &reqs[9]);
    
    //send and recieve front
    MPI_Isend(frontOut, cf.ngi*cf.ngj*cf.ng*cf.nv, MPI_DOUBLE, cf.zPlus,           0, cf.comm, &reqs[10]);
    MPI_Irecv(frontIn,  cf.ngi*cf.ngj*cf.ng*cf.nv, MPI_DOUBLE, cf.zPlus, MPI_ANY_TAG, cf.comm, &reqs[11]);
    
    MPI_Waitall(12,reqs, MPI_STATUSES_IGNORE);

    if (cf.xMinus != MPI_PROC_NULL){
        for (int v=0; v<cf.nv; ++v){
            for (int k=0; k<cf.ngk; ++k){
                for (int j=0; j<cf.ngj; ++j){
                    for (int i=cf.ng; i<cf.ng+cf.ng; ++i){
                        idx = v*cf.ng*cf.ngj*cf.ngk + k*cf.ng*cf.ngj + j*cf.ng + (i-cf.ng);
                        hostV(i,j,k,v) = leftIn[idx];
                    }
                }
            }
        }
    }
        
    if (cf.xPlus != MPI_PROC_NULL){
        for (int v=0; v<cf.nv; ++v){
            for (int k=0; k<cf.ngk; ++k){
                for (int j=cf.ng; j<cf.ng+cf.ng; ++j){
                    for (int i=0; i<cf.ngi; ++i){
                        idx = v*cf.ngi*cf.ng*cf.ngk + k*cf.ngi*cf.ng + (j-cf.ng)*cf.ngi + i;
                        hostV(i,j,k,v) = bottomIn[idx];
                    }
                }
            }
        }
    }
    
    if (cf.yMinus != MPI_PROC_NULL){
        for (int v=0; v<cf.nv; ++v){
            for (int k=cf.ng; k<cf.ng+cf.ng; ++k){
                for (int j=0; j<cf.ngj; ++j){
                    for (int i=0; i<cf.ngi; ++i){
                        idx = v*cf.ngi*cf.ngj*cf.ng + (k-cf.ng)*cf.ngi*cf.ngj + j*cf.ngi + i;
                        hostV(i,j,k,v) = backIn[idx];
                    }
                }
            }
        }
    }

    if (cf.yPlus != MPI_PROC_NULL){
        for (int v=0; v<cf.nv; ++v){
            for (int k=0; k<cf.ngk; ++k){
                for (int j=0; j<cf.ngj; ++j){
                    for (int i=cf.ng; i<cf.ng+cf.ng; ++i){
                        idx = v*cf.ng*cf.ngj*cf.ngk + k*cf.ng*cf.ngj + j*cf.ng + (i-cf.ng);
                        hostV(i+cf.nci-cf.ng,j,k,v) = rightIn[idx];
                    }
                }
            }
        }
    }
    
    if (cf.zMinus != MPI_PROC_NULL){
        for (int v=0; v<cf.nv; ++v){
            for (int k=0; k<cf.ngk; ++k){
                for (int j=cf.ng; j<cf.ng+cf.ng; ++j){
                    for (int i=0; i<cf.ngi; ++i){
                        idx = v*cf.ngi*cf.ng*cf.ngk + k*cf.ngi*cf.ng + (j-cf.ng)*cf.ngi + i;
                        hostV(i+cf.nci-cf.ng,j,k,v) = topIn[idx];
                    }
                }
            }
        }
    }
    
    if (cf.zPlus != MPI_PROC_NULL){
        for (int v=0; v<cf.nv; ++v){
            for (int k=cf.ng; k<cf.ng+cf.ng; ++k){
                for (int j=0; j<cf.ngj; ++j){
                    for (int i=0; i<cf.ngi; ++i){
                        idx = v*cf.ngi*cf.ngj*cf.ng + (k-cf.ng)*cf.ngi*cf.ngj + j*cf.ngi + i;
                        hostV(i+cf.nci-cf.ng,j,k,v) = frontIn[idx];
                    }
                }
            }
        }
    }

    Kokkos::deep_copy(deviceV,hostV);

    free(leftOut);
    free(leftIn);
    free(rightOut);
    free(rightIn);
    free(bottomOut);
    free(bottomIn);
    free(topOut);
    free(topIn);
    free(backOut);
    free(backIn);
    free(frontOut);
    free(frontIn);
}
