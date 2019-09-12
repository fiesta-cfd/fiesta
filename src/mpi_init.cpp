#include "mpi_init.hpp"

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
