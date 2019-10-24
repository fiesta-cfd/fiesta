#include "mpi.hpp"
#include "lsdebug.hpp"

struct inputConfig mpi_init(struct inputConfig cf){

    int rem;

    int dims[3] = {cf.xProcs,cf.yProcs,cf.zProcs};
    //int periods[3] = {0,0,0};
    //printf("Periods = {%d,%d,%d}\n",cf.xPer,cf.yPer,cf.zPer);
    int periods[3] = {cf.xPer,cf.yPer,cf.zPer};
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

void haloExchange(struct inputConfig cf, Kokkos::View<double****> &deviceV){

    Kokkos::View<double****> leftSend("leftSend",cf.ng,cf.ngj,cf.ngk,cf.nv+4);
    Kokkos::View<double****> leftRecv("leftRecv",cf.ng,cf.ngj,cf.ngk,cf.nv+4);
    Kokkos::View<double****> rightSend("rightSend",cf.ng,cf.ngj,cf.ngk,cf.nv+4);
    Kokkos::View<double****> rightRecv("rightRecv",cf.ng,cf.ngj,cf.ngk,cf.nv+4);

    Kokkos::View<double****> bottomSend("bottomSend",cf.ngi,cf.ng,cf.ngk,cf.nv+4);
    Kokkos::View<double****> bottomRecv("bottomRecv",cf.ngi,cf.ng,cf.ngk,cf.nv+4);
    Kokkos::View<double****> topSend("topSend",cf.ngi,cf.ng,cf.ngk,cf.nv+4);
    Kokkos::View<double****> topRecv("topRecv",cf.ngi,cf.ng,cf.ngk,cf.nv+4);

    Kokkos::View<double****> backSend("backSend",cf.ngi,cf.ngj,cf.ng,cf.nv+4);
    Kokkos::View<double****> backRecv("backRecv",cf.ngi,cf.ngj,cf.ng,cf.nv+4);
    Kokkos::View<double****> frontSend("frontSend",cf.ngi,cf.ngj,cf.ng,cf.nv+4);
    Kokkos::View<double****> frontRecv("frontRecv",cf.ngi,cf.ngj,cf.ng,cf.nv+4);

    typename Kokkos::View<double****>::HostMirror leftSend_H = Kokkos::create_mirror_view(leftSend);
    typename Kokkos::View<double****>::HostMirror leftRecv_H = Kokkos::create_mirror_view(leftRecv);
    typename Kokkos::View<double****>::HostMirror rightSend_H = Kokkos::create_mirror_view(rightSend);
    typename Kokkos::View<double****>::HostMirror rightRecv_H = Kokkos::create_mirror_view(rightRecv);
    typename Kokkos::View<double****>::HostMirror bottomSend_H = Kokkos::create_mirror_view(bottomSend);
    typename Kokkos::View<double****>::HostMirror bottomRecv_H = Kokkos::create_mirror_view(bottomRecv);
    typename Kokkos::View<double****>::HostMirror topSend_H = Kokkos::create_mirror_view(topSend);
    typename Kokkos::View<double****>::HostMirror topRecv_H = Kokkos::create_mirror_view(topRecv);
    typename Kokkos::View<double****>::HostMirror backSend_H = Kokkos::create_mirror_view(backSend);
    typename Kokkos::View<double****>::HostMirror backRecv_H = Kokkos::create_mirror_view(backRecv);
    typename Kokkos::View<double****>::HostMirror frontSend_H = Kokkos::create_mirror_view(frontSend);
    typename Kokkos::View<double****>::HostMirror frontRecv_H = Kokkos::create_mirror_view(frontRecv);

    MPI_Request reqs[12];

    typedef Kokkos::MDRangePolicy<Kokkos::Rank<4>> policy_rind;
    policy_rind xPol = policy_rind({0,0,0,0},{cf.ng ,cf.ngj,cf.ngk,cf.nv+4});
    policy_rind yPol = policy_rind({0,0,0,0},{cf.ngi,cf.ng ,cf.ngk,cf.nv+4});
    policy_rind zPol = policy_rind({0,0,0,0},{cf.ngi,cf.ngj,cf.ng ,cf.nv+4});
    
    Kokkos::parallel_for( xPol, KOKKOS_LAMBDA  (const int i, const int j, const int k, const int v){
        leftSend(i,j,k,v) = deviceV(cf.ng+i,j,k,v);
        rightSend(i,j,k,v) = deviceV(i+cf.nci,j,k,v);
    });
    Kokkos::parallel_for( yPol, KOKKOS_LAMBDA  (const int i, const int j, const int k, const int v){
        bottomSend(i,j,k,v) = deviceV(i,cf.ng+j,k,v);
        topSend(i,j,k,v) = deviceV(i,j+cf.ncj,k,v);
    });
    Kokkos::parallel_for( zPol, KOKKOS_LAMBDA  (const int i, const int j, const int k, const int v){
        backSend(i,j,k,v) = deviceV(i,j,cf.ng+k,v);
        frontSend(i,j,k,v) = deviceV(i,j,k+cf.nck,v);
    });
    
    Kokkos::deep_copy(leftSend_H,  leftSend  );
    Kokkos::deep_copy(rightSend_H, rightSend );
    Kokkos::deep_copy(bottomSend_H,bottomSend);
    Kokkos::deep_copy(topSend_H,   topSend   );
    Kokkos::deep_copy(backSend_H,  backSend  );
    Kokkos::deep_copy(frontSend_H, frontSend );

    //send and recieve left
    MPI_Isend(leftSend_H.data(), cf.ng*cf.ngj*cf.ngk*(cf.nv+4), MPI_DOUBLE, cf.xMinus,           0, cf.comm, &reqs[0]);
    MPI_Irecv(leftRecv_H.data(), cf.ng*cf.ngj*cf.ngk*(cf.nv+4), MPI_DOUBLE, cf.xMinus, MPI_ANY_TAG, cf.comm, &reqs[1]);

    //send and recieve right
    MPI_Isend(rightSend_H.data(), cf.ng*cf.ngj*cf.ngk*(cf.nv+4), MPI_DOUBLE, cf.xPlus,           0, cf.comm, &reqs[2]);
    MPI_Irecv(rightRecv_H.data(), cf.ng*cf.ngj*cf.ngk*(cf.nv+4), MPI_DOUBLE, cf.xPlus, MPI_ANY_TAG, cf.comm, &reqs[3]);
    
    //send and recieve bottom
    MPI_Isend(bottomSend_H.data(), cf.ngi*cf.ng*cf.ngk*(cf.nv+4), MPI_DOUBLE, cf.yMinus,           0, cf.comm, &reqs[4]);
    MPI_Irecv(bottomRecv_H.data(), cf.ngi*cf.ng*cf.ngk*(cf.nv+4), MPI_DOUBLE, cf.yMinus, MPI_ANY_TAG, cf.comm, &reqs[5]);
    
    //send and recieve top
    MPI_Isend(topSend_H.data(), cf.ngi*cf.ng*cf.ngk*(cf.nv+4), MPI_DOUBLE, cf.yPlus,           0, cf.comm, &reqs[6]);
    MPI_Irecv(topRecv_H.data(), cf.ngi*cf.ng*cf.ngk*(cf.nv+4), MPI_DOUBLE, cf.yPlus, MPI_ANY_TAG, cf.comm, &reqs[7]);
    
    //send and recieve back
    MPI_Isend(backSend_H.data(), cf.ngi*cf.ngj*cf.ng*(cf.nv+4), MPI_DOUBLE, cf.zMinus,           0, cf.comm, &reqs[8]);
    MPI_Irecv(backRecv_H.data(), cf.ngi*cf.ngj*cf.ng*(cf.nv+4), MPI_DOUBLE, cf.zMinus, MPI_ANY_TAG, cf.comm, &reqs[9]);
    
    //send and recieve front
    MPI_Isend(frontSend_H.data(), cf.ngi*cf.ngj*cf.ng*(cf.nv+4), MPI_DOUBLE, cf.zPlus,           0, cf.comm, &reqs[10]);
    MPI_Irecv(frontRecv_H.data(), cf.ngi*cf.ngj*cf.ng*(cf.nv+4), MPI_DOUBLE, cf.zPlus, MPI_ANY_TAG, cf.comm, &reqs[11]);
    
    MPI_Waitall(12,reqs, MPI_STATUSES_IGNORE);

    Kokkos::deep_copy(leftRecv,  leftRecv_H  );
    Kokkos::deep_copy(rightRecv, rightRecv_H );
    Kokkos::deep_copy(bottomRecv,bottomRecv_H);
    Kokkos::deep_copy(topRecv,   topRecv_H   );
    Kokkos::deep_copy(backRecv,  backRecv_H  );
    Kokkos::deep_copy(frontRecv, frontRecv_H );

    Kokkos::parallel_for( xPol, KOKKOS_LAMBDA  (const int i, const int j, const int k, const int v){
        deviceV(i,j,k,v) = leftRecv(i,j,k,v);
        deviceV(cf.ngi-cf.ng+i,j,k,v) = rightRecv(i,j,k,v);
    });
    Kokkos::parallel_for( yPol, KOKKOS_LAMBDA  (const int i, const int j, const int k, const int v){
        deviceV(i,j,k,v) = bottomRecv(i,j,k,v);
        deviceV(i,cf.ngj-cf.ng+j,k,v) = topRecv(i,j,k,v);
    });
    Kokkos::parallel_for( zPol, KOKKOS_LAMBDA  (const int i, const int j, const int k, const int v){
        deviceV(i,j,k,v) = backRecv(i,j,k,v);
        deviceV(i,j,cf.ngk-cf.ng+k,v) = frontRecv(i,j,k,v);
    });
}
