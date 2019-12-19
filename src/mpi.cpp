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

    if (cf.ndim == 2){
        cf.ngk = 1;
        cf.kStart = 0;
        cf.kEnd = 1;
        cf.nck = 1;
        cf.nk = 2;
    }else{
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
    }


    return cf;
}

void haloExchange(struct inputConfig cf, Kokkos::View<double**> &deviceV){
    int cv = 0;
    if (cf.ceq == 1){
        cv = 5;
    }

    Kokkos::View<double**> leftSend("leftSend",cf.ng*cf.ngj*cf.ngk,cf.nv+cv);
    Kokkos::View<double**> leftRecv("leftRecv",cf.ng*cf.ngj*cf.ngk,cf.nv+cv);
    Kokkos::View<double**> rightSend("rightSend",cf.ng*cf.ngj*cf.ngk,cf.nv+cv);
    Kokkos::View<double**> rightRecv("rightRecv",cf.ng*cf.ngj*cf.ngk,cf.nv+cv);

    Kokkos::View<double**> bottomSend("bottomSend",cf.ngi*cf.ng*cf.ngk,cf.nv+cv);
    Kokkos::View<double**> bottomRecv("bottomRecv",cf.ngi*cf.ng*cf.ngk,cf.nv+cv);
    Kokkos::View<double**> topSend("topSend",cf.ngi*cf.ng*cf.ngk,cf.nv+cv);
    Kokkos::View<double**> topRecv("topRecv",cf.ngi*cf.ng*cf.ngk,cf.nv+cv);

    typename Kokkos::View<double**>::HostMirror leftSend_H = Kokkos::create_mirror_view(leftSend);
    typename Kokkos::View<double**>::HostMirror leftRecv_H = Kokkos::create_mirror_view(leftRecv);
    typename Kokkos::View<double**>::HostMirror rightSend_H = Kokkos::create_mirror_view(rightSend);
    typename Kokkos::View<double**>::HostMirror rightRecv_H = Kokkos::create_mirror_view(rightRecv);
    typename Kokkos::View<double**>::HostMirror bottomSend_H = Kokkos::create_mirror_view(bottomSend);
    typename Kokkos::View<double**>::HostMirror bottomRecv_H = Kokkos::create_mirror_view(bottomRecv);
    typename Kokkos::View<double**>::HostMirror topSend_H = Kokkos::create_mirror_view(topSend);
    typename Kokkos::View<double**>::HostMirror topRecv_H = Kokkos::create_mirror_view(topRecv);
    //if (cf.ndim == 3){
        Kokkos::View<double**> backSend("backSend",cf.ngi*cf.ngj*cf.ng,cf.nv+cv);
        Kokkos::View<double**> backRecv("backRecv",cf.ngi*cf.ngj*cf.ng,cf.nv+cv);
        Kokkos::View<double**> frontSend("frontSend",cf.ngi*cf.ngj*cf.ng,cf.nv+cv);
        Kokkos::View<double**> frontRecv("frontRecv",cf.ngi*cf.ngj*cf.ng,cf.nv+cv);
        typename Kokkos::View<double**>::HostMirror backSend_H = Kokkos::create_mirror_view(backSend);
        typename Kokkos::View<double**>::HostMirror backRecv_H = Kokkos::create_mirror_view(backRecv);
        typename Kokkos::View<double**>::HostMirror frontSend_H = Kokkos::create_mirror_view(frontSend);
        typename Kokkos::View<double**>::HostMirror frontRecv_H = Kokkos::create_mirror_view(frontRecv);
    //}


    typedef Kokkos::MDRangePolicy<Kokkos::Rank<4>> policy_rind;
    policy_rind xPol = policy_rind({0,0,0,0},{cf.ng ,cf.ngj,cf.ngk,cf.nv+cv});
    policy_rind yPol = policy_rind({0,0,0,0},{cf.ngi,cf.ng ,cf.ngk,cf.nv+cv});
    policy_rind zPol = policy_rind({0,0,0,0},{cf.ngi,cf.ngj,cf.ng ,cf.nv+cv});
    
    Kokkos::parallel_for( xPol, KOKKOS_LAMBDA  (const int i, const int j, const int k, const int v){
        int idxl   = cf.ng + i + j * cf.ngi + k * cf.ngi * cf.ngj;
        int idxr   = i + cf.nci + j * cf.ngi + k * cf.ngi * cf.ngj;
        int idx_bc = i + j * cf.ng + k * cf.ng * cf.ngj;
        leftSend(idx_bc,v) = deviceV(idxl,v);
        rightSend(idx_bc,v) = deviceV(idxr,v);
    });
    Kokkos::parallel_for( yPol, KOKKOS_LAMBDA  (const int i, const int j, const int k, const int v){
        int idxb   = i + (cf.ng + j) * cf.ngi + k * cf.ngi * cf.ngj;
        int idxt   = i + (j + cf.ncj) * cf.ngi + k * cf.ngi * cf.ngj;
        int idx_bc = i + j * cf.ngi + k * cf.ngi * cf.ng;
        bottomSend(idx_bc,v) = deviceV(idxb,v);
        topSend(idx_bc,v) = deviceV(idxt,v);
    });
    if (cf.ndim == 3){
        Kokkos::parallel_for( zPol, KOKKOS_LAMBDA  (const int i, const int j, const int k, const int v){
            int idxb   = i + j * cf.ngi + (cf.ng + k) * cf.ngi * cf.ngj;
            int idxf   = i + j * cf.ngi + (k + cf.nck) * cf.ngi * cf.ngj;
            int idx_bc = i + j * cf.ngi + k * cf.ngi * cf.ngj;
            backSend(idx_bc,v) = deviceV(idxb,v);
            frontSend(idx_bc,v) = deviceV(idxf,v);
        });
    }
    
    Kokkos::deep_copy(leftSend_H,  leftSend  );
    Kokkos::deep_copy(rightSend_H, rightSend );
    Kokkos::deep_copy(bottomSend_H,bottomSend);
    Kokkos::deep_copy(topSend_H,   topSend   );
    if (cf.ndim == 3){
        Kokkos::deep_copy(backSend_H,  backSend  );
        Kokkos::deep_copy(frontSend_H, frontSend );
    }

    int wait_count = 8;
    if (cf.ndim == 3)
        wait_count = 12;
    MPI_Request reqs[wait_count];
    

    //send and recieve left
    MPI_Isend(leftSend_H.data(), cf.ng*cf.ngj*cf.ngk*(cf.nv+cv), MPI_DOUBLE, cf.xMinus,           0, cf.comm, &reqs[0]);
    MPI_Irecv(leftRecv_H.data(), cf.ng*cf.ngj*cf.ngk*(cf.nv+cv), MPI_DOUBLE, cf.xMinus, MPI_ANY_TAG, cf.comm, &reqs[1]);

    //send and recieve right
    MPI_Isend(rightSend_H.data(), cf.ng*cf.ngj*cf.ngk*(cf.nv+cv), MPI_DOUBLE, cf.xPlus,           0, cf.comm, &reqs[2]);
    MPI_Irecv(rightRecv_H.data(), cf.ng*cf.ngj*cf.ngk*(cf.nv+cv), MPI_DOUBLE, cf.xPlus, MPI_ANY_TAG, cf.comm, &reqs[3]);
    
    //send and recieve bottom
    MPI_Isend(bottomSend_H.data(), cf.ngi*cf.ng*cf.ngk*(cf.nv+cv), MPI_DOUBLE, cf.yMinus,           0, cf.comm, &reqs[4]);
    MPI_Irecv(bottomRecv_H.data(), cf.ngi*cf.ng*cf.ngk*(cf.nv+cv), MPI_DOUBLE, cf.yMinus, MPI_ANY_TAG, cf.comm, &reqs[5]);
    
    //send and recieve top
    MPI_Isend(topSend_H.data(), cf.ngi*cf.ng*cf.ngk*(cf.nv+cv), MPI_DOUBLE, cf.yPlus,           0, cf.comm, &reqs[6]);
    MPI_Irecv(topRecv_H.data(), cf.ngi*cf.ng*cf.ngk*(cf.nv+cv), MPI_DOUBLE, cf.yPlus, MPI_ANY_TAG, cf.comm, &reqs[7]);
    
    if (cf.ndim == 3){
        //send and recieve back
        MPI_Isend(backSend_H.data(), cf.ngi*cf.ngj*cf.ng*(cf.nv+cv), MPI_DOUBLE, cf.zMinus,           0, cf.comm, &reqs[8]);
        MPI_Irecv(backRecv_H.data(), cf.ngi*cf.ngj*cf.ng*(cf.nv+cv), MPI_DOUBLE, cf.zMinus, MPI_ANY_TAG, cf.comm, &reqs[9]);
        
        //send and recieve front
        MPI_Isend(frontSend_H.data(), cf.ngi*cf.ngj*cf.ng*(cf.nv+cv), MPI_DOUBLE, cf.zPlus,           0, cf.comm, &reqs[10]);
        MPI_Irecv(frontRecv_H.data(), cf.ngi*cf.ngj*cf.ng*(cf.nv+cv), MPI_DOUBLE, cf.zPlus, MPI_ANY_TAG, cf.comm, &reqs[11]);
    }
    
    MPI_Waitall(wait_count,reqs, MPI_STATUSES_IGNORE);

    Kokkos::deep_copy(leftRecv,  leftRecv_H  );
    Kokkos::deep_copy(rightRecv, rightRecv_H );
    Kokkos::deep_copy(bottomRecv,bottomRecv_H);
    Kokkos::deep_copy(topRecv,   topRecv_H   );
    if (cf.ndim == 3){
        Kokkos::deep_copy(backRecv,  backRecv_H  );
        Kokkos::deep_copy(frontRecv, frontRecv_H );
    }

    Kokkos::parallel_for( xPol, KOKKOS_LAMBDA  (const int i, const int j, const int k, const int v){
        int idxl   = i + j * cf.ngi + k * cf.ngi * cf.ngj;
        int idxr   = cf.ngi - cf.ng + i + j * cf.ngi + k * cf.ngi * cf.ngj;
        int idx_bc = i + j * cf.ng + k * cf.ng * cf.ngj;
        deviceV(idxl,v) = leftRecv(idx_bc,v);
        deviceV(idxr,v) = rightRecv(idx_bc,v);
    });
    Kokkos::parallel_for( yPol, KOKKOS_LAMBDA  (const int i, const int j, const int k, const int v){
        int idxb   = i + j * cf.ngi + k * cf.ngi * cf.ngj;
        int idxt   = i + (cf.ngj - cf.ng + j) * cf.ngi + k * cf.ngi * cf.ngj;
        int idx_bc = i + j * cf.ng + k * cf.ng * cf.ngj;
        deviceV(idxb,v) = bottomRecv(idx_bc,v);
        deviceV(idxt,v) = topRecv(idx_bc,v);
    });
    if (cf.ndim == 3){
        Kokkos::parallel_for( zPol, KOKKOS_LAMBDA  (const int i, const int j, const int k, const int v){
            int idxb   = i + j * cf.ngi + k * cf.ngi * cf.ngj;
            int idxf   = i + j * cf.ngi + (cf.ngk - cf.ng + k) * cf.ngi * cf.ngj;
            int idx_bc = i + j * cf.ng + k * cf.ng * cf.ngj;
            deviceV(idxb,v) = backRecv(idx_bc,v);
            deviceV(idxf,v) = frontRecv(idx_bc,v);
        });
    }
}
