#include "mpi.hpp"
#include "debug.hpp"

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

mpiBuffers::mpiBuffers(struct inputConfig cf){
    int cv = 0;
    if (cf.ceq == 1){
        cv = 5;
    }

    leftSend   = Kokkos::View<double****,FS_LAYOUT>("leftSend",cf.ng,cf.ngj,cf.ngk,cf.nv+cv);
    leftRecv   = Kokkos::View<double****,FS_LAYOUT>("leftRecv",cf.ng,cf.ngj,cf.ngk,cf.nv+cv);
    rightSend  = Kokkos::View<double****,FS_LAYOUT>("rightSend",cf.ng,cf.ngj,cf.ngk,cf.nv+cv);
    rightRecv  = Kokkos::View<double****,FS_LAYOUT>("rightRecv",cf.ng,cf.ngj,cf.ngk,cf.nv+cv);
    bottomSend = Kokkos::View<double****,FS_LAYOUT>("bottomSend",cf.ngi,cf.ng,cf.ngk,cf.nv+cv);
    bottomRecv = Kokkos::View<double****,FS_LAYOUT>("bottomRecv",cf.ngi,cf.ng,cf.ngk,cf.nv+cv);
    topSend    = Kokkos::View<double****,FS_LAYOUT>("topSend",cf.ngi,cf.ng,cf.ngk,cf.nv+cv);
    topRecv    = Kokkos::View<double****,FS_LAYOUT>("topRecv",cf.ngi,cf.ng,cf.ngk,cf.nv+cv);
    backSend   = Kokkos::View<double****,FS_LAYOUT>("backSend",cf.ngi,cf.ngj,cf.ng,cf.nv+cv);
    backRecv   = Kokkos::View<double****,FS_LAYOUT>("backRecv",cf.ngi,cf.ngj,cf.ng,cf.nv+cv);
    frontSend  = Kokkos::View<double****,FS_LAYOUT>("frontSend",cf.ngi,cf.ngj,cf.ng,cf.nv+cv);
    frontRecv  = Kokkos::View<double****,FS_LAYOUT>("frontRecv",cf.ngi,cf.ngj,cf.ng,cf.nv+cv);
    leftSend_H = Kokkos::create_mirror_view(leftSend);
    leftRecv_H = Kokkos::create_mirror_view(leftRecv);
    rightSend_H = Kokkos::create_mirror_view(rightSend);
    rightRecv_H = Kokkos::create_mirror_view(rightRecv);
    bottomSend_H = Kokkos::create_mirror_view(bottomSend);
    bottomRecv_H = Kokkos::create_mirror_view(bottomRecv);
    topSend_H = Kokkos::create_mirror_view(topSend);
    topRecv_H = Kokkos::create_mirror_view(topRecv);
    backSend_H = Kokkos::create_mirror_view(backSend);
    backRecv_H = Kokkos::create_mirror_view(backRecv);
    frontSend_H = Kokkos::create_mirror_view(frontSend);
    frontRecv_H = Kokkos::create_mirror_view(frontRecv);


}

void haloExchange(struct inputConfig cf, FS4D &deviceV, class mpiBuffers &m){
    int cv = 0;
    if (cf.ceq == 1){
        cv = 5;
    }

    typedef Kokkos::MDRangePolicy<Kokkos::Rank<4>> policy_rind;
    policy_rind xPol = policy_rind({0,0,0,0},{cf.ng ,cf.ngj,cf.ngk,cf.nv+cv});
    policy_rind yPol = policy_rind({0,0,0,0},{cf.ngi,cf.ng ,cf.ngk,cf.nv+cv});
    policy_rind zPol = policy_rind({0,0,0,0},{cf.ngi,cf.ngj,cf.ng ,cf.nv+cv});
    
    Kokkos::parallel_for( xPol, KOKKOS_LAMBDA  (const int i, const int j, const int k, const int v){
        m.leftSend(i,j,k,v) = deviceV(cf.ng+i,j,k,v);
        m.rightSend(i,j,k,v) = deviceV(i+cf.nci,j,k,v);
    });
    Kokkos::parallel_for( yPol, KOKKOS_LAMBDA  (const int i, const int j, const int k, const int v){
        m.bottomSend(i,j,k,v) = deviceV(i,cf.ng+j,k,v);
        m.topSend(i,j,k,v) = deviceV(i,j+cf.ncj,k,v);
    });
    if (cf.ndim == 3){
        Kokkos::parallel_for( zPol, KOKKOS_LAMBDA  (const int i, const int j, const int k, const int v){
            m.backSend(i,j,k,v) = deviceV(i,j,cf.ng+k,v);
            m.frontSend(i,j,k,v) = deviceV(i,j,k+cf.nck,v);
        });
    }
    
    Kokkos::deep_copy(m.leftSend_H,  m.leftSend  );
    Kokkos::deep_copy(m.rightSend_H, m.rightSend );
    Kokkos::deep_copy(m.bottomSend_H,m.bottomSend);
    Kokkos::deep_copy(m.topSend_H,   m.topSend   );
    if (cf.ndim == 3){
        Kokkos::deep_copy(m.backSend_H,  m.backSend  );
        Kokkos::deep_copy(m.frontSend_H, m.frontSend );
    }

    int wait_count = 8;
    if (cf.ndim == 3)
        wait_count = 12;
    MPI_Request reqs[wait_count];
    

    //send and recieve left
    MPI_Isend(m.leftSend_H.data(), cf.ng*cf.ngj*cf.ngk*(cf.nv+cv), MPI_DOUBLE, cf.xMinus,           0, cf.comm, &reqs[0]);
    MPI_Irecv(m.leftRecv_H.data(), cf.ng*cf.ngj*cf.ngk*(cf.nv+cv), MPI_DOUBLE, cf.xMinus, MPI_ANY_TAG, cf.comm, &reqs[1]);

    //send and recieve right
    MPI_Isend(m.rightSend_H.data(), cf.ng*cf.ngj*cf.ngk*(cf.nv+cv), MPI_DOUBLE, cf.xPlus,           0, cf.comm, &reqs[2]);
    MPI_Irecv(m.rightRecv_H.data(), cf.ng*cf.ngj*cf.ngk*(cf.nv+cv), MPI_DOUBLE, cf.xPlus, MPI_ANY_TAG, cf.comm, &reqs[3]);
    
    //send and recieve bottom
    MPI_Isend(m.bottomSend_H.data(), cf.ngi*cf.ng*cf.ngk*(cf.nv+cv), MPI_DOUBLE, cf.yMinus,           0, cf.comm, &reqs[4]);
    MPI_Irecv(m.bottomRecv_H.data(), cf.ngi*cf.ng*cf.ngk*(cf.nv+cv), MPI_DOUBLE, cf.yMinus, MPI_ANY_TAG, cf.comm, &reqs[5]);
    
    //send and recieve top
    MPI_Isend(m.topSend_H.data(), cf.ngi*cf.ng*cf.ngk*(cf.nv+cv), MPI_DOUBLE, cf.yPlus,           0, cf.comm, &reqs[6]);
    MPI_Irecv(m.topRecv_H.data(), cf.ngi*cf.ng*cf.ngk*(cf.nv+cv), MPI_DOUBLE, cf.yPlus, MPI_ANY_TAG, cf.comm, &reqs[7]);
    
    if (cf.ndim == 3){
        //send and recieve back
        MPI_Isend(m.backSend_H.data(), cf.ngi*cf.ngj*cf.ng*(cf.nv+cv), MPI_DOUBLE, cf.zMinus,           0, cf.comm, &reqs[8]);
        MPI_Irecv(m.backRecv_H.data(), cf.ngi*cf.ngj*cf.ng*(cf.nv+cv), MPI_DOUBLE, cf.zMinus, MPI_ANY_TAG, cf.comm, &reqs[9]);
        
        //send and recieve front
        MPI_Isend(m.frontSend_H.data(), cf.ngi*cf.ngj*cf.ng*(cf.nv+cv), MPI_DOUBLE, cf.zPlus,           0, cf.comm, &reqs[10]);
        MPI_Irecv(m.frontRecv_H.data(), cf.ngi*cf.ngj*cf.ng*(cf.nv+cv), MPI_DOUBLE, cf.zPlus, MPI_ANY_TAG, cf.comm, &reqs[11]);
    }
    
    MPI_Waitall(wait_count,reqs, MPI_STATUSES_IGNORE);

    Kokkos::deep_copy(m.leftRecv,  m.leftRecv_H  );
    Kokkos::deep_copy(m.rightRecv, m.rightRecv_H );
    Kokkos::deep_copy(m.bottomRecv,m.bottomRecv_H);
    Kokkos::deep_copy(m.topRecv,   m.topRecv_H   );
    if (cf.ndim == 3){
        Kokkos::deep_copy(m.backRecv,  m.backRecv_H  );
        Kokkos::deep_copy(m.frontRecv, m.frontRecv_H );
    }

    Kokkos::parallel_for( xPol, KOKKOS_LAMBDA  (const int i, const int j, const int k, const int v){
        deviceV(i,j,k,v) = m.leftRecv(i,j,k,v);
        deviceV(cf.ngi-cf.ng+i,j,k,v) = m.rightRecv(i,j,k,v);
    });
    Kokkos::parallel_for( yPol, KOKKOS_LAMBDA  (const int i, const int j, const int k, const int v){
        deviceV(i,j,k,v) = m.bottomRecv(i,j,k,v);
        deviceV(i,cf.ngj-cf.ng+j,k,v) = m.topRecv(i,j,k,v);
    });
    if (cf.ndim == 3){
        Kokkos::parallel_for( zPol, KOKKOS_LAMBDA  (const int i, const int j, const int k, const int v){
            deviceV(i,j,k,v) = m.backRecv(i,j,k,v);
            deviceV(i,j,cf.ngk-cf.ng+k,v) = m.frontRecv(i,j,k,v);
        });
    }
}
