/*
  Copyright 2019-2021 The University of New Mexico

  This file is part of FIESTA.
  
  FIESTA is free software: you can redistribute it and/or modify it under the
  terms of the GNU Lesser General Public License as published by the Free
  Software Foundation, either version 3 of the License, or (at your option) any
  later version.
  
  FIESTA is distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
  A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
  details.
  
  You should have received a copy of the GNU Lesser General Public License
  along with FIESTA.  If not, see <https://www.gnu.org/licenses/>.
*/

#include "mpi.hpp"
#include "debug.hpp"
#include <type_traits>
// #include "mpipcl.h"


void mpi_init(struct inputConfig &cf)
{
  int rem;

  int dims[3] = {cf.xProcs, cf.yProcs, cf.zProcs};
  int periods[3] = {cf.xPer, cf.yPer, cf.zPer};
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
  cf.nci = floor(cf.glbl_nci / cf.xProcs);
  if (coords[0] < rem) {
    cf.nci = cf.nci + 1;
    cf.iStart = cf.nci * coords[0];
  } else {
    cf.iStart = rem * (cf.nci + 1) + (coords[0] - rem) * cf.nci;
  }
  cf.iEnd = cf.iStart + cf.nci;
  cf.ni = cf.nci + 1;
  cf.ngi = cf.nci + 2 * cf.ng;

  rem = cf.glbl_ncj % cf.yProcs;
  cf.ncj = floor(cf.glbl_ncj / cf.yProcs);
  if (coords[1] < rem) {
    cf.ncj = cf.ncj + 1;
    cf.jStart = cf.ncj * coords[1];
  } else {
    cf.jStart = rem * (cf.ncj + 1) + (coords[1] - rem) * cf.ncj;
  }
  cf.jEnd = cf.jStart + cf.ncj;
  cf.nj = cf.ncj + 1;
  cf.ngj = cf.ncj + 2 * cf.ng;

  if (cf.ndim == 2) {
    cf.ngk = 1;
    cf.kStart = 0;
    cf.kEnd = 1;
    cf.nck = 1;
    cf.nk = 1;
  } else {
    rem = cf.glbl_nck % cf.zProcs;
    cf.nck = floor(cf.glbl_nck / cf.zProcs);
    if (coords[2] < rem) {
      cf.nck = cf.nck + 1;
      cf.kStart = cf.nck * coords[2];
    } else {
      cf.kStart = rem * (cf.nck + 1) + (coords[2] - rem) * cf.nck;
    }
    cf.kEnd = cf.kStart + cf.nck;
    cf.nk = cf.nck + 1;
    cf.ngk = cf.nck + 2 * cf.ng;
  }
  cf.xmp = coords[0];
  cf.ymp = coords[0];
  cf.zmp = coords[0];
  cf.globalGridDims[0] = cf.glbl_ni;
  cf.globalGridDims[1] = cf.glbl_nj;
  cf.globalGridDims[2] = cf.glbl_nk;
  cf.globalCellDims[0] = cf.glbl_nci;
  cf.globalCellDims[1] = cf.glbl_ncj;
  cf.globalCellDims[2] = cf.glbl_nck;
  cf.localGridDims[0] = cf.ni;
  cf.localGridDims[1] = cf.nj;
  cf.localGridDims[2] = cf.nk;
  cf.localCellDims[0] = cf.nci;
  cf.localCellDims[1] = cf.ncj;
  cf.localCellDims[2] = cf.nck;
  cf.subdomainOffset[0] = cf.iStart;
  cf.subdomainOffset[1] = cf.jStart;
  cf.subdomainOffset[2] = cf.kStart;
}

void mpiHaloExchange::haloExchange() {

  Kokkos::Profiling::pushRegion("mpi::haloExchange");
  
  MPI_Request reqs[12];

  for (int i = 0; i < 12; i++) reqs[i] = MPI_REQUEST_NULL;

  Kokkos::Profiling::pushRegion("mpi::haloExchange::receiveHalo");
  receiveHalo(reqs);
  Kokkos::Profiling::popRegion(); // receiveHalo
  Kokkos::Profiling::pushRegion("mpi::haloExchange::SendHalo");
  sendHalo(reqs + 6);
  Kokkos::Profiling::popRegion(); // sendHalo

  // Wait for the sends and receives to finish
  Kokkos::Profiling::pushRegion("mpi::haloExchange::waitall");
  MPI_Waitall(12, reqs, MPI_STATUSES_IGNORE);
  Kokkos::Profiling::popRegion(); // mpi::haloExchange::waitall

  Kokkos::Profiling::pushRegion("mpi::haloExchange::unpackHalo");
  unpackHalo();
  Kokkos::Profiling::popRegion(); // unpackHalo

  Kokkos::Profiling::popRegion(); // mpi::haloExchange
}

// Create MPI datatypes for the subarray borders that can be used by 
// MPI routines to send/receive in place. Note that many MPIs have 
// really bad performance doing this!
directHaloExchange::directHaloExchange(struct inputConfig &c, FS4D &v)
  : mpiHaloExchange(c, v) {
  int order;
  if (std::is_same<FS_LAYOUT, Kokkos::LayoutLeft>::value) {
    order = MPI_ORDER_FORTRAN;
  } else if (std::is_same<FS_LAYOUT, Kokkos::LayoutRight>::value) {
    order = MPI_ORDER_C;
  } else {
    cerr << "Invalid array order in mpiBuffers.\n";
    exit(-1);
  }

  int bigsizes[4] = {cf.ngi, cf.ngj, cf.ngk, cf.nvt};  

  int xsubsizes[4] = {cf.ng, cf.ngj, cf.ngk, cf.nvt};
  int leftRecvStarts[4] = {0, 0, 0, 0};
  MPI_Type_create_subarray(4, bigsizes, xsubsizes, leftRecvStarts, order, MPI_DOUBLE, &leftRecvSubArray);
  MPI_Type_commit(&leftRecvSubArray);

  int leftSendStarts[4] = {cf.ng, 0, 0, 0};
  MPI_Type_create_subarray(4, bigsizes, xsubsizes, leftSendStarts, order, MPI_DOUBLE, &leftSendSubArray);
  MPI_Type_commit(&leftSendSubArray);

  int rightRecvStarts[4] = {cf.ngi - cf.ng, 0, 0, 0};
  MPI_Type_create_subarray(4, bigsizes, xsubsizes, rightRecvStarts, order, MPI_DOUBLE, &rightRecvSubArray);
  MPI_Type_commit(&rightRecvSubArray);

  int rightSendStarts[4] = {cf.nci, 0, 0, 0};
  MPI_Type_create_subarray(4, bigsizes, xsubsizes, rightSendStarts, order, MPI_DOUBLE, &rightSendSubArray);
  MPI_Type_commit(&rightSendSubArray);


  int ysubsizes[4] = {cf.ngi, cf.ng, cf.ngk, cf.nvt};

  int bottomRecvStarts[4] = {0, 0, 0, 0};
  MPI_Type_create_subarray(4, bigsizes, ysubsizes, bottomRecvStarts, order, MPI_DOUBLE, &bottomRecvSubArray);
  MPI_Type_commit(&bottomRecvSubArray);

  int bottomSendStarts[4] = {0, cf.ng, 0, 0};
  MPI_Type_create_subarray(4, bigsizes, ysubsizes, bottomSendStarts, order, MPI_DOUBLE, &bottomSendSubArray);
  MPI_Type_commit(&bottomSendSubArray);

  int topRecvStarts[4] = {0, cf.ngj - cf.ng, 0, 0};
  MPI_Type_create_subarray(4, bigsizes, ysubsizes, topRecvStarts, order, MPI_DOUBLE, &topRecvSubArray);
  MPI_Type_commit(&topRecvSubArray);

  int topSendStarts[4] = {0, cf.ncj, 0, 0};
  MPI_Type_create_subarray(4, bigsizes, ysubsizes, topSendStarts, order, MPI_DOUBLE, &topSendSubArray);
  MPI_Type_commit(&topSendSubArray);

  if (cf.ndim == 3) {
    int zsubsizes[4] = {cf.ngi, cf.ngj, cf.ng, cf.nvt};
    
    int backRecvStarts[4] = {0, 0, 0, 0};
    MPI_Type_create_subarray(4, bigsizes, zsubsizes, backRecvStarts, order, MPI_DOUBLE, &backRecvSubArray);
    MPI_Type_commit(&backRecvSubArray);

    int backSendStarts[4] = {0, 0, cf.ng, 0};
    MPI_Type_create_subarray(4, bigsizes, zsubsizes, backSendStarts, order, MPI_DOUBLE, &backSendSubArray);
    MPI_Type_commit(&backSendSubArray);

    int frontRecvStarts[4] = {0, 0, cf.ngk - cf.ng, 0};
    MPI_Type_create_subarray(4, bigsizes, zsubsizes, frontRecvStarts, order, MPI_DOUBLE, &frontRecvSubArray);
    MPI_Type_commit(&frontRecvSubArray);

    int frontSendStarts[4] = {0, 0, cf.nck, 0};
    MPI_Type_create_subarray(4, bigsizes, zsubsizes, frontSendStarts, order, MPI_DOUBLE, &frontSendSubArray);
    MPI_Type_commit(&frontSendSubArray);
  }
}

#define FIESTA_HALO_TAG 7223
 
void directHaloExchange::receiveHalo(MPI_Request reqs[6]) {

  int wait_count = 0;
///////////////////////////////////////////////
// Post all halo exchange receives
/////////////////////////////////////////////// 
  MPI_Irecv(deviceV.data(), 1, leftRecvSubArray, cf.xMinus, FIESTA_HALO_TAG, cf.comm, &reqs[wait_count++]);
  MPI_Irecv(deviceV.data(), 1, rightRecvSubArray, cf.xPlus, FIESTA_HALO_TAG, cf.comm, &reqs[wait_count++]);

  MPI_Irecv(deviceV.data(), 1, bottomRecvSubArray, cf.yMinus, FIESTA_HALO_TAG, cf.comm, &reqs[wait_count++]);
  MPI_Irecv(deviceV.data(), 1, topRecvSubArray, cf.yPlus, FIESTA_HALO_TAG, cf.comm, &reqs[wait_count++]);
  
  if (cf.ndim == 3) {
    MPI_Irecv(deviceV.data(), 1, backRecvSubArray, cf.zMinus, FIESTA_HALO_TAG, cf.comm, &reqs[wait_count++]);
    MPI_Irecv(deviceV.data(), 1, frontRecvSubArray, cf.zPlus, FIESTA_HALO_TAG, cf.comm, &reqs[wait_count++]); 
  }
}

void directHaloExchange::sendHalo(MPI_Request reqs[6]) {
///////////////////////////////////////////////
// Post all halo exchange sends
/////////////////////////////////////////////// 
  int wait_count = 0;
  // We fence here to make sure any computation on the host is finished.
  Kokkos::fence();

  // Now that we know the data on the host is sane, we queue up all the sends. MPI *should* be able to 
  // pipeline these, though it doesn't know they're all coming. A neighbor collective might give it
  // more flexibility here.
  MPI_Isend(deviceV.data(), 1, leftSendSubArray, cf.xMinus, FIESTA_HALO_TAG, cf.comm, &reqs[wait_count++]);
  MPI_Isend(deviceV.data(), 1, rightSendSubArray, cf.xPlus, FIESTA_HALO_TAG, cf.comm, &reqs[wait_count++]);
  MPI_Isend(deviceV.data(), 1, bottomSendSubArray, cf.yMinus, FIESTA_HALO_TAG, cf.comm, &reqs[wait_count++]);
  MPI_Isend(deviceV.data(), 1, topSendSubArray, cf.yPlus, FIESTA_HALO_TAG, cf.comm, &reqs[wait_count++]);
  if (cf.ndim == 3) {
    MPI_Isend(deviceV.data(), 1, backSendSubArray, cf.zMinus, FIESTA_HALO_TAG, cf.comm, &reqs[wait_count++]);
    MPI_Isend(deviceV.data(), 1, frontSendSubArray, cf.zPlus, FIESTA_HALO_TAG, cf.comm, &reqs[wait_count++]);
  }
}

packedHaloExchange::packedHaloExchange(struct inputConfig &c, FS4D &v) 
  : mpiHaloExchange(c, v) {

  leftSend   = Kokkos::View<double****,FS_LAYOUT>("leftSend",cf.ng,cf.ngj,cf.ngk,cf.nvt);
  leftRecv   = Kokkos::View<double****,FS_LAYOUT>("leftRecv",cf.ng,cf.ngj,cf.ngk,cf.nvt);
  rightSend  = Kokkos::View<double****,FS_LAYOUT>("rightSend",cf.ng,cf.ngj,cf.ngk,cf.nvt);
  rightRecv  = Kokkos::View<double****,FS_LAYOUT>("rightRecv",cf.ng,cf.ngj,cf.ngk,cf.nvt);
  bottomSend = Kokkos::View<double****,FS_LAYOUT>("bottomSend",cf.ngi,cf.ng,cf.ngk,cf.nvt);
  bottomRecv = Kokkos::View<double****,FS_LAYOUT>("bottomRecv",cf.ngi,cf.ng,cf.ngk,cf.nvt);
  topSend    = Kokkos::View<double****,FS_LAYOUT>("topSend",cf.ngi,cf.ng,cf.ngk,cf.nvt);
  topRecv    = Kokkos::View<double****,FS_LAYOUT>("topRecv",cf.ngi,cf.ng,cf.ngk,cf.nvt);
  backSend   = Kokkos::View<double****,FS_LAYOUT>("backSend",cf.ngi,cf.ngj,cf.ng,cf.nvt);
  backRecv   = Kokkos::View<double****,FS_LAYOUT>("backRecv",cf.ngi,cf.ngj,cf.ng,cf.nvt);
  frontSend  = Kokkos::View<double****,FS_LAYOUT>("frontSend",cf.ngi,cf.ngj,cf.ng,cf.nvt);
  frontRecv  = Kokkos::View<double****,FS_LAYOUT>("frontRecv",cf.ngi,cf.ngj,cf.ng,cf.nvt);
}


// This code requires potentially fences to handle parallelism between MPI and the GPU. MPI
// may implicilty fence if it uses the GPU for the copy to the host, but we have no wauy of
// knowing that, so according to the spec we have to fence pessimisticly here.
void packedHaloExchange::sendHalo(MPI_Request reqs[6]) {
  auto xPol = Kokkos::MDRangePolicy<Kokkos::Rank<4>>({0, 0, 0, 0}, {cf.ng, cf.ngj, cf.ngk, cf.nvt});
  auto yPol = Kokkos::MDRangePolicy<Kokkos::Rank<4>>({0, 0, 0, 0}, {cf.ngi, cf.ng, cf.ngk, cf.nvt});
  auto zPol = Kokkos::MDRangePolicy<Kokkos::Rank<4>>({0, 0, 0, 0}, {cf.ngi, cf.ngj, cf.ng, cf.nvt});

  // XXX Should we do all packs, then fence one, then all sends, or do repeated 
  // pack/fence/send/pack/fence/send? The tradeoff here concurrency between the sends and packs 
  // versus the overheads of the fences, as well as potential memoruy bandwidth limits
  // on the card. For now we do the later, but need to model this to be able to choose
  // appropiately.
  int mng = cf.ng;
  int mnci = cf.nci;
  auto mdeviceV = deviceV;
  Kokkos::parallel_for(
      xPol,
      KOKKOS_CLASS_LAMBDA(const int i, const int j, const int k, const int v)
      {
        leftSend(i, j, k, v) = mdeviceV(mng + i, j, k, v);
        rightSend(i, j, k, v) = mdeviceV(i + mnci, j, k, v);
      });

  Kokkos::fence();
  MPI_Isend(leftSend.data(), cf.ng*cf.ngj*cf.ngk*(cf.nvt), MPI_DOUBLE, cf.xMinus, FIESTA_HALO_TAG, cf.comm, &reqs[0]);
  MPI_Isend(rightSend.data(), cf.ng*cf.ngj*cf.ngk*(cf.nvt), MPI_DOUBLE, cf.xPlus, FIESTA_HALO_TAG, cf.comm, &reqs[1]);

  int mncj = cf.ncj;
  Kokkos::parallel_for( yPol, 
   KOKKOS_CLASS_LAMBDA(const int i, const int j, const int k, const int v) {
      bottomSend(i,j,k,v) = mdeviceV(i,mng+j,k,v);
      topSend(i,j,k,v) = mdeviceV(i,j+mncj,k,v);
    } );

  Kokkos::fence();
  MPI_Isend(bottomSend.data(), cf.ngi*cf.ng*cf.ngk*(cf.nvt), MPI_DOUBLE, cf.yMinus, FIESTA_HALO_TAG, cf.comm, &reqs[2]);
  MPI_Isend(topSend.data(), cf.ngi*cf.ng*cf.ngk*(cf.nvt), MPI_DOUBLE, cf.yPlus, FIESTA_HALO_TAG, cf.comm, &reqs[3]);

  if (cf.ndim == 3){
    int mnck = cf.nck;
    Kokkos::parallel_for( zPol, 
      KOKKOS_CLASS_LAMBDA(const int i, const int j, const int k, const int v) {
        backSend(i,j,k,v) = mdeviceV(i,j,mng+k,v);
        frontSend(i,j,k,v) = mdeviceV(i,j,k+mnck,v);
      } );

    Kokkos::fence();
    MPI_Isend(backSend.data(), cf.ngi*cf.ngj*cf.ng*(cf.nvt), MPI_DOUBLE, cf.zMinus, FIESTA_HALO_TAG, cf.comm, &reqs[4]);
    MPI_Isend(frontSend.data(), cf.ngi*cf.ngj*cf.ng*(cf.nvt), MPI_DOUBLE, cf.zPlus, FIESTA_HALO_TAG, cf.comm, &reqs[5]);
  }
}

void packedHaloExchange::receiveHalo(MPI_Request reqs[6]) {
  int wait_count = 0;

  MPI_Irecv(leftRecv.data(), cf.ng * cf.ngj * cf.ngk * (cf.nvt),
            MPI_DOUBLE, cf.xMinus, FIESTA_HALO_TAG, cf.comm, &reqs[wait_count++]);
  MPI_Irecv(rightRecv.data(), cf.ng * cf.ngj * cf.ngk * (cf.nvt),
            MPI_DOUBLE, cf.xPlus, FIESTA_HALO_TAG, cf.comm, &reqs[wait_count++]);
  MPI_Irecv(bottomRecv.data(), cf.ngi * cf.ng * cf.ngk * (cf.nvt),
            MPI_DOUBLE, cf.yMinus, FIESTA_HALO_TAG, cf.comm, &reqs[wait_count++]);
  MPI_Irecv(topRecv.data(), cf.ngi * cf.ng * cf.ngk * (cf.nvt),
            MPI_DOUBLE, cf.yPlus, FIESTA_HALO_TAG, cf.comm, &reqs[wait_count++]);
  if (cf.ndim == 3) {
    MPI_Irecv(backRecv.data(), cf.ngi * cf.ngj * cf.ng * (cf.nvt),
              MPI_DOUBLE, cf.zMinus, FIESTA_HALO_TAG, cf.comm, &reqs[wait_count++]);
    MPI_Irecv(frontRecv.data(), cf.ngi * cf.ngj * cf.ng * (cf.nvt),
              MPI_DOUBLE, cf.zPlus, FIESTA_HALO_TAG, cf.comm, &reqs[wait_count++]);
  }

}

void packedHaloExchange::unpackHalo()
{
  auto xPol = Kokkos::MDRangePolicy<Kokkos::Rank<4>>({0, 0, 0, 0}, {cf.ng, cf.ngj, cf.ngk, cf.nvt});
  auto yPol = Kokkos::MDRangePolicy<Kokkos::Rank<4>>({0, 0, 0, 0}, {cf.ngi, cf.ng, cf.ngk, cf.nvt});
  auto zPol = Kokkos::MDRangePolicy<Kokkos::Rank<4>>({0, 0, 0, 0}, {cf.ngi, cf.ngj, cf.ng, cf.nvt});
 
  int mng = cf.ng;
  int mngi = cf.ngi;
  auto mdeviceV = deviceV;
  Kokkos::parallel_for( xPol,
    KOKKOS_CLASS_LAMBDA(const int i, const int j, const int k, const int v) {
      mdeviceV(i,j,k,v) = leftRecv(i,j,k,v);
      mdeviceV(mngi-mng+i,j,k,v) = rightRecv(i,j,k,v);
    });

  int mngj = cf.ngj;
  Kokkos::parallel_for( yPol, 
    KOKKOS_CLASS_LAMBDA(const int i, const int j, const int k, const int v) {
      mdeviceV(i,j,k,v) = bottomRecv(i,j,k,v);
      mdeviceV(i,mngj-mng+j,k,v) = topRecv(i,j,k,v);
    });

  if (cf.ndim == 3){
    int mngk = cf.ngk;
    Kokkos::parallel_for( zPol, 
      KOKKOS_CLASS_LAMBDA(const int i, const int j, const int k, const int v) {
        mdeviceV(i,j,k,v) = backRecv(i,j,k,v);
        mdeviceV(i,j,mngk-mng+k,v) = frontRecv(i,j,k,v);
      });
  }
}

void copyHaloExchange::sendHalo(MPI_Request reqs[6]) {
  auto xPol = Kokkos::MDRangePolicy<Kokkos::Rank<4>>({0, 0, 0, 0}, {cf.ng, cf.ngj, cf.ngk, cf.nvt});
  auto yPol = Kokkos::MDRangePolicy<Kokkos::Rank<4>>({0, 0, 0, 0}, {cf.ngi, cf.ng, cf.ngk, cf.nvt});
  auto zPol = Kokkos::MDRangePolicy<Kokkos::Rank<4>>({0, 0, 0, 0}, {cf.ngi, cf.ngj, cf.ng, cf.nvt});

  // x direction pack, copy, and send
  int mng = cf.ng;
  int mnci = cf.nci;
  auto mdeviceV = deviceV;
  Kokkos::parallel_for( xPol, 
    KOKKOS_CLASS_LAMBDA(const int i, const int j, const int k, const int v) { 
      leftSend(i,j,k,v) = mdeviceV(mng+i,j,k,v);
      rightSend(i,j,k,v) = mdeviceV(i+mnci,j,k,v);
    } );

  Kokkos::deep_copy(leftSend_H, leftSend);
  MPI_Isend(leftSend_H.data(), cf.ng*cf.ngj*cf.ngk*(cf.nvt), MPI_DOUBLE, cf.xMinus, FIESTA_HALO_TAG, cf.comm, &reqs[0]);
  Kokkos::deep_copy(rightSend_H, rightSend);
  MPI_Isend(rightSend_H.data(), cf.ng*cf.ngj*cf.ngk*(cf.nvt), MPI_DOUBLE, cf.xPlus, FIESTA_HALO_TAG, cf.comm, &reqs[1]);

  // y direction pack, copy, and send
  int mncj = cf.ncj;
  Kokkos::parallel_for( yPol, 
   KOKKOS_CLASS_LAMBDA(const int i, const int j, const int k, const int v) {
      bottomSend(i,j,k,v) = mdeviceV(i,mng+j,k,v);
      topSend(i,j,k,v) = mdeviceV(i,j+mncj,k,v);
    } );
  Kokkos::deep_copy(bottomSend_H, bottomSend);
  MPI_Isend(bottomSend_H.data(), cf.ngi*cf.ng*cf.ngk*(cf.nvt), MPI_DOUBLE, cf.yMinus, FIESTA_HALO_TAG, cf.comm, &reqs[2]);
  Kokkos::deep_copy(topSend_H, topSend);
  MPI_Isend(topSend_H.data(), cf.ngi*cf.ng*cf.ngk*(cf.nvt), MPI_DOUBLE, cf.yPlus, FIESTA_HALO_TAG, cf.comm, &reqs[3]);

  // z direction pack, copy, and send
  if (cf.ndim == 3){
    int mnck = cf.nck;
    Kokkos::parallel_for( zPol, 
      KOKKOS_CLASS_LAMBDA(const int i, const int j, const int k, const int v) {
        backSend(i,j,k,v) = mdeviceV(i,j,mng+k,v);
        frontSend(i,j,k,v) = mdeviceV(i,j,k+mnck,v);
      } );

    Kokkos::deep_copy(backSend_H, backSend);
    MPI_Isend(backSend_H.data(), cf.ngi*cf.ngj*cf.ng*(cf.nvt), MPI_DOUBLE, cf.zMinus, FIESTA_HALO_TAG, cf.comm, &reqs[4]);
    Kokkos::deep_copy(frontSend_H, frontSend);
    MPI_Isend(frontSend_H.data(), cf.ngi*cf.ngj*cf.ng*(cf.nvt), MPI_DOUBLE, cf.zPlus, FIESTA_HALO_TAG, cf.comm, &reqs[5]);
  }
}

void copyHaloExchange::receiveHalo(MPI_Request reqs[6]) {
  int wait_count = 0;

  MPI_Irecv(leftRecv_H.data(), cf.ng * cf.ngj * cf.ngk * (cf.nvt),
            MPI_DOUBLE, cf.xMinus, FIESTA_HALO_TAG, cf.comm, &reqs[wait_count++]);
  MPI_Irecv(rightRecv_H.data(), cf.ng * cf.ngj * cf.ngk * (cf.nvt),
            MPI_DOUBLE, cf.xPlus, FIESTA_HALO_TAG, cf.comm, &reqs[wait_count++]);
  MPI_Irecv(bottomRecv_H.data(), cf.ngi * cf.ng * cf.ngk * (cf.nvt),
            MPI_DOUBLE, cf.yMinus, FIESTA_HALO_TAG, cf.comm, &reqs[wait_count++]);
  MPI_Irecv(topRecv_H.data(), cf.ngi * cf.ng * cf.ngk * (cf.nvt),
            MPI_DOUBLE, cf.yPlus, FIESTA_HALO_TAG, cf.comm, &reqs[wait_count++]);
  if (cf.ndim == 3) {
    MPI_Irecv(backRecv_H.data(), cf.ngi * cf.ngj * cf.ng * (cf.nvt),
              MPI_DOUBLE, cf.zMinus, FIESTA_HALO_TAG, cf.comm, &reqs[wait_count++]);
    MPI_Irecv(frontRecv_H.data(), cf.ngi * cf.ngj * cf.ng * (cf.nvt),
              MPI_DOUBLE, cf.zPlus, FIESTA_HALO_TAG, cf.comm, &reqs[wait_count++]);
  }
}


void copyHaloExchange::unpackHalo()
{
  // I think these all fence before and after anyway, so no point in
  // trying to overlap them with the unpacks
  Kokkos::deep_copy(leftRecv, leftRecv_H);
  Kokkos::deep_copy(rightRecv, rightRecv_H);
  
  Kokkos::deep_copy(bottomRecv, bottomRecv_H);
  Kokkos::deep_copy(topRecv, topRecv_H);

  if (cf.ndim == 3){
    Kokkos::deep_copy(backRecv, backRecv_H);
    Kokkos::deep_copy(frontRecv, frontRecv_H);
  }
  
  // After everything is in the device-side arrays, just call superclass unpack.
  packedHaloExchange::unpackHalo();
}


copyHaloExchange::copyHaloExchange(struct inputConfig &c, FS4D &v) 
  : packedHaloExchange(c, v) {
  leftSend_H   = Kokkos::create_mirror_view(leftSend);
  leftRecv_H   = Kokkos::create_mirror_view(leftRecv);
  rightSend_H  = Kokkos::create_mirror_view(rightSend);
  rightRecv_H  = Kokkos::create_mirror_view(rightRecv);
  bottomSend_H = Kokkos::create_mirror_view(bottomSend);
  bottomRecv_H = Kokkos::create_mirror_view(bottomRecv);
  topSend_H    = Kokkos::create_mirror_view(topSend);
  topRecv_H    = Kokkos::create_mirror_view(topRecv);
  backSend_H   = Kokkos::create_mirror_view(backSend);
  backRecv_H   = Kokkos::create_mirror_view(backRecv);
  frontSend_H  = Kokkos::create_mirror_view(frontSend);
  frontRecv_H  = Kokkos::create_mirror_view(frontRecv);
}
