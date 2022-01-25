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
#include "log2.hpp"
#include "kokkosTypes.hpp"

#define FIESTA_FORWARD_TAG 1
#define FIESTA_BACKWARD_TAG 2
#define FIESTA_HALO_TAG 7723

void mpi_init(struct inputConfig &cf){
  int rem;
  bool chunkable;

  int dims[3] = {cf.xProcs, cf.yProcs, cf.zProcs};
  int periods[3] = {cf.xPer, cf.yPer, cf.zPer};
  int coords[3];

  /* Get basic MPI parameters */
  MPI_Comm_size(MPI_COMM_WORLD, &cf.numProcs);
  MPI_Comm_rank(MPI_COMM_WORLD, &cf.rank);
  Log::debug("MPI_INIT A");

  /* Create cartesian topology and get rank dimensions and neighbors */
  MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods, 0, &cf.comm);
  MPI_Comm_rank(cf.comm, &cf.rank);
  MPI_Cart_coords(cf.comm, cf.rank, 3, coords);
  MPI_Cart_shift(cf.comm, 0, 1, &cf.xMinus, &cf.xPlus);
  MPI_Cart_shift(cf.comm, 1, 1, &cf.yMinus, &cf.yPlus);
  MPI_Cart_shift(cf.comm, 2, 1, &cf.zMinus, &cf.zPlus);
  int idx=0;
  int coordsN[3];
  int rankN;  
  bool validCoord;

  Log::debug("MPI_INIT B");

  for(int ih=-1;ih<2;++ih){
    for(int jh=-1;jh<2;++jh){
      for(int kh=-1;kh<2;++kh){
        if(ih==0 && jh==0 && kh == 0) continue;
        coordsN[0]=coords[0]+ih;
        coordsN[1]=coords[1]+jh;
        coordsN[2]=coords[2]+kh;

        validCoord=true;
        for (int d=0; d<3; ++d){
          //if (coordsN[d] >= 0 && coordsN[d] < dims[d])
          //  validCoord=true;
          if (!periods[d] && (coordsN[d] >= dims[d] || coordsN[d] < 0))
            validCoord=false;
        }

        if(validCoord){
          MPI_Cart_rank(cf.comm,coordsN,&rankN);
        }else
          rankN = MPI_PROC_NULL;

        cf.proc[idx]=rankN;

        //Log::debug("Coords ({},{},{}): {} {} ",ih,jh,kh,idx,cf.proc[idx]);

        idx+=1;
      }
    }
  }

  /* Distribute grid cells to mpi ranks including uneven remainders */
  rem = cf.glbl_nci % cf.xProcs;
  chunkable=false;
  if(rem==0)
    chunkable=true;
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
  if(rem==0 && chunkable)
    chunkable=true;
  else
    chunkable=false;
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
    if(rem==0 && chunkable)
      chunkable=true;
    else
      chunkable=false;
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
  cf.globalGridDims.push_back(cf.glbl_ni);
  cf.globalGridDims.push_back(cf.glbl_nj);
  cf.globalCellDims.push_back(cf.glbl_nci);
  cf.globalCellDims.push_back(cf.glbl_ncj);
  cf.localGridDims.push_back(cf.ni);
  cf.localGridDims.push_back(cf.nj);
  cf.localCellDims.push_back(cf.nci);
  cf.localCellDims.push_back(cf.ncj);
  cf.subdomainOffset.push_back(cf.iStart);
  cf.subdomainOffset.push_back(cf.jStart);
  if(cf.ndim==3){
    cf.globalGridDims.push_back(cf.glbl_nk);
    cf.globalCellDims.push_back(cf.glbl_nck);
    cf.localGridDims.push_back(cf.nk);
    cf.localCellDims.push_back(cf.nck);
    cf.subdomainOffset.push_back(cf.kStart);
  }
  if (cf.chunkable && !chunkable){
    Log::warning("Grid and procs do not support chunking, disabling. (Grid must be evenly divisible by procs to enable chunking.)");
  }
  cf.chunkable=cf.chunkable && chunkable;
}

void mpiHaloExchange::haloExchange(){
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
  MPI_Type_create_subarray(4, bigsizes, xsubsizes, leftRecvStarts, order, MPI_FSCAL, &leftRecvSubArray);
  MPI_Type_commit(&leftRecvSubArray);

  int leftSendStarts[4] = {cf.ng, 0, 0, 0};
  MPI_Type_create_subarray(4, bigsizes, xsubsizes, leftSendStarts, order, MPI_FSCAL, &leftSendSubArray);
  MPI_Type_commit(&leftSendSubArray);

  int rightRecvStarts[4] = {cf.ngi - cf.ng, 0, 0, 0};
  MPI_Type_create_subarray(4, bigsizes, xsubsizes, rightRecvStarts, order, MPI_FSCAL, &rightRecvSubArray);
  MPI_Type_commit(&rightRecvSubArray);

  int rightSendStarts[4] = {cf.nci, 0, 0, 0};
  MPI_Type_create_subarray(4, bigsizes, xsubsizes, rightSendStarts, order, MPI_FSCAL, &rightSendSubArray);
  MPI_Type_commit(&rightSendSubArray);


  int ysubsizes[4] = {cf.ngi, cf.ng, cf.ngk, cf.nvt};

  int bottomRecvStarts[4] = {0, 0, 0, 0};
  MPI_Type_create_subarray(4, bigsizes, ysubsizes, bottomRecvStarts, order, MPI_FSCAL, &bottomRecvSubArray);
  MPI_Type_commit(&bottomRecvSubArray);

  int bottomSendStarts[4] = {0, cf.ng, 0, 0};
  MPI_Type_create_subarray(4, bigsizes, ysubsizes, bottomSendStarts, order, MPI_FSCAL, &bottomSendSubArray);
  MPI_Type_commit(&bottomSendSubArray);

  int topRecvStarts[4] = {0, cf.ngj - cf.ng, 0, 0};
  MPI_Type_create_subarray(4, bigsizes, ysubsizes, topRecvStarts, order, MPI_FSCAL, &topRecvSubArray);
  MPI_Type_commit(&topRecvSubArray);

  int topSendStarts[4] = {0, cf.ncj, 0, 0};
  MPI_Type_create_subarray(4, bigsizes, ysubsizes, topSendStarts, order, MPI_FSCAL, &topSendSubArray);
  MPI_Type_commit(&topSendSubArray);

  if (cf.ndim == 3) {
    int zsubsizes[4] = {cf.ngi, cf.ngj, cf.ng, cf.nvt};
    
    int backRecvStarts[4] = {0, 0, 0, 0};
    MPI_Type_create_subarray(4, bigsizes, zsubsizes, backRecvStarts, order, MPI_FSCAL, &backRecvSubArray);
    MPI_Type_commit(&backRecvSubArray);

    int backSendStarts[4] = {0, 0, cf.ng, 0};
    MPI_Type_create_subarray(4, bigsizes, zsubsizes, backSendStarts, order, MPI_FSCAL, &backSendSubArray);
    MPI_Type_commit(&backSendSubArray);

    int frontRecvStarts[4] = {0, 0, cf.ngk - cf.ng, 0};
    MPI_Type_create_subarray(4, bigsizes, zsubsizes, frontRecvStarts, order, MPI_FSCAL, &frontRecvSubArray);
    MPI_Type_commit(&frontRecvSubArray);

    int frontSendStarts[4] = {0, 0, cf.nck, 0};
    MPI_Type_create_subarray(4, bigsizes, zsubsizes, frontSendStarts, order, MPI_FSCAL, &frontSendSubArray);
    MPI_Type_commit(&frontSendSubArray);
  }
}

 
void directHaloExchange::receiveHalo(MPI_Request reqs[6]) {
///////////////////////////////////////////////
// Post all halo exchange receives
/////////////////////////////////////////////// 
  int wait_count = 0;
  MPI_Irecv(deviceV.data(), 1, leftRecvSubArray, cf.xMinus, FIESTA_FORWARD_TAG, cf.comm, &reqs[wait_count++]);
  MPI_Irecv(deviceV.data(), 1, rightRecvSubArray, cf.xPlus, FIESTA_BACKWARD_TAG, cf.comm, &reqs[wait_count++]);

  MPI_Irecv(deviceV.data(), 1, bottomRecvSubArray, cf.yMinus, FIESTA_FORWARD_TAG , cf.comm, &reqs[wait_count++]);
  MPI_Irecv(deviceV.data(), 1, topRecvSubArray, cf.yPlus, FIESTA_BACKWARD_TAG, cf.comm, &reqs[wait_count++]);
  
  if (cf.ndim == 3) {
    MPI_Irecv(deviceV.data(), 1, backRecvSubArray, cf.zMinus, FIESTA_FORWARD_TAG , cf.comm, &reqs[wait_count++]);
    MPI_Irecv(deviceV.data(), 1, frontRecvSubArray, cf.zPlus, FIESTA_BACKWARD_TAG, cf.comm, &reqs[wait_count++]); 
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
  MPI_Isend(deviceV.data(), 1, leftSendSubArray, cf.xMinus, FIESTA_BACKWARD_TAG, cf.comm, &reqs[wait_count++]);
  MPI_Isend(deviceV.data(), 1, rightSendSubArray, cf.xPlus, FIESTA_FORWARD_TAG, cf.comm, &reqs[wait_count++]);
  MPI_Isend(deviceV.data(), 1, bottomSendSubArray, cf.yMinus, FIESTA_BACKWARD_TAG, cf.comm, &reqs[wait_count++]);
  MPI_Isend(deviceV.data(), 1, topSendSubArray, cf.yPlus, FIESTA_FORWARD_TAG, cf.comm, &reqs[wait_count++]);
  if (cf.ndim == 3) {
    MPI_Isend(deviceV.data(), 1, backSendSubArray, cf.zMinus, FIESTA_BACKWARD_TAG, cf.comm, &reqs[wait_count++]);
    MPI_Isend(deviceV.data(), 1, frontSendSubArray, cf.zPlus, FIESTA_FORWARD_TAG, cf.comm, &reqs[wait_count++]);
  }
}

packedHaloExchange::packedHaloExchange(struct inputConfig &c, FS4D &v) 
  : mpiHaloExchange(c, v) {

  leftSend   = FS4D("leftSend",cf.ng,cf.ngj,cf.ngk,cf.nvt);
  leftRecv   = FS4D("leftRecv",cf.ng,cf.ngj,cf.ngk,cf.nvt);
  rightSend  = FS4D("rightSend",cf.ng,cf.ngj,cf.ngk,cf.nvt);
  rightRecv  = FS4D("rightRecv",cf.ng,cf.ngj,cf.ngk,cf.nvt);
  bottomSend = FS4D("bottomSend",cf.ngi,cf.ng,cf.ngk,cf.nvt);
  bottomRecv = FS4D("bottomRecv",cf.ngi,cf.ng,cf.ngk,cf.nvt);
  topSend    = FS4D("topSend",cf.ngi,cf.ng,cf.ngk,cf.nvt);
  topRecv    = FS4D("topRecv",cf.ngi,cf.ng,cf.ngk,cf.nvt);
  backSend   = FS4D("backSend",cf.ngi,cf.ngj,cf.ng,cf.nvt);
  backRecv   = FS4D("backRecv",cf.ngi,cf.ngj,cf.ng,cf.nvt);
  frontSend  = FS4D("frontSend",cf.ngi,cf.ngj,cf.ng,cf.nvt);
  frontRecv  = FS4D("frontRecv",cf.ngi,cf.ngj,cf.ng,cf.nvt);
}


// This code requires potentially fences to handle parallelism between MPI and the GPU. MPI
// may implicilty fence if it uses the GPU for the copy to the host, but we have no wauy of
// knowing that, so according to the spec we have to fence pessimisticly here.
void packedHaloExchange::sendHalo(MPI_Request reqs[6]){
  // XXX Should we do all packs, then fence one, then all sends, or do repeated 
  // pack/fence/send/pack/fence/send? The tradeoff here concurrency between the sends and packs 
  // versus the overheads of the fences, as well as potential memoruy bandwidth limits
  // on the card. For now we do the later, but need to model this to be able to choose
  // appropiately.
  size_t bufferLength = 0;

  bufferLength = cf.ng*cf.ngj*cf.ngk*cf.nvt;
  pack({-1,0,0},deviceV,leftSend);
  pack({+1,0,0},deviceV,rightSend);
  Kokkos::fence();
  MPI_Isend(leftSend.data(),  bufferLength, MPI_FSCAL, cf.xMinus, FIESTA_BACKWARD_TAG, cf.comm, &reqs[0]);
  MPI_Isend(rightSend.data(), bufferLength, MPI_FSCAL, cf.xPlus,  FIESTA_FORWARD_TAG,  cf.comm, &reqs[1]);

  bufferLength = cf.ngi*cf.ng*cf.ngk*cf.nvt;
  pack({0,-1,0},deviceV,bottomSend);
  pack({0,+1,0},deviceV,topSend);
  Kokkos::fence();
  MPI_Isend(bottomSend.data(), bufferLength, MPI_FSCAL, cf.yMinus, FIESTA_BACKWARD_TAG, cf.comm, &reqs[2]);
  MPI_Isend(topSend.data(),    bufferLength, MPI_FSCAL, cf.yPlus,  FIESTA_FORWARD_TAG,  cf.comm, &reqs[3]);

  if (cf.ndim == 3){
    bufferLength = cf.ngi*cf.ngj*cf.ng*cf.nvt;
    pack({0,0,-1},deviceV,backSend);
    pack({0,0,+1},deviceV,frontSend);
    Kokkos::fence();
    MPI_Isend(backSend.data(),  bufferLength, MPI_FSCAL, cf.zMinus, FIESTA_BACKWARD_TAG, cf.comm, &reqs[4]);
    MPI_Isend(frontSend.data(), bufferLength, MPI_FSCAL, cf.zPlus,  FIESTA_FORWARD_TAG,  cf.comm, &reqs[5]);
  }
}

void packedHaloExchange::receiveHalo(MPI_Request reqs[6]){
  size_t bufferLength = 0;

  bufferLength = cf.ng*cf.ngj*cf.ngk*cf.nvt;
  MPI_Irecv(leftRecv.data(),  bufferLength, MPI_FSCAL, cf.xMinus, FIESTA_FORWARD_TAG, cf.comm, &reqs[0]);
  MPI_Irecv(rightRecv.data(), bufferLength, MPI_FSCAL, cf.xPlus, FIESTA_BACKWARD_TAG, cf.comm, &reqs[1]);

  bufferLength = cf.ngi*cf.ng*cf.ngk*cf.nvt;
  MPI_Irecv(bottomRecv.data(), bufferLength, MPI_FSCAL, cf.yMinus, FIESTA_FORWARD_TAG, cf.comm, &reqs[2]);
  MPI_Irecv(topRecv.data(),    bufferLength, MPI_FSCAL, cf.yPlus, FIESTA_BACKWARD_TAG, cf.comm, &reqs[3]);

  if (cf.ndim == 3) {
    bufferLength = cf.ngi*cf.ngj*cf.ng*cf.nvt;
    MPI_Irecv(backRecv.data(),  bufferLength, MPI_FSCAL, cf.zMinus, FIESTA_FORWARD_TAG, cf.comm, &reqs[4]);
    MPI_Irecv(frontRecv.data(), bufferLength, MPI_FSCAL, cf.zPlus, FIESTA_BACKWARD_TAG, cf.comm, &reqs[5]);
  }
}

void packedHaloExchange::unpackHalo(){
  unpack({-1,0,0},deviceV,leftRecv);
  unpack({+1,0,0},deviceV,rightRecv);

  unpack({0,-1,0},deviceV,bottomRecv);
  unpack({0,+1,0},deviceV,topRecv);

  if (cf.ndim == 3){
    unpack({0,0,-1},deviceV,backRecv);
    unpack({0,0,+1},deviceV,frontRecv);
  }
}

void copyHaloExchange::sendHalo(MPI_Request reqs[6]){
  size_t bufferLength = 0;

  // x direction pack, copy, and send
  bufferLength = cf.ng*cf.ngj*cf.ngk*cf.nvt;
  pack({-1,0,0},deviceV,leftSend);
  pack({+1,0,0},deviceV,rightSend);

  Kokkos::deep_copy(leftSend_H, leftSend);
  MPI_Isend(leftSend_H.data(),  bufferLength, MPI_FSCAL, cf.xMinus, FIESTA_BACKWARD_TAG, cf.comm, &reqs[0]);
  Kokkos::deep_copy(rightSend_H, rightSend);
  MPI_Isend(rightSend_H.data(), bufferLength, MPI_FSCAL, cf.xPlus,  FIESTA_FORWARD_TAG,  cf.comm, &reqs[1]);

  // y direction pack, copy, and send
  bufferLength = cf.ngi*cf.ng*cf.ngk*cf.nvt;
  pack({0,-1,0},deviceV,bottomSend);
  pack({0,+1,0},deviceV,topSend);

  Kokkos::deep_copy(bottomSend_H, bottomSend);
  MPI_Isend(bottomSend_H.data(), bufferLength, MPI_FSCAL, cf.yMinus, FIESTA_BACKWARD_TAG, cf.comm, &reqs[2]);
  Kokkos::deep_copy(topSend_H, topSend);
  MPI_Isend(topSend_H.data(),    bufferLength, MPI_FSCAL, cf.yPlus,  FIESTA_FORWARD_TAG,  cf.comm, &reqs[3]);

  // z direction pack, copy, and send
  if (cf.ndim == 3){
    bufferLength = cf.ngi*cf.ngj*cf.ng*cf.nvt;
    pack({0,0,-1},deviceV,backSend);
    pack({0,0,+1},deviceV,frontSend);

    Kokkos::deep_copy(backSend_H, backSend);
    MPI_Isend(backSend_H.data(),  bufferLength, MPI_FSCAL, cf.zMinus, FIESTA_BACKWARD_TAG, cf.comm, &reqs[4]);
    Kokkos::deep_copy(frontSend_H, frontSend);
    MPI_Isend(frontSend_H.data(), bufferLength, MPI_FSCAL, cf.zPlus,  FIESTA_FORWARD_TAG,  cf.comm, &reqs[5]);
  }
}

void copyHaloExchange::receiveHalo(MPI_Request reqs[6]){
  size_t bufferLength = 0;

  bufferLength = cf.ng*cf.ngj*cf.ngk*cf.nvt;
  MPI_Irecv(leftRecv_H.data(),  bufferLength, MPI_FSCAL, cf.xMinus, FIESTA_FORWARD_TAG, cf.comm, &reqs[0]);
  MPI_Irecv(rightRecv_H.data(), bufferLength, MPI_FSCAL, cf.xPlus, FIESTA_BACKWARD_TAG, cf.comm, &reqs[1]);

  bufferLength = cf.ngi*cf.ng*cf.ngk*cf.nvt;
  MPI_Irecv(bottomRecv_H.data(), bufferLength, MPI_FSCAL, cf.yMinus, FIESTA_FORWARD_TAG, cf.comm, &reqs[2]);
  MPI_Irecv(topRecv_H.data(),    bufferLength, MPI_FSCAL, cf.yPlus, FIESTA_BACKWARD_TAG, cf.comm, &reqs[3]);

  if (cf.ndim == 3) {
    bufferLength = cf.ngi*cf.ngj*cf.ng*cf.nvt;
    MPI_Irecv(backRecv_H.data(),  bufferLength, MPI_FSCAL, cf.zMinus, FIESTA_FORWARD_TAG, cf.comm, &reqs[4]);
    MPI_Irecv(frontRecv_H.data(), bufferLength, MPI_FSCAL, cf.zPlus, FIESTA_BACKWARD_TAG, cf.comm, &reqs[5]);
  }
  //int wait_count = 0;

  //MPI_Irecv(leftRecv_H.data(), cf.ng * cf.ngj * cf.ngk * (cf.nvt),
  //          MPI_FSCAL, cf.xMinus, ANTI_CLOCKWISE, cf.comm, &reqs[wait_count++]);
  //MPI_Irecv(rightRecv_H.data(), cf.ng * cf.ngj * cf.ngk * (cf.nvt),
  //          MPI_FSCAL, cf.xPlus, CLOCKWISE, cf.comm, &reqs[wait_count++]);
  //MPI_Irecv(bottomRecv_H.data(), cf.ngi * cf.ng * cf.ngk * (cf.nvt),
  //          MPI_FSCAL, cf.yMinus, ANTI_CLOCKWISE, cf.comm, &reqs[wait_count++]);
  //MPI_Irecv(topRecv_H.data(), cf.ngi * cf.ng * cf.ngk * (cf.nvt),
  //          MPI_FSCAL, cf.yPlus, CLOCKWISE, cf.comm, &reqs[wait_count++]);
  //if (cf.ndim == 3) {
  //  MPI_Irecv(backRecv_H.data(), cf.ngi * cf.ngj * cf.ng * (cf.nvt),
  //            MPI_FSCAL, cf.zMinus, ANTI_CLOCKWISE, cf.comm, &reqs[wait_count++]);
  //  MPI_Irecv(frontRecv_H.data(), cf.ngi * cf.ngj * cf.ng * (cf.nvt),
  //            MPI_FSCAL, cf.zPlus, CLOCKWISE, cf.comm, &reqs[wait_count++]);
  //}
}

void copyHaloExchange::unpackHalo(){
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

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
orderedHaloExchange::orderedHaloExchange(struct inputConfig &c, FS4D &v) 
  : mpiHaloExchange(c, v) {

  leftSend   = Kokkos::View<FSCAL****,FS_LAYOUT>("leftSend",cf.ng,cf.ngj,cf.ngk,cf.nvt);
  leftRecv   = Kokkos::View<FSCAL****,FS_LAYOUT>("leftRecv",cf.ng,cf.ngj,cf.ngk,cf.nvt);
  rightSend  = Kokkos::View<FSCAL****,FS_LAYOUT>("rightSend",cf.ng,cf.ngj,cf.ngk,cf.nvt);
  rightRecv  = Kokkos::View<FSCAL****,FS_LAYOUT>("rightRecv",cf.ng,cf.ngj,cf.ngk,cf.nvt);
  bottomSend = Kokkos::View<FSCAL****,FS_LAYOUT>("bottomSend",cf.ngi,cf.ng,cf.ngk,cf.nvt);
  bottomRecv = Kokkos::View<FSCAL****,FS_LAYOUT>("bottomRecv",cf.ngi,cf.ng,cf.ngk,cf.nvt);
  topSend    = Kokkos::View<FSCAL****,FS_LAYOUT>("topSend",cf.ngi,cf.ng,cf.ngk,cf.nvt);
  topRecv    = Kokkos::View<FSCAL****,FS_LAYOUT>("topRecv",cf.ngi,cf.ng,cf.ngk,cf.nvt);
  backSend   = Kokkos::View<FSCAL****,FS_LAYOUT>("backSend",cf.ngi,cf.ngj,cf.ng,cf.nvt);
  backRecv   = Kokkos::View<FSCAL****,FS_LAYOUT>("backRecv",cf.ngi,cf.ngj,cf.ng,cf.nvt);
  frontSend  = Kokkos::View<FSCAL****,FS_LAYOUT>("frontSend",cf.ngi,cf.ngj,cf.ng,cf.nvt);
  frontRecv  = Kokkos::View<FSCAL****,FS_LAYOUT>("frontRecv",cf.ngi,cf.ngj,cf.ng,cf.nvt);
}

void orderedHaloExchange::haloExchange(){
  int waitCount = 4;
  int buffSize = 0;
  MPI_Request reqs[waitCount];

  //////// X
  packFace({-1,0,0},deviceV,leftSend);  // pack left face
  packFace({+1,0,0},deviceV,rightSend); // pack right face
  Kokkos::fence();

  waitCount=0;
  buffSize=cf.ng*cf.ngj*cf.ngk*cf.nvt;
  MPI_Irecv(leftRecv.data(),  buffSize, MPI_FSCAL, cf.xMinus, FIESTA_FORWARD_TAG,  cf.comm, &reqs[waitCount++]);
  MPI_Irecv(rightRecv.data(), buffSize, MPI_FSCAL, cf.xPlus,  FIESTA_BACKWARD_TAG, cf.comm, &reqs[waitCount++]);
  MPI_Isend(leftSend.data(),  buffSize, MPI_FSCAL, cf.xMinus, FIESTA_BACKWARD_TAG, cf.comm, &reqs[waitCount++]);
  MPI_Isend(rightSend.data(), buffSize, MPI_FSCAL, cf.xPlus,  FIESTA_FORWARD_TAG,  cf.comm, &reqs[waitCount++]);
  MPI_Waitall(waitCount, reqs, MPI_STATUSES_IGNORE);

  unpackFace({-1,0,0},deviceV,leftRecv);
  unpackFace({+1,0,0},deviceV,rightRecv);
  Kokkos::fence();

  //////// Y
  packFace({0,-1,0},deviceV,bottomSend);
  packFace({0,+1,0},deviceV,topSend);
  Kokkos::fence();

  waitCount=0;
  buffSize=cf.ngi*cf.ng*cf.ngk*cf.nvt;
  MPI_Irecv(bottomRecv.data(), buffSize, MPI_FSCAL, cf.yMinus, FIESTA_FORWARD_TAG,  cf.comm, &reqs[waitCount++]);
  MPI_Irecv(topRecv.data(),    buffSize, MPI_FSCAL, cf.yPlus,  FIESTA_BACKWARD_TAG, cf.comm, &reqs[waitCount++]);
  MPI_Isend(bottomSend.data(), buffSize, MPI_FSCAL, cf.yMinus, FIESTA_BACKWARD_TAG, cf.comm, &reqs[waitCount++]);
  MPI_Isend(topSend.data(),    buffSize, MPI_FSCAL, cf.yPlus,  FIESTA_FORWARD_TAG,  cf.comm, &reqs[waitCount++]);
  MPI_Waitall(waitCount, reqs, MPI_STATUSES_IGNORE);

  unpackFace({0,-1,0},deviceV,bottomRecv);
  unpackFace({0,+1,0},deviceV,topRecv);
  Kokkos::fence();

  //////// Z
  if (cf.ndim == 3){
    packFace({0,0,-1},deviceV,backSend);
    packFace({0,0,+1},deviceV,frontSend);
    Kokkos::fence();

    waitCount=0;
    buffSize=cf.ngi*cf.ngj*cf.ng*cf.nvt;
    MPI_Irecv(backRecv.data(),  buffSize, MPI_FSCAL, cf.zMinus, FIESTA_FORWARD_TAG,  cf.comm, &reqs[waitCount++]);
    MPI_Irecv(frontRecv.data(), buffSize, MPI_FSCAL, cf.zPlus,  FIESTA_BACKWARD_TAG, cf.comm, &reqs[waitCount++]);
    MPI_Isend(backSend.data(),  buffSize, MPI_FSCAL, cf.zMinus, FIESTA_BACKWARD_TAG, cf.comm, &reqs[waitCount++]);
    MPI_Isend(frontSend.data(), buffSize, MPI_FSCAL, cf.zPlus,  FIESTA_FORWARD_TAG,  cf.comm, &reqs[waitCount++]);
    MPI_Waitall(waitCount, reqs, MPI_STATUSES_IGNORE);

    unpackFace({0,0,-1},deviceV,backRecv);
    unpackFace({0,0,+1},deviceV,frontRecv);
    Kokkos::fence();
  }
}

unorderedHaloExchange::unorderedHaloExchange(struct inputConfig &c, FS4D &v) 
  : mpiHaloExchange(c, v) {

  int si,sj,sk;
  int idx=0;

  for(int ih=-1;ih<2;++ih){
    for(int jh=-1;jh<2;++jh){
      for(int kh=-1;kh<2;++kh){
        if(ih==0 && jh==0 && kh == 0) continue;
        si=(ih!=0)*cf.ng+(ih==0)*cf.nci;
        sj=(jh!=0)*cf.ng+(jh==0)*cf.ncj;
        sk=(kh!=0)*cf.ng+(kh==0)*cf.nck;

        sendBuffers[idx] = Kokkos::View<FSCAL****,FS_LAYOUT>("send",si,sj,sk,cf.nvt);
        recvBuffers[idx] = Kokkos::View<FSCAL****,FS_LAYOUT>("recv",si,sj,sk,cf.nvt);
        buffSize[idx] = si*sj*sk*cf.nvt;

        idx += 1;
      }
    }
  }
}

void unorderedHaloExchange::haloExchange(){
  int waitCount = 0;
  int idx = 0;
  int tag = 0;
  MPI_Request reqs[52];
  for (int i = 0; i < 52; i++) reqs[i] = MPI_REQUEST_NULL;

  idx=0;
  for(int ih=-1;ih<2;++ih){
    for(int jh=-1;jh<2;++jh){
      for(int kh=-1;kh<2;++kh){
        if(ih==0 && jh==0 && kh == 0) continue;

        tag=100*(ih+1)+10*(jh+1)+(kh+1);

        //Log::debug("Recieve ({},{},{}): {}, {} {} - {}",ih,jh,kh,idx,waitCount,tag,buffSize[idx]);
        MPI_Irecv(recvBuffers[idx].data(), buffSize[idx], MPI_FSCAL, cf.proc[idx], tag,  cf.comm, &reqs[waitCount]);
        waitCount += 1;
        idx += 1;
      }
    }
  }

  idx=0;
  for(int ih=-1;ih<2;++ih){
    for(int jh=-1;jh<2;++jh){
      for(int kh=-1;kh<2;++kh){
        if(ih==0 && jh==0 && kh == 0) continue;

        tag=100*(-ih+1)+10*(-jh+1)+(-kh+1);
        //tag=100*(ih+1)+10*(jh+1)+(kh+1);

        pack({ih,jh,kh},deviceV,sendBuffers[idx]);
        Kokkos::fence();
        //Log::debug("Send ({},{},{}): {} {} {} - {}",ih,jh,kh,idx,waitCount,tag,buffSize[idx]);
        MPI_Isend(sendBuffers[idx].data(), buffSize[idx], MPI_FSCAL, cf.proc[idx], tag,  cf.comm, &reqs[waitCount]);
        waitCount += 1;
        idx += 1;
      }
    }
  }

  MPI_Waitall(waitCount, reqs, MPI_STATUSES_IGNORE);
  //MPI_Status stat[52];
  //MPI_Waitall(waitCount, reqs, stat);
  //for (int i=0; i<52; ++i)
  //  Log::debugAll("{}",stat[i].MPI_ERROR);

  idx=0;
  for(int ih=-1;ih<2;++ih){
    for(int jh=-1;jh<2;++jh){
      for(int kh=-1;kh<2;++kh){
        if(ih==0 && jh==0 && kh == 0) continue;

        unpack({ih,jh,kh},deviceV,recvBuffers[idx]);
        idx += 1;
      }
    }
  }
  Kokkos::fence();
}

orderedHostHaloExchange::orderedHostHaloExchange(struct inputConfig &c, FS4D &v) 
  : mpiHaloExchange(c, v) {

  leftSend   = Kokkos::View<FSCAL****,FS_LAYOUT>("leftSend",cf.ng,cf.ngj,cf.ngk,cf.nvt);
  leftRecv   = Kokkos::View<FSCAL****,FS_LAYOUT>("leftRecv",cf.ng,cf.ngj,cf.ngk,cf.nvt);
  rightSend  = Kokkos::View<FSCAL****,FS_LAYOUT>("rightSend",cf.ng,cf.ngj,cf.ngk,cf.nvt);
  rightRecv  = Kokkos::View<FSCAL****,FS_LAYOUT>("rightRecv",cf.ng,cf.ngj,cf.ngk,cf.nvt);
  bottomSend = Kokkos::View<FSCAL****,FS_LAYOUT>("bottomSend",cf.ngi,cf.ng,cf.ngk,cf.nvt);
  bottomRecv = Kokkos::View<FSCAL****,FS_LAYOUT>("bottomRecv",cf.ngi,cf.ng,cf.ngk,cf.nvt);
  topSend    = Kokkos::View<FSCAL****,FS_LAYOUT>("topSend",cf.ngi,cf.ng,cf.ngk,cf.nvt);
  topRecv    = Kokkos::View<FSCAL****,FS_LAYOUT>("topRecv",cf.ngi,cf.ng,cf.ngk,cf.nvt);
  backSend   = Kokkos::View<FSCAL****,FS_LAYOUT>("backSend",cf.ngi,cf.ngj,cf.ng,cf.nvt);
  backRecv   = Kokkos::View<FSCAL****,FS_LAYOUT>("backRecv",cf.ngi,cf.ngj,cf.ng,cf.nvt);
  frontSend  = Kokkos::View<FSCAL****,FS_LAYOUT>("frontSend",cf.ngi,cf.ngj,cf.ng,cf.nvt);
  frontRecv  = Kokkos::View<FSCAL****,FS_LAYOUT>("frontRecv",cf.ngi,cf.ngj,cf.ng,cf.nvt);
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

void orderedHostHaloExchange::haloExchange(){
  int waitCount = 4;
  int buffSize = 0;
  MPI_Request reqs[waitCount];

  //////// X
  packFace({-1,0,0},deviceV,leftSend);  // pack left face
  packFace({+1,0,0},deviceV,rightSend); // pack right face
  Kokkos::fence();
  Kokkos::deep_copy(leftSend_H, leftSend);
  Kokkos::deep_copy(rightSend_H, rightSend);

  waitCount=0;
  buffSize=cf.ng*cf.ngj*cf.ngk*cf.nvt;
  MPI_Irecv(leftRecv_H.data(),  buffSize, MPI_FSCAL, cf.xMinus, FIESTA_FORWARD_TAG,  cf.comm, &reqs[waitCount++]);
  MPI_Irecv(rightRecv_H.data(), buffSize, MPI_FSCAL, cf.xPlus,  FIESTA_BACKWARD_TAG, cf.comm, &reqs[waitCount++]);
  MPI_Isend(leftSend_H.data(),  buffSize, MPI_FSCAL, cf.xMinus, FIESTA_BACKWARD_TAG, cf.comm, &reqs[waitCount++]);
  MPI_Isend(rightSend_H.data(), buffSize, MPI_FSCAL, cf.xPlus,  FIESTA_FORWARD_TAG,  cf.comm, &reqs[waitCount++]);
  MPI_Waitall(waitCount, reqs, MPI_STATUSES_IGNORE);

  Kokkos::deep_copy(leftRecv, leftRecv_H);
  Kokkos::deep_copy(rightRecv, rightRecv_H);
  unpackFace({-1,0,0},deviceV,leftRecv);
  unpackFace({+1,0,0},deviceV,rightRecv);
  Kokkos::fence();

  //////// Y
  packFace({0,-1,0},deviceV,bottomSend);
  packFace({0,+1,0},deviceV,topSend);
  Kokkos::fence();
  Kokkos::deep_copy(bottomSend_H, bottomSend);
  Kokkos::deep_copy(topSend_H, topSend);

  waitCount=0;
  buffSize=cf.ngi*cf.ng*cf.ngk*cf.nvt;
  MPI_Irecv(bottomRecv_H.data(), buffSize, MPI_FSCAL, cf.yMinus, FIESTA_FORWARD_TAG,  cf.comm, &reqs[waitCount++]);
  MPI_Irecv(topRecv_H.data(),    buffSize, MPI_FSCAL, cf.yPlus,  FIESTA_BACKWARD_TAG, cf.comm, &reqs[waitCount++]);
  MPI_Isend(bottomSend_H.data(), buffSize, MPI_FSCAL, cf.yMinus, FIESTA_BACKWARD_TAG, cf.comm, &reqs[waitCount++]);
  MPI_Isend(topSend_H.data(),    buffSize, MPI_FSCAL, cf.yPlus,  FIESTA_FORWARD_TAG,  cf.comm, &reqs[waitCount++]);
  MPI_Waitall(waitCount, reqs, MPI_STATUSES_IGNORE);

  Kokkos::deep_copy(bottomRecv, bottomRecv_H);
  Kokkos::deep_copy(topRecv, topRecv_H);
  unpackFace({0,-1,0},deviceV,bottomRecv);
  unpackFace({0,+1,0},deviceV,topRecv);
  Kokkos::fence();

  //////// Z
  if (cf.ndim == 3){
    packFace({0,0,-1},deviceV,backSend);
    packFace({0,0,+1},deviceV,frontSend);
    Kokkos::fence();
    Kokkos::deep_copy(backSend_H, backSend);
    Kokkos::deep_copy(frontSend_H, frontSend);

    waitCount=0;
    buffSize=cf.ngi*cf.ngj*cf.ng*cf.nvt;
    MPI_Irecv(backRecv_H.data(),  buffSize, MPI_FSCAL, cf.zMinus, FIESTA_FORWARD_TAG,  cf.comm, &reqs[waitCount++]);
    MPI_Irecv(frontRecv_H.data(), buffSize, MPI_FSCAL, cf.zPlus,  FIESTA_BACKWARD_TAG, cf.comm, &reqs[waitCount++]);
    MPI_Isend(backSend_H.data(),  buffSize, MPI_FSCAL, cf.zMinus, FIESTA_BACKWARD_TAG, cf.comm, &reqs[waitCount++]);
    MPI_Isend(frontSend_H.data(), buffSize, MPI_FSCAL, cf.zPlus,  FIESTA_FORWARD_TAG,  cf.comm, &reqs[waitCount++]);
    MPI_Waitall(waitCount, reqs, MPI_STATUSES_IGNORE);

    Kokkos::deep_copy(backRecv, backRecv_H);
    Kokkos::deep_copy(frontRecv, frontRecv_H);
    unpackFace({0,0,-1},deviceV,backRecv);
    unpackFace({0,0,+1},deviceV,frontRecv);
    Kokkos::fence();
  }
}

void mpiHaloExchange::pack(std::vector<int> ion, FS4D &var, FS4D &buff){
  int si = (ion[0]!=0)*cf.ng + (ion[0]==0)*cf.nci; // (s)ize of send buffer
  int sj = (ion[1]!=0)*cf.ng + (ion[1]==0)*cf.ncj;
  int sk = (ion[2]!=0)*cf.ng + (ion[2]==0)*cf.nck;

  int bi=(ion[0]<0)*cf.ng + (ion[0]>0)*cf.nci + (ion[0]==0)*cf.ng;  //(b)eginning indexes for send buffer
  int bj=(ion[1]<0)*cf.ng + (ion[1]>0)*cf.ncj + (ion[1]==0)*cf.ng;
  int bk=(ion[2]<0)*cf.ng + (ion[2]>0)*cf.nck + (ion[2]==0)*cf.ng;

  Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<4>>({0,0,0,0}, {si,sj,sk,cf.nvt}),
      KOKKOS_LAMBDA(const int i, const int j, const int k, const int v){
        buff(i, j, k, v) = var(bi+i, bj+j, bk+k, v);
      });
}

void mpiHaloExchange::unpack(std::vector<int> ion, FS4D &var, FS4D &buff){
  int si = (ion[0]!=0)*cf.ng + (ion[0]==0)*cf.nci; // (s)ize of send buffer
  int sj = (ion[1]!=0)*cf.ng + (ion[1]==0)*cf.ncj;
  int sk = (ion[2]!=0)*cf.ng + (ion[2]==0)*cf.nck;

  int bi=(ion[0]>0)*(cf.ngi-cf.ng)+(ion[0]==0)*cf.ng;  //(b)eggining indexes for unpack location
  int bj=(ion[1]>0)*(cf.ngj-cf.ng)+(ion[1]==0)*cf.ng;
  int bk=(ion[2]>0)*(cf.ngk-cf.ng)+(ion[2]==0)*cf.ng;

  Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<4>>({0,0,0,0}, {si,sj,sk,cf.nvt}),
    KOKKOS_LAMBDA(const int i, const int j, const int k, const int v){
      var(bi+i,bj+j,bk+k,v) = buff(i,j,k,v);
    });
}

void mpiHaloExchange::packFace(std::vector<int> ion, FS4D &var, FS4D &buff){
  int si = (ion[0]!=0)*cf.ng + (ion[0]==0)*cf.ngi; // (s)ize of send buffer
  int sj = (ion[1]!=0)*cf.ng + (ion[1]==0)*cf.ngj;
  int sk = (ion[2]!=0)*cf.ng + (ion[2]==0)*cf.ngk;

  int bi=(ion[0]<0)*cf.ng + (ion[0]>0)*cf.nci;  //(b)eginning indexes for send buffer
  int bj=(ion[1]<0)*cf.ng + (ion[1]>0)*cf.ncj;
  int bk=(ion[2]<0)*cf.ng + (ion[2]>0)*cf.nck;

  Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<4>>({0,0,0,0}, {si,sj,sk,cf.nvt}),
      KOKKOS_LAMBDA(const int i, const int j, const int k, const int v){
        buff(i, j, k, v) = var(bi+i, bj+j, bk+k, v);
      });
}

void mpiHaloExchange::unpackFace(std::vector<int> ion, FS4D &var, FS4D &buff){
  int si = (ion[0]!=0)*cf.ng + (ion[0]==0)*cf.ngi; // (s)ize of send buffer
  int sj = (ion[1]!=0)*cf.ng + (ion[1]==0)*cf.ngj;
  int sk = (ion[2]!=0)*cf.ng + (ion[2]==0)*cf.ngk;

  int bi=(ion[0]>0)*(cf.ngi-cf.ng);  //(b)eggining indexes for unpack location
  int bj=(ion[1]>0)*(cf.ngj-cf.ng);
  int bk=(ion[2]>0)*(cf.ngk-cf.ng);

  Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<4>>({0,0,0,0}, {si,sj,sk,cf.nvt}),
    KOKKOS_LAMBDA(const int i, const int j, const int k, const int v){
      var(bi+i,bj+j,bk+k,v) = buff(i,j,k,v);
    });
}

void orderedHaloExchange::sendHalo(MPI_Request reqs[6]) {}
void orderedHaloExchange::receiveHalo(MPI_Request reqs[6]) {}
void orderedHaloExchange::unpackHalo(){}
void unorderedHaloExchange::sendHalo(MPI_Request reqs[6]) {}
void unorderedHaloExchange::receiveHalo(MPI_Request reqs[6]) {}
void unorderedHaloExchange::unpackHalo(){}
void orderedHostHaloExchange::sendHalo(MPI_Request reqs[6]) {}
void orderedHostHaloExchange::receiveHalo(MPI_Request reqs[6]) {}
void orderedHostHaloExchange::unpackHalo(){}
