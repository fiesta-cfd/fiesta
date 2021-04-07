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
#include "mpipcl.h"

void mpi_init(struct inputConfig &cf) {

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

mpiBuffers::mpiBuffers(struct inputConfig cf) {
  allSend = Kokkos::View<double****, FS_LAYOUT>("allSend", cf.ngi, cf.ngj, cf.ngk, cf.nvt);

  leftSend = Kokkos::View<double ****, FS_LAYOUT>("leftSend", cf.ng, cf.ngj,
                                                  cf.ngk, cf.nvt);
  leftRecv = Kokkos::View<double ****, FS_LAYOUT>("leftRecv", cf.ng, cf.ngj,
                                                  cf.ngk, cf.nvt);
  rightSend = Kokkos::View<double ****, FS_LAYOUT>("rightSend", cf.ng, cf.ngj,
                                                   cf.ngk, cf.nvt);
  rightRecv = Kokkos::View<double ****, FS_LAYOUT>("rightRecv", cf.ng, cf.ngj,
                                                   cf.ngk, cf.nvt);
  bottomSend = Kokkos::View<double ****, FS_LAYOUT>("bottomSend", cf.ngi, cf.ng,
                                                    cf.ngk, cf.nvt);
  bottomRecv = Kokkos::View<double ****, FS_LAYOUT>("bottomRecv", cf.ngi, cf.ng,
                                                    cf.ngk, cf.nvt);
  topSend = Kokkos::View<double ****, FS_LAYOUT>("topSend", cf.ngi, cf.ng,
                                                 cf.ngk, cf.nvt);
  topRecv = Kokkos::View<double ****, FS_LAYOUT>("topRecv", cf.ngi, cf.ng,
                                                 cf.ngk, cf.nvt);
  if (cf.ndim == 3){
    backSend = Kokkos::View<double ****, FS_LAYOUT>("backSend", cf.ngi, cf.ngj,
                                                    cf.ng, cf.nvt);
    backRecv = Kokkos::View<double ****, FS_LAYOUT>("backRecv", cf.ngi, cf.ngj,
                                                    cf.ng, cf.nvt);
    frontSend = Kokkos::View<double ****, FS_LAYOUT>("frontSend", cf.ngi, cf.ngj,
                                                     cf.ng, cf.nvt);
    frontRecv = Kokkos::View<double ****, FS_LAYOUT>("frontRecv", cf.ngi, cf.ngj,
                                                     cf.ng, cf.nvt);
  }else{
    backSend = Kokkos::View<double ****, FS_LAYOUT>("backSend",1,1,1,1);
    backRecv = Kokkos::View<double ****, FS_LAYOUT>("backRecv",1,1,1,1);
    frontSend = Kokkos::View<double ****, FS_LAYOUT>("frontSend",1,1,1,1);
    frontRecv = Kokkos::View<double ****, FS_LAYOUT>("frontRecv",1,1,1,1);
  }
  allSend_H = Kokkos::create_mirror_view(allSend);
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

#define FIESTA_HALO_TAG 7223

void haloExchange(struct inputConfig cf, FS4D &deviceV, class mpiBuffers &m) {

  Kokkos::Profiling::pushRegion("mpi::haloExchange");

  typedef Kokkos::MDRangePolicy<Kokkos::Rank<4>> policy_rind;
  policy_rind xPol = policy_rind({0, 0, 0, 0}, {cf.ng, cf.ngj, cf.ngk, cf.nvt});
  policy_rind yPol = policy_rind({0, 0, 0, 0}, {cf.ngi, cf.ng, cf.ngk, cf.nvt});
  policy_rind zPol = policy_rind({0, 0, 0, 0}, {cf.ngi, cf.ngj, cf.ng, cf.nvt});

  int wait_count = 0;
  MPI_Request reqs[12];


  int mngi = cf.ngi;
  int mngj = cf.ngj;
  int mngk = cf.ngk;
  int mnci = cf.nci;
  int mncj = cf.ncj;
  int mnck = cf.nck;
  int mng = cf.ng;

/*
// MPI Datatypes implementation
///////////////////////////////////////////////
// Post all halo exchange receives
/////////////////////////////////////////////// 
  MPI_Datatype leftSubArray;
  MPI_Datatype rightSubArray;
  MPI_Datatype bottomSubArray;
  MPI_Datatype topSubArray;
  MPI_Datatype backSubArray;
  MPI_Datatype frontSubArray;

  int bigsizes[4] = {cf.ngi, cf.ngj, cf.ngk, cf.nvt};  

  int xsubsizes[4] = {cf.ng, cf.ngj, cf.ngk, cf.nvt};
  
  int leftStarts[4] = {cf.ng, 0, 0, 0};
  MPI_Type_create_subarray(4, bigsizes, xsubsizes, leftStarts, MPI_ORDER_C, MPI_DOUBLE, &leftSubArray);
  MPI_Type_commit(&leftSubArray);

  int rightStarts[4] = {cf.nci, 0, 0, 0};
  MPI_Type_create_subarray(4, bigsizes, xsubsizes, rightStarts, MPI_ORDER_C, MPI_DOUBLE, &rightSubArray);
  MPI_Type_commit(&rightSubArray);


  int ysubsizes[4] = {cf.ngi, cf.ng, cf.ngk, cf.nvt};

  int bottomStarts[4] = {0, cf.ng, 0, 0};
  MPI_Type_create_subarray(4, bigsizes, ysubsizes, bottomStarts, MPI_ORDER_C, MPI_DOUBLE, &bottomSubArray);
  MPI_Type_commit(&bottomSubArray);

  int topStarts[4] = {0, cf.ncj, 0, 0};
  MPI_Type_create_subarray(4, bigsizes, ysubsizes, topStarts, MPI_ORDER_C, MPI_DOUBLE, &topSubArray);
  MPI_Type_commit(&topSubArray);

  if (cf.ndim == 3) {
    int zsubsizes[4] = {cf.ngi, cf.ngj, cf.ng, cf.nvt};
    
    int backStarts[4] = {0, 0,cf.ng, 0};
    MPI_Type_create_subarray(4, bigsizes, zsubsizes, backStarts, MPI_ORDER_C, MPI_DOUBLE, &backSubArray);
    MPI_Type_commit(&backSubArray);

    int frontStarts[4] = {0, 0, cf.nck, 0};
    MPI_Type_create_subarray(4, bigsizes, zsubsizes, frontStarts, MPI_ORDER_C, MPI_DOUBLE, &frontSubArray);
    MPI_Type_commit(&frontSubArray);
  }


  Kokkos::Profiling::pushRegion("mpi::haloExchange::postRecv");
  MPI_Irecv(m.leftRecv_H.data(), cf.ng * cf.ngj * cf.ngk * (cf.nvt), 
            MPI_DOUBLE, cf.xMinus, FIESTA_HALO_TAG, cf.comm, &reqs[wait_count++]);
  MPI_Irecv(m.rightRecv_H.data(), cf.ng * cf.ngj * cf.ngk * (cf.nvt),
            MPI_DOUBLE, cf.xPlus, FIESTA_HALO_TAG, cf.comm, &reqs[wait_count++]);
  MPI_Irecv(m.bottomRecv_H.data(), cf.ngi * cf.ng * cf.ngk * (cf.nvt),
            MPI_DOUBLE, cf.yMinus, FIESTA_HALO_TAG, cf.comm, &reqs[wait_count++]);
  MPI_Irecv(m.topRecv_H.data(), cf.ngi * cf.ng * cf.ngk * (cf.nvt), 
            MPI_DOUBLE, cf.yPlus, FIESTA_HALO_TAG, cf.comm, &reqs[wait_count++]);
  if (cf.ndim == 3) {
    MPI_Irecv(m.backRecv_H.data(), cf.ngi * cf.ngj * cf.ng * (cf.nvt),
              MPI_DOUBLE, cf.zMinus, FIESTA_HALO_TAG, cf.comm, &reqs[wait_count++]);
    MPI_Irecv(m.frontRecv_H.data(), cf.ngi * cf.ngj * cf.ng * (cf.nvt),
              MPI_DOUBLE, cf.zPlus, FIESTA_HALO_TAG, cf.comm, &reqs[wait_count++]);
  }
  Kokkos::Profiling::popRegion(); // postRecv


  Kokkos::deep_copy(deviceV, m.allSend); // Remove for CUDA-Aware

///////////////////////////////////////////////
// Post all halo exchange sends
/////////////////////////////////////////////// 
  Kokkos::Profiling::pushRegion("mpi::haloExchange::copy-send");

  MPI_Isend(m.allSend.data(), 1, leftSubArray, cf.xMinus, FIESTA_HALO_TAG, cf.comm, &reqs[wait_count++]);
  MPI_Isend(m.allSend.data(), 1, rightSubArray, cf.xPlus, FIESTA_HALO_TAG, cf.comm, &reqs[wait_count++]);
  MPI_Isend(m.allSend.data(), 1, bottomSubArray, cf.yMinus, FIESTA_HALO_TAG, cf.comm, &reqs[wait_count++]);
  MPI_Isend(m.allSend.data(), 1, topSubArray, cf.yPlus, FIESTA_HALO_TAG, cf.comm, &reqs[wait_count++]);
  if (cf.ndim == 3) {
    MPI_Isend(m.allSend.data(), 1, backSubArray, cf.zMinus, FIESTA_HALO_TAG, cf.comm, &reqs[wait_count++]);
    MPI_Isend(m.allSend.data(), 1, frontSubArray, cf.zPlus, FIESTA_HALO_TAG, cf.comm, &reqs[wait_count++]);
  }

  Kokkos::Profiling::popRegion(); // copy-send

  // Wait for the sends and receives to finish
  Kokkos::Profiling::pushRegion("mpi::haloExchange::waitall");
  MPI_Waitall(wait_count, reqs, MPI_STATUSES_IGNORE);
  Kokkos::Profiling::popRegion(); // waitall
*/



///////////////////////////////////////////////
// Post all halo exchange receives
///////////////////////////////////////////////
  Kokkos::Profiling::pushRegion("mpi::haloExchange::postRecv");
  MPI_Irecv(m.leftRecv_H.data(), cf.ng * cf.ngj * cf.ngk * (cf.nvt), 
            MPI_DOUBLE, cf.xMinus, FIESTA_HALO_TAG, cf.comm, &reqs[wait_count++]);
  MPI_Irecv(m.rightRecv_H.data(), cf.ng * cf.ngj * cf.ngk * (cf.nvt),
            MPI_DOUBLE, cf.xPlus, FIESTA_HALO_TAG, cf.comm, &reqs[wait_count++]);
  MPI_Irecv(m.bottomRecv_H.data(), cf.ngi * cf.ng * cf.ngk * (cf.nvt),
            MPI_DOUBLE, cf.yMinus, FIESTA_HALO_TAG, cf.comm, &reqs[wait_count++]);
  MPI_Irecv(m.topRecv_H.data(), cf.ngi * cf.ng * cf.ngk * (cf.nvt), 
            MPI_DOUBLE, cf.yPlus, FIESTA_HALO_TAG, cf.comm, &reqs[wait_count++]);
  if (cf.ndim == 3) {
    MPI_Irecv(m.backRecv_H.data(), cf.ngi * cf.ngj * cf.ng * (cf.nvt),
              MPI_DOUBLE, cf.zMinus, FIESTA_HALO_TAG, cf.comm, &reqs[wait_count++]);
    MPI_Irecv(m.frontRecv_H.data(), cf.ngi * cf.ngj * cf.ng * (cf.nvt),
              MPI_DOUBLE, cf.zPlus, FIESTA_HALO_TAG, cf.comm, &reqs[wait_count++]);
  }
  Kokkos::Profiling::popRegion(); // postRecv

///////////////////////////////////////////////
// Pack and send in x direction
///////////////////////////////////////////////
  Kokkos::Profiling::pushRegion("mpi::haloExchange::copy-send");

  Kokkos::Profiling::pushRegion("mpi::haloExchange::copy-send::x-direction");
  Kokkos::Profiling::pushRegion("mpi::haloExchange::copy-send::x-direction::copyout");
  Kokkos::parallel_for(xPol, KOKKOS_LAMBDA(const int i, const int j, const int k, const int v) {
        m.leftSend(i, j, k, v) = deviceV(mng + i, j, k, v);
        m.rightSend(i, j, k, v) = deviceV(i + mnci, j, k, v);
      });
  Kokkos::deep_copy(m.leftSend_H, m.leftSend);
  Kokkos::deep_copy(m.rightSend_H, m.rightSend);
  Kokkos::Profiling::popRegion(); // x-direction-copyout

  Kokkos::Profiling::pushRegion("mpi::haloExchange::copy-send::x-direction::send");
  MPI_Isend(m.leftSend_H.data(), cf.ng * cf.ngj * cf.ngk * (cf.nvt), MPI_DOUBLE,
            cf.xMinus, FIESTA_HALO_TAG, cf.comm, &reqs[wait_count++]);
  MPI_Isend(m.rightSend_H.data(), cf.ng * cf.ngj * cf.ngk * (cf.nvt),
            MPI_DOUBLE, cf.xPlus, FIESTA_HALO_TAG, cf.comm, &reqs[wait_count++]);
  Kokkos::Profiling::popRegion(); // x-direction-postSend

  Kokkos::Profiling::popRegion(); // copy-send-x-direction

///////////////////////////////////////////////
// Pack and send in y direction
///////////////////////////////////////////////
  Kokkos::Profiling::pushRegion("mpi::haloExchange::copy-send::y-direction");

  Kokkos::Profiling::pushRegion("mpi::haloExchange::copy-send::y-direction::copyout");
  Kokkos::parallel_for(yPol, KOKKOS_LAMBDA(const int i, const int j, const int k, const int v) {
        m.bottomSend(i, j, k, v) = deviceV(i, mng + j, k, v);
        m.topSend(i, j, k, v) = deviceV(i, j + mncj, k, v);
      });
  Kokkos::deep_copy(m.bottomSend_H, m.bottomSend);
  Kokkos::deep_copy(m.topSend_H, m.topSend);
  Kokkos::Profiling::popRegion(); // y-direction-copyout

  Kokkos::Profiling::pushRegion("mpi::haloExchange::copy-send::y-direction::send");
  MPI_Isend(m.bottomSend_H.data(), cf.ngi * cf.ng * cf.ngk * (cf.nvt),
            MPI_DOUBLE, cf.yMinus, FIESTA_HALO_TAG, cf.comm, &reqs[wait_count++]);
  MPI_Isend(m.topSend_H.data(), cf.ngi * cf.ng * cf.ngk * (cf.nvt), MPI_DOUBLE,
            cf.yPlus, FIESTA_HALO_TAG, cf.comm, &reqs[wait_count++]);
  Kokkos::Profiling::popRegion(); // y-direction-postSend

  Kokkos::Profiling::popRegion(); // copy-send-y-direction

///////////////////////////////////////////////
// Pack and send exchange in z direction
///////////////////////////////////////////////
  if (cf.ndim == 3) {
    Kokkos::Profiling::pushRegion("mpi::haloExchange::copy-send::z-direction");

    Kokkos::Profiling::pushRegion("mpi::haloExchange::copy-send::z-direction::copyout");
    Kokkos::parallel_for( zPol,
        KOKKOS_LAMBDA(const int i, const int j, const int k, const int v) {
          m.backSend(i, j, k, v) = deviceV(i, j, mng + k, v);
          m.frontSend(i, j, k, v) = deviceV(i, j, k + mnck, v);
        });
    Kokkos::deep_copy(m.backSend_H, m.backSend);
    Kokkos::deep_copy(m.frontSend_H, m.frontSend);
    Kokkos::Profiling::popRegion(); // z-direction-copyout

    Kokkos::Profiling::pushRegion("mpi::haloExchange::copy-send::z-direction::send");
    MPI_Isend(m.backSend_H.data(), cf.ngi * cf.ngj * cf.ng * (cf.nvt),
              MPI_DOUBLE, cf.zMinus, FIESTA_HALO_TAG, cf.comm, &reqs[wait_count++]);
    MPI_Isend(m.frontSend_H.data(), cf.ngi * cf.ngj * cf.ng * (cf.nvt),
              MPI_DOUBLE, cf.zPlus, FIESTA_HALO_TAG, cf.comm, &reqs[wait_count++]);
    Kokkos::Profiling::popRegion(); // z-direction-postSend

    Kokkos::Profiling::popRegion(); // copy-send-z-direction
  }
  Kokkos::Profiling::popRegion(); // copy-send

  // Wait for the sends and receives to finish
  Kokkos::Profiling::pushRegion("mpi::haloExchange::waitall");
  MPI_Waitall(wait_count, reqs, MPI_STATUSES_IGNORE);
  Kokkos::Profiling::popRegion(); // waitall

  // Now that the sends and receives are done, copy data back to where it belongs
  // In theory we could wait for individual receives to finish with waitany or waitsome
  // and do this as they complete, but that's much more complicated. In the long term,
  // switching to CUDA-aware MPI and letting it do that is best way to achieve that.
  Kokkos::Profiling::pushRegion("mpi::haloExchange::copyin");

  Kokkos::Profiling::pushRegion("mpi::haloExchange::copyin::x-direction");
  Kokkos::deep_copy(m.leftRecv, m.leftRecv_H);
  Kokkos::deep_copy(m.rightRecv, m.rightRecv_H);
  Kokkos::parallel_for( xPol, KOKKOS_LAMBDA(const int i, const int j, const int k, const int v) {
        deviceV(i, j, k, v) = m.leftRecv(i, j, k, v);
        deviceV(mngi - mng + i, j, k, v) = m.rightRecv(i, j, k, v);
      });
  Kokkos::Profiling::popRegion(); // x-direction-copyin

  Kokkos::Profiling::pushRegion("mpi::haloExchange::copyin::y-direction");
  Kokkos::deep_copy(m.bottomRecv, m.bottomRecv_H);
  Kokkos::deep_copy(m.topRecv, m.topRecv_H);
  Kokkos::parallel_for( yPol, KOKKOS_LAMBDA(const int i, const int j, const int k, const int v) {
        deviceV(i, j, k, v) = m.bottomRecv(i, j, k, v);
        deviceV(i, mngj - mng + j, k, v) = m.topRecv(i, j, k, v);
      });
  Kokkos::Profiling::popRegion(); // y-direction-copyin

  if (cf.ndim == 3) {
    Kokkos::Profiling::pushRegion("mpi::haloExchange::copyin::z-direction");
    Kokkos::deep_copy(m.backRecv, m.backRecv_H);
    Kokkos::deep_copy(m.frontRecv, m.frontRecv_H);

    Kokkos::parallel_for( zPol,
        KOKKOS_LAMBDA(const int i, const int j, const int k, const int v) {
          deviceV(i, j, k, v) = m.backRecv(i, j, k, v);
          deviceV(i, j, mngk - mng + k, v) = m.frontRecv(i, j, k, v);
        });
    Kokkos::Profiling::popRegion(); // z-direction-copyin
  }
  Kokkos::Profiling::popRegion(); // mpi::haloExchange::copyin

  Kokkos::Profiling::popRegion(); // mpi::haloExchange
}
