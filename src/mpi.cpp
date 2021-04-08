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
  //all = Kokkos::View<double****, FS_LAYOUT>("all", cf.ngi, cf.ngj, cf.ngk, cf.nvt);
  //all_H = Kokkos::create_mirror_view(all);

  // MPI Datatypes implementation
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

void haloExchange(struct inputConfig cf, FS4D &deviceV, class mpiBuffers &m) {

  Kokkos::Profiling::pushRegion("mpi::haloExchange");

  typedef Kokkos::MDRangePolicy<Kokkos::Rank<4>> policy_rind;
  policy_rind xPol = policy_rind({0, 0, 0, 0}, {cf.ng, cf.ngj, cf.ngk, cf.nvt});
  policy_rind yPol = policy_rind({0, 0, 0, 0}, {cf.ngi, cf.ng, cf.ngk, cf.nvt});
  policy_rind zPol = policy_rind({0, 0, 0, 0}, {cf.ngi, cf.ngj, cf.ng, cf.nvt});

  int wait_count = 0;
  MPI_Request reqs[12];

  // GPU -> CPU
  //Kokkos::Profiling::pushRegion("mpi::haloExchange::deepcopy1");
  //Kokkos::deep_copy(deviceV, m.all); // Remove for CUDA-Aware
  //Kokkos::Profiling::popRegion(); // mpi::haloExchange::deepcopy1
 
///////////////////////////////////////////////
// Post all halo exchange receives
/////////////////////////////////////////////// 
  Kokkos::Profiling::pushRegion("mpi::haloExchange::postRecv");
  MPI_Irecv(deviceV.data(), 1, m.leftRecvSubArray, cf.xMinus, FIESTA_HALO_TAG, cf.comm, &reqs[wait_count++]);
  MPI_Irecv(deviceV.data(), 1, m.rightRecvSubArray, cf.xPlus, FIESTA_HALO_TAG, cf.comm, &reqs[wait_count++]);
  MPI_Irecv(deviceV.data(), 1, m.bottomRecvSubArray, cf.yMinus, FIESTA_HALO_TAG, cf.comm, &reqs[wait_count++]);
  MPI_Irecv(deviceV.data(), 1, m.topRecvSubArray, cf.yPlus, FIESTA_HALO_TAG, cf.comm, &reqs[wait_count++]);
  
  if (cf.ndim == 3) {
    MPI_Irecv(deviceV.data(), 1, m.backRecvSubArray, cf.zMinus, FIESTA_HALO_TAG, cf.comm, &reqs[wait_count++]);
    MPI_Irecv(deviceV.data(), 1, m.frontRecvSubArray, cf.zPlus, FIESTA_HALO_TAG, cf.comm, &reqs[wait_count++]); 
  }
  Kokkos::Profiling::popRegion(); // mpi::haloExchange::postRecv

///////////////////////////////////////////////
// Post all halo exchange sends
/////////////////////////////////////////////// 
  Kokkos::Profiling::pushRegion("mpi::haloExchange::postSend");
  MPI_Isend(deviceV.data(), 1, m.leftSendSubArray, cf.xMinus, FIESTA_HALO_TAG, cf.comm, &reqs[wait_count++]);
  MPI_Isend(deviceV.data(), 1, m.rightSendSubArray, cf.xPlus, FIESTA_HALO_TAG, cf.comm, &reqs[wait_count++]);
  MPI_Isend(deviceV.data(), 1, m.bottomSendSubArray, cf.yMinus, FIESTA_HALO_TAG, cf.comm, &reqs[wait_count++]);
  MPI_Isend(deviceV.data(), 1, m.topSendSubArray, cf.yPlus, FIESTA_HALO_TAG, cf.comm, &reqs[wait_count++]);
  if (cf.ndim == 3) {
    MPI_Isend(deviceV.data(), 1, m.backSendSubArray, cf.zMinus, FIESTA_HALO_TAG, cf.comm, &reqs[wait_count++]);
    MPI_Isend(deviceV.data(), 1, m.frontSendSubArray, cf.zPlus, FIESTA_HALO_TAG, cf.comm, &reqs[wait_count++]);
  }
  Kokkos::Profiling::popRegion(); // mpi::haloExchange::postSend

  // Wait for the sends and receives to finish
  Kokkos::Profiling::pushRegion("mpi::haloExchange::waitall");
  MPI_Waitall(wait_count, reqs, MPI_STATUSES_IGNORE);
  Kokkos::Profiling::popRegion(); // mpi::haloExchange::waitall

  // CPU -> GPU
  //Kokkos::Profiling::pushRegion("mpi::haloExchange::deepcopy2");
  //Kokkos::deep_copy(m.all, deviceV); // Remove for CUDA-Aware
  //Kokkos::Profiling::popRegion(); // mpi::haloExchange::deepcopy2

  Kokkos::Profiling::popRegion(); // mpi::haloExchange
}
