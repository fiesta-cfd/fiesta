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

#ifndef MPI_INIT_H
#define MPI_INIT_H

#include "Kokkos_Core.hpp"
#include "kokkosTypes.hpp"
#include "input.hpp"
#include "mpi.h"

void mpi_init(struct inputConfig &cf);

class mpiBuffers {

public:
  mpiBuffers(struct inputConfig cf);
  FS4D all;  
  FS4DH all_H;
  
  MPI_Datatype leftRecvSubArray;
  MPI_Datatype rightRecvSubArray;
  MPI_Datatype bottomRecvSubArray;
  MPI_Datatype topRecvSubArray;
  MPI_Datatype backRecvSubArray;
  MPI_Datatype frontRecvSubArray;

  MPI_Datatype leftSendSubArray;
  MPI_Datatype rightSendSubArray;
  MPI_Datatype bottomSendSubArray;
  MPI_Datatype topSendSubArray;
  MPI_Datatype backSendSubArray;
  MPI_Datatype frontSendSubArray;
  };

void haloExchange(struct inputConfig cf, FS4D &deviceV, class mpiBuffers &m);

#endif
