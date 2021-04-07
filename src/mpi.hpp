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

  FS4D leftSend;
  FS4D leftRecv;
  FS4D rightSend;
  FS4D rightRecv;
  FS4D bottomSend;
  FS4D bottomRecv;
  FS4D topSend;
  FS4D topRecv;
  FS4D backSend;
  FS4D backRecv;
  FS4D frontSend;
  FS4D frontRecv;

  FS4DH leftSend_H;
  FS4DH leftRecv_H;
  FS4DH rightSend_H;
  FS4DH rightRecv_H;
  FS4DH bottomSend_H;
  FS4DH bottomRecv_H;
  FS4DH topSend_H;
  FS4DH topRecv_H;
  FS4DH backSend_H;
  FS4DH backRecv_H;
  FS4DH frontSend_H;
  FS4DH frontRecv_H;
};

void haloExchange(struct inputConfig cf, FS4D &deviceV, class mpiBuffers &m);

#endif
