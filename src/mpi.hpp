#ifndef MPI_INIT_H
#define MPI_INIT_H

#include "Kokkos_Core.hpp"
#include "fiesta.hpp"
#include "input.hpp"
#include <mpi.h>

struct inputConfig mpi_init(struct inputConfig cf);

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
