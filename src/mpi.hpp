#ifndef MPI_INIT_H
#define MPI_INIT_H

#include "input.hpp"
#include <mpi.h>
#include "Kokkos_Core.hpp"

struct inputConfig  mpi_init(struct inputConfig cf);

class mpiBuffers {

public:

    mpiBuffers(struct inputConfig cf);

    Kokkos::View<double****> leftSend;
    Kokkos::View<double****> leftRecv;
    Kokkos::View<double****> rightSend;
    Kokkos::View<double****> rightRecv;
    Kokkos::View<double****> bottomSend;
    Kokkos::View<double****> bottomRecv;
    Kokkos::View<double****> topSend;
    Kokkos::View<double****> topRecv;
    Kokkos::View<double****> backSend;
    Kokkos::View<double****> backRecv;
    Kokkos::View<double****> frontSend;
    Kokkos::View<double****> frontRecv;

    typename Kokkos::View<double****>::HostMirror leftSend_H;
    typename Kokkos::View<double****>::HostMirror leftRecv_H;
    typename Kokkos::View<double****>::HostMirror rightSend_H;
    typename Kokkos::View<double****>::HostMirror rightRecv_H;
    typename Kokkos::View<double****>::HostMirror bottomSend_H;
    typename Kokkos::View<double****>::HostMirror bottomRecv_H;
    typename Kokkos::View<double****>::HostMirror topSend_H;
    typename Kokkos::View<double****>::HostMirror topRecv_H;
    typename Kokkos::View<double****>::HostMirror backSend_H;
    typename Kokkos::View<double****>::HostMirror backRecv_H;
    typename Kokkos::View<double****>::HostMirror frontSend_H;
    typename Kokkos::View<double****>::HostMirror frontRecv_H;
};

void haloExchange(struct inputConfig cf, Kokkos::View<double****> &deviceV, class mpiBuffers &m);

#endif
