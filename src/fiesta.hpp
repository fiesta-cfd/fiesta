#ifndef FIESTA_HPP
#define FIESTA_HPP

#include "Kokkos_Core.hpp"
//typedef typename Kokkos::LayoutLeft FS_LAYOUT
//typedef typename Kokkos::View<double******, FS_LAYOUT> FS6D;
//typedef typename Kokkos::View<double*****, FS_LAYOUT> FS5D;
//typedef typename Kokkos::View<double****, FS_LAYOUT> FS4D;
//typedef typename Kokkos::View<double***, FS_LAYOUT> FS3D;
//typedef typename Kokkos::View<double**, FS_LAYOUT> FS2D;
//typedef typename Kokkos::View<double*, FS_LAYOUT> FS1D;

//namespace Fiesta {
    typedef typename Kokkos::View<double******> FS6D;
    typedef typename Kokkos::View<double*****> FS5D;
    typedef typename Kokkos::View<double****> FS4D;
    typedef typename Kokkos::View<double***> FS3D;
    typedef typename Kokkos::View<double**> FS2D;
    typedef typename Kokkos::View<double*> FS1D;

    typedef typename Kokkos::View<double******>::HostMirror FS6DH;
    typedef typename Kokkos::View<double*****>::HostMirror FS5DH;
    typedef typename Kokkos::View<double****>::HostMirror FS4DH;
    typedef typename Kokkos::View<double***>::HostMirror FS3DH;
    typedef typename Kokkos::View<double**>::HostMirror FS2DH;
    typedef typename Kokkos::View<double*>::HostMirror FS1DH;
//}

#endif
