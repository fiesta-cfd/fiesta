#ifndef FIESTA_HPP
#define FIESTA_HPP

#include "Kokkos_Core.hpp"

#define FS_LAYOUT  Kokkos::LayoutLeft

typedef typename Kokkos::View<double******, FS_LAYOUT> FS6D;
typedef typename Kokkos::View<double*****,  FS_LAYOUT> FS5D;
typedef typename Kokkos::View<double****,   FS_LAYOUT> FS4D;
typedef typename Kokkos::View<double***,    FS_LAYOUT> FS3D;
typedef typename Kokkos::View<double**,     FS_LAYOUT> FS2D;
typedef typename Kokkos::View<double*,      FS_LAYOUT> FS1D;

typedef typename Kokkos::View<double******, FS_LAYOUT>::HostMirror FS6DH;
typedef typename Kokkos::View<double*****,  FS_LAYOUT>::HostMirror FS5DH;
typedef typename Kokkos::View<double****,   FS_LAYOUT>::HostMirror FS4DH;
typedef typename Kokkos::View<double***,    FS_LAYOUT>::HostMirror FS3DH;
typedef typename Kokkos::View<double**,     FS_LAYOUT>::HostMirror FS2DH;
typedef typename Kokkos::View<double*,      FS_LAYOUT>::HostMirror FS1DH;

typedef typename Kokkos::MDRangePolicy<Kokkos::Rank<2>> policy_f;
typedef typename Kokkos::MDRangePolicy<Kokkos::Rank<3>> policy_f3;
typedef typename Kokkos::MDRangePolicy<Kokkos::Rank<4>> policy_f4;

#endif
