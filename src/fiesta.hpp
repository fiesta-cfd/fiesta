#ifndef FIESTA_HPP
#define FIESTA_HPP

#include "Kokkos_Core.hpp"
#include "debug.hpp"
#include <fstream>
#include <iomanip>
#include <ostream>
#include <sstream>

//#define FS_LAYOUT  Kokkos::LayoutRight
//#define FS_LAYOUT  Kokkos::LayoutLeft
#define FS_LAYOUT Kokkos::DefaultExecutionSpace::array_layout

// double view types
typedef typename Kokkos::View<double ******, FS_LAYOUT> FS6D;
typedef typename Kokkos::View<double *****, FS_LAYOUT> FS5D;
typedef typename Kokkos::View<double ****, FS_LAYOUT> FS4D;
typedef typename Kokkos::View<double ***, FS_LAYOUT> FS3D;
typedef typename Kokkos::View<double **, FS_LAYOUT> FS2D;
typedef typename Kokkos::View<double *, FS_LAYOUT> FS1D;

typedef typename Kokkos::View<double ******, FS_LAYOUT>::HostMirror FS6DH;
typedef typename Kokkos::View<double *****, FS_LAYOUT>::HostMirror FS5DH;
typedef typename Kokkos::View<double ****, FS_LAYOUT>::HostMirror FS4DH;
typedef typename Kokkos::View<double ***, FS_LAYOUT>::HostMirror FS3DH;
typedef typename Kokkos::View<double **, FS_LAYOUT>::HostMirror FS2DH;
typedef typename Kokkos::View<double *, FS_LAYOUT>::HostMirror FS1DH;

// int view types
typedef typename Kokkos::View<int ******, FS_LAYOUT> FS6D_I;
typedef typename Kokkos::View<int *****, FS_LAYOUT> FS5D_I;
typedef typename Kokkos::View<int ****, FS_LAYOUT> FS4D_I;
typedef typename Kokkos::View<int ***, FS_LAYOUT> FS3D_I;
typedef typename Kokkos::View<int **, FS_LAYOUT> FS2D_I;
typedef typename Kokkos::View<int *, FS_LAYOUT> FS1D_I;

typedef typename Kokkos::View<int ******, FS_LAYOUT>::HostMirror FS6DH_I;
typedef typename Kokkos::View<int *****, FS_LAYOUT>::HostMirror FS5DH_I;
typedef typename Kokkos::View<int ****, FS_LAYOUT>::HostMirror FS4DH_I;
typedef typename Kokkos::View<int ***, FS_LAYOUT>::HostMirror FS3DH_I;
typedef typename Kokkos::View<int **, FS_LAYOUT>::HostMirror FS2DH_I;
typedef typename Kokkos::View<int *, FS_LAYOUT>::HostMirror FS1DH_I;

typedef typename Kokkos::MDRangePolicy<Kokkos::Rank<2>> policy_f;
typedef typename Kokkos::MDRangePolicy<Kokkos::Rank<3>> policy_f3;
typedef typename Kokkos::MDRangePolicy<Kokkos::Rank<4>> policy_f4;
typedef typename Kokkos::MDRangePolicy<Kokkos::Rank<5>> policy_f5;

#endif
