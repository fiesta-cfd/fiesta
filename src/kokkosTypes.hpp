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

#ifndef FIESTA_HPP
#define FIESTA_HPP

#ifdef HAVE_SINGLE
#define FSCAL float
#define MPI_FSCAL MPI_FLOAT
#else
#define FSCAL double
#define MPI_FSCAL MPI_DOUBLE
#endif

#include "Kokkos_Core.hpp"
#include "debug.hpp"
#include <fstream>
#include <iomanip>
#include <ostream>
#include <sstream>

//#define FS_LAYOUT  Kokkos::LayoutRight
//#define FS_LAYOUT  Kokkos::LayoutLeft
#define FS_LAYOUT Kokkos::DefaultExecutionSpace::array_layout

// FSCAL view types
typedef typename Kokkos::View<FSCAL ******, FS_LAYOUT> FS6D;
typedef typename Kokkos::View<FSCAL *****, FS_LAYOUT> FS5D;
typedef typename Kokkos::View<FSCAL ****, FS_LAYOUT> FS4D;
typedef typename Kokkos::View<FSCAL ***, FS_LAYOUT> FS3D;
typedef typename Kokkos::View<FSCAL **, FS_LAYOUT> FS2D;
typedef typename Kokkos::View<FSCAL *, FS_LAYOUT> FS1D;

typedef typename Kokkos::View<FSCAL ******, FS_LAYOUT>::HostMirror FS6DH;
typedef typename Kokkos::View<FSCAL *****, FS_LAYOUT>::HostMirror FS5DH;
typedef typename Kokkos::View<FSCAL ****, FS_LAYOUT>::HostMirror FS4DH;
typedef typename Kokkos::View<FSCAL ***, FS_LAYOUT>::HostMirror FS3DH;
typedef typename Kokkos::View<FSCAL **, FS_LAYOUT>::HostMirror FS2DH;
typedef typename Kokkos::View<FSCAL *, FS_LAYOUT>::HostMirror FS1DH;

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
