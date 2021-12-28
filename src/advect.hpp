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

#ifndef ADVECT_H
#define ADVECT_H

struct advect2D {
  FS4D dvar;
  FS2D fluxx, fluxy;
  FS1D cd;
  int v;

  advect2D(FS4D d_, FS2D fx_, FS2D fy_, FS1D cd_, int v_)
      : dvar(d_), fluxx(fx_), fluxy(fy_), cd(cd_), v(v_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, const int j) const {
    dvar(i, j, 0, v) =
        -((fluxx(i, j) - fluxx(i - 1, j)) + (fluxy(i, j) - fluxy(i, j - 1)));
  }
};

#endif
