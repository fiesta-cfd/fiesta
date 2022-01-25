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

#ifndef presgrad_H
#define presgrad_H
struct applyPressureGradient2D {
  FS4D dvar;
  FS2D p;
  FS1D cd;

  applyPressureGradient2D(FS4D dvar_, FS2D p_, FS1D cd_)
      : dvar(dvar_), p(p_), cd(cd_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, const int j) const {
    FSCAL dx = cd(1);
    FSCAL dy = cd(2);
    // calculate pressure gradient across cell in each direction using 4th order
    // central difference
    FSCAL dxp = (p(i-2,j) - 8.0*p(i-1,j) + 8.0*p(i+1,j) - p(i+2,j)) /(12.0*dx);
    FSCAL dyp = (p(i,j-2) - 8.0*p(i,j-1) + 8.0*p(i,j+1) - p(i,j+2)) /(12.0*dy);

    dvar(i, j, 0, 0) -= dxp; 
    dvar(i, j, 0, 1) -= dyp;
  }
};

struct applyPressureGradient3D {
  FS4D dvar,varx;
  FS3D p;
  FS1D cd;

  applyPressureGradient3D(FS4D dvar_, FS4D varx_, FS3D p_, FS1D cd_)
      : dvar(dvar_), varx(varx_), p(p_), cd(cd_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, const int j, const int k) const {
    FSCAL dx = cd(1);
    FSCAL dy = cd(2);
    FSCAL dz = cd(2);

    // calculate pressure gradient across cell in each direction using 4th order
    // central difference
    FSCAL dxp = (p(i-2,j,k) - 8.0*p(i-1,j,k) + 8.0*p(i+1,j,k) - p(i+2,j,k))/(12.0*dx);
    FSCAL dyp = (p(i,j-2,k) - 8.0*p(i,j-1,k) + 8.0*p(i,j+1,k) - p(i,j+2,k))/(12.0*dy);
    FSCAL dzp = (p(i,j,k-2) - 8.0*p(i,j,k-1) + 8.0*p(i,j,k+1) - p(i,j,k+2))/(12.0*dz);

    dvar(i, j, k, 0) = dvar(i, j, k, 0) - dxp;
    dvar(i, j, k, 1) = dvar(i, j, k, 1) - dyp;
    dvar(i, j, k, 2) = dvar(i, j, k, 2) - dzp;

    //varx(i,j,k,20) = -dyp;

  }
};
struct applyGenPressureGradient3D {

  FS4D dvar;
  FS3D p;
  FS5D m;
  FS1D cd;

  applyGenPressureGradient3D(FS4D dvar_, FS5D m_, FS3D p_, FS1D cd_)
      : dvar(dvar_), m(m_), p(p_), cd(cd_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, const int j, const int k) const {

    FSCAL dxipxix;
    FSCAL detpetx;
    FSCAL dztpztx;
    FSCAL dxipxiy;
    FSCAL detpety;
    FSCAL dztpzty;
    FSCAL dxipxiz;
    FSCAL detpetz;
    FSCAL dztpztz;

    dxipxix = (p(i - 2, j, k) * m(i - 2, j, k, 0, 0) -
               8.0 * p(i - 1, j, k) * m(i - 1, j, k, 0, 0) +
               8.0 * p(i + 1, j, k) * m(i + 1, j, k, 0, 0) -
               p(i + 2, j, k) * m(i + 2, j, k, 0, 0)) /
              (12.0);
    detpetx = (p(i, j - 2, k) * m(i, j - 2, k, 1, 0) -
               8.0 * p(i, j - 1, k) * m(i, j - 1, k, 1, 0) +
               8.0 * p(i, j + 1, k) * m(i, j + 1, k, 1, 0) -
               p(i, j + 2, k) * m(i, j + 2, k, 1, 0)) /
              (12.0);
    dztpztx = (p(i, j, k - 2) * m(i, j, k - 2, 2, 0) -
               8.0 * p(i, j, k - 1) * m(i, j, k - 1, 2, 0) +
               8.0 * p(i, j, k + 1) * m(i, j, k + 1, 2, 0) -
               p(i, j, k + 2) * m(i, j, k + 2, 2, 0)) /
              (12.0);

    dxipxiy = (p(i - 2, j, k) * m(i - 2, j, k, 0, 1) -
               8.0 * p(i - 1, j, k) * m(i - 1, j, k, 0, 1) +
               8.0 * p(i + 1, j, k) * m(i + 1, j, k, 0, 1) -
               p(i + 2, j, k) * m(i + 2, j, k, 0, 1)) /
              (12.0);
    detpety = (p(i, j - 2, k) * m(i, j - 2, k, 1, 1) -
               8.0 * p(i, j - 1, k) * m(i, j - 1, k, 1, 1) +
               8.0 * p(i, j + 1, k) * m(i, j + 1, k, 1, 1) -
               p(i, j + 2, k) * m(i, j + 2, k, 1, 1)) /
              (12.0);
    dztpzty = (p(i, j, k - 2) * m(i, j, k - 2, 2, 1) -
               8.0 * p(i, j, k - 1) * m(i, j, k - 1, 2, 1) +
               8.0 * p(i, j, k + 1) * m(i, j, k + 1, 2, 1) -
               p(i, j, k + 2) * m(i, j, k + 2, 2, 1)) /
              (12.0);

    dxipxiz = (p(i - 2, j, k) * m(i - 2, j, k, 0, 2) -
               8.0 * p(i - 1, j, k) * m(i - 1, j, k, 0, 2) +
               8.0 * p(i + 1, j, k) * m(i + 1, j, k, 0, 2) -
               p(i + 2, j, k) * m(i + 2, j, k, 0, 2)) /
              (12.0);
    detpetz = (p(i, j - 2, k) * m(i, j - 2, k, 1, 2) -
               8.0 * p(i, j - 1, k) * m(i, j - 1, k, 1, 2) +
               8.0 * p(i, j + 1, k) * m(i, j + 1, k, 1, 2) -
               p(i, j + 2, k) * m(i, j + 2, k, 1, 2)) /
              (12.0);
    dztpztz = (p(i, j, k - 2) * m(i, j, k - 2, 2, 2) -
               8.0 * p(i, j, k - 1) * m(i, j, k - 1, 2, 2) +
               8.0 * p(i, j, k + 1) * m(i, j, k + 1, 2, 2) -
               p(i, j, k + 2) * m(i, j, k + 2, 2, 2)) /
              (12.0);

    dvar(i, j, k, 0) = (dvar(i, j, k, 0) - (dxipxix + detpetx + dztpztx));
    dvar(i, j, k, 1) = (dvar(i, j, k, 1) - (dxipxiy + detpety + dztpzty));
    dvar(i, j, k, 2) = (dvar(i, j, k, 2) - (dxipxiz + detpetz + dztpztz));
  }
};
struct applyGenPressureGradient2D {

  FS4D dvar;
  FS2D p;
  FS4D m;
  FS1D cd;

  applyGenPressureGradient2D(FS4D dvar_, FS4D m_, FS2D p_, FS1D cd_)
      : dvar(dvar_), m(m_), p(p_), cd(cd_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, const int j) const {

    FSCAL dxipxix;
    FSCAL detpetx;
    FSCAL dxipxiy;
    FSCAL detpety;

    dxipxix = (p(i - 2, j) * m(i - 2, j, 0, 0) -
               8.0 * p(i - 1, j) * m(i - 1, j, 0, 0) +
               8.0 * p(i + 1, j) * m(i + 1, j, 0, 0) -
               p(i + 2, j) * m(i + 2, j, 0, 0)) /
              (12.0);
    detpetx = (p(i, j - 2) * m(i, j - 2, 1, 0) -
               8.0 * p(i, j - 1) * m(i, j - 1, 1, 0) +
               8.0 * p(i, j + 1) * m(i, j + 1, 1, 0) -
               p(i, j + 2) * m(i, j + 2, 1, 0)) /
              (12.0);
    dxipxiy = (p(i - 2, j) * m(i - 2, j, 0, 1) -
               8.0 * p(i - 1, j) * m(i - 1, j, 0, 1) +
               8.0 * p(i + 1, j) * m(i + 1, j, 0, 1) -
               p(i + 2, j) * m(i + 2, j, 0, 1)) /
              (12.0);
    detpety = (p(i, j - 2) * m(i, j - 2, 1, 1) -
               8.0 * p(i, j - 1) * m(i, j - 1, 1, 1) +
               8.0 * p(i, j + 1) * m(i, j + 1, 1, 1) -
               p(i, j + 2) * m(i, j + 2, 1, 1)) /
              (12.0);

    dvar(i, j, 0, 0) = (dvar(i, j, 0, 0) - (dxipxix + detpetx));
    dvar(i, j, 0, 1) = (dvar(i, j, 0, 1) - (dxipxiy + detpety));
  }
};

struct applyPressure {

  FS4D dvar;
  FS3D p;
  Kokkos::View<FSCAL *> cd;

  applyPressure(FS4D dvar_, FS3D p_, Kokkos::View<FSCAL *> cd_)
      : dvar(dvar_), p(p_), cd(cd_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, const int j, const int k) const {
    FSCAL dxp = (p(i - 2, j, k) - 8.0 * p(i - 1, j, k) + 8.0 * p(i + 1, j, k) -
                  p(i + 2, j, k)) /
                 (12.0 * cd(1));
    FSCAL dyp = (p(i, j - 2, k) - 8.0 * p(i, j - 1, k) + 8.0 * p(i, j + 1, k) -
                  p(i, j + 2, k)) /
                 (12.0 * cd(2));
    FSCAL dzp = (p(i, j, k - 2) - 8.0 * p(i, j, k - 1) + 8.0 * p(i, j, k + 1) -
                  p(i, j, k + 2)) /
                 (12.0 * cd(3));

    dvar(i, j, k, 0) = dvar(i, j, k, 0) - dxp;
    dvar(i, j, k, 1) = dvar(i, j, k, 1) - dyp;
    dvar(i, j, k, 2) = dvar(i, j, k, 2) - dzp;
  }
};
#endif
