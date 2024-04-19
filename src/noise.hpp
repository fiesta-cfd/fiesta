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

#ifndef NOISE_HPP
#define NOISE_HPP
#include "log2.hpp"

// 3D

struct detectNoise3D {
  FS4D var,varx;
  FS3D_I noise;
  int v;
  FSCAL dh;
  FSCAL coff;
  Kokkos::View<FSCAL *> cd;

  detectNoise3D(FS4D var_, FS4D varx_, FS3D_I n_, FSCAL dh_, FSCAL c_,
                Kokkos::View<FSCAL *> cd_, int v_)
      : var(var_), varx(varx_), noise(n_), dh(dh_), coff(c_), cd(cd_), v(v_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int ii, const int jj, const int kk) const {

    int nv = (int)cd(0) + 4;
    int ng = (int)cd(5);

    FSCAL dx = cd(1);
    FSCAL dy = cd(2);
    FSCAL dz = cd(3);

    int i = 2 * ii + ng;
    int j = 2 * jj + ng;
    int k = 2 * kk + ng;

    int idx=0;
    int jdx=0;
    int kdx=0;

    FSCAL a = sqrt(165.0 * dx * dy * dz)/15840.0;

    FSCAL c = -a*(
               //  2.0*var(i-1,j+1,k+1,v) +   3.0*var(i,j+1,k+1,v) + 2.0*var(i+1,j+1,k+1,v)
               //+ 3.0*var(i-1,j  ,k+1,v) +   6.0*var(i,j  ,k+1,v) + 3.0*var(i+1,j  ,k+1,v)
               //+ 2.0*var(i-1,j-1,k+1,v) +   3.0*var(i,j-1,k+1,v) + 2.0*var(i+1,j-1,k+1,v)

               //+ 3.0*var(i-1,j+1,k  ,v) +   6.0*var(i,j+1,k  ,v) + 3.0*var(i+1,j+1,k  ,v)
               //+ 6.0*var(i-1,j  ,k  ,v) + -88.0*var(i,j  ,k  ,v) + 6.0*var(i+1,j  ,k  ,v)
               //+ 3.0*var(i-1,j-1,k  ,v) +   6.0*var(i,j-1,k  ,v) + 3.0*var(i+1,j-1,k  ,v)

               //+ 2.0*var(i-1,j+1,k-1,v) +   3.0*var(i,j+1,k-1,v) + 2.0*var(i+1,j+1,k-1,v)
               //+ 3.0*var(i-1,j  ,k-1,v) +   6.0*var(i,j  ,k-1,v) + 3.0*var(i+1,j  ,k-1,v)
               //+ 2.0*var(i-1,j-1,k-1,v) +   3.0*var(i,j-1,k-1,v) + 2.0*var(i+1,j-1,k-1,v)
               //);
                 2.0*var(i-1,j+1,k+1,v) +   2.0*var(i,j+1,k+1,v) + 2.0*var(i+1,j+1,k+1,v)
               + 2.0*var(i-1,j  ,k+1,v) +   8.0*var(i,j  ,k+1,v) + 2.0*var(i+1,j  ,k+1,v)
               + 2.0*var(i-1,j-1,k+1,v) +   2.0*var(i,j-1,k+1,v) + 2.0*var(i+1,j-1,k+1,v)

               + 2.0*var(i-1,j+1,k  ,v) +   8.0*var(i,j+1,k  ,v) + 2.0*var(i+1,j+1,k  ,v)
               + 8.0*var(i-1,j  ,k  ,v) + -88.0*var(i,j  ,k  ,v) + 8.0*var(i+1,j  ,k  ,v)
               + 2.0*var(i-1,j-1,k  ,v) +   8.0*var(i,j-1,k  ,v) + 2.0*var(i+1,j-1,k  ,v)

               + 2.0*var(i-1,j+1,k-1,v) +   2.0*var(i,j+1,k-1,v) + 2.0*var(i+1,j+1,k-1,v)
               + 2.0*var(i-1,j  ,k-1,v) +   8.0*var(i,j  ,k-1,v) + 2.0*var(i+1,j  ,k-1,v)
               + 2.0*var(i-1,j-1,k-1,v) +   2.0*var(i,j-1,k-1,v) + 2.0*var(i+1,j-1,k-1,v)
               );

    FSCAL cref = dh * sqrt(dx*dy*dz/12.0);

    for (idx=-1;idx<2;++idx){
      for (jdx=-1;jdx<2;++jdx){
        for (kdx=-1;kdx<2;++kdx){
          noise(i+idx,j+jdx,k+kdx)=0;
        }
      }
    }

    if (abs(c) >= cref) {

      if (coff == 0.0) {
        for (idx=-1;idx<2;++idx){
          for (jdx=-1;jdx<2;++jdx){
            for (kdx=-1;kdx<2;++kdx){
              noise(i+idx,j+jdx,k+kdx)=1;
            }
          }
        }

      }else{

        for (idx=-1;idx<2;++idx){
          for (jdx=-1;jdx<2;++jdx){
            for (kdx=-1;kdx<2;++kdx){
              if(abs(var(i+idx,j+jdx,k+kdx,nv+1)) < coff){
                noise(i+idx,j+jdx,k+kdx)=1;
              }
            }
          }
        }

      }
      
    } // end noise detected

    //if(v==1)
    //for (idx=-1;idx<2;++idx){
    //  for (jdx=-1;jdx<2;++jdx){
    //    for (kdx=-1;kdx<2;++kdx){
    //      varx(i+idx,j+jdx,k+kdx,6)=c;
    //    }
    //  }
    //}

  }
};

struct removeNoise3D {
  FS4D dvar;
  FS4D var;
  FS4D varx;
  FS3D_I noise;
  FSCAL dt;
  Kokkos::View<FSCAL *> cd;
  int v;

  removeNoise3D(FS4D dvar_, FS4D var_, FS4D varx_, FS3D_I n_, FSCAL dt_,
                Kokkos::View<FSCAL *> cd_, int v_)
      : dvar(dvar_), var(var_), varx(varx_), noise(n_), dt(dt_), cd(cd_), v(v_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, const int j, const int k) const {

    int ng = (int)cd(5);

    FSCAL dx = cd(1);
    FSCAL dy = cd(2);
    FSCAL dz = cd(3);
    FSCAL lap,dnoise;

    lap = (var(i-1,j,k,v)-2*var(i,j,k,v)+var(i+1,j,k,v))/(dx*dx)
        + (var(i,j-1,k,v)-2*var(i,j,k,v)+var(i,j+1,k,v))/(dy*dy)
        + (var(i,j,k-1,v)-2*var(i,j,k,v)+var(i,j,k+1,v))/(dz*dz);

    dnoise = dt*(dx*dx+dy*dy+dz*dz)*noise(i,j,k)*lap;
    var(i,j,k,v) += dnoise;

    //if(v==1){
    //  varx(i,j,k,7) = noise(i,j,k);
    //  varx(i,j,k,8) = dnoise;
    //}
  }
};

// 2D

struct detectNoise2D {
  FS4D var,varx;
  FS2D_I noise;
  int v;
  FSCAL dh;
  FSCAL coff;
  Kokkos::View<FSCAL *> cd;

  detectNoise2D(FS4D var_, FS4D varx_, FS2D_I n_, FSCAL dh_, FSCAL c_,
                Kokkos::View<FSCAL *> cd_, int v_)
      : var(var_), varx(varx_), noise(n_), dh(dh_), coff(c_), cd(cd_), v(v_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int ii, const int jj) const {

    int nv = (int)cd(0) + 3;
    int ng = (int)cd(5);

    FSCAL dx = cd(1);
    FSCAL dy = cd(2);

    int i = 2 * ii + ng;
    int j = 2 * jj + ng;

    FSCAL a = sqrt(6.0 * dx * dy);

    FSCAL c =
        -(a / 192.0) * (var(i - 1, j + 1, 0, v) + 2 * var(i, j + 1, 0, v) +
                        var(i + 1, j + 1, 0, v) + 2 * var(i - 1, j, 0, v) -
                        12 * var(i, j, 0, v) + 2 * var(i + 1, j, 0, v) +
                        var(i - 1, j - 1, 0, v) + 2 * var(i, j - 1, 0, v) +
                        var(i + 1, j - 1, 0, v));

    FSCAL cref = dh * a / 16.0;

    for (int idx=-1;idx<2;++idx)
      for (int jdx=-1;jdx<2;++jdx)
        noise(i+idx,j+jdx)=0;

    if (abs(c) >= cref) {
      if (coff == 0.0) {
        for (int idx=-1;idx<2;++idx)
          for (int jdx=-1;jdx<2;++jdx)
            noise(i+idx,j+jdx)=1;
      } else {
        for (int idx=-1;idx<2;++idx)
          for (int jdx=-1;jdx<2;++jdx)
            if(var(i+idx,j+jdx,0,nv+1) < coff)
              noise(i+idx,j+jdx)=1;
      }
    }
    //for (int idx=-1;idx<2;++idx)
    //  for (int jdx=-1;jdx<2;++jdx){
    //    varx(i+idx,j+jdx,0,7)=noise(i+idx,j+jdx);
    //    varx(i+idx,j+jdx,0,8)=c;
    //  }
  }
};

struct removeNoise2D {
  FS4D dvar;
  FS4D var,varx;
  FS2D_I noise;
  int v;
  FSCAL dt;
  Kokkos::View<FSCAL *> cd;

  removeNoise2D(FS4D dvar_, FS4D var_, FS4D varx_, FS2D_I n_, FSCAL dt_,
                Kokkos::View<FSCAL *> cd_, int v_)
      : dvar(dvar_), var(var_), varx(varx_), noise(n_), dt(dt_), cd(cd_), v(v_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, const int j) const {

    int ng = (int)cd(5);

    FSCAL dx = cd(1);
    FSCAL dy = cd(2);

    FSCAL lgrad = (var(i, j, 0, v) - var(i - 1, j, 0, v)) / dx;
    FSCAL rgrad = (var(i + 1, j, 0, v) - var(i, j, 0, v)) / dx;
    FSCAL bgrad = (var(i, j, 0, v) - var(i, j - 1, 0, v)) / dy;
    FSCAL tgrad = (var(i, j + 1, 0, v) - var(i, j, 0, v)) / dy;

    FSCAL dvar = dt * (dx * dx + dy * dy) * noise(i, j) *
                       ((rgrad - lgrad) / dx + (tgrad - bgrad) / dy);
    var(i,j,0,v) += dvar;
    //varx(i,j,0,9) = dvar;
  }
};
#endif
