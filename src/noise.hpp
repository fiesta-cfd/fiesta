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

struct detectNoise2D {
  FS4D var;
  FS2D_I noise;
  int v;
  double dh;
  double coff;
  Kokkos::View<double *> cd;

  detectNoise2D(FS4D var_, FS2D_I n_, double dh_, double c_,
                Kokkos::View<double *> cd_, int v_)
      : var(var_), noise(n_), dh(dh_), coff(c_), cd(cd_), v(v_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int ii, const int jj) const {

    int nv = (int)cd(0) + 3;
    int ng = (int)cd(5);

    double dx = cd(1);
    double dy = cd(2);

    int i = 2 * ii + ng;
    int j = 2 * jj + ng;

    double a = sqrt(6.0 * dx * dy);

    double c =
        -(a / 192.0) * (var(i - 1, j + 1, 0, v) + 2 * var(i, j + 1, 0, v) +
                        var(i + 1, j + 1, 0, v) + 2 * var(i - 1, j, 0, v) -
                        12 * var(i, j, 0, v) + 2 * var(i + 1, j, 0, v) +
                        var(i - 1, j - 1, 0, v) + 2 * var(i, j - 1, 0, v) +
                        var(i + 1, j - 1, 0, v));

    double cref = dh * a / 16.0;

    //        var(i-1,j+1,0,11) = 0.0;
    //        var(i  ,j+1,0,11) = 0.0;
    //        var(i+1,j+1,0,11) = 0.0;
    //        var(i-1,j  ,0,11) = 0.0;
    //        var(i  ,j  ,0,11) = 0.0;
    //        var(i+1,j  ,0,11) = 0.0;
    //        var(i-1,j-1,0,11) = 0.0;
    //        var(i  ,j-1,0,11) = 0.0;
    //        var(i+1,j-1,0,11) = 0.0;

    noise(i - 1, j + 1) = 0;
    noise(i, j + 1) = 0;
    noise(i + 1, j + 1) = 0;
    noise(i - 1, j) = 0;
    noise(i, j) = 0;
    noise(i + 1, j) = 0;
    noise(i - 1, j - 1) = 0;
    noise(i, j - 1) = 0;
    noise(i + 1, j - 1) = 0;

    if (abs(c) >= cref) {
      if (coff == 0.0) {
        noise(i - 1, j + 1) = 1;
        noise(i, j + 1) = 1;
        noise(i + 1, j + 1) = 1;
        noise(i - 1, j) = 1;
        noise(i, j) = 1;
        noise(i + 1, j) = 1;
        noise(i - 1, j - 1) = 1;
        noise(i, j - 1) = 1;
        noise(i + 1, j - 1) = 1;
      } else {
        if (var(i - 1, j + 1, 0, nv + 1) < coff)
          noise(i - 1, j + 1) = 1;
        if (var(i, j + 1, 0, nv + 1) < coff)
          noise(i, j + 1) = 1;
        if (var(i + 1, j + 1, 0, nv + 1) < coff)
          noise(i + 1, j + 1) = 1;
        if (var(i - 1, j, 0, nv + 1) < coff)
          noise(i - 1, j) = 1;
        if (var(i, j, 0, nv + 1) < coff)
          noise(i, j) = 1;
        if (var(i + 1, j, 0, nv + 1) < coff)
          noise(i + 1, j) = 1;
        if (var(i - 1, j - 1, 0, nv + 1) < coff)
          noise(i - 1, j - 1) = 1;
        if (var(i, j - 1, 0, nv + 1) < coff)
          noise(i, j - 1) = 1;
        if (var(i + 1, j - 1, 0, nv + 1) < coff)
          noise(i + 1, j - 1) = 1;
      }
    }

    //        var(i-1,j+1,0,11) = noise(i-1,j+1);
    //        var(i  ,j+1,0,11) = noise(i  ,j+1);
    //        var(i+1,j+1,0,11) = noise(i+1,j+1);
    //        var(i-1,j  ,0,11) = noise(i-1,j  );
    //        var(i  ,j  ,0,11) = noise(i  ,j  );
    //        var(i+1,j  ,0,11) = noise(i+1,j  );
    //        var(i-1,j-1,0,11) = noise(i-1,j-1);
    //        var(i  ,j-1,0,11) = noise(i  ,j-1);
    //        var(i+1,j-1,0,11) = noise(i+1,j-1);
    //
    //        var(i-1,j+1,0,10) = c;
    //        var(i  ,j+1,0,10) = c;
    //        var(i+1,j+1,0,10) = c;
    //        var(i-1,j  ,0,10) = c;
    //        var(i  ,j  ,0,10) = c;
    //        var(i+1,j  ,0,10) = c;
    //        var(i-1,j-1,0,10) = c;
    //        var(i  ,j-1,0,10) = c;
    //        var(i+1,j-1,0,10) = c;
  }
};

struct removeNoise2D {
  FS4D dvar;
  FS4D var;
  FS2D_I noise;
  int v;
  double dt;
  Kokkos::View<double *> cd;

  removeNoise2D(FS4D dvar_, FS4D var_, FS2D_I n_, double dt_,
                Kokkos::View<double *> cd_, int v_)
      : dvar(dvar_), var(var_), noise(n_), dt(dt_), cd(cd_), v(v_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, const int j) const {

    int ng = (int)cd(5);

    double dx = cd(1);
    double dy = cd(2);

    double lgrad = (var(i, j, 0, v) - var(i - 1, j, 0, v)) / dx;
    double rgrad = (var(i + 1, j, 0, v) - var(i, j, 0, v)) / dx;
    double bgrad = (var(i, j, 0, v) - var(i, j - 1, 0, v)) / dy;
    double tgrad = (var(i, j + 1, 0, v) - var(i, j, 0, v)) / dy;

    dvar(i, j, 0, v) = dt * (dx * dx + dy * dy) * noise(i, j) *
                       ((rgrad - lgrad) / dx + (tgrad - bgrad) / dy);
    //        var(i,j,0,12) = dvar(i,j,0,v);
  }
};

struct updateNoise2D {
  FS4D dvar;
  FS4D var;
  int v;

  updateNoise2D(FS4D dvar_, FS4D var_, int v_)
      : dvar(dvar_), var(var_), v(v_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, const int j) const {

    var(i, j, 0, v) += dvar(i, j, 0, v);
  }
};
#endif
