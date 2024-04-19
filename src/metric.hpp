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

#ifndef METRICS_HPP
#define METRICS_HPP
struct computeMetrics2D {

  FS4D metrics;
  FS4D grid;

  computeMetrics2D(FS4D m_, FS4D g_) : metrics(m_), grid(g_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, const int j) const {

    FSCAL x_xi, y_xi, x_et, y_et, jac;
    int ii, jj;

    ii = i - 3;
    jj = j - 3;

    x_xi = (grid(ii + 1, jj, 0, 0) + grid(ii + 1, jj + 1, 0, 0)) / 2.0 -
           (grid(ii, jj, 0, 0) + grid(ii, jj + 1, 0, 0)) / 2.0;
    y_xi = (grid(ii + 1, jj, 0, 1) + grid(ii + 1, jj + 1, 0, 1)) / 2.0 -
           (grid(ii, jj, 0, 1) + grid(ii, jj + 1, 0, 1)) / 2.0;

    x_et = (grid(ii, jj + 1, 0, 0) + grid(ii + 1, jj + 1, 0, 0)) / 2.0 -
           (grid(ii, jj, 0, 0) + grid(ii + 1, jj, 0, 0)) / 2.0;
    y_et = (grid(ii, jj + 1, 0, 1) + grid(ii + 1, jj + 1, 0, 1)) / 2.0 -
           (grid(ii, jj, 0, 1) + grid(ii + 1, jj, 0, 1)) / 2.0;

    jac = 1.0 / (x_xi * y_et - x_et * y_xi);

    metrics(i, j, 0, 0) = jac * y_et;
    metrics(i, j, 0, 1) = -jac * x_et;
    metrics(i, j, 1, 0) = -jac * y_xi;
    metrics(i, j, 1, 1) = jac * x_xi;
  }
};

struct computeMetrics3D {

  FS5D metrics;
  FS4D grid;

  computeMetrics3D(FS5D m_, FS4D g_) : metrics(m_), grid(g_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, const int j, const int k) const {

    int ii = i - 3;
    int jj = j - 3;
    int kk = k - 3;

    FSCAL xmx = (grid(ii, jj, kk, 0) + grid(ii, jj + 1, kk, 0) +
                  grid(ii, jj, kk + 1, 0) + grid(ii, jj + 1, kk + 1, 0)) /
                 4.0;
    FSCAL xpx =
        (grid(ii + 1, jj, kk, 0) + grid(ii + 1, jj + 1, kk, 0) +
         grid(ii + 1, jj, kk + 1, 0) + grid(ii + 1, jj + 1, kk + 1, 0)) /
        4.0;
    FSCAL xmy = (grid(ii, jj, kk, 1) + grid(ii, jj + 1, kk, 1) +
                  grid(ii, jj, kk + 1, 1) + grid(ii, jj + 1, kk + 1, 1)) /
                 4.0;
    FSCAL xpy =
        (grid(ii + 1, jj, kk, 1) + grid(ii + 1, jj + 1, kk, 1) +
         grid(ii + 1, jj, kk + 1, 1) + grid(ii + 1, jj + 1, kk + 1, 1)) /
        4.0;
    FSCAL xmz = (grid(ii, jj, kk, 2) + grid(ii, jj + 1, kk, 2) +
                  grid(ii, jj, kk + 1, 2) + grid(ii, jj + 1, kk + 1, 2)) /
                 4.0;
    FSCAL xpz =
        (grid(ii + 1, jj, kk, 2) + grid(ii + 1, jj + 1, kk, 2) +
         grid(ii + 1, jj, kk + 1, 2) + grid(ii + 1, jj + 1, kk + 1, 2)) /
        4.0;

    FSCAL zmx = (grid(ii, jj, kk, 0) + grid(ii + 1, jj, kk, 0) +
                  grid(ii, jj + 1, kk, 0) + grid(ii + 1, jj + 1, kk, 0)) /
                 4.0;
    FSCAL zpx =
        (grid(ii, jj, kk + 1, 0) + grid(ii + 1, jj, kk + 1, 0) +
         grid(ii, jj + 1, kk + 1, 0) + grid(ii + 1, jj + 1, kk + 1, 0)) /
        4.0;
    FSCAL zmy = (grid(ii, jj, kk, 1) + grid(ii + 1, jj, kk, 1) +
                  grid(ii, jj + 1, kk, 1) + grid(ii + 1, jj + 1, kk, 1)) /
                 4.0;
    FSCAL zpy =
        (grid(ii, jj, kk + 1, 1) + grid(ii + 1, jj, kk + 1, 1) +
         grid(ii, jj + 1, kk + 1, 1) + grid(ii + 1, jj + 1, kk + 1, 1)) /
        4.0;
    FSCAL zmz = (grid(ii, jj, kk, 2) + grid(ii + 1, jj, kk, 2) +
                  grid(ii, jj + 1, kk, 2) + grid(ii + 1, jj + 1, kk, 2)) /
                 4.0;
    FSCAL zpz =
        (grid(ii, jj, kk + 1, 2) + grid(ii + 1, jj, kk + 1, 2) +
         grid(ii, jj + 1, kk + 1, 2) + grid(ii + 1, jj + 1, kk + 1, 2)) /
        4.0;

    FSCAL ymx = (grid(ii, jj, kk, 0) + grid(ii + 1, jj, kk, 0) +
                  grid(ii, jj, kk + 1, 0) + grid(ii + 1, jj, kk + 1, 0)) /
                 4.0;
    FSCAL ypx =
        (grid(ii, jj + 1, kk, 0) + grid(ii + 1, jj + 1, kk, 0) +
         grid(ii, jj + 1, kk + 1, 0) + grid(ii + 1, jj + 1, kk + 1, 0)) /
        4.0;
    FSCAL ymy = (grid(ii, jj, kk, 1) + grid(ii + 1, jj, kk, 1) +
                  grid(ii, jj, kk + 1, 1) + grid(ii + 1, jj, kk + 1, 1)) /
                 4.0;
    FSCAL ypy =
        (grid(ii, jj + 1, kk, 1) + grid(ii + 1, jj + 1, kk, 1) +
         grid(ii, jj + 1, kk + 1, 1) + grid(ii + 1, jj + 1, kk + 1, 1)) /
        4.0;
    FSCAL ymz = (grid(ii, jj, kk, 2) + grid(ii + 1, jj, kk, 2) +
                  grid(ii, jj, kk + 1, 2) + grid(ii + 1, jj, kk + 1, 2)) /
                 4.0;
    FSCAL ypz =
        (grid(ii, jj + 1, kk, 2) + grid(ii + 1, jj + 1, kk, 2) +
         grid(ii, jj + 1, kk + 1, 2) + grid(ii + 1, jj + 1, kk + 1, 2)) /
        4.0;

    FSCAL xxi = xpx - xmx;
    FSCAL xet = ypx - ymx;
    FSCAL xzt = zpx - zmx;

    FSCAL yxi = xpy - xmy;
    FSCAL yet = ypy - ymy;
    FSCAL yzt = zpy - zmy;

    FSCAL zxi = xpz - xmz;
    FSCAL zet = ypz - ymz;
    FSCAL zzt = zpz - zmz;

    FSCAL jac =
        1.0 / (xxi * (yet * zzt - yzt * zet) + xet * (yzt * zxi - yxi * zzt) +
               xzt * (yxi * zet - yet * zxi));

    metrics(i, j, k, 0, 0) = jac * (yet * zzt - yzt * zet);
    metrics(i, j, k, 0, 1) = jac * (xzt * zet - xet * zzt);
    metrics(i, j, k, 0, 2) = jac * (xet * yzt - xzt * yet);

    metrics(i, j, k, 1, 0) = jac * (yzt * zxi - yxi * zzt);
    metrics(i, j, k, 1, 1) = jac * (xxi * zzt - xzt * zxi);
    metrics(i, j, k, 1, 2) = jac * (xzt * yxi - xxi * yzt);

    metrics(i, j, k, 2, 0) = jac * (yxi * zet - yet * zxi);
    metrics(i, j, k, 2, 1) = jac * (xet * zxi - xxi * zet);
    metrics(i, j, k, 2, 2) = jac * (xxi * yet - xet * yxi);
  }
};
#endif
