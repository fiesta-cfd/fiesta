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

#include "kokkosTypes.hpp"
#ifndef FST_H
#define FST_H

#include "input.hpp"
#include "writer.hpp"

class hdfWriter : public writer {

public:
  hdfWriter(struct inputConfig, rk_func *f);
  ~hdfWriter();
  //void writeGrid(struct inputConfig cf, const FS4D gridD, const char *fname);
  //void writeSPGrid(struct inputConfig cf, const FS4D gridD, const char *fname);

  void writeSolution(struct inputConfig cf, rk_func *f, int tdx, double time);
  void writeRestart(struct inputConfig cf, rk_func *f, int tdx, double time);

  void readSolution(struct inputConfig cf, FS4D &deviceG, FS4D &deviceV);
  void readTerrain(struct inputConfig cf, FS4D &deviceG);

  template<typename T>
  void writeHDF(struct inputConfig cf, rk_func *f, int tdx,
                              double time, T* x, T* var, string name);

private:
  double *xdp;
  double *ydp;
  double *zdp;
  float *xsp;
  float *ysp;
  float *zsp;

  float* psp;
  double* pdp;
  int* pi;

  double *vdp;
  float *vsp;
  double *readV;
  FS4DH gridH;
  FS4DH varH;
  FS4DH varxH;
};

#endif
