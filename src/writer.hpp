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

#ifndef WRITER_H
#define WRITER_H

#include "Kokkos_Core.hpp"
#include "kokkosTypes.hpp"
#include "input.hpp"
#include "rkfunction.hpp"
#include "timer.hpp"
#include <map>
#include <string>

class writer {

public:
  // writer(struct inputConfig &cf_, FS4D gridD, FS4D varD);
  writer();
  ~writer();

  //virtual void writeGrid(struct inputConfig cf, const FS4D gridD,
  //                       const char *fname) = 0;
  //virtual void writeSPGrid(struct inputConfig cf, const FS4D gridD,
  //                         const char *fname) = 0;

  virtual void writeSolution(struct inputConfig cf, class rk_func *f, int tdx,
                             double time) = 0;
  virtual void writeRestart(struct inputConfig cf, class rk_func *f, int tdx,
                            double time) = 0;

  virtual void readSolution(struct inputConfig cf, FS4D &gridD, FS4D &varD) = 0;
  virtual void readTerrain(struct inputConfig cf, FS4D &gridD) = 0;

  // protected:
  //    double *xdp;
  //    double *ydp;
  //    double *zdp;
  //    float *xsp;
  //    float *ysp;
  //    float *zsp;
  //
  //    double *v;
  //    float *vsp;
  //    double *readV;
  //
  //    FS4DH gridH;
  //    FS4DH varH;
};

#endif
