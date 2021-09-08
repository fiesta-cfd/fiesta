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

#include "reader.hpp"
#include "kokkosTypes.hpp"
#include "h5.hpp"
#include "rkfunction.hpp"
#include <sstream>
#include <string>
#include <vector>
#include "log2.hpp"
#include "debug.hpp"

void readRestart(struct inputConfig cf, rk_func *f) {
  h5Writer<double> writer;
  writer.openRead(cf.restartName);

  size_t idx;

  auto readV = (double *)malloc(cf.ni * cf.nj * cf.nk * sizeof(double));
  auto gridH = Kokkos::create_mirror_view(f->grid);
  auto varH = Kokkos::create_mirror_view(f->var);

  // read grid
  for (int v=0; v<cf.ndim; ++v){
    std::stringstream pathString;
    pathString << "/Grid/Dimension" << v;
    writer.read(pathString.str(), cf.ndim, cf.globalGridDims, cf.localGridDims, cf.subdomainOffset, readV);
    for (int i = 0; i < cf.ni; ++i) {
      for (int j = 0; j < cf.nj; ++j) {
        for (int k = 0; k < cf.nk; ++k) {
          idx = (cf.ni * cf.nj) * k + cf.ni * j + i;
          gridH(i, j, k, v) = readV[idx];
        }
      }
    }
  }

  // read cell data
  int koffset, ii, jj, kk;
  if (cf.ndim == 3) koffset = cf.ng;
  else koffset = 0;
  for (int v=0; v<cf.nvt; ++v){
    std::stringstream pathString;
    pathString << "/Solution/Variable" << setw(2) << setfill('0') << v;
    writer.read(pathString.str(), cf.ndim, cf.globalCellDims, cf.localCellDims, cf.subdomainOffset, readV);
    for (int k = koffset; k < cf.nck + koffset; ++k) {
      for (int j = cf.ng; j < cf.ncj + cf.ng; ++j) {
        for (int i = cf.ng; i < cf.nci + cf.ng; ++i) {
          ii = i - cf.ng;
          jj = j - cf.ng;
          kk = k - koffset;
          idx = (cf.nci * cf.ncj) * kk + cf.nci * jj + ii;
          varH(i, j, k, v) = readV[idx];
        }
      }
    }
  }

  //readAttribute("t",cf.t);
  //readAttribute("time",cf.time);
  
  writer.close();
  Kokkos::deep_copy(f->grid,gridH);
  Kokkos::deep_copy(f->var,varH);
}

void readTerrain(struct inputConfig cf, rk_func *f) {
  h5Writer<double> writer;
  writer.openRead(cf.restartName);

  auto readV = (double *)malloc(cf.ni * cf.nj * sizeof(double));
  auto gridH = Kokkos::create_mirror_view(f->grid);

  std::vector<size_t> gridDims,offset,gridCount;
  gridDims.push_back(cf.globalGridDims[0]);
  gridDims.push_back(cf.globalGridDims[1]);
  offset.push_back(cf.subdomainOffset[0]);
  offset.push_back(cf.subdomainOffset[1]);
  gridCount.push_back(cf.localGridDims[0]);
  gridCount.push_back(cf.localGridDims[1]);

  size_t idx,ii,jj,kk;
  double h, dz;

  // read grid
  writer.read("Height", 2, cf.globalCellDims, cf.localCellDims, cf.subdomainOffset, readV);
  for (int i = 0; i < cf.ni; ++i) {
    for (int j = 0; j < cf.nj; ++j) {
      idx = cf.ni * j + i;
      ii = cf.iStart + i;
      jj = cf.jStart + j;

      h = readV[idx];
      dz = (cf.h-h)/cf.glbl_nck;
      for (int k = 0; k < cf.nk; ++k) {
        kk = cf.kStart + k;
        gridH(i, j, k, 0) = cf.tdx*ii;
        gridH(i, j, k, 1) = cf.tdy*jj;
        gridH(i, j, k, 2) = h + dz*kk;
      }
    }
  }

  writer.close();
  Kokkos::deep_copy(f->grid,gridH);
}
