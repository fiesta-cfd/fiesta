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
#include <filesystem>
#include "fiesta.hpp"
#include <fmt/core.h>

void readRestart(struct inputConfig &cf, std::unique_ptr<class rk_func>&f) {
  if (cf.rank==0){
    if (!std::filesystem::exists(cf.restartName)){
      Log::error("Restart file '{}' does not exist.",cf.restartName);
      exit(EXIT_FAILURE);
    }
  }

  h5Writer<FSCAL> writer;
  writer.openRead(cf.restartName);

  size_t idx;
  std::string path;

  auto readV = (FSCAL *)malloc(cf.ni * cf.nj * cf.nk * sizeof(FSCAL));
  auto varH = Kokkos::create_mirror_view(f->var);

  // read grid
  if (cf.grid==1){
    auto gridH = Kokkos::create_mirror_view(f->grid);
    for (int v=0; v<cf.ndim; ++v){
      path = fmt::format("/Grid/Dimension{}",v);
      writer.read(path, cf.ndim, cf.globalGridDims, cf.localGridDims, cf.subdomainOffset, readV);
      for (int i = 0; i < cf.ni; ++i) {
        for (int j = 0; j < cf.nj; ++j) {
          for (int k = 0; k < cf.nk; ++k) {
            idx = (cf.ni * cf.nj) * k + cf.ni * j + i;
            gridH(i, j, k, v) = readV[idx];
          }
        }
      }
    }
    Kokkos::deep_copy(f->grid,gridH);
  }

  // read cell data
  int koffset, ii, jj, kk;
  if (cf.ndim == 3) koffset = cf.ng;
  else koffset = 0;
  for (int v=0; v<cf.nvt; ++v){
    path = fmt::format("/Solution/{}",f->varNames[v]);
    writer.read(path, cf.ndim, cf.globalCellDims, cf.localCellDims, cf.subdomainOffset, readV);
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
  Kokkos::deep_copy(f->var,varH);

  std::string temp_title;
  int restart_version;

  writer.readAttribute("restart_file_version",restart_version);
  if(restart_version != FIESTA_RESTART_VERSION){
    Log::error("Cannot read restart file '{}' with version {}. Only version {} files are readable by Fiesta version {}",
        cf.restartName,restart_version,FIESTA_RESTART_VERSION,FIESTA_VERSION);
    exit(EXIT_FAILURE);
  }


  writer.readAttribute("title",temp_title);

  if(cf.restartReset==false){
    writer.readAttribute("time_index",cf.tstart);
    writer.readAttribute("time",cf.time);
  }else{
    Log::message("Restart reset enabled, using index and time from input file.");
  }

  if(cf.tinterval){
    cf.tend=cf.tstart+cf.nt;
  }else{
    cf.nt=cf.tend-cf.tstart;
  }
  cf.t=cf.tstart;
  Log::message("Restart Title: '{}'",temp_title);
  Log::message("Restart Properties: t={} time={:.2g}",cf.tstart,cf.time);
  
  writer.close();
}

void readTerrain(struct inputConfig &cf, std::unique_ptr<class rk_func>&f) {
  h5Writer<FSCAL> writer;
  writer.openRead(cf.restartName);

  auto readV = (FSCAL *)malloc(cf.ni * cf.nj * sizeof(FSCAL));
  auto gridH = Kokkos::create_mirror_view(f->grid);

  std::vector<size_t> gridDims,offset,gridCount;
  gridDims.push_back(cf.globalGridDims[0]);
  gridDims.push_back(cf.globalGridDims[1]);
  offset.push_back(cf.subdomainOffset[0]);
  offset.push_back(cf.subdomainOffset[1]);
  gridCount.push_back(cf.localGridDims[0]);
  gridCount.push_back(cf.localGridDims[1]);

  size_t idx,ii,jj,kk;
  FSCAL h, dz;

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
