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

#include "hdf.hpp"
#include "kokkosTypes.hpp"
#include "hdf5.h"
#include "output.hpp"
#include "rkfunction.hpp"
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include "xdmf.hpp"
//#ifndef NOMPI
#include "mpi.h"
#include <vector>
#include "timer.hpp"
#include <filesystem>
#include "fmt/core.h"
//#endif

using namespace std;
using Kokkos::ALL;
using Kokkos::subview;

template <typename T, typename C>
void invertArray(int ndim, T* out, C* in){
  for (int i=0; i<ndim; ++i){
    out[i] = in[ndim-1-i];
  }
}

 void write_xmf(string fname, string hname, double time,
                struct inputConfig &cf, vector<string> vNames, vector<string> vxNames ){}

// open and hdf5 file for writing
hid_t openHDF5ForWrite(MPI_Comm comm, MPI_Info info, string fname){
  hid_t fid,pid;

  pid = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(pid, comm, info);
  fid = H5Fcreate(fname.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, pid);
  H5Pclose(pid);

  return fid;
}

// open an hdf5 file for reading only
hid_t openHDF5ForRead(string fname){
  hid_t fid;
  fid = H5Fopen(fname.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  return fid;
}

// close an hdf5 file
void close_h5(hid_t fid){
  H5Fclose(fid);
}

// function to return hdf5 type id given c++ type
template <typename H> hid_t getH5Type2();
template <> hid_t getH5Type2<float>(){ return H5T_NATIVE_FLOAT; }
template <> hid_t getH5Type2<double>(){ return H5T_NATIVE_DOUBLE; }
template <> hid_t getH5Type2<int>(){ return H5T_NATIVE_INT; }

 template <typename S>
 void write_h5(hid_t group_id, string dname, int ndim,
                    hsize_t* dims_global, hsize_t* dims_local, hsize_t* offset, S* data){}


// // writer function for hdf5
 template <typename S>
 void write_h52(hid_t group_id, string dname, int ndim,
                    size_t* in_dims_global, size_t* in_dims_local, size_t* in_offset, S* data){}

hdfWriter::hdfWriter(struct inputConfig cf, rk_func *f) {
  xdp = (double *)malloc(cf.ni * cf.nj * cf.nk * sizeof(double));
  xsp = (float *)malloc(cf.ni * cf.nj * cf.nk * sizeof(float));

  vsp = (float *)malloc(cf.nci * cf.ncj * cf.nck * sizeof(float));
  vdp = (double *)malloc(cf.nci * cf.ncj * cf.nck * sizeof(double));

  readV = (double *)malloc(cf.nci * cf.ncj * cf.nck * cf.nv * sizeof(double));

  gridH = Kokkos::create_mirror_view(f->grid);
  varH = Kokkos::create_mirror_view(f->var);
  varxH = Kokkos::create_mirror_view(f->varx);
}


 template<typename T>
 void hdfWriter::writeHDF(struct inputConfig cf, rk_func *f, int tdx,
                               double time, T* x, T* var, string name) {}

// check dimensions of restart file against expected values
void checkDataDimensions(hid_t filespace, int ndim, hsize_t* dims){

  int rank;
  //herr_t status;
  hsize_t dimsg[ndim];

  // get and check rank
  rank      = H5Sget_simple_extent_ndims(filespace);
  if (ndim != rank){
    printf("Number of dimensions of restart file different than expected\n");
    exit(EXIT_FAILURE);
  }

  // get dimensions
  H5Sget_simple_extent_dims(filespace, dimsg, NULL);
  //status    = H5Sget_simple_extent_dims(filespace, dimsg, NULL);

  // check dimensions
  for (int d=0; d<ndim; ++d){
    if (dims[d] != dimsg[d]){
      printf("Dimension extents of restart file different than expected\n");
      exit (EXIT_FAILURE);
    }
  }

}

// reader function for hdf5
template <typename T>
void read_H5(hid_t fid, string path, int ndim, hsize_t* dims, hsize_t* offset, hsize_t* count, T* data){

  // identifiers
  hid_t dset_id, filespace, memspace, dtype_id;
  //herr_t status;

  // get hd5 datatype
  dtype_id = getH5Type2<T>();

  // open dataset
  dset_id = H5Dopen(fid,path.c_str(),H5P_DEFAULT);

  // get filespace
  filespace = H5Dget_space(dset_id);    /* Get filespace handle first. */

  // check dimensions
  checkDataDimensions(filespace, ndim, dims);

  // specify memory space for this mpi rank
  memspace =  H5Screate_simple(ndim, count, NULL);

  // create hyperslab and read data
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);
  H5Dread(dset_id, dtype_id, memspace, filespace, H5P_DEFAULT, data);
  //status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);
  //status = H5Dread(dset_id, dtype_id, memspace, filespace, H5P_DEFAULT, data);

  // close identifiers
  H5Sclose(memspace);
  H5Sclose(filespace);
  H5Dclose(dset_id);
}

// Read Cell Data from Restart File
void hdfWriter::readSolution(struct inputConfig cf, FS4D &gridD, FS4D &varD) {

  hsize_t offset[cf.ndim], gridCount[cf.ndim], cellCount[cf.ndim];
  hsize_t gridDims[cf.ndim], cellDims[cf.ndim];
  invertArray(cf.ndim,offset,cf.subdomainOffset);
  invertArray(cf.ndim,gridCount,cf.localGridDims);
  invertArray(cf.ndim,cellCount,cf.localCellDims);
  invertArray(cf.ndim,gridDims,cf.globalGridDims);
  invertArray(cf.ndim,cellDims,cf.globalCellDims);

  // open restart file for reading
  hid_t fid;
  fid = openHDF5ForRead(cf.restartName);
  int idx;

  // read grid
  for (int v=0; v<cf.ndim; ++v){
    stringstream pathString;
    pathString << "/Grid/Dimension" << v;
    read_H5(fid, pathString.str(), cf.ndim, gridDims, offset, gridCount, xdp);
    for (int i = 0; i < cf.ni; ++i) {
      for (int j = 0; j < cf.nj; ++j) {
        for (int k = 0; k < cf.nk; ++k) {
          idx = (cf.ni * cf.nj) * k + cf.ni * j + i;
          gridH(i, j, k, v) = xdp[idx];
        }
      }
    }
  }

  // read cell data
  int koffset, ii, jj, kk;
  if (cf.ndim == 3) koffset = cf.ng;
  else koffset = 0;
  for (int v=0; v<cf.nvt; ++v){
    stringstream pathString;
    pathString << "/Solution/Variable" << setw(2) << setfill('0') << v;
    read_H5(fid, pathString.str(), cf.ndim, cellDims, offset, cellCount, vdp);
    for (int k = koffset; k < cf.nck + koffset; ++k) {
      for (int j = cf.ng; j < cf.ncj + cf.ng; ++j) {
        for (int i = cf.ng; i < cf.nci + cf.ng; ++i) {
          ii = i - cf.ng;
          jj = j - cf.ng;
          kk = k - koffset;
          idx = (cf.nci * cf.ncj) * kk + cf.nci * jj + ii;
          varH(i, j, k, v) = vdp[idx];
        }
      }
    }
  }
  
  H5Fclose(fid);
  Kokkos::deep_copy(gridD,gridH);
  Kokkos::deep_copy(varD,varH);
}

// Read Terrain data into grid array
void hdfWriter::readTerrain(struct inputConfig cf, FS4D &gridD) {

  hsize_t offset[2], gridCount[2], gridDims[2];
  invertArray(2,offset,cf.subdomainOffset);
  invertArray(2,gridCount,cf.localGridDims);
  invertArray(2,gridDims,cf.globalGridDims);

  // open restart file for reading
  hid_t fid;
  fid = openHDF5ForRead(cf.terrainName);
  int idx,ii,jj,kk;
  double h, dz;

  printf("#### %d: %d, %d\n",cf.rank,cf.kStart,cf.nk);

  // read grid
  read_H5(fid, "Height", 2, gridDims, offset, gridCount, xdp);
  for (int i = 0; i < cf.ni; ++i) {
    for (int j = 0; j < cf.nj; ++j) {
      idx = cf.ni * j + i;
      ii = cf.iStart + i;
      jj = cf.jStart + j;

      h = xdp[idx];
      dz = (cf.h-h)/cf.glbl_nck;
      for (int k = 0; k < cf.nk; ++k) {
        kk = cf.kStart + k;
        gridH(i, j, k, 0) = cf.tdx*ii;
        gridH(i, j, k, 1) = cf.tdy*jj;
        gridH(i, j, k, 2) = h + dz*kk;
      }
    }
  }

  H5Fclose(fid);
  Kokkos::deep_copy(gridD,gridH);
}

// write solution file
void hdfWriter::writeSolution(struct inputConfig cf, rk_func *f, int tdx,
                              double time) {
    writeHDF(cf, f, tdx, time, xsp, vsp, "sol");
}

// write restart file
void hdfWriter::writeRestart(struct inputConfig cf, rk_func *f, int tdx,
                             double time) {
    writeHDF(cf, f, tdx, time, xdp, vdp, "restart");
}
