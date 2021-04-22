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
#include "block.hpp"
#include "kokkosTypes.hpp"
#include "hdf5.h"
#include "output.hpp"
#include "rkfunction.hpp"
#include "lua.hpp"
#include "luaReader.hpp"
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <string>
#include <vector>
#include "xdmf.hpp"
#include "timer.hpp"
#include "fmt/core.h"
#include <filesystem>
//#ifndef NOMPI
#include "mpi.h"
//#endif

using namespace std;
using fmt::format;

template <typename T, typename C>
void blockWriter::reverseArray(int ndim, T* out, C* in){
  for (int i=0; i<ndim; ++i){
    out[i] = in[ndim-1-i];
  }
}

// open and hdf5 file for writing
hid_t blockWriter::openHDF5ForWrite(MPI_Comm comm, MPI_Info info, string fname){
  hid_t fid,pid;

  pid = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(pid, comm, info);
  fid = H5Fcreate(fname.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, pid);
  H5Pclose(pid);

  return fid;
}

// close an hdf5 file
void blockWriter::close_h5(hid_t fid){
  H5Fclose(fid);
}

size_t blockWriter::frq() { return freq; }

// function to return hdf5 type id given c++ type
template <> hid_t blockWriter::getH5Type<float>(){ return H5T_NATIVE_FLOAT; }
template <> hid_t blockWriter::getH5Type<double>(){ return H5T_NATIVE_DOUBLE; }
template <> hid_t blockWriter::getH5Type<int>(){ return H5T_NATIVE_INT; }

blockWriter::blockWriter(){}
blockWriter::blockWriter(struct inputConfig& cf, rk_func* f, string name_, string path_, bool avg_, size_t frq_,
  vector<size_t> start_, vector<size_t> end_, vector<size_t> stride_):
  name(name_),path(path_),avg(avg_),freq(frq_),gStart(start_),gEnd(end_),stride(stride_){

  int strideDelta;
  size_t gMin;

  pad = (int)log10(cf.tend) + 1;

  // initialize parameters
  lElems=1;
  lElemsG=1;
  myColor=MPI_UNDEFINED;
  slicePresent=true;

  for (int i=0; i<cf.ndim; ++i){
    //adjust block global size if it does not line up with strides
    strideDelta=(gEnd[i]-gStart[i])%stride[i];
    gEnd[i]-=strideDelta;
    gExt.push_back((gEnd[i]-gStart[i])/stride[i]+1);
    gExtG.push_back(gExt[i]+1);
  }

  for (int i=0; i<cf.ndim; ++i){
    lStart.push_back(0);
    lEnd.push_back(cf.localCellDims[i]-1);
    lEndG.push_back(cf.localCellDims[i]-1);
    lOffset.push_back(0);
  }

  for (int i=0; i<cf.ndim; ++i){
    lExt.push_back(lEnd[i]-lStart[i]+1);
    lExtG.push_back(lEnd[i]-lStart[i]+1);
  }

  // check if any part of the block slice is in this subdomain
  for (int i=0; i<cf.ndim; ++i){
    slicePresent = slicePresent 
                  && ( (cf.subdomainOffset[i]+cf.localCellDims[i]) >= gStart[i]
                  &&    cf.subdomainOffset[i] <= gEnd[i]);
  }

  // set color and create subcommunicator
  if (slicePresent){
    myColor=1;
  }
  // create slice communicator
  MPI_Comm_split(cf.comm,myColor,cf.rank,&sliceComm);

  if (slicePresent){
    for (int i=0; i<cf.ndim; ++i){
      offsetDelta=gStart[i]-cf.subdomainOffset[i]; //location of left edge of slice wrt left edge of subdomain
      if(offsetDelta > 0) lStart[i]=offsetDelta;   //if slice starts inside subdomain, set local start to slice edge
      else lOffset[i] = -(offsetDelta);            //if subdomain is within slice, adjust offset

      offsetDelta=(cf.subdomainOffset[i]+cf.localCellDims[i])-gEnd[i]; //location of right edge of slice wrt right edge of domain
      if(offsetDelta > 0){
        lEnd[i]=cf.localCellDims[i]-offsetDelta;   // if slice ends in subdomain, set local end to slide edge
      }

      //adjust start and offset if it does not land on a stride
      strideDelta = (lOffset[i]-gStart[i])%stride[i];
      lStart[i]+=(stride[i]-strideDelta)%stride[i];
      lOffset[i]+=(stride[i]-strideDelta)%stride[i];
      lOffset[i]/=stride[i];

      //adjust local end index if it does not align with a stride
      strideDelta=(lEnd[i]-lStart[i])%stride[i];
      lEnd[i]-=strideDelta;

      // reduce global dimensions if stride falls off the end of the domain
      if ((cf.subdomainOffset[i]+cf.localCellDims[i]) > gEnd[i]){
        if (((cf.localCellDims[i]-1)-lEnd[i]) < (stride[i]-1)){
          lEnd[i]-=stride[i];
          gEnd[i]-=stride[i];
          gExt[i]-=1;
          gExtG[i]-=1;
        }
        lEndG[i]=lEnd[i]+stride[i];
      }

      // compute local data and grid extend based on computed start and end indexex
      lExt[i]=(lEnd[i]-lStart[i])/stride[i]+1;
      lExtG[i]=(lEndG[i]-lStart[i])/stride[i]+1;

      // compute length of 1d buffer
      lElems*=lExt[i];
      lElemsG*=lExtG[i];
    }

    // processes must agree on global dimensions
    // last "righ-most" process may have had to reduce
    // global dimensions to accomadate the stride
    // so the minimum extent and end index are used
    gMin=0;
    for (int i=0; i<cf.ndim; ++i){
      MPI_Allreduce(&gEnd[i],&gMin,1,MPI_INT,MPI_MIN,sliceComm);
      gEnd[i]=gMin;

      MPI_Allreduce(&gExt[i],&gMin,1,MPI_INT,MPI_MIN,sliceComm);
      gExt[i]=gMin;

      MPI_Allreduce(&gExtG[i],&gMin,1,MPI_INT,MPI_MIN,sliceComm);
      gExtG[i]=gMin;
    }
  }
  
  // allocate pack buffers
  fill_n(back_inserter(varData),lElems,0.0);
  fill_n(back_inserter(gridData),lElemsG,0.0);

  // create kokkos host mirrors
  gridH = Kokkos::create_mirror_view(f->grid);
  varH = Kokkos::create_mirror_view(f->var);
  varxH = Kokkos::create_mirror_view(f->varx);
}

void blockWriter::write(struct inputConfig cf, rk_func *f, int tdx, double time) {
  string baseFormat = format("{{}}-{{:0{}d}}",pad);
  string blockBase = format(baseFormat,name,tdx);
  string hdfName   = format("{}.h5",blockBase);
  string hdfPath   = format("{}/{}",path,hdfName);
  string xmfPath   = format("{}/{}.xmf",path,blockBase);
  
  cf.log->message(format("[{}] Writing '{}'",cf.t,hdfPath));
  fiestaTimer writeTimer = fiestaTimer();

  if (myColor==1){
    Kokkos::deep_copy(gridH,f->grid);
    Kokkos::deep_copy(varH,f->var);
    Kokkos::deep_copy(varxH,f->varx);

    hid_t file_id = openHDF5ForWrite(sliceComm, MPI_INFO_NULL, hdfPath);

    hid_t group_id = H5Gcreate(file_id, "/Grid", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    for (int vn=0; vn<cf.ndim; ++vn){
      dataPack(cf.ndim, 0, lStart, lEndG, lExtG, stride, gridData, gridH,vn,false);
      write_h5<float>(group_id, format("Dimension{}",vn), cf.ndim, gExtG, lExtG, lOffset, gridData); 
    }
    H5Gclose(group_id);

    group_id = H5Gcreate(file_id, "/Solution", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    for (int vn=0; vn<cf.nvt; ++vn){
      dataPack(cf.ndim, cf.ng, lStart, lEnd, lExt, stride, varData, varH,vn,true);
      write_h5<float>(group_id, format("Variable{:02d}",vn), cf.ndim, gExt, lExt, lOffset, varData); 
    }
    for (size_t vn = 0; vn < f->varxNames.size(); ++vn) {
      dataPack(cf.ndim, cf.ng, lStart, lEnd, lExt, stride, varData, varxH,vn,true);
      write_h5<float>(group_id, format("Variable{:02d}",vn+cf.nvt), cf.ndim, gExt, lExt, lOffset, varData); 
    }
    H5Gclose(group_id);

    close_h5(file_id);

    filesystem::path hdfFilePath{hdfPath};
    auto fsize = filesystem::file_size(hdfFilePath);
    double ftime = writeTimer.check();
    double frate = (fsize/1048576.0)/ftime;

    if (fsize > 1073741824)
      cf.log->message(format("[{}] '{}': {:.2f} GiB in {:.2f}s ({:.2f} MiB/s)",cf.t,hdfPath,fsize/1073741824.0,ftime,frate));  //divide by bytes per GiB (1024*1024*1024)
    else
      cf.log->message(format("[{}] '{}': {:.2f} MiB in {:.2f}s ({:.2f} MiB/s)",cf.t,hdfPath,fsize/1048576.0,ftime,frate));     //divide by bytes per MiB (1024*1024)
  }

  cf.log->message(format("[{}] Writing '{}'",cf.t,xmfPath));
  if (myColor==1){
    writeXMF(xmfPath, hdfName, time, cf.ndim, gExt.data(), cf.nvt, f->varNames,f->varxNames);
  }
}

// writer function for hdf5
template <typename S>
void blockWriter::write_h5(hid_t group_id, string dname, int ndim,
                   vector<size_t> in_dims_global, vector<size_t> in_dims_local, vector<size_t> in_offset, vector<S>& data){

  vector<hsize_t> gdims(ndim,0);
  vector<hsize_t> ldims(ndim,0);
  vector<hsize_t> offset(ndim,0);

  // reverse order of array indexes to "c"
  std::reverse_copy(in_dims_global.begin(),in_dims_global.end(),gdims.begin());
  std::reverse_copy(in_dims_local.begin(),in_dims_local.end(),ldims.begin());
  std::reverse_copy(in_offset.begin(),in_offset.end(),offset.begin());

  // identifiers
  hid_t filespace, memspace, dset_id, plist_id, dtype_id;

  // get type id
  dtype_id = getH5Type<S>();

  // create global filespace and local memoryspace
  filespace = H5Screate_simple(ndim, gdims.data(), NULL);
  memspace  = H5Screate_simple(ndim, ldims.data(), NULL);

  // create dataset
  dset_id = H5Dcreate(group_id, dname.c_str(), dtype_id, filespace,
                                        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  H5Sclose(filespace);

  // select hyperslab in global filespace
  filespace = H5Dget_space(dset_id);
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset.data(), NULL, ldims.data(), NULL);
  //status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, dims_local, NULL);

  // create property list for collective dataset write
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

  // write data
  H5Dwrite(dset_id, dtype_id, memspace, filespace, plist_id, data.data());
  //status = H5Dwrite(dset_id, dtype_id, memspace, filespace, plist_id, data);

  // close sets, spaces and lists
  H5Dclose(dset_id);
  H5Sclose(filespace);
  H5Sclose(memspace);
  H5Pclose(plist_id);
}

template<typename T>
void blockWriter::dataPack(int ndim, int ng, const vector<size_t>& start, const vector<size_t>& end, const vector<size_t>& extent, const vector<size_t>& stride, vector<T>& dest, const FS4DH& source, const int vn, const bool average){
  int ii,jj,kk;
  int idx;
  double mean;

  if (ndim==3){
    for (int i = start[0]; i < end[0]+1; i+=stride[0]) {
      for (int j = start[1]; j < end[1]+1; j+=stride[1]) {
        for (int k = start[2]; k < end[2]+1; k+=stride[2]) {
          ii=(i-start[0])/stride[0];
          jj=(j-start[1])/stride[1];
          kk=(k-start[2])/stride[2];
          idx = (extent[0] * extent[1]) * kk + extent[0] * jj + ii;
          if (average){
            mean=0.0;
            for (int ia=0; ia<stride[0]; ++ia){
              for (int ja=0; ja<stride[1]; ++ja){
                for (int ka=0; ka<stride[2]; ++ka){
                  mean += source (i+ia+ng,j+ja+ng,k+ka+ng,vn);
                }
              }
            }
            dest[idx] = mean/(stride[0]*stride[1]*stride[2]);
          }else{
            dest[idx] = source(i+ng, j+ng, k+ng,vn);
          }
        }
      }
    }
  }else{
    for (int i = start[0]; i < end[0]+1; i+=stride[0]) {
      for (int j = start[1]; j < end[1]+1; j+=stride[1]) {
        ii=(i-start[0])/stride[0];
        jj=(j-start[1])/stride[1];
        idx = extent[0] * jj + ii;
        if (average){
          mean=0.0;
          for (int ia=0; ia<stride[0]; ++ia){
            for (int ja=0; ja<stride[1]; ++ja){
              mean += source (i+ia+ng,j+ja+ng,0,vn);
            }
          }
          dest[idx] = mean/(stride[0]*stride[1]);
        }else{
          dest[idx] = source(i+ng, j+ng, 0,vn);
        }
      }
    }
  }
}

blockWriter::~blockWriter(){
  //free(varData);
  //free(gridData);
}
