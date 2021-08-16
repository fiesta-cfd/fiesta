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
#include "h5.hpp"
//#endif

using namespace std;
using fmt::format;

template <typename T>
size_t blockWriter<T>::frq() { return freq; }

template <typename T>
blockWriter<T>::blockWriter(){}

template <typename T>
blockWriter<T>::blockWriter(struct inputConfig& cf, rk_func* f, string name_, string path_, bool avg_, size_t frq_):
  name(name_),path(path_),avg(avg_),freq(frq_){

  if (cf.rank==0){
    if (!std::filesystem::exists(path)){
      Fiesta::Log::message("Creating directory: '{}'",path);
      std::filesystem::create_directories(path);
    }
  }

  pad = (int)log10(cf.tend) + 1;

  writeVarx=false;

  // initialize parameters
  lElems=1;
  lElemsG=1;
  myColor=1;
  slicePresent=true;
  MPI_Comm_split(cf.comm,myColor,cf.rank,&sliceComm);

  for (int i=0; i<cf.ndim; ++i){
    gStart.push_back(0);
    gEnd.push_back(cf.globalCellDims[i]-1);
    gExt.push_back(cf.globalCellDims[i]);
    gExtG.push_back(cf.globalCellDims[i]+1);
    stride.push_back(1);

    lStart.push_back(0);
    lEnd.push_back(cf.localCellDims[i]-1);
    lEndG.push_back(cf.localCellDims[i]);
    lExt.push_back(cf.localCellDims[i]);
    lExtG.push_back(cf.localCellDims[i]+1);

    lOffset.push_back(cf.subdomainOffset[i]);

    lElems*=lExt[i];
    lElemsG*=lExtG[i];
  }

  // allocate pack buffers
  fill_n(back_inserter(varData),lElems,0.0);
  fill_n(back_inserter(gridData),lElemsG,0.0);

  // create kokkos host mirrors
  gridH = Kokkos::create_mirror_view(f->grid);
  varH = Kokkos::create_mirror_view(f->var);
  varxH = Kokkos::create_mirror_view(f->varx);

  //Fiesta::Log::debug("{}: gStart={} gEnd={} stride={}",name,gStart,gEnd,stride);
  //MPI_Barrier(cf.comm);
  //Fiesta::Log::debugAll("{}: lStart={} lEnd={} lExt={}",name,gStart,gEnd,stride);

}

template <typename T>
blockWriter<T>::blockWriter(struct inputConfig& cf, rk_func* f, string name_, string path_, bool avg_, size_t frq_,
  vector<size_t> start_, vector<size_t> end_, vector<size_t> stride_):
  name(name_),path(path_),avg(avg_),freq(frq_),gStart(start_),gEnd(end_),stride(stride_){

  int strideDelta;
  size_t gMin;

  if (cf.rank==0){
    if (!std::filesystem::exists(path)){
      Fiesta::Log::message("Creating directory: '{}'",path);
      std::filesystem::create_directories(path);
    }
  }

  pad = (int)log10(cf.tend) + 1;

  // initialize parameters
  lElems=1;
  lElemsG=1;
  myColor=MPI_UNDEFINED;
  slicePresent=true;
  writeVarx=true;

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
    }

      //MPI_Barrier(sliceComm);
      //Fiesta::Log::debugAll("PRE {}: gStart={} lStart={} offset={}",name,gStart,lStart,lOffset);
      //MPI_Barrier(sliceComm);

    for (int i=0; i<cf.ndim; ++i){
      //adjust start and offset if it does not land on a stride if the slice does not start in this domain
      offsetDelta = lOffset[i]-gStart[i];
        strideDelta = (lOffset[i])%stride[i];
      if( strideDelta > 0 ){
        lStart[i]+=stride[i]-strideDelta;
        lOffset[i]+=stride[i]-strideDelta;
      }
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

      // compute local data and grid extent based on computed start and end indexes
      lExt[i]=(lEnd[i]-lStart[i])/stride[i]+1;
      lExtG[i]=(lEndG[i]-lStart[i])/stride[i]+1;

      // compute length of 1d buffer
      lElems*=lExt[i];
      lElemsG*=lExtG[i];
    }

    // Processes must agree on global dimensions
    // last "righ-most" process may have had to reduce
    // global dimensions to accomadate the stride.
    // Processes will use the minimum extent and end index.
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

  if (slicePresent){
    //Fiesta::Log::debug("{}: gStart={} gEnd={} gExt={} gExtG={} stride={}",name,gStart,gEnd,gExt,gExtG,stride);
    //MPI_Barrier(sliceComm);
    //if (myColor==1)
      //Fiesta::Log::debugAll("{}: lStart={} lEnd={} lExt={} offset={}",name,lStart,lEnd,lExt,lOffset);
    //MPI_Barrier(sliceComm);
  }
}

template<typename T>
void blockWriter<T>::write(struct inputConfig cf, rk_func *f, int tdx, double time) {
  string baseFormat = format("{{}}-{{:0{}d}}",pad);
  string blockBase = format(baseFormat,name,tdx);
  string hdfName   = format("{}.h5",blockBase);
  string hdfPath   = format("{}/{}",path,hdfName);
  string xmfPath   = format("{}/{}.xmf",path,blockBase);
  
  Fiesta::Log::message("[{}] Writing '{}'",cf.t,hdfPath);
  fiestaTimer writeTimer = fiestaTimer();

  if (myColor==1){

    h5Writer<T> writer;
    writer.open(sliceComm, MPI_INFO_NULL, hdfPath);

    writer.openGroup("/Grid");
    Kokkos::deep_copy(gridH,f->grid);
    for (int vn=0; vn<cf.ndim; ++vn){
      dataPack(cf.ndim, 0, lStart, lEndG, lExtG, stride, gridData, gridH,vn,false);
      writer.write(format("Dimension{}",vn), cf.ndim, gExtG, lExtG, lOffset, gridData); 
    }
    writer.closeGroup();

    writer.openGroup("/Solution");
    Kokkos::deep_copy(varH,f->var);
    for (int vn=0; vn<cf.nvt; ++vn){
      dataPack(cf.ndim, cf.ng, lStart, lEnd, lExt, stride, varData, varH,vn,avg);
      writer.write(format("Variable{:02d}",vn), cf.ndim, gExt, lExt, lOffset, varData); 
    }
    if (writeVarx){
      Kokkos::deep_copy(varxH,f->varx);
      for (size_t vn = 0; vn < f->varxNames.size(); ++vn) {
        dataPack(cf.ndim, cf.ng, lStart, lEnd, lExt, stride, varData, varxH,vn,avg);
        writer.write(format("Variable{:02d}",vn+cf.nvt), cf.ndim, gExt, lExt, lOffset, varData); 
      }
    }
    writer.closeGroup();

    writer.close();
  }

  MPI_Barrier(cf.comm);

  if (cf.rank==0){
    filesystem::path hdfFilePath{hdfPath};
    auto fsize = filesystem::file_size(hdfFilePath);
    double ftime = writeTimer.check();
    double frate = (fsize/1048576.0)/ftime;

    if (fsize > 1073741824)
      Fiesta::Log::message("[{}] '{}': {:.2f} GiB in {:.2f}s ({:.2f} MiB/s)",cf.t,hdfPath,fsize/1073741824.0,ftime,frate);  //divide by bytes per GiB (1024*1024*1024)
    else
      Fiesta::Log::message("[{}] '{}': {:.2f} MiB in {:.2f}s ({:.2f} MiB/s)",cf.t,hdfPath,fsize/1048576.0,ftime,frate);     //divide by bytes per MiB (1024*1024)
  }

  Fiesta::Log::message("[{}] Writing '{}'",cf.t,xmfPath);
  if (myColor==1){
    writeXMF(xmfPath, hdfName, time, cf.ndim, gExt.data(), cf.nvt, writeVarx,f->varNames,f->varxNames);
  }
}

template<typename T>
void blockWriter<T>::dataPack(int ndim, int ng, const vector<size_t>& start, const vector<size_t>& end, const vector<size_t>& extent, const vector<size_t>& stride, vector<T>& dest, const FS4DH& source, const int vn, const bool average){
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

template class blockWriter<float>;
template class blockWriter<double>;
