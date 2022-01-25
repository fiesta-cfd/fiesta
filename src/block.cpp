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
#include <numeric>
#include "xdmf.hpp"
#include "timer.hpp"
#include "fmt/core.h"
#include <filesystem>
#ifdef HAVE_MPI
#include "mpi.h"
#endif
#include "h5.hpp"
#include "log2.hpp"
#include "fiesta.hpp"

using namespace std;
using fmt::format;

template <typename T>
size_t blockWriter<T>::frq() { return freq; }

template <typename T>
blockWriter<T>::blockWriter(){}

template <typename T>
blockWriter<T>::blockWriter(struct inputConfig& cf, std::unique_ptr<class rk_func>& f, string name_, string path_, bool avg_, size_t frq_, bool appStep_):
  name(name_), path(path_), avg(avg_), freq(frq_), appStep(appStep_),
  gStart(std::vector<size_t>(cf.ndim,0)),
  gExt(cf.globalCellDims),
  gEnd(cf.globalCellDims),
  gExtG(cf.globalCellDims),
  lStart(std::vector<size_t>(cf.ndim,0)),
  lExt(cf.localCellDims),
  lEndG(cf.localCellDims),
  lEnd(cf.localCellDims),
  lExtG(cf.localCellDims),
  lOffset(cf.subdomainOffset),
  stride(std::vector<size_t>(cf.ndim,1)),
  gOrigin(cf.dxvec),
  iodx(cf.dxvec)
  {

  // initialize parameters
  myColor=1;
  slicePresent=true;
  writeVarx=false;
  pad = (int)log10(cf.tend) + 1;
  chunkable = cf.chunkable;
#ifdef HAVE_MPI
  sliceComm = cf.comm;
  sliceRank = cf.rank;
  sliceSize = cf.numProcs;
#endif

  if (cf.rank==0){
    if (!std::filesystem::exists(path)){
      Log::message("Creating directory: '{}'",path);
      std::filesystem::create_directories(path);
    }
  }

  for (int i=0; i<cf.ndim; ++i){
    gEnd[i] -= 1;
    gExtG[i] += 1;
    lEnd[i] -= 1;
    lExtG[i] += 1;
    gOrigin[i] *= gStart[i];
  }

  lElems  = std::accumulate(lExt.begin(), lExt.end(), 1, std::multiplies<int>());
  lElemsG = std::accumulate(lExtG.begin(), lExtG.end(), 1, std::multiplies<int>());

  // allocate pack buffers
  fill_n(back_inserter(varData),lElems,0.0);
  fill_n(back_inserter(gridData),lElemsG,0.0);

  // create kokkos host mirrors
  gridH = Kokkos::create_mirror_view(f->grid);
  varH = Kokkos::create_mirror_view(f->var);
  varxH = Kokkos::create_mirror_view(f->varx);

  // if(sliceRank==0)
  //   Log::debugAll("IO VIEW: name={}, rank={}, size={}, chunkable={} compressible={}",name,sliceRank,sliceSize,cf.chunkable,cf.compressible);

}

template <typename T>
blockWriter<T>::blockWriter(struct inputConfig& cf, std::unique_ptr<class rk_func>& f, string name_, string path_, bool avg_, size_t frq_,
  vector<size_t> start_, vector<size_t> end_, vector<size_t> stride_, bool appStep_):
  name(name_),path(path_),avg(avg_),freq(frq_),
  gStart(start_),
  gEnd(end_),
  stride(stride_),
  appStep(appStep_){

#ifdef HAVE_MPI
  sliceRank=MPI_PROC_NULL;
  sliceSize=0;
#endif

  for (int i=0; i<cf.ndim; ++i){
    if (gStart[i] > cf.globalCellDims[i]){
      Log::error("Starting index of ioview '{}' exceeds domain bounds.",name);
      exit(EXIT_FAILURE);
    }

    if (gEnd[i] > cf.globalCellDims[i]){
      Log::error("Ending index of ioview '{}' exceeds domain bounds.",name);
      exit(EXIT_FAILURE);
    }

    if (gStart[i] > gEnd[i]){
      Log::error("Sarting index of ioview '{}' exceeds ending index.",name);
      exit(EXIT_FAILURE);
    }

    if (stride[i] > 4){
      Log::error("Stride of ioview '{}' exceeds 4.",name);
      exit(EXIT_FAILURE);
    }

    if (gEnd[i] > gStart[i] && gEnd[i]-gStart[i] < stride[i]){
      Log::error("Stride of ioview '{}' exceeds view dimensions.",name);
      exit(EXIT_FAILURE);
    }

  }

  int strideDelta;
  size_t gMin;

  if (cf.rank==0){
    if (!std::filesystem::exists(path)){
      Log::message("Creating directory: '{}'",path);
      std::filesystem::create_directories(path);
    }
  }
  pad = (int)log10(cf.tend) + 1;

  // initialize parameters
  lElems=1;
  lElemsG=1;
  myColor=1;
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
#ifdef HAVE_MPI
  myColor=MPI_UNDEFINED;
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
#endif

  if (slicePresent){
#ifdef HAVE_MPI
    MPI_Comm_rank(sliceComm, &sliceRank);
    MPI_Comm_size(sliceComm, &sliceSize);
#endif
    for (int i=0; i<cf.ndim; ++i){
      offsetDelta=gStart[i]-cf.subdomainOffset[i]; //location of left edge of slice wrt left edge of subdomain
      if(offsetDelta > 0) lStart[i]=offsetDelta;   //if slice starts inside subdomain, set local start to slice edge
      else lOffset[i] = -(offsetDelta);            //if subdomain is within slice, adjust offset

      offsetDelta=(cf.subdomainOffset[i]+cf.localCellDims[i])-gEnd[i]; //location of right edge of slice wrt right edge of domain
      if(offsetDelta > 0){
        lEnd[i]=cf.localCellDims[i]-offsetDelta;   // if slice ends in subdomain, set local end to slide edge
      }
    }

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
#ifdef HAVE_MPI
    for (int i=0; i<cf.ndim; ++i){
      MPI_Allreduce(&gEnd[i],&gMin,1,MPI_INT,MPI_MIN,sliceComm);
      gEnd[i]=gMin;

      MPI_Allreduce(&gExt[i],&gMin,1,MPI_INT,MPI_MIN,sliceComm);
      gExt[i]=gMin;

      MPI_Allreduce(&gExtG[i],&gMin,1,MPI_INT,MPI_MIN,sliceComm);
      gExtG[i]=gMin;

      gOrigin.push_back(gStart[i]*cf.dxvec[i]);
      iodx.push_back(cf.dxvec[i]*stride[i]);
    }
#endif
  }
 
  // allocate pack buffers
  fill_n(back_inserter(varData),lElems,0.0);
  fill_n(back_inserter(gridData),lElemsG,0.0);

  // create kokkos host mirrors
  gridH = Kokkos::create_mirror_view(f->grid);
  varH = Kokkos::create_mirror_view(f->var);
  varxH = Kokkos::create_mirror_view(f->varx);
  /* if(sliceRank==0) */
  /*   Log::debugAll("IO VIEW: name={}, rank={}, size={}, start={}, end={}, stride={}, extent={}",name,sliceRank,sliceSize,gStart,gEnd,stride,gExt); */

  chunkable=false;
  Log::warning("Chunking it not supported for IO views with strides.  Chunking will be disabled for {}.",name);
}

template<typename T>
void blockWriter<T>::write(struct inputConfig cf, std::unique_ptr<class rk_func>&f, int tdx, FSCAL time) {
  string baseFormat,blockBase;
  if(appStep){
    baseFormat = format("{{}}-{{:0{}d}}",pad);
    blockBase = format(baseFormat,name,tdx);
  }else{
    baseFormat = format("{{}}");
    blockBase = format(baseFormat,name);
  }
  string hdfName   = format("{}.h5",blockBase);
  string hdfPath   = format("{}/{}",path,hdfName);
  string xmfPath   = format("{}/{}.xmf",path,blockBase);
 
  Log::message("[{}] Writing '{}'",cf.t,hdfPath);
  Timer::fiestaTimer writeTimer = Timer::fiestaTimer();

  //#ifdef HAVE_MPI
  //  MPI_Barrier(cf.comm);
  //#endif

  if (myColor==1){
    h5Writer<T> writer;
#ifdef HAVE_MPI
    writer.open(sliceComm, MPI_INFO_NULL, hdfPath);
#else
    writer.open(hdfPath);
#endif

    if (cf.grid > 0){
      writer.openGroup("/Grid");
      Kokkos::deep_copy(gridH,f->grid);
      for (int vn=0; vn<cf.ndim; ++vn){
        dataPack(cf.ndim, 0, lStart, lEndG, lExtG, stride, gridData, gridH,vn,false);
        writer.write(format("Dimension{}",vn), cf.ndim, gExtG, lExtG, lOffset, gridData, chunkable, cf.compressible); 
      }
      writer.closeGroup();
    }

    writer.openGroup("/Solution");
    Kokkos::deep_copy(varH,f->var);
    for (int vn=0; vn<cf.nvt; ++vn){
      dataPack(cf.ndim, cf.ng, lStart, lEnd, lExt, stride, varData, varH,vn,avg);
      writer.write(f->varNames[vn], cf.ndim, gExt, lExt, lOffset, varData, chunkable, cf.compressible); 
    }
    if (writeVarx){
      Kokkos::deep_copy(varxH,f->varx);
      for (size_t vn = 0; vn < f->varxNames.size(); ++vn) {
        dataPack(cf.ndim, cf.ng, lStart, lEnd, lExt, stride, varData, varxH,vn,avg);
        writer.write(f->varxNames[vn], cf.ndim, gExt, lExt, lOffset, varData, chunkable, cf.compressible); 
      }
    }
    writer.closeGroup();

    writer.openGroup("/Properties");
    writer.writeAttribute("time_index",tdx);
    writer.writeAttribute("time",time);
    writer.writeAttribute("title",cf.title);
    writer.writeAttribute("fiesta_version",std::string(FIESTA_VERSION));
    writer.writeAttribute("fiesta_options",std::string(FIESTA_OPTIONS));
    writer.writeAttribute("fiesta_btime",std::string(FIESTA_BTIME));
    writer.writeAttribute("restart_file_version",FIESTA_RESTART_VERSION);
    writer.closeGroup();

    writer.close();
    //#ifdef HAVE_MPI
    //    MPI_Barrier(sliceComm);
    //#endif
  }

#ifdef HAVE_MPI
  MPI_Barrier(cf.comm);
#endif

  if (cf.rank==0){
    filesystem::path hdfFilePath{hdfPath};
    auto fsize = filesystem::file_size(hdfFilePath);
    FSCAL ftime = writeTimer.check();
    FSCAL frate = (fsize/1048576.0)/ftime;

    if (fsize > 1073741824)
      Log::message("[{}] '{}': {:.2f} GiB in {:.2f}s ({:.2f} MiB/s)",cf.t,hdfPath,fsize/1073741824.0,ftime,frate);  //divide by bytes per GiB (1024*1024*1024)
    else
      Log::message("[{}] '{}': {:.2f} MiB in {:.2f}s ({:.2f} MiB/s)",cf.t,hdfPath,fsize/1048576.0,ftime,frate);     //divide by bytes per MiB (1024*1024)
  }

  //Log::message("[{}] Writing '{}'",cf.t,xmfPath);
  if (myColor==1){
    writeXMF(xmfPath, hdfName, cf.grid, time, cf.ndim, gExt.data(),gOrigin,iodx, cf.nvt, writeVarx,f->varNames,f->varxNames);
  }
}

template<typename T>
void blockWriter<T>::dataPack(int ndim, int ng, const vector<size_t>& start, const vector<size_t>& end, const vector<size_t>& extent, const vector<size_t>& stride, vector<T>& dest, const FS4DH& source, const int vn, const bool average){
  int ii,jj,kk;
  int idx;
  FSCAL mean;

  if (ndim==3){
    for (size_t i = start[0]; i < end[0]+1; i+=stride[0]) {
      for (size_t j = start[1]; j < end[1]+1; j+=stride[1]) {
        for (size_t k = start[2]; k < end[2]+1; k+=stride[2]) {
          ii=(i-start[0])/stride[0];
          jj=(j-start[1])/stride[1];
          kk=(k-start[2])/stride[2];
          idx = (extent[0] * extent[1]) * kk + extent[0] * jj + ii;
          if (average){
            mean=0.0;
            for (size_t ia=0; ia<stride[0]; ++ia){
              for (size_t ja=0; ja<stride[1]; ++ja){
                for (size_t ka=0; ka<stride[2]; ++ka){
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
    for (size_t i = start[0]; i < end[0]+1; i+=stride[0]) {
      for (size_t j = start[1]; j < end[1]+1; j+=stride[1]) {
        ii=(i-start[0])/stride[0];
        jj=(j-start[1])/stride[1];
        idx = extent[0] * jj + ii;
        if (average){
          mean=0.0;
          for (size_t ia=0; ia<stride[0]; ++ia){
            for (size_t ja=0; ja<stride[1]; ++ja){
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
