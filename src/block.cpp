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
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include "xdmf.hpp"
//#ifndef NOMPI
#include "mpi.h"
//#endif

using namespace std;

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
//template <typename H> blockWriter::getH5Type();
template <> hid_t blockWriter::getH5Type<float>(){ return H5T_NATIVE_FLOAT; }
template <> hid_t blockWriter::getH5Type<double>(){ return H5T_NATIVE_DOUBLE; }
template <> hid_t blockWriter::getH5Type<int>(){ return H5T_NATIVE_INT; }


//blockWriter::blockWriter(struct inputConfig& cf_, rk_func *f_) : cf(cf_),f(f_) {
blockWriter::blockWriter(){}
blockWriter::blockWriter(struct inputConfig& cf, rk_func *f){
  lStart  = new size_t[cf.ndim];
  lEnd    = new size_t[cf.ndim];
  lEndG   = new size_t[cf.ndim];
  lExt    = new size_t[cf.ndim];
  lExtG   = new size_t[cf.ndim];
  lOffset = new size_t[cf.ndim];
  gStart  = new size_t[cf.ndim];
  gEnd    = new size_t[cf.ndim];
  gExt    = new size_t[cf.ndim];
  gExtG   = new size_t[cf.ndim];
  stride  = new size_t[cf.ndim];
  lElems=1;
  lElemsG=1;

  myColor=MPI_UNDEFINED;
  slicePresent=true;

  int strideDelta;

  luaReader L(cf.inputFname);
  L.getIOBlock(cf.ndim,name,path,freq,gStart,gEnd);
  L.close();
  cf.log->debug("ioblock: ",name,", ",path,",",freq,", (",gStart[0],",",gStart[1],",",gStart[2],")","(",gEnd[0],",",gEnd[1],",",gEnd[2],")");
  stride[0]=3;
  stride[1]=3;
  stride[2]=3;

  for (int i=0; i<cf.ndim; ++i){
    //adjust block global size if it does not line up with strides
    strideDelta=(gEnd[i]-gStart[i])%stride[i];
    gEnd[i]-=strideDelta;

    gExt[i]=(gEnd[i]-gStart[i])/stride[i]+1;
    gExtG[i]=gExt[i]+1;
    //if (i==2) cout << cf.rank << " A : " << strideDelta << ", " << gStart[i] << "," << gEnd[i] << ", " << gExt[i] << ", " << gExtG[i] << endl;
  }

  for (int i=0; i<cf.ndim; ++i){
    lStart[i]=0;
    lEnd[i]=cf.localCellDims[i]-1;
    lEndG[i]=cf.localCellDims[i]-1; //none
    lOffset[i]=0;
  }

  for (int i=0; i<cf.ndim; ++i){
    lExt[i]=lEnd[i]-lStart[i]+1;
    lExtG[i]=lEnd[i]-lStart[i]+1; //+2
  }

  for (int i=0; i<cf.ndim; ++i){
    slicePresent = slicePresent && ((cf.subdomainOffset[i]+cf.localCellDims[i]) >= gStart[i] && cf.subdomainOffset[i] <= gEnd[i]);
  }

  if (slicePresent){
    myColor=1;
  }


  // create slice communicator
  MPI_Comm_split(cf.comm,myColor,cf.rank,&sliceComm);
  if (slicePresent){
    int sliceRank;
    int sliceWorld;
    MPI_Comm_rank(sliceComm,&sliceRank);
    MPI_Comm_size(sliceComm,&sliceWorld);
    //size_t gSize[cf.numProcs];
    size_t gMin;


    for (int i=0; i<cf.ndim; ++i){

      //lEndG[i]=lEnd[i];//+1

      offsetDelta=gStart[i]-cf.subdomainOffset[i]; //location of left edge of slice wrt left edge of subdomain
      if(offsetDelta > 0) lStart[i]=offsetDelta;   //if slice starts inside subdomain, set local start to slice edge
      else lOffset[i] = -(offsetDelta);            //if subdomain is within slice, adjust offset

      offsetDelta=(cf.subdomainOffset[i]+cf.localCellDims[i])-gEnd[i]; //location of right edge of slice wrt right edge of domain
      if(offsetDelta > 0){
        lEnd[i]=cf.localCellDims[i]-offsetDelta;   // if slice ends in subdomain, set local end to slide edge
        //lEndG[i]=lEnd[i]+1;
      }


      //if (i==2) cout << cf.rank << " B : " << lOffset[i] << ", " << lStart[i] << endl;
      //if (i==2) cout << cf.rank << " B.2 : " << lStart[i] << ", " << gStart[i] << "," << stride[i] << endl;

      //adjust start and offset if it does not land on a stride
      strideDelta = (lOffset[i]-gStart[i])%stride[i];
      lStart[i]+=(stride[i]-strideDelta)%stride[i];
      lOffset[i]+=(stride[i]-strideDelta)%stride[i];
      lOffset[i]/=stride[i];

      //if (i==2) cout << cf.rank << " C : " << lOffset[i] << ", " << lStart[i] << ", " << strideDelta << endl;
      //if (i==2) cout << cf.rank << " A : " << strideDelta << ", " << gEnd[i] << ", " << gExt[i] << ", " << gExtG[i] << endl;

      //adjust local end index if it does not align with a stride
      //if (i==2) cout << cf.rank << " D : " << lEnd[i] << ", " << lEndG[i] << endl;
      strideDelta=(lEnd[i]-lStart[i])%stride[i];
      lEnd[i]-=strideDelta;
      //lEndG[i]=lEnd[i];
      //if (i==2) cout << cf.rank << " E : " << lEnd[i] << ", " << lEndG[i] << endl;

      //cf.log->debugAll(sliceRank,"A ",i,": ",cf.localCellDims[i],", ",lEnd[i],", ",stride[i]);
      if ((cf.subdomainOffset[i]+cf.localCellDims[i]) > gEnd[i]){
        if (((cf.localCellDims[i]-1)-lEnd[i]) < (stride[i]-1)){
          lEnd[i]-=stride[i];
          gEnd[i]-=stride[i];
          gExt[i]-=1;
          gExtG[i]-=1;
          //cout << sliceRank << ", " << i << ": " << gEnd[i] << ", " << gExt[i] << "\n";
          //cf.log->debugAll(sliceRank," ",i,": ",gEnd[i],", ",gExt[i]);
        }
        lEndG[i]=lEnd[i]+stride[i];
      }
      //if (i==2) cout << cf.rank << " F : " << gEnd[i] << ", " << gExt[i] << ", " << gExtG[i] << ", " << lEnd[i] << ", " << lEndG[i] << "\n";
      //processes must come to agreement about gEnd, gExt and gExtG
      //MPI_Allgather(&gEnd[i],1,MPI_INT,gSize,1,MPI_INT,sliceComm);
      
      //lEndG[i]=lEnd[i]+stride[i];
      lExt[i]=(lEnd[i]-lStart[i])/stride[i]+1;
      lExtG[i]=(lEndG[i]-lStart[i])/stride[i]+1;
      //if (i==2) cout << cf.rank << " G : " << lExt[i] << ", " << lExtG[i] << endl;

      lElems*=lExt[i];
      lElemsG*=lExtG[i];


    }
    //cout << cf.rank << ": (" << lExt[0]   << "," << lExt[1]  << "," << lExt[2]  << "), "
    //                <<   "(" << lExtG[0]  << "," << lExtG[1] << "," << lExtG[2] << ")" << endl;
    //cout << cf.rank << ": (" << lElems << ", " << lElemsG << ")\n";
    for (int i=0; i<cf.ndim; ++i){
      MPI_Allreduce(&gEnd[i],&gMin,1,MPI_INT,MPI_MIN,sliceComm);
      gEnd[i]=gMin;

      MPI_Allreduce(&gExt[i],&gMin,1,MPI_INT,MPI_MIN,sliceComm);
      gExt[i]=gMin;

      MPI_Allreduce(&gExtG[i],&gMin,1,MPI_INT,MPI_MIN,sliceComm);
      gExtG[i]=gMin;

      if(slicePresent){
        //if (i==2) cout << cf.rank << " H : (" << i << ") " << lOffset[i] << ", " << lStart[i] << ", " << lEnd[i] << ", "  << lExt[i] << ", " << lEndG[i] << ", " << lExtG[i] << "\n";
        //if (i==2) cout << cf.rank << " H : (" << i << ") " << gStart[i] << ", " << gEnd[i] << ", "  << gExt[i] << ", " << gExtG[i] << "\n";
      }
      //MPI_Barrier(sliceComm);
    }
  }
  

  varData = (float*)malloc(lElems*sizeof(float*));
  gridData = (float*)malloc(lElemsG*sizeof(float*));

  gridH = Kokkos::create_mirror_view(f->grid);
  varH = Kokkos::create_mirror_view(f->var);
  varxH = Kokkos::create_mirror_view(f->varx);
}

//template<typename T>
//void blockWriter::writeHDF(struct inputConfig cf, rk_func *f, int tdx, double time, T* x, T* var, string name) {
void blockWriter::write(struct inputConfig cf, rk_func *f, int tdx, double time) {


  // calcualte string width for time index
  pad = (int)log10(cf.nt) + 1;
  stringstream lineName,lineBase,lineHDF,lineXMF;
  lineBase << name << "-" << setw(pad) << setfill('0') << tdx;
  lineName << path << "/" << lineBase.str() << ".h5";
  lineXMF << path << "/" << lineBase.str() << ".xmf";
  lineHDF << lineBase.str() << ".h5";
  cf.log->message("[",cf.t,"] ","Writing ",lineName.str());

  if (myColor==1){
    Kokkos::deep_copy(gridH,f->grid);
    Kokkos::deep_copy(varH,f->var);

    //int sliceRank;
    //int sliceWorld;
    //MPI_Comm_rank(sliceComm,&sliceRank);
    //MPI_Comm_size(sliceComm,&sliceWorld);

    hid_t linefile_id = openHDF5ForWrite(sliceComm, MPI_INFO_NULL, lineName.str());

    hid_t linegroup_id = H5Gcreate(linefile_id, "/Grid", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    for (int vn=0; vn<cf.ndim; ++vn){
      dataPack(cf.ndim, 0, lStart, lEndG, lExtG, stride, gridData, gridH,vn);

      stringstream lineStr;
      lineStr << "Dimension" << vn;
      write_h5<float>(linegroup_id, lineStr.str(), cf.ndim, gExtG, lExtG, lOffset, gridData); 
    }
    H5Gclose(linegroup_id);

    linegroup_id = H5Gcreate(linefile_id, "/Solution", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    for (int vn=0; vn<cf.nvt; ++vn){
      dataPack(cf.ndim, cf.ng, lStart, lEnd, lExt, stride, varData, varH,vn);
       
      stringstream lineStr;
      lineStr << "Variable" << setw(2) << setfill('0') << vn;
      write_h5<float>(linegroup_id, lineStr.str(), cf.ndim, gExt, lExt, lOffset, varData); 
    }
    for (size_t vn = 0; vn < f->varxNames.size(); ++vn) {
      dataPack(cf.ndim, cf.ng, lStart, lEnd, lExt, stride, varData, varxH,vn);
       
      stringstream lineStr;
      lineStr << "Variable" << setw(2) << setfill('0') << vn+cf.nvt;
      write_h5<float>(linegroup_id, lineStr.str(), cf.ndim, gExt, lExt, lOffset, varData); 
    }
    H5Gclose(linegroup_id);

    close_h5(linefile_id);

  }

  cf.log->message("[",cf.t,"] ","Writing ",lineXMF.str());
  if (myColor==1){
    writeXMF(lineXMF.str(), lineHDF.str(), time, cf.ndim, gExt, cf.nvt, f->varNames,f->varxNames);
  }
}

// writer function for hdf5
template <typename S>
void blockWriter::write_h5(hid_t group_id, string dname, int ndim,
                   size_t* in_dims_global, size_t* in_dims_local, size_t* in_offset, S* data){

  hsize_t dims_global[ndim];
  hsize_t dims_local[ndim];
  hsize_t offset[ndim];

  for (int i=0; i<ndim; ++i){
    dims_global[i]=in_dims_global[ndim-1-i];
    dims_local[i] =in_dims_local[ndim-1-i];
    offset[i]     =in_offset[ndim-1-i];
  }

  //invertArray(ndim,dims_global,in_dims_global);
  //invertArray(ndim,dims_local,in_dims_local);
  //invertArray(ndim,offset,in_offset);

  // identifiers
  hid_t filespace, memspace, dset_id, plist_id, dtype_id;
  //herr_t status;

  // get type id
  dtype_id = getH5Type<S>();

  // create global filespace and local memoryspace
  filespace = H5Screate_simple(ndim, dims_global, NULL);
  memspace  = H5Screate_simple(ndim, dims_local, NULL);

  // create dataset
  dset_id = H5Dcreate(group_id, dname.c_str(), dtype_id, filespace,
                                        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  H5Sclose(filespace);

  // select hyperslab in global filespace
  filespace = H5Dget_space(dset_id);
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, dims_local, NULL);
  //status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, dims_local, NULL);

  // create property list for collective dataset write
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

  // write data
  H5Dwrite(dset_id, dtype_id, memspace, filespace, plist_id, data);
  //status = H5Dwrite(dset_id, dtype_id, memspace, filespace, plist_id, data);

  // close sets, spaces and lists
  H5Dclose(dset_id);
  H5Sclose(filespace);
  H5Sclose(memspace);
  H5Pclose(plist_id);
}

template<typename T>
void blockWriter::dataPack(int ndim, int ng, size_t* start, size_t* end, size_t* extent, size_t* stride, T* dest, FS4DH& source,int vn){
  int ii,jj,kk;
  int idx;


  if (ndim==3){
    for (int i = start[0]; i < end[0]+1; i+=stride[0]) {
      for (int j = start[1]; j < end[1]+1; j+=stride[1]) {
        for (int k = start[2]; k < end[2]+1; k+=stride[2]) {
          ii=(i-start[0])/stride[0];
          jj=(j-start[1])/stride[1];
          kk=(k-start[2])/stride[2];
          idx = (extent[0] * extent[1]) * kk + extent[0] * jj + ii;
          //cout << "######## (" << i << ", " << j << ", " << k << ") " << idx << endl;
          dest[idx] = source(i+ng, j+ng, k+ng,vn);
        }
      }
    }
  }else{
    for (int i = start[0]; i < end[0]+1; i+=stride[0]) {
      for (int j = start[1]; j < end[1]+1; j+=stride[1]) {
        ii=(i-start[0])/stride[0];
        jj=(j-start[1])/stride[1];
        idx = extent[0] * jj + ii;
        dest[idx] = source(i+ng, j+ng, 0,vn);
      }
    }
  }
}

blockWriter::~blockWriter(){
  free(varData);
  free(gridData);
}
