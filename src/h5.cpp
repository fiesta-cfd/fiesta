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

#include <string>
#include <vector>
#include "hdf5.h"
#ifdef HAVE_MPI
#include "mpi.h"
#endif
#include "h5.hpp"
#include <algorithm>
#include "log2.hpp"
#include "debug.hpp"
#include "kokkosTypes.hpp"

// open and hdf5 file for writing
template <typename T>
h5Writer<T>::h5Writer(){}

#ifdef HAVE_MPI
template <typename T>
void h5Writer<T>::open(MPI_Comm comm, MPI_Info info, std::string fname){
  hid_t pid;

  pid = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(pid, comm, info);
  file_id = H5Fcreate(fname.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, pid);
  H5Pclose(pid);
  MPI_Barrier(comm);
}
#else
template <typename T>
void h5Writer<T>::open(std::string fname){
  hid_t pid;

  pid = H5Pcreate(H5P_FILE_ACCESS);
  file_id = H5Fcreate(fname.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, pid);
  H5Pclose(pid);
}
#endif

template <typename T>
void h5Writer<T>::openRead(std::string fname){
  file_id = H5Fopen(fname.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
}

// close an hdf5 file
template <typename T>
void h5Writer<T>::close(){
  H5Fclose(file_id);
}

// writer function for hdf5
template <typename T>
void h5Writer<T>::write(std::string dname, int ndim,
  std::vector<size_t> in_dims_global, std::vector<size_t> in_dims_local, std::vector<size_t> in_offset,
  std::vector<T>& data,
  bool chunkable,
  bool compress){
  

  std::vector<hsize_t> gdims(ndim,0);
  std::vector<hsize_t> ldims(ndim,0);
  std::vector<hsize_t> offset(ndim,0);

  // reverse order of array indexes to "c"
  std::reverse_copy(in_dims_global.begin(),in_dims_global.end(),gdims.begin());
  std::reverse_copy(in_dims_local.begin(),in_dims_local.end(),ldims.begin());
  std::reverse_copy(in_offset.begin(),in_offset.end(),offset.begin());

  // identifiers
  hid_t filespace, memspace, dset_id, plist_id, dtype_id, dcpl_id;
  unsigned szip_options_mask, szip_pixels_per_block;

  // get type id
  if (std::is_same<T,FSCAL>::value) dtype_id = H5T_NATIVE_DOUBLE;
  if (std::is_same<T,float>::value) dtype_id = H5T_NATIVE_FLOAT;
  if (std::is_same<T,int>::value) dtype_id = H5T_NATIVE_INT;

  if (chunkable){
    dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(dcpl_id,ndim,ldims.data());
    if(compress){
      szip_options_mask=H5_SZIP_NN_OPTION_MASK;
      szip_pixels_per_block=32;
      H5Pset_szip (dcpl_id, szip_options_mask, szip_pixels_per_block);
    }
  }else{
    dcpl_id = H5P_DEFAULT;
  }

  // create global filespace and local memoryspace
  filespace = H5Screate_simple(ndim, gdims.data(), NULL);
  memspace  = H5Screate_simple(ndim, ldims.data(), NULL);

  // create dataset
  dset_id = H5Dcreate(group_id, dname.c_str(), dtype_id, filespace,
                                        H5P_DEFAULT, dcpl_id, H5P_DEFAULT);

  H5Sclose(filespace);

  // select hyperslab in global filespace
  filespace = H5Dget_space(dset_id);
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset.data(), NULL, ldims.data(), NULL);
  //status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, dims_local, NULL);

  // create property list for collective dataset write
  plist_id = H5Pcreate(H5P_DATASET_XFER);
#ifdef HAVE_MPI
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
#endif

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
template<typename S>
void h5Writer<T>::writeAttribute(std::string name, S data){
  hid_t dspace_id, att_id, dtype_id;
  //bool stringAttribute;

  if (std::is_same<S,FSCAL>::value) dtype_id = H5T_NATIVE_DOUBLE;
  if (std::is_same<S,float>::value) dtype_id = H5T_NATIVE_FLOAT;
  if (std::is_same<S,int>::value) dtype_id = H5T_NATIVE_INT;
  //if (std::is_same<S,std::string>::value) stringAttribute=true;

  //if (stringAttribute){
  //  dspace_id = H5Screate(H5S_SCALAR);
  //  dtype_id = H5Tcopy(H5T_C_S1);
  //  H5Tset_size(dtype_id, data.length());
  //  H5Tset_strpad(dtype_id,H5T_STR_NULLTERM);
  //  att_id = H5Acreate2(group_id, name.c_str(), dtype_id, dspace_id, H5P_DEFAULT, H5P_DEFAULT);
  //  H5Awrite(att_id,dtype_id,data.c_str());
  //  H5Aclose(att_id);
  //  H5Sclose(dspace_id);
  //}else{
    dspace_id = H5Screate(H5S_SCALAR);
    att_id = H5Acreate2(group_id,name.c_str(),dtype_id,dspace_id, H5P_DEFAULT,H5P_DEFAULT);
    H5Awrite(att_id,dtype_id,&data);

    H5Aclose(att_id);
    H5Sclose(dspace_id);
  //}

}

template<>
template<>
void h5Writer<double>::writeAttribute(std::string name, std::string data){
  hid_t dspace_id, att_id, dtype_id;

  dspace_id = H5Screate(H5S_SCALAR);
  dtype_id = H5Tcopy(H5T_C_S1);
  H5Tset_size(dtype_id, data.length());
  H5Tset_strpad(dtype_id,H5T_STR_NULLTERM);
  att_id = H5Acreate2(group_id, name.c_str(), dtype_id, dspace_id, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(att_id,dtype_id,data.c_str());
  H5Aclose(att_id);
  H5Sclose(dspace_id);
}
template<>
template<>
void h5Writer<float>::writeAttribute(std::string name, std::string data){
  hid_t dspace_id, att_id, dtype_id;

  dspace_id = H5Screate(H5S_SCALAR);
  dtype_id = H5Tcopy(H5T_C_S1);
  H5Tset_size(dtype_id, data.length());
  H5Tset_strpad(dtype_id,H5T_STR_NULLTERM);
  att_id = H5Acreate2(group_id, name.c_str(), dtype_id, dspace_id, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(att_id,dtype_id,data.c_str());
  H5Aclose(att_id);
  H5Sclose(dspace_id);
}

template <typename T>
void h5Writer<T>::read(std::string path, int ndim, std::vector<size_t> in_dims_global, std::vector<size_t> in_dims_local, std::vector<size_t> in_offset, T* data){
  std::vector<hsize_t> gdims(ndim,0);
  std::vector<hsize_t> ldims(ndim,0);
  std::vector<hsize_t> offset(ndim,0);

  // reverse order of array indexes to "c"
  std::reverse_copy(in_dims_global.begin(),in_dims_global.end(),gdims.begin());
  std::reverse_copy(in_dims_local.begin(),in_dims_local.end(),ldims.begin());
  std::reverse_copy(in_offset.begin(),in_offset.end(),offset.begin());

  //Log::debugAll("{} {} {}",gdims,ldims,offset);

  // identifiers
  hid_t dset_id, filespace, memspace, dtype_id;
  //herr_t status;

  // get hd5 datatype
  if (std::is_same<T,FSCAL>::value) dtype_id = H5T_NATIVE_DOUBLE;
  if (std::is_same<T,float>::value) dtype_id = H5T_NATIVE_FLOAT;
  if (std::is_same<T,int>::value) dtype_id = H5T_NATIVE_INT;

  // open dataset
  dset_id = H5Dopen(file_id,path.c_str(),H5P_DEFAULT);

  // get filespace
  filespace = H5Dget_space(dset_id);    /* Get filespace handle first. */

  // check dimensions
  checkDataDimensions(filespace, ndim, gdims);

  // specify memory space for this mpi rank
  memspace =  H5Screate_simple(ndim, ldims.data(), NULL);

  // create hyperslab and read data
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset.data(), NULL, ldims.data(), NULL);
  H5Dread(dset_id, dtype_id, memspace, filespace, H5P_DEFAULT, data);

  // close identifiers
  H5Sclose(memspace);
  H5Sclose(filespace);
  H5Dclose(dset_id);
}

template <typename T>
template <typename S>
void h5Writer<T>::readAttribute(std::string name, S& data){
  hid_t gpid,att_id,dtype_id;

  if (std::is_same<S,FSCAL>::value) dtype_id = H5T_NATIVE_DOUBLE;
  if (std::is_same<S,float>::value) dtype_id = H5T_NATIVE_FLOAT;
  if (std::is_same<S,int>::value) dtype_id = H5T_NATIVE_INT;

  //char tmp[96];
  //att_id = H5Aopen(gpid,"Name",H5P_DEFAULT);
  //dtype_id = H5Aget_type(att_id);
  //dtypemem_id = H5Tget_native_type(dtype_id,H5T_DIR_ASCEND);
  //H5Aread(att_id,dtypemem_id,tmp);
  //H5Aclose(att_id);

  gpid = H5Gopen(file_id,"Properties",H5P_DEFAULT);
  att_id = H5Aopen(gpid,name.c_str(),H5P_DEFAULT);
  H5Aread(att_id,dtype_id,&data);

  H5Aclose(att_id);
  H5Gclose(gpid);
}

template <>
template <>
void h5Writer<float>::readAttribute(std::string name, std::string& data){
  hid_t gpid,att_id,dtype_id,dtypemem_id;

  gpid = H5Gopen(file_id,"Properties",H5P_DEFAULT);
  char tmp[128];
  att_id = H5Aopen(gpid,name.c_str(),H5P_DEFAULT);
  dtype_id = H5Aget_type(att_id);
  dtypemem_id = H5Tget_native_type(dtype_id,H5T_DIR_ASCEND);
  H5Aread(att_id,dtypemem_id,tmp);
  data = tmp;
  
  H5Aclose(att_id);
  H5Gclose(gpid);
}

template <>
template <>
void h5Writer<double>::readAttribute(std::string name, std::string& data){
  hid_t gpid,att_id,dtype_id,dtypemem_id;

  gpid = H5Gopen(file_id,"Properties",H5P_DEFAULT);
  char tmp[128];
  att_id = H5Aopen(gpid,name.c_str(),H5P_DEFAULT);
  dtype_id = H5Aget_type(att_id);
  dtypemem_id = H5Tget_native_type(dtype_id,H5T_DIR_ASCEND);
  H5Aread(att_id,dtypemem_id,tmp);
  data = tmp;
  
  H5Aclose(att_id);
  H5Gclose(gpid);
}



template <typename T>
void h5Writer<T>::openGroup(std::string name){
  group_id = H5Gcreate(file_id, name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
}

template <typename T>
void h5Writer<T>::closeGroup(){
  H5Gclose(group_id);
}

template <typename T>
void h5Writer<T>::checkDataDimensions(hid_t filespace, int ndim, std::vector<hsize_t> dims){
  int rank;
  //herr_t status;
  hsize_t dimsg[ndim];

  // get and check rank
  rank      = H5Sget_simple_extent_ndims(filespace);
  if (ndim != rank){
    Log::error("Expected {}D restart but found {}D restart.",ndim,rank);
    exit(EXIT_FAILURE);
  }

  // get dimensions
  H5Sget_simple_extent_dims(filespace, dimsg, NULL);

  // check dimensions
  for (int d=0; d<ndim; ++d){
    if (dims[d] != dimsg[d]){
      Log::error("Extents of restart file different than expected. (Got {} but expected {} in direction {}.)\n",dimsg[d],dims[d],d);
      exit (EXIT_FAILURE);
    }
  }
}

template class h5Writer<float>;
template class h5Writer<double>;

template void h5Writer<float>::writeAttribute<int>(std::string,int);
template void h5Writer<float>::writeAttribute<FSCAL>(std::string,FSCAL);
//template void h5Writer<float>::writeAttribute<std::string>(std::string,std::string);

template void h5Writer<double>::writeAttribute<int>(std::string,int);
template void h5Writer<double>::writeAttribute<FSCAL>(std::string,FSCAL);
//template void h5Writer<FSCAL>::writeAttribute<std::string>(std::string,std::string);

template void h5Writer<float>::readAttribute<int>(std::string,int&);
template void h5Writer<float>::readAttribute<FSCAL>(std::string,FSCAL&);
//template void h5Writer<float>::readAttribute<std::string>(std::string,std::string&);

template void h5Writer<double>::readAttribute<int>(std::string,int&);
template void h5Writer<double>::readAttribute<FSCAL>(std::string,FSCAL&);
//template void h5Writer<FSCAL>::readAttribute<std::string>(std::string,std::string&);
