
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

// close an hdf5 file
template <typename T>
void h5Writer<T>::close(){
  H5Fclose(file_id);
}

// writer function for hdf5
template <typename T>
void h5Writer<T>::write(std::string dname, int ndim,
  std::vector<size_t> in_dims_global, std::vector<size_t> in_dims_local, std::vector<size_t> in_offset, std::vector<T>& data){

  std::vector<hsize_t> gdims(ndim,0);
  std::vector<hsize_t> ldims(ndim,0);
  std::vector<hsize_t> offset(ndim,0);

  // reverse order of array indexes to "c"
  std::reverse_copy(in_dims_global.begin(),in_dims_global.end(),gdims.begin());
  std::reverse_copy(in_dims_local.begin(),in_dims_local.end(),ldims.begin());
  std::reverse_copy(in_offset.begin(),in_offset.end(),offset.begin());

  // identifiers
  hid_t filespace, memspace, dset_id, plist_id, dtype_id;

  // get type id
  if (std::is_same<T,double>::value) dtype_id = H5T_NATIVE_DOUBLE;
  if (std::is_same<T,float>::value) dtype_id = H5T_NATIVE_FLOAT;
  if (std::is_same<T,int>::value) dtype_id = H5T_NATIVE_INT;

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

template <typename T>
void h5Writer<T>::openGroup(std::string name){
  group_id = H5Gcreate(file_id, name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
}

template <typename T>
void h5Writer<T>::closeGroup(){
  H5Gclose(group_id);
}

template class h5Writer<float>;
template class h5Writer<double>;
