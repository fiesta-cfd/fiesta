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
//#endif

using namespace std;

template <typename T, typename C>
void invertArray(int ndim, T* out, C* in){
  for (int i=0; i<ndim; ++i){
    out[i] = in[ndim-1-i];
  }
}

void writeDataItem(FILE* xmf, string path, int ndim, int* dims){
  if (ndim == 1)
    fprintf(xmf, "       <DataItem Dimensions=\"%d\" NumberType=\"Float\" " "Precision=\"4\" Format=\"HDF\">\n",dims[0]);
  else if (ndim == 2)
    fprintf(xmf, "       <DataItem Dimensions=\"%d %d\" NumberType=\"Float\" " "Precision=\"4\" Format=\"HDF\">\n",dims[0],dims[1]);
  else
    fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" " "Precision=\"4\" Format=\"HDF\">\n",dims[0],dims[1],dims[2]);
  fprintf(xmf, "        %s\n", path.c_str());
  fprintf(xmf, "       </DataItem>\n");
}
void write_xmf(string fname, string hname, double time,
               struct inputConfig &cf, vector<string> vNames, vector<string> vxNames ){

    int dims[cf.ndim];
    invertArray(cf.ndim,dims,cf.globalCellDims);
    int ndim = cf.ndim;
    //int nv = cf.nv;
    int nvt = cf.nvt;

    int gdims[ndim];
    for (int i=0; i<ndim; ++i) gdims[i] = dims[i]+1;
    stringstream path;

    FILE *xmf = 0;
    xmf = fopen(fname.c_str(), "w");

    // Header
    fprintf(xmf, "<?xml version=\"1.0\" ?>\n");
    fprintf(xmf, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
    fprintf(xmf, "<Xdmf Version=\"3.0\">\n");
    fprintf(xmf, " <Domain>\n");


    // Grid Header
    fprintf(xmf, "   <Grid Name=\"mesh1\" GridType=\"Uniform\">\n");
    fprintf(xmf, "     <Time Value=\"%e\" />\n",time);
    if (ndim == 2){
      fprintf(xmf, "     <Topology TopologyType=\"2DSMesh\" NumberOfElements=\"%d %d\"/>\n", dims[0]+1, dims[1]+1);
      fprintf(xmf, "     <Geometry GeometryType=\"X_Y\">\n");
    }else{
      fprintf(xmf, "     <Topology TopologyType=\"3DSMesh\" NumberOfElements=\"%d %d %d\"/>\n", dims[0]+1, dims[1]+1, dims[2]+1);
      fprintf(xmf, "     <Geometry GeometryType=\"X_Y_Z\">\n");
    }

    // Grid Coordinate Arrays
    for (int d=0; d<ndim; ++d){
      path.str("");
      path << hname << ":/Grid/Dimension" << d;
      writeDataItem(xmf, path.str(), ndim, gdims);
    }
    fprintf(xmf, "     </Geometry>\n");

    // Momentum Vector
    fprintf(xmf, "     <Attribute Name=\"Momentum\" AttributeType=\"Vector\" " "Center=\"Cell\">\n");
    if (ndim == 2)
      fprintf(xmf, "      <DataItem Dimensions=\"%d %d 2\" Function=\"JOIN($0,$1)\" " "ItemType=\"Function\">\n",dims[0],dims[1]);
    else
      fprintf(xmf, "      <DataItem Dimensions=\"%d %d %d 3\" Function=\"JOIN($0,$1,$2)\" " "ItemType=\"Function\">\n", dims[0], dims[1],dims[2]);

    for (int d=0; d<ndim; ++d){
      path.str("");
      path << hname << ":/Solution/Variable" << setw(2) << setfill('0') << d;
      writeDataItem(xmf, path.str(), ndim, dims);
    }
    fprintf(xmf, "      </DataItem>\n");
    fprintf(xmf, "     </Attribute>\n");

    // Other Variables
    fprintf(xmf, "     <Attribute Name=\"Energy\" AttributeType=\"Scalar\" " "Center=\"Cell\">\n");
    path.str("");
    path << hname << ":/Solution/Variable" << setw(2) << setfill('0') << ndim;
    writeDataItem(xmf, path.str(), ndim, dims);
    fprintf(xmf, "     </Attribute>\n");

    for (int var = ndim+1; var < nvt; ++var) {
      fprintf(xmf, "     <Attribute Name=\"%s\" AttributeType=\"Scalar\" " "Center=\"Cell\">\n",vNames[var].c_str());
      path.str("");
      path << hname << ":/Solution/Variable" << setw(2) << setfill('0') << var;
      writeDataItem(xmf, path.str(), ndim, dims);
      fprintf(xmf, "     </Attribute>\n");
    }

    // Velocity Vector
    fprintf(xmf, "     <Attribute Name=\"Velocity\" AttributeType=\"Vector\" " "Center=\"Cell\">\n");
    if (ndim == 2)
      fprintf(xmf, "      <DataItem Dimensions=\"%d %d 2\" Function=\"JOIN($0,$1)\" " "ItemType=\"Function\">\n",dims[0],dims[1]);
    else
      fprintf(xmf, "      <DataItem Dimensions=\"%d %d %d 3\" Function=\"JOIN($0,$1,$2)\" " "ItemType=\"Function\">\n", dims[0], dims[1],dims[2]);

    for (int d=0; d<ndim; ++d){
      path.str("");
      path << hname << ":/Solution/Variable" << setw(2) << setfill('0') << nvt+d;
      writeDataItem(xmf, path.str(), ndim, dims);
    }
    fprintf(xmf, "      </DataItem>\n");
    fprintf(xmf, "     </Attribute>\n");

    // Other Extra Variables
    for (size_t var = ndim; var < vxNames.size(); ++var) {
      fprintf(xmf, "     <Attribute Name=\"%s\" AttributeType=\"Scalar\" " "Center=\"Cell\">\n",vxNames[var].c_str());
      path.str("");
      path << hname << ":/Solution/Variable" << setw(2) << setfill('0') << var+nvt;
      writeDataItem(xmf, path.str(), ndim, dims);
      fprintf(xmf, "     </Attribute>\n");
    }

    // End Field Variables
    fprintf(xmf, "   </Grid>\n");

    fprintf(xmf, " </Domain>\n");
    fprintf(xmf, "</Xdmf>\n");
    fclose(xmf);
}

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
template <typename H> hid_t getH5Type();
template <> hid_t getH5Type<float>(){ return H5T_NATIVE_FLOAT; }
template <> hid_t getH5Type<double>(){ return H5T_NATIVE_DOUBLE; }
template <> hid_t getH5Type<int>(){ return H5T_NATIVE_INT; }


// writer function for hdf5
template <typename S>
void write_h5(hid_t group_id, string dname, int ndim,
                   hsize_t* dims_global, hsize_t* dims_local, hsize_t* offset, S* data){

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
                              double time, T* x, T* var, string name) {

  // calcualte string width for time index
  int pad = (int)log10(cf.nt) + 1;
  int idx;

  // invert dimension arrays to c order
  hsize_t offset[cf.ndim], gridCount[cf.ndim], cellCount[cf.ndim];
  hsize_t gridDims[cf.ndim], cellDims[cf.ndim];
  invertArray(cf.ndim,offset,cf.subdomainOffset);
  invertArray(cf.ndim,gridCount,cf.localGridDims);
  invertArray(cf.ndim,cellCount,cf.localCellDims);
  invertArray(cf.ndim,gridDims,cf.globalGridDims);
  invertArray(cf.ndim,cellDims,cf.globalCellDims);

  // format hdf5 and xdmf filenames
  stringstream baseName, xmfName, hdfName, hdfBaseName;
  baseName << name << "-" << setw(pad) << setfill('0') << tdx;
  hdfBaseName  << baseName.str() << ".h5";
  hdfName << cf.pathName << "/" << hdfBaseName.str();
  xmfName << cf.pathName << "/" << baseName.str() << ".xmf";


  // ############################################################## \\

  int ist=40;
  int jst=40;
  int kst=0;

  int ind=50;
  int jnd=60;
  int knd=cf.glbl_nck-1;

  int isz=ind-ist+1;
  int jsz=jnd-jst+1;
  int ksz=knd-kst+1;

  int ilo=0;
  int jlo=0;
  int klo=0;

  int ihi=cf.nci-1;
  int jhi=cf.ncj-1;
  int khi=cf.nck-1;

  int iex=1;
  int jex=1;
  int kex=1;

  int lo[3];
  int hi[3];
  int dms[3];
  int ofs[3];

  vector<int> actDims;
  string actDimNames;
  //if (ksz>1){ actDims.push_back(2); actDimNames.append(to_string(2));}
  //if (jsz>1){ actDims.push_back(1); actDimNames.append(to_string(1));}
  //if (isz>1){ actDims.push_back(0); actDimNames.append(to_string(0));}

  actDims.push_back(2); actDimNames.append(to_string(2));
  actDims.push_back(1); actDimNames.append(to_string(1));
  actDims.push_back(0); actDimNames.append(to_string(0));


  cf.log->debug("Active Dimensions: ",actDims.size(),"(",actDimNames,")");

  int myColor=MPI_UNDEFINED;
  if ( (cf.iEnd >=ist && cf.iStart <= ind)
    && (cf.jEnd >=jst && cf.jStart <= jnd)
    && (cf.kEnd >=kst && cf.kStart <= knd) ){
    myColor=1;
    if (ist > cf.iStart) ilo=ist-cf.iStart;
    if (jst > cf.jStart) jlo=jst-cf.jStart;
    if (kst > cf.kStart) klo=kst-cf.kStart;

    if (ind < cf.iEnd) ihi=cf.nci -(cf.iEnd-ind);
    if (jnd < cf.jEnd) jhi=cf.ncj -(cf.jEnd-jnd);
    if (knd < cf.kEnd) khi=cf.nck -(cf.kEnd-knd);

    iex=ihi-ilo+1;
    jex=jhi-jlo+1;
    kex=khi-klo+1;

    lo[0]=ilo;
    lo[1]=jlo;
    lo[2]=klo;

    hi[0]=ihi;
    hi[1]=jhi;
    hi[2]=khi;

    dms[0]=isz;
    dms[1]=jsz;
    dms[2]=ksz;

    ofs[0]=0;
    ofs[1]=0;
    ofs[2]=0;

    if (ist<cf.iStart) ofs[0] = cf.iStart-ist;
    if (jst<cf.jStart) ofs[1] = cf.jStart-jst;
    if (kst<cf.kStart) ofs[2] = cf.kStart-kst;
  }
  //if ( (ti >=cf.iStart && ti< cf.iEnd) && (tj >=cf.jStart && tj< cf.jEnd) ) myColor=1;
  MPI_Comm sliceComm;
  MPI_Comm_split(cf.comm,myColor,cf.rank,&sliceComm);
  if (myColor==1){
    hsize_t slcgDims[actDims.size()];
    hsize_t slcgCount[actDims.size()];

    hsize_t slcDims[actDims.size()];
    hsize_t slcCount[actDims.size()];
    hsize_t slcOffset[actDims.size()];

    Kokkos::deep_copy(gridH,f->grid);
    Kokkos::deep_copy(varH,f->var);

    string countStr(": ");
    int hdx=0;
    for (int d : actDims){
      slcCount[hdx]=hi[d]-lo[d]+1;
      slcDims[hdx]=dms[d];

      slcgCount[hdx]=hi[d]-lo[d]+2;
      slcgDims[hdx]=dms[d]+1;

      slcOffset[hdx]=ofs[d];

      countStr.append(to_string(slcCount[hdx]));
      countStr.append(",");
      countStr.append(to_string(slcDims[hdx]));
      countStr.append(",");
      countStr.append(to_string(slcOffset[hdx]));
      countStr.append(" : ");
      hdx++;
    }

    int slcLen = (iex*jex*kex);
    T* slcData = (T*)malloc(slcLen*sizeof(T*));

    slcLen = (ihi-ilo+2)*(jhi-jlo+2)*(khi-klo+2);
    T* slcGrid = (T*)malloc(slcLen*sizeof(T*));

    int ii,jj,kk;
    int koffset;

    stringstream lineName,lineBase,lineHDF,lineXMF;
    lineBase << "line-" << setw(pad) << setfill('0') << tdx;
    lineName << cf.pathName << "/line/" << lineBase.str() << ".h5";
    lineXMF << cf.pathName << "/line/" << lineBase.str() << ".xmf";
    lineHDF << lineBase.str() << ".h5";

    cf.log->debug("[",cf.t,"] ","Writing ",lineName.str());

    hid_t linefile_id = openHDF5ForWrite(sliceComm, MPI_INFO_NULL, lineName.str());

    hid_t linegroup_id = H5Gcreate(linefile_id, "/Grid", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    
    for (int vn=0; vn<cf.ndim; ++vn){
      for (int i = ilo; i < ihi+2; ++i) {
        for (int j = jlo; j < jhi+2; ++j) {
          for (int k = klo; k < khi+2; ++k) {
            ii=i-ilo;
            jj=j-jlo;
            kk=k-klo;
            idx = ((iex+1) * (jex+1)) * kk + (iex+1) * jj + ii;
            slcGrid[idx] = gridH(i, j, k, vn);
          }
        }
      }
      stringstream lineStr;
      lineStr << "Dimension" << vn;
      int slcDim = actDims.size();
      cf.log->debug("Dataset name: ",lineStr.str());
      write_h5<T>(linegroup_id, lineStr.str(), slcDim, slcgDims, slcgCount, slcOffset, slcGrid); 
    }

    H5Gclose(linegroup_id);
    linegroup_id = H5Gcreate(linefile_id, "/Solution", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    if (cf.ndim == 3) koffset = cf.ng;
    else koffset = 0;
    for (int vn=0; vn<cf.nvt; ++vn){
      for (int i = ilo; i < ihi+1; ++i) {
        for (int j = jlo; j < jhi+1; ++j) {
          for (int k = klo; k < khi+1; ++k) {
            ii=i-ilo;
            jj=j-jlo;
            kk=k-klo;
            idx = (iex * jex) * kk + iex * jj + ii;
            slcData[idx] = varH(i+cf.ng, j+cf.ng, k+koffset, vn);
          }
        }
      }
      stringstream lineStr;
      lineStr << "Variable" << setw(2) << setfill('0') << vn;
      int slcDim = actDims.size();
      cf.log->debug("Dataset name: ",lineStr.str());
      write_h5<T>(linegroup_id, lineStr.str(), slcDim, slcDims, slcCount, slcOffset, slcData); 
    }

    if (cf.ndim == 3) koffset = cf.ng;
    else koffset = 0;
    for (size_t vn = 0; vn < f->varxNames.size(); ++vn) {
      for (int i = ilo; i < ihi+1; ++i) {
        for (int j = jlo; j < jhi+1; ++j) {
          for (int k = klo; k < khi+1; ++k) {
            ii=i-ilo;
            jj=j-jlo;
            kk=k-klo;
            idx = (iex * jex) * kk + iex * jj + ii;
            slcData[idx] = varxH(i+cf.ng, j+cf.ng, k+koffset, vn);
          }
        }
      }
      stringstream lineStr;
      lineStr << "Variable" << setw(2) << setfill('0') << vn+cf.nvt;
      int slcDim = actDims.size();
      cf.log->debug("Dataset name: ",lineStr.str());
      write_h5<T>(linegroup_id, lineStr.str(), slcDim, slcDims, slcCount, slcOffset, slcData); 
    }

    cf.log->debug("[",cf.t,"] ","Writing ",lineXMF.str());
    vector<int> xmfDims;
    for (auto d : slcDims) xmfDims.push_back(d);
    //xmfDims.push_back(iex);
    //xmfDims.push_back(jiex);
    //xmfDims.push_back(kex);

    writeXMF(lineXMF.str(), lineHDF.str(), time, actDims.size(), actDims, xmfDims, cf.nvt, f->varNames,f->varxNames);
    //writeXMF(xmfName.str(), hdfBaseName.str(), time, cf, f->varNames,f->varxNames);

    //int test;
    int sliceRank;
    int sliceWorld;
    MPI_Comm_rank(sliceComm,&sliceRank);
    MPI_Comm_size(sliceComm,&sliceWorld);

    cout << "@#@#@#@#@# " << sliceRank << countStr << endl;
    
    cout << "#~#~#~#~#~#~#~#~#~#~" << cf.rank << "-" << sliceRank << "/" << sliceWorld << ": "
         << ilo << "/" << ihi << " "
         << jlo << "/" << jhi << " "
         << klo << "/" << khi << endl;

    //cout << "########################" << cf.rank << " : " << myColor << endl;

    //T* linex = (T*)malloc(cf.nk*sizeof(T));

    ////stringstream lineName;
    ////lineName << cf.pathName << "/line/" << "line-" << setw(pad) << setfill('0') << tdx << ".h5";



    //int sliceDim=1;
    //hsize_t sliceDims = cf.glbl_nk;
    //hsize_t sliceCount=cf.nk;
    //hsize_t sliceOffset=cf.kStart;




    //int activeDims[] = {2};

    //for (int vn : activeDims){
    //  for (int k = 0; k<cf.nk; ++k){
    //    //idx = (cf.ni * cf.nj) * k + cf.ni * tj + ti;
    //    linex[k]=gridH(ti,tj,k,vn);
    //  }
    //  string lineStr("Dimension",vn);

    //  write_h5<T>(linegroup_id, lineStr, sliceDim, &sliceDims, &sliceCount, &sliceOffset, linex); 
    //}
    H5Gclose(linegroup_id);

    close_h5(linefile_id);
  }

  // ############################################################## \\

  // identifiers
  hid_t file_id, group_id;

  // write a message
  if (cf.rank == 0) {
    cout << c(0, YEL) << left << setw(22)
         << "    Writing File: " << c(0, NON) << c(0, CYA) << left
         << "'" + hdfName.str() + "'" << c(0, NON) << endl;
  }

  //{
    //stringstream message;
    //message << "Writing " << hdfName.str();
    //cf.log->message(message.str());
  //}
  cf.log->message("[",cf.t,"] ","Writing ",hdfName.str());


  // open file
  file_id = openHDF5ForWrite(cf.comm, MPI_INFO_NULL, hdfName.str());

  // create grid group and write grid data
  group_id = H5Gcreate(file_id, "/Grid", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  Kokkos::deep_copy(gridH,f->grid);
  for (int vn = 0; vn < cf.ndim; ++vn) {
    stringstream vname;
    vname << "Dimension" << setw(1) << setfill('0') << vn;
    for (int i = 0; i < cf.ni; ++i) {
      for (int j = 0; j < cf.nj; ++j) {
        for (int k = 0; k < cf.nk; ++k) {
          idx = (cf.ni * cf.nj) * k + cf.ni * j + i;
          x[idx] = gridH(i, j, k, vn);
        }
      }
    }
    write_h5<T>(group_id, vname.str(), cf.ndim, gridDims, gridCount, offset, x); 
  }
  H5Gclose(group_id);

  // create solution group and write solution data
  group_id = H5Gcreate(file_id, "/Solution", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  Kokkos::deep_copy(varH,f->var);
  for (int vn = 0; vn < cf.nvt; ++vn) {
    // Format Dataset Name
    stringstream vname;
    vname << "Variable" << setw(2) << setfill('0') << vn;

    int koffset;
    if (cf.ndim == 3) koffset = cf.ng;
    else koffset = 0;
    for (int k = koffset; k < cf.nck + koffset; ++k) {
      for (int j = cf.ng; j < cf.ncj + cf.ng; ++j) {
        for (int i = cf.ng; i < cf.nci + cf.ng; ++i) {
          int ii = i - cf.ng;
          int jj = j - cf.ng;
          int kk = k - koffset;
          idx = (cf.nci * cf.ncj) * kk + cf.nci * jj + ii;
          var[idx] = varH(i, j, k, vn);
        }
      }
    }
    write_h5<T>(group_id, vname.str(), cf.ndim, cellDims, cellCount, offset, var); 
  }
  Kokkos::deep_copy(varxH,f->varx);
  for (size_t vn = 0; vn < f->varxNames.size(); ++vn) {
    // Format Dataset Name
    stringstream vname;
    vname << "Variable" << setw(2) << setfill('0') << vn+cf.nvt;

    int koffset;
    if (cf.ndim == 3) koffset = cf.ng;
    else koffset = 0;
    for (int k = koffset; k < cf.nck + koffset; ++k) {
      for (int j = cf.ng; j < cf.ncj + cf.ng; ++j) {
        for (int i = cf.ng; i < cf.nci + cf.ng; ++i) {
          int ii = i - cf.ng;
          int jj = j - cf.ng;
          int kk = k - koffset;
          idx = (cf.nci * cf.ncj) * kk + cf.nci * jj + ii;
          var[idx] = varxH(i, j, k, vn);
        }
      }
    }
    write_h5<T>(group_id, vname.str(), cf.ndim, cellDims, cellCount, offset, var); 
  }
  H5Gclose(group_id);

  //{
  //  stringstream message;
  //  message << "Writing " << xmfName.str();
  //  cf.log->message(message.str());
  //}
  cf.log->message("[",cf.t,"] ","Writing ",xmfName.str());
  write_xmf(xmfName.str(), hdfBaseName.str(), time, cf, f->varNames,f->varxNames);

  close_h5(file_id);
}

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
  dtype_id = getH5Type<T>();

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

hdfWriter::~hdfWriter(){
  free(xdp);
  free(xsp);
  free(vsp);
  free(vdp);
  free(psp);
  free(pdp);
  free(pi );
  free(readV);
       
  //delete gridH;
  //delete varH;
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

// deprecated 
//void hdfWriter::writeSPGrid(struct inputConfig cf, const FS4D gridD,
//                            const char *fname) {}

// deprecated 
//void hdfWriter::writeGrid(struct inputConfig cf, const FS4D gridD,
//                          const char *fname) {}
