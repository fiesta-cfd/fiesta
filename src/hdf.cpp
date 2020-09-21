#include "hdf.hpp"
#include "kokkosTypes.hpp"
#include "hdf5.h"
#include "output.hpp"
#include "rkfunction.hpp"
#include "particle.hpp"
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
//#ifndef NOMPI
//#include "mpi.h"
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
    fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" " "Precision=\"4\" Format=\"HDF\">\n",dims[0],dims[2],dims[2]);
  fprintf(xmf, "        %s\n", path.c_str());
  fprintf(xmf, "       </DataItem>\n");
}
//void write_xmf(string fname, string hname, double time, int ndim, int *dims, int nv, int nvt){
void write_xmf(string fname, string hname, double time, struct inputConfig &cf, int np, vector<string> vNames, vector<string> vxNames ){

    int dims[cf.ndim];
    invertArray(cf.ndim,dims,cf.globalCellDims);
    int ndim = cf.ndim;
    int nv = cf.nv;
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
    for (int var = ndim; var < vxNames.size(); ++var) {
      fprintf(xmf, "     <Attribute Name=\"%s\" AttributeType=\"Scalar\" " "Center=\"Cell\">\n",vxNames[var].c_str());
      path.str("");
      path << hname << ":/Solution/Variable" << setw(2) << setfill('0') << var+nvt;
      writeDataItem(xmf, path.str(), ndim, dims);
      fprintf(xmf, "     </Attribute>\n");
    }

    // End Field Variables
    fprintf(xmf, "   </Grid>\n");

    //Particles
    if (cf.particle == 1){
      fprintf(xmf, "   <Grid Name=\"particles\" >\n");
        fprintf(xmf, "     <Topology TopologyType=\"Polyvertex\" NumberOfElements=\"%d\"/>\n", np);
      if (ndim == 2){
        fprintf(xmf, "     <Geometry GeometryType=\"X_Y\">\n");
      }else{
        fprintf(xmf, "     <Geometry GeometryType=\"X_Y_Z\">\n");
      }
      for (int d=0; d<ndim; ++d){
        path.str("");
        path << hname << ":/Particles/Dim" << d;
        writeDataItem(xmf, path.str(), 1, &np);
      }
      fprintf(xmf, "     </Geometry>\n");

      fprintf(xmf, "     <Attribute Name=\"State\" " "AttributeType=\"Scalar\" Center=\"Node\">\n");
      path.str("");
      path << hname << ":/Particles/state";
      writeDataItem(xmf, path.str(), 1, &np);
      fprintf(xmf, "     </Attribute>\n");
 
      fprintf(xmf, "     <Attribute Name=\"u\" " "AttributeType=\"Scalar\" Center=\"Node\">\n");
      path.str("");
      path << hname << ":/Particles/xvel";
      writeDataItem(xmf, path.str(), 1, &np);
      fprintf(xmf, "     </Attribute>\n");
 
      fprintf(xmf, "     <Attribute Name=\"r\" " "AttributeType=\"Scalar\" Center=\"Node\">\n");
      path.str("");
      path << hname << ":/Particles/r";
      writeDataItem(xmf, path.str(), 1, &np);
      fprintf(xmf, "     </Attribute>\n");

      fprintf(xmf, "     <Attribute Name=\"m\" " "AttributeType=\"Scalar\" Center=\"Node\">\n");
      path.str("");
      path << hname << ":/Particles/m";
      writeDataItem(xmf, path.str(), 1, &np);
      fprintf(xmf, "     </Attribute>\n");

      fprintf(xmf, "   </Grid>\n");
    }
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
  herr_t status;

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
  status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, dims_local, NULL);

  // create property list for collective dataset write
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

  // write data
  status = H5Dwrite(dset_id, dtype_id, memspace, filespace, plist_id, data);

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

  psp = (float*)malloc(cf.p_np*sizeof(float));
  pdp = (double*)malloc(cf.p_np*sizeof(double));
  pi = (int*)malloc(cf.p_np*sizeof(int));

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
  stringstream baseName, xmfName, hdfName;
  baseName << name << "-" << setw(pad) << setfill('0') << tdx;
  hdfName << baseName.str() << ".h5";
  xmfName << baseName.str() << ".xmf";


  // identifiers
  hid_t file_id, group_id;

  // write a message
  if (cf.rank == 0) {
    cout << c(0, YEL) << left << setw(22)
         << "    Writing File: " << c(0, NON) << c(0, CYA) << left
         << "'" + hdfName.str() + "'" << c(0, NON) << endl;
  }

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
  for (int vn = 0; vn < cf.nvt; ++vn) {
    // Format Dataset Name
    stringstream vname;
    vname << "Variable" << setw(2) << setfill('0') << vn;

    Kokkos::deep_copy(varH,f->var);
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
  for (int vn = 0; vn < f->varxNames.size(); ++vn) {
    // Format Dataset Name
    stringstream vname;
    vname << "Variable" << setw(2) << setfill('0') << vn+cf.nvt;

    Kokkos::deep_copy(varxH,f->varx);
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

  if (cf.particle == 1){

    // get global particle count and offsets
    int nParticles[cf.numProcs];
    MPI_Allgather(&cf.p_np,1,MPI_INT,nParticles,1,MPI_INT,cf.comm);
    hsize_t pOffset = 0;
    for (int i=0; i<cf.rank; ++i)
      pOffset += nParticles[i];
    hsize_t globalPartDim =0;
    for (int i=0; i<cf.numProcs; ++i)
      globalPartDim += nParticles[i];
    hsize_t localPartDims = cf.p_np;
    
    // temporary arrays for particle writing
    float * partX = (float*)malloc(cf.p_np*sizeof(float));
    int * partS = (int*)malloc(cf.p_np*sizeof(int));
    Kokkos::View<particleStruct2D *>::HostMirror parH =
        Kokkos::create_mirror_view(f->particles);
    Kokkos::deep_copy(parH,f->particles);

    // create particle group and write partivle data
    group_id = H5Gcreate(file_id, "/Particles", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    {
      string vname("Dim0");
      for (int p=0; p<cf.p_np; ++p) partX[p] = parH(p).x;
      write_h5(group_id, vname, 1, &globalPartDim, &localPartDims, &pOffset, partX); 

      vname = string("Dim1");
      for (int p=0; p<cf.p_np; ++p) partX[p] = parH(p).y;
      write_h5(group_id, vname, 1, &globalPartDim, &localPartDims, &pOffset, partX); 

      vname =string("state");
      for (int p=0; p<cf.p_np; ++p) partS[p] = parH(p).state;
      write_h5(group_id, vname, 1, &globalPartDim, &localPartDims, &pOffset, partS); 

      vname =string("xvel");
      for (int p=0; p<cf.p_np; ++p) partX[p] = parH(p).u;
      write_h5(group_id, vname, 1, &globalPartDim, &localPartDims, &pOffset, partX); 

      vname =string("yvel");
      for (int p=0; p<cf.p_np; ++p) partX[p] = parH(p).v;
      write_h5(group_id, vname, 1, &globalPartDim, &localPartDims, &pOffset, partX); 

      vname =string("r");
      for (int p=0; p<cf.p_np; ++p) partX[p] = parH(p).r;
      write_h5(group_id, vname, 1, &globalPartDim, &localPartDims, &pOffset, partX); 

      vname =string("m");
      for (int p=0; p<cf.p_np; ++p) partX[p] = parH(p).m;
      write_h5(group_id, vname, 1, &globalPartDim, &localPartDims, &pOffset, partX); 
    }
    H5Gclose(group_id);
    write_xmf(xmfName.str(), hdfName.str(), time, cf, globalPartDim,f->varNames,f->varxNames);
  } else {
    write_xmf(xmfName.str(), hdfName.str(), time, cf, 0,f->varNames,f->varxNames);
  }

  close_h5(file_id);
}

// check dimensions of restart file against expected values
void checkDataDimensions(hid_t filespace, int ndim, hsize_t* dims){

  int rank;
  herr_t status;
  hsize_t dimsg[ndim];

  // get and check rank
  rank      = H5Sget_simple_extent_ndims(filespace);
  if (ndim != rank){
    printf("Number of dimensions of restart file different than expected\n");
    exit(EXIT_FAILURE);
  }

  // get dimensions
  status    = H5Sget_simple_extent_dims(filespace, dimsg, NULL);

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
  herr_t status;

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
  status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);
  status = H5Dread(dset_id, dtype_id, memspace, filespace, H5P_DEFAULT, data);

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
void hdfWriter::writeSPGrid(struct inputConfig cf, const FS4D gridD,
                            const char *fname) {}

// deprecated 
void hdfWriter::writeGrid(struct inputConfig cf, const FS4D gridD,
                          const char *fname) {}
