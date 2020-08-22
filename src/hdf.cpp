#include "hdf.hpp"
#include "fiesta.hpp"
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

using namespace std;

void writeDataItem(FILE* xmf, string path, int ndim, int* dims){
  if (ndim == 2)
    fprintf(xmf, "       <DataItem Dimensions=\"%d %d\" NumberType=\"Float\" " "Precision=\"4\" Format=\"HDF\">\n",dims[0],dims[1]);
  else
    fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" " "Precision=\"4\" Format=\"HDF\">\n",dims[0],dims[2],dims[2]);
  fprintf(xmf, "        %s\n", path.c_str());
  fprintf(xmf, "       </DataItem>\n");
}
void write_xmf(string fname, string hname, int ndim, int *dims, int nv, int nvt){
    int gdims[ndim];
    for (int i=0; i<ndim; ++i) gdims[i] = dims[i]+1;
    stringstream path;

    FILE *xmf = 0;
    xmf = fopen(fname.c_str(), "w");
    fprintf(xmf, "<?xml version=\"1.0\" ?>\n");
    fprintf(xmf, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
    fprintf(xmf, "<Xdmf Version=\"3.0\">\n");
    fprintf(xmf, " <Domain>\n");

    fprintf(xmf, "   <Grid Name=\"mesh1\" GridType=\"Uniform\">\n");
    if (ndim == 2){
      fprintf(xmf, "     <Topology TopologyType=\"2DSMesh\" NumberOfElements=\"%d %d\"/>\n", dims[0]+1, dims[1]+1);
      fprintf(xmf, "     <Geometry GeometryType=\"X_Y\">\n");
    }else{
      fprintf(xmf, "     <Topology TopologyType=\"3DSMesh\" NumberOfElements=\"%d %d %d\"/>\n", dims[0]+1, dims[1]+1, dims[2]+1);
      fprintf(xmf, "     <Geometry GeometryType=\"X_Y_Z\">\n");
    }

    for (int d=0; d<ndim; ++d){
      path.str("");
      path << hname << ":/Grid/Dimension" << d;
      writeDataItem(xmf, path.str(), ndim, gdims);
    }
    fprintf(xmf, "     </Geometry>\n");

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

    fprintf(xmf, "     <Attribute Name=\"Energy\" AttributeType=\"Scalar\" " "Center=\"Cell\">\n");
    path.str("");
    path << hname << ":/Solution/Variable" << setw(2) << setfill('0') << ndim;
    writeDataItem(xmf, path.str(), ndim, dims);
    fprintf(xmf, "     </Attribute>\n");

    for (int var = ndim+1; var < nv; ++var) {
      fprintf(xmf, "     <Attribute Name=\"Density%02d\" AttributeType=\"Scalar\" " "Center=\"Cell\">\n", var - ndim);
      path.str("");
      path << hname << ":/Solution/Variable" << setw(2) << setfill('0') << var;
      writeDataItem(xmf, path.str(), ndim, dims);
      fprintf(xmf, "     </Attribute>\n");
    }
    for (int var = nv; var < nvt; ++var) {
      fprintf(xmf, "     <Attribute Name=\"Variable%02d\" " "AttributeType=\"Scalar\" Center=\"Cell\">\n", var - nv);
      path.str("");
      path << hname << ":/Solution/Variable" << setw(2) << setfill('0') << var;
      writeDataItem(xmf, path.str(), ndim, dims);
      fprintf(xmf, "     </Attribute>\n");
    }
    fprintf(xmf, "   </Grid>\n");
    fprintf(xmf, " </Domain>\n");
    fprintf(xmf, "</Xdmf>\n");
    fclose(xmf);
}

hid_t open_h5(MPI_Comm comm, MPI_Info info, string fname){
  hid_t fid,pid;

  pid = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(pid, comm, info);
  fid = H5Fcreate(fname.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, pid);
  H5Pclose(pid);

  return fid;
}

void close_h5(hid_t fid){
  H5Fclose(fid);
}

template <typename H> hid_t getH5Type();

template <> hid_t getH5Type<float>(){ return H5T_NATIVE_FLOAT; }
template <> hid_t getH5Type<double>(){ return H5T_NATIVE_DOUBLE; }
template <> hid_t getH5Type<int>(){ return H5T_NATIVE_INT; }

template <typename S>
void write_h5(hid_t group_id, string dname, S* data, int ndim,
                   hsize_t* dims_global, hsize_t* dims_local, hsize_t* offset){

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

fstWriter::fstWriter(struct inputConfig cf, rk_func *f) {
  xdp = (double *)malloc(cf.ni * cf.nj * cf.nk * sizeof(double));
  xsp = (float *)malloc(cf.ni * cf.nj * cf.nk * sizeof(float));

  vsp = (float *)malloc(cf.nci * cf.ncj * cf.nck * sizeof(float));
  v = (double *)malloc(cf.nci * cf.ncj * cf.nck * sizeof(double));

  readV = (double *)malloc(cf.nci * cf.ncj * cf.nck * cf.nv * sizeof(double));

  gridH = Kokkos::create_mirror_view(f->grid);
  varH = Kokkos::create_mirror_view(f->var);
}

void fstWriter::writeSPGrid(struct inputConfig cf, const FS4D gridD,
                            const char *fname) {}
void fstWriter::writeGrid(struct inputConfig cf, const FS4D gridD,
                          const char *fname) {}

template<typename T>
void fstWriter::writeHDF(struct inputConfig cf, rk_func *f, int tdx,
                              double time, T* x, T* var, string name) {

  int pad = (int)log10(cf.nt) + 1;

  // format hdf5 and xdmf filenames
  stringstream baseName, xmfName, hdfName;
  baseName << name << "-" << setw(pad) << setfill('0') << tdx;
  hdfName << baseName.str() << ".h5";
  xmfName << baseName.str() << ".xmf";

  hsize_t dimsg[cf.ndim];                 // dataset dimensions
  hsize_t dims[cf.ndim];                  // chunk dimensions
  hsize_t offset[cf.ndim];                  // chunk dimensions
  int idx;

  // get offsets
  if (cf.ndim == 2){
    offset[0] = cf.jStart;
    offset[1] = cf.iStart;
  }else{
    offset[0] = cf.kStart;
    offset[1] = cf.jStart;
    offset[2] = cf.iStart;
  }

  if (cf.ndim == 2){
    dimsg[0] = cf.glbl_nj;
    dimsg[1] = cf.glbl_ni;
    dims[0] = cf.nj;
    dims[1] = cf.ni;
  }else{
    dimsg[0] = cf.glbl_nk;
    dimsg[1] = cf.glbl_nj;
    dimsg[2] = cf.glbl_ni;
    dims[0] = cf.nk;
    dims[1] = cf.nj;
    dims[2] = cf.ni;
  }

  // identifiers
  hid_t file_id, group_id;

  if (cf.rank == 0) {
    cout << c(0, YEL) << left << setw(22)
         << "    Writing File: " << c(0, NON) << c(0, CYA) << left
         << "'" + hdfName.str() + "'" << c(0, NON) << endl;
  }

  // open file
  file_id = open_h5(cf.comm, MPI_INFO_NULL, hdfName.str());

  // create grid group
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
    write_h5<T>(group_id, vname.str(), x, cf.ndim, dimsg, dims, offset); 
  }
  H5Gclose(group_id);

  // create solution group
  if (cf.ndim == 2){
    dimsg[0] = cf.glbl_ncj;
    dimsg[1] = cf.glbl_nci;
    dims[0] = cf.ncj;
    dims[1] = cf.nci;
  }else{
    dimsg[0] = cf.glbl_nck;
    dimsg[1] = cf.glbl_ncj;
    dimsg[2] = cf.glbl_nci;
    dims[0] = cf.nck;
    dims[1] = cf.ncj;
    dims[2] = cf.nci;
  }

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
    write_h5<T>(group_id, vname.str(), var, cf.ndim, dimsg, dims, offset); 
  }
  H5Gclose(group_id);

  close_h5(file_id);

  int cdims[cf.ndim];
  if (cf.ndim == 2){
    cdims[0] = cf.ncj;
    cdims[1] = cf.nci;
  }else{
    cdims[0] = cf.nck;
    cdims[1] = cf.ncj;
    cdims[2] = cf.nci;
  }

  write_xmf(xmfName.str(), hdfName.str(), cf.ndim, cdims, cf.nv, cf.nvt);
}

void fstWriter::writeSolution(struct inputConfig cf, rk_func *f, int tdx,
                              double time) {
    writeHDF(cf, f, tdx, time, xsp, vsp, "sol");
}

void fstWriter::writeRestart(struct inputConfig cf, rk_func *f, int tdx,
                             double time) {
    writeHDF(cf, f, tdx, time, xdp, v, "restart");
}

void fstWriter::readSolution(struct inputConfig cf, FS4D &gridD, FS4D &varD) {}
