#include "fiesta.hpp"
#include "hdf.hpp"
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "output.hpp"
#include "hdf5.h"
#include <cstdlib>
#include <cstdio>
#include "rkfunction.hpp"
//#include "mpi.hpp"

using namespace std;

fstWriter::fstWriter(struct inputConfig cf, rk_func *f){
    xdp = (double*)malloc(cf.ni*cf.nj*cf.nk*sizeof(double));
    ydp = (double*)malloc(cf.ni*cf.nj*cf.nk*sizeof(double));
    zdp = (double*)malloc(cf.ni*cf.nj*cf.nk*sizeof(double));

    xsp = (float*)malloc(cf.ni*cf.nj*cf.nk*sizeof(float));
    ysp = (float*)malloc(cf.ni*cf.nj*cf.nk*sizeof(float));
    zsp = (float*)malloc(cf.ni*cf.nj*cf.nk*sizeof(float));

    vsp = (float*)malloc(cf.nci*cf.ncj*cf.nck*sizeof(float));
    v = (double*)malloc(cf.nci*cf.ncj*cf.nck*sizeof(double));

    readV = (double *)malloc(cf.nci*cf.ncj*cf.nck*cf.nv*sizeof(double));
    
    gridH = Kokkos::create_mirror_view(f->grid);
    varH = Kokkos::create_mirror_view(f->var);
}

void fstWriter::writeSPGrid(struct inputConfig cf, const FS4D gridD, const char * fname){
}
void fstWriter::writeGrid(struct inputConfig cf, const FS4D gridD ,const char * fname){
}

void fstWriter::writeSolution(struct inputConfig cf, rk_func *f, int tdx, double time){

    if (cf.ndim == 2){
        stringstream filenameStream;
        filenameStream << "solution-" << setw(7) << setfill('0') << tdx << ".h5";
        string filename = filenameStream.str();

        /*
         * HDF5 APIs definitions
         */     
        hid_t       file_id, group_id, dset_id;   // file, group and dataset identifiers
        hid_t       filespace, memspace;          // file and memory dataspace identifiers
        hsize_t     dimsg[2];                     // dataset dimensions
        hsize_t     dims[2];                      // chunk dimensions
        hsize_t count[2];                         // hyperslab selection parameters
        hsize_t stride[2];
        hsize_t block[2];
        hsize_t offset[2];
        hid_t   plist_id;                         // property list identifier
        herr_t  status;

        // MPI variables
        int mpi_size, mpi_rank;
        MPI_Comm comm  = cf.comm;
        MPI_Info info  = MPI_INFO_NULL;

        // Set up file access property list with parallel I/O access
        plist_id = H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_fapl_mpio(plist_id, cf.comm, info);

        // Create a new file collectively and release property list identifier.
        file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
        H5Pclose(plist_id);


        if (cf.rank == 0){
            cout << c(0,YEL) << left << setw(22) << "    Writing Solution: " << c(0,NON)
                 << c(0,CYA) << left << "'" + filename + "'" << c(0,NON) << endl;
        }

        Kokkos::deep_copy(gridH,f->grid);

        int idx;
        group_id = H5Gcreate(file_id, "/Grid", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        // Cycle through dimensions to write grid
        for (int vn=0; vn<2; ++vn){
            stringstream vnamestream;
            vnamestream << "Dimension" << setw(1) << setfill('0') << vn;
            string vname = vnamestream.str();
            for (int i=0; i<cf.ni; ++i){
                for (int j=0; j<cf.nj; ++j){
                    //for (int k=0; k<cf.nk; ++k){
                        idx = cf.ni*j+i;
                        xsp[idx] = (float)gridH(i,j,0,vn);
                        //ysp[idx] = gridH(i,j,k,1);
                        //zsp[idx] = gridH(i,j,k,2);
                    //}
                }
            }

            dimsg[1] = cf.glbl_ni;
            dimsg[0] = cf.glbl_nj;
            dims[1] = cf.ni;   
            dims[0] = cf.nj;   

            filespace = H5Screate_simple(cf.ndim, dimsg, NULL); 
            memspace  = H5Screate_simple(cf.ndim, dims, NULL); 


            // Create Dataset
            dset_id = H5Dcreate(group_id, vname.c_str(), H5T_NATIVE_FLOAT, filespace,
                    H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            H5Sclose(filespace);

            /* 
             * Each process defines dataset in memory and writes it to the hyperslab
             * in the file.
             */
            offset[1] = cf.iStart;
            offset[0] = cf.jStart;

            // Select hyperslab in the file.
            filespace = H5Dget_space(dset_id);
            status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, dims, NULL);

            // Create property list for collective dataset write.
            plist_id = H5Pcreate(H5P_DATASET_XFER);
            H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
            
            status = H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace, filespace,
                      plist_id, xsp);

            // Close/release resources.
            H5Dclose(dset_id);
            H5Sclose(filespace);
            H5Sclose(memspace);
            H5Pclose(plist_id);
        }
        H5Gclose(group_id);

        // Create group for solution datasets
        group_id = H5Gcreate(file_id, "/Solution", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        Kokkos::deep_copy(varH,f->var);
        for (int vn=0; vn<cf.nvt; ++vn){
            // Format Dataset Name
            stringstream vnamestream;
            vnamestream << "Variable" << setw(2) << setfill('0') << vn;
            string vname = vnamestream.str();

            if (cf.ndim == 3){
                for (int k=cf.ng; k<cf.nck+cf.ng; ++k){
                    for (int j=cf.ng; j<cf.ncj+cf.ng; ++j){
                        for (int i=cf.ng; i<cf.nci+cf.ng; ++i){
                            int ii = i - cf.ng;
                            int jj = j - cf.ng;
                            int kk = k - cf.ng;
                            idx = (cf.nci*cf.ncj)*kk+cf.nci*jj+ii;
                            vsp[idx] = varH(i,j,k,vn);
                        }
                    }
                }
            }else{
                for (int j=cf.ng; j<cf.ncj+cf.ng; ++j){
                    for (int i=cf.ng; i<cf.nci+cf.ng; ++i){
                        int ii = i - cf.ng;
                        int jj = j - cf.ng;
                        idx = cf.nci*jj+ii;
                        vsp[idx] = varH(i,j,0,vn);
                    }
                }
            }

            // Create the dataspace for the dataset.
            dimsg[1] = cf.glbl_nci;
            dimsg[0] = cf.glbl_ncj;
            dims[1] = cf.nci;   
            dims[0] = cf.ncj;   
            filespace = H5Screate_simple(cf.ndim, dimsg, NULL); 
            memspace  = H5Screate_simple(cf.ndim, dims, NULL); 

            // Create Dataset
            dset_id = H5Dcreate(group_id, vname.c_str(), H5T_NATIVE_FLOAT, filespace,
                    H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            H5Sclose(filespace);

            /* 
             * Each process defines dataset in memory and writes it to the hyperslab
             * in the file.
             */
            count[1] = cf.nci;
            count[0] = cf.ncj;
            offset[1] = cf.iStart;
            offset[0] = cf.jStart;

            // Select hyperslab in the file.
            filespace = H5Dget_space(dset_id);
            status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);

            // Create property list for collective dataset write.
            plist_id = H5Pcreate(H5P_DATASET_XFER);
            H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
            
            status = H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace, filespace,
                      plist_id, vsp);

            // Close/release resources.
            H5Dclose(dset_id);
            H5Sclose(filespace);
            H5Sclose(memspace);
            H5Pclose(plist_id);
        }
        H5Gclose(group_id);
        H5Fclose(file_id);

        // Write XDMF Descriptor file 
        if (cf.rank == 0){
            stringstream xilenameStream;
            xilenameStream << "solution-" << setw(7) << setfill('0') << tdx << ".xmf";
            string xilename = xilenameStream.str();

            int NX = cf.glbl_nci;
            int NY = cf.glbl_ncj;
            FILE *xmf = 0;
            xmf = fopen(xilename.c_str(), "w");
            fprintf(xmf, "<?xml version=\"1.0\" ?>\n");
            fprintf(xmf, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
            fprintf(xmf, "<Xdmf Version=\"3.0\">\n");
            fprintf(xmf, " <Domain>\n");
            fprintf(xmf, "   <Grid Name=\"mesh1\" GridType=\"Uniform\">\n");
            fprintf(xmf, "     <Topology TopologyType=\"2DSMesh\" NumberOfElements=\"%d %d\"/>\n", NY+1, NX+1);
            fprintf(xmf, "     <Geometry GeometryType=\"X_Y\">\n");
            fprintf(xmf, "       <DataItem Dimensions=\"%d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", (NY+1), (NX+1));
            fprintf(xmf, "        %s:/Grid/Dimension0\n",filename.c_str());
            fprintf(xmf, "       </DataItem>\n");
            fprintf(xmf, "       <DataItem Dimensions=\"%d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", (NY+1), (NX+1));
            fprintf(xmf, "        %s:/Grid/Dimension1\n",filename.c_str());
            fprintf(xmf, "       </DataItem>\n");
            fprintf(xmf, "     </Geometry>\n");
            fprintf(xmf, "     <Attribute Name=\"Momentum\" AttributeType=\"Vector\" Center=\"Cell\">\n");
            fprintf(xmf, "      <DataItem Dimensions=\"%d %d 2\" Function=\"JOIN($0,$1)\" ItemType=\"Function\">\n", NY, NX);
            fprintf(xmf, "       <DataItem Dimensions=\"%d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", NY, NX);
            fprintf(xmf, "        %s:/Solution/Variable00\n",filename.c_str());
            fprintf(xmf, "       </DataItem>\n");
            fprintf(xmf, "       <DataItem Dimensions=\"%d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", NY, NX);
            fprintf(xmf, "        %s:/Solution/Variable01\n",filename.c_str());
            fprintf(xmf, "       </DataItem>\n");
            fprintf(xmf, "      </DataItem>\n");
            fprintf(xmf, "     </Attribute>\n");
            fprintf(xmf, "     <Attribute Name=\"Energy\" AttributeType=\"Scalar\" Center=\"Cell\">\n");
            fprintf(xmf, "       <DataItem Dimensions=\"%d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", NY, NX);
            fprintf(xmf, "        %s:/Solution/Variable02\n",filename.c_str());
            fprintf(xmf, "       </DataItem>\n");
            fprintf(xmf, "     </Attribute>\n");
            for (int var=3; var<cf.nv; ++var){
                fprintf(xmf, "     <Attribute Name=\"Density%02d\" AttributeType=\"Scalar\" Center=\"Cell\">\n",var-2);
                fprintf(xmf, "       <DataItem Dimensions=\"%d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", NY, NX);
                fprintf(xmf, "        %s:/Solution/Variable%02d\n",filename.c_str(),var);
                fprintf(xmf, "       </DataItem>\n");
                fprintf(xmf, "     </Attribute>\n");
            }
            for (int var=cf.nv; var<cf.nvt; ++var){
                fprintf(xmf, "     <Attribute Name=\"Variable%02d\" AttributeType=\"Scalar\" Center=\"Cell\">\n",var-cf.nv);
                fprintf(xmf, "       <DataItem Dimensions=\"%d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", NY, NX);
                fprintf(xmf, "        %s:/Solution/Variable%02d\n",filename.c_str(),var);
                fprintf(xmf, "       </DataItem>\n");
                fprintf(xmf, "     </Attribute>\n");
            }
            fprintf(xmf, "   </Grid>\n");
            fprintf(xmf, " </Domain>\n");
            fprintf(xmf, "</Xdmf>\n");
            fclose(xmf);
        }
    }else{
        stringstream filenameStream;
        filenameStream << "solution-" << setw(7) << setfill('0') << tdx << ".h5";
        string filename = filenameStream.str();

        /*
         * HDF5 APIs definitions
         */     
        hid_t       file_id, group_id, dset_id;   // file, group and dataset identifiers
        hid_t       filespace, memspace;          // file and memory dataspace identifiers
        hsize_t     dimsg[3];                     // dataset dimensions
        hsize_t     dims[3];                      // chunk dimensions
        hsize_t count[3];                         // hyperslab selection parameters
        hsize_t stride[3];
        hsize_t block[3];
        hsize_t offset[3];
        hid_t   plist_id;                         // property list identifier
        herr_t  status;

        // MPI variables
        int mpi_size, mpi_rank;
        MPI_Comm comm  = cf.comm;
        MPI_Info info  = MPI_INFO_NULL;

        // Set up file access property list with parallel I/O access
        plist_id = H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_fapl_mpio(plist_id, cf.comm, info);

        // Create a new file collectively and release property list identifier.
        file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
        H5Pclose(plist_id);


        if (cf.rank == 0){
            cout << c(0,YEL) << left << setw(22) << "    Writing Solution: " << c(0,NON)
                 << c(0,CYA) << left << "'" + filename + "'" << c(0,NON) << endl;
        }

        Kokkos::deep_copy(gridH,f->grid);

        int idx;
        group_id = H5Gcreate(file_id, "/Grid", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        // Cycle through dimensions to write grid
        for (int vn=0; vn<3; ++vn){
            stringstream vnamestream;
            vnamestream << "Dimension" << setw(1) << setfill('0') << vn;
            string vname = vnamestream.str();
            for (int i=0; i<cf.ni; ++i){
                for (int j=0; j<cf.nj; ++j){
                    for (int k=0; k<cf.nk; ++k){
                        idx = (cf.ni*cf.nj)*k+cf.ni*j+i;
                        xsp[idx] = (float)gridH(i,j,k,vn);
                    }
                }
            }

            dimsg[2] = cf.glbl_ni;
            dimsg[1] = cf.glbl_nj;
            dimsg[0] = cf.glbl_nk;
            dims[2] = cf.ni;   
            dims[1] = cf.nj;   
            dims[0] = cf.nk;   

            filespace = H5Screate_simple(cf.ndim, dimsg, NULL); 
            memspace  = H5Screate_simple(cf.ndim, dims, NULL); 


            // Create Dataset
            dset_id = H5Dcreate(group_id, vname.c_str(), H5T_NATIVE_FLOAT, filespace,
                    H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            H5Sclose(filespace);

            /* 
             * Each process defines dataset in memory and writes it to the hyperslab
             * in the file.
             */
            offset[2] = cf.iStart;
            offset[1] = cf.jStart;
            offset[0] = cf.kStart;

            // Select hyperslab in the file.
            filespace = H5Dget_space(dset_id);
            status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, dims, NULL);

            // Create property list for collective dataset write.
            plist_id = H5Pcreate(H5P_DATASET_XFER);
            H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
            
            status = H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace, filespace,
                      plist_id, xsp);

            // Close/release resources.
            H5Dclose(dset_id);
            H5Sclose(filespace);
            H5Sclose(memspace);
            H5Pclose(plist_id);
        }
        H5Gclose(group_id);

        // Create group for solution datasets
        group_id = H5Gcreate(file_id, "/Solution", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        Kokkos::deep_copy(varH,f->var);
        for (int vn=0; vn<cf.nvt; ++vn){
            // Format Dataset Name
            stringstream vnamestream;
            vnamestream << "Variable" << setw(2) << setfill('0') << vn;
            string vname = vnamestream.str();

            for (int k=cf.ng; k<cf.nck+cf.ng; ++k){
                for (int j=cf.ng; j<cf.ncj+cf.ng; ++j){
                    for (int i=cf.ng; i<cf.nci+cf.ng; ++i){
                        int ii = i - cf.ng;
                        int jj = j - cf.ng;
                        int kk = k - cf.ng;
                        idx = (cf.nci*cf.ncj)*kk+cf.nci*jj+ii;
                        vsp[idx] = varH(i,j,k,vn);
                    }
                }
            }

            // Create the dataspace for the dataset.
            dimsg[2] = cf.glbl_nci;
            dimsg[1] = cf.glbl_ncj;
            dimsg[0] = cf.glbl_nck;
            dims[2] = cf.nci;   
            dims[1] = cf.ncj;   
            dims[0] = cf.nck;   
            filespace = H5Screate_simple(cf.ndim, dimsg, NULL); 
            memspace  = H5Screate_simple(cf.ndim, dims, NULL); 

            // Create Dataset
            dset_id = H5Dcreate(group_id, vname.c_str(), H5T_NATIVE_FLOAT, filespace,
                    H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            H5Sclose(filespace);

            /* 
             * Each process defines dataset in memory and writes it to the hyperslab
             * in the file.
             */
            count[2] = cf.nci;
            count[1] = cf.ncj;
            count[0] = cf.nck;
            offset[2] = cf.iStart;
            offset[1] = cf.jStart;
            offset[0] = cf.kStart;

            // Select hyperslab in the file.
            filespace = H5Dget_space(dset_id);
            status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);

            // Create property list for collective dataset write.
            plist_id = H5Pcreate(H5P_DATASET_XFER);
            H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
            
            status = H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace, filespace,
                      plist_id, vsp);

            // Close/release resources.
            H5Dclose(dset_id);
            H5Sclose(filespace);
            H5Sclose(memspace);
            H5Pclose(plist_id);
        }
        H5Gclose(group_id);
        H5Fclose(file_id);

        // Write XDMF Descriptor file 
        if (cf.rank == 0){
            stringstream xilenameStream;
            xilenameStream << "solution-" << setw(7) << setfill('0') << tdx << ".xmf";
            string xilename = xilenameStream.str();

            int NX = cf.glbl_nci;
            int NY = cf.glbl_ncj;
            int NZ = cf.glbl_nck;
            FILE *xmf = 0;
            xmf = fopen(xilename.c_str(), "w");
            fprintf(xmf, "<?xml version=\"1.0\" ?>\n");
            fprintf(xmf, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
            fprintf(xmf, "<Xdmf Version=\"3.0\">\n");
            fprintf(xmf, " <Domain>\n");
            fprintf(xmf, "   <Grid Name=\"mesh1\" GridType=\"Uniform\">\n");
            fprintf(xmf, "     <Topology TopologyType=\"3DSMesh\" NumberOfElements=\"%d %d %d\"/>\n", NZ+1, NY+1, NX+1);
            fprintf(xmf, "     <Geometry GeometryType=\"X_Y_Z\">\n");
            fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", (NZ+1), (NY+1), (NX+1));
            fprintf(xmf, "        %s:/Grid/Dimension0\n",filename.c_str());
            fprintf(xmf, "       </DataItem>\n");
            fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", (NZ+1), (NY+1), (NX+1));
            fprintf(xmf, "        %s:/Grid/Dimension1\n",filename.c_str());
            fprintf(xmf, "       </DataItem>\n");
            fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", (NZ+1), (NY+1), (NX+1));
            fprintf(xmf, "        %s:/Grid/Dimension2\n",filename.c_str());
            fprintf(xmf, "       </DataItem>\n");
            fprintf(xmf, "     </Geometry>\n");
            fprintf(xmf, "     <Attribute Name=\"Momentum\" AttributeType=\"Vector\" Center=\"Cell\">\n");
            fprintf(xmf, "      <DataItem Dimensions=\"%d %d %d 3\" Function=\"JOIN($0,$1,$2)\" ItemType=\"Function\">\n", NZ, NY, NX);
            fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", NZ, NY, NX);
            fprintf(xmf, "        %s:/Solution/Variable00\n",filename.c_str());
            fprintf(xmf, "       </DataItem>\n");
            fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", NZ, NY, NX);
            fprintf(xmf, "        %s:/Solution/Variable01\n",filename.c_str());
            fprintf(xmf, "       </DataItem>\n");
            fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", NZ, NY, NX);
            fprintf(xmf, "        %s:/Solution/Variable02\n",filename.c_str());
            fprintf(xmf, "       </DataItem>\n");
            fprintf(xmf, "      </DataItem>\n");
            fprintf(xmf, "     </Attribute>\n");
            fprintf(xmf, "     <Attribute Name=\"Energy\" AttributeType=\"Scalar\" Center=\"Cell\">\n");
            fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", NZ, NY, NX);
            fprintf(xmf, "        %s:/Solution/Variable03\n",filename.c_str());
            fprintf(xmf, "       </DataItem>\n");
            fprintf(xmf, "     </Attribute>\n");
            for (int var=4; var<cf.nv; ++var){
                fprintf(xmf, "     <Attribute Name=\"Density%02d\" AttributeType=\"Scalar\" Center=\"Cell\">\n",var-3);
                fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", NZ, NY, NX);
                fprintf(xmf, "        %s:/Solution/Variable%02d\n",filename.c_str(),var);
                fprintf(xmf, "       </DataItem>\n");
                fprintf(xmf, "     </Attribute>\n");
            }
            fprintf(xmf, "   </Grid>\n");
            fprintf(xmf, " </Domain>\n");
            fprintf(xmf, "</Xdmf>\n");
            fclose(xmf);
        }
        
    }

}
void fstWriter::writeRestart(struct inputConfig cf, rk_func *f, int tdx, double time){
}


void fstWriter::readSolution(struct inputConfig cf, FS4D &gridD, FS4D &varD){
}
