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

#include "vtk.hpp"
#include "kokkosTypes.hpp"
//#include "Kokkos_Core.hpp"
#include "debug.hpp"
#include "output.hpp"
#include "rkfunction.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

using namespace std;

serialVTKWriter::serialVTKWriter(struct inputConfig cf, FS4D gridD, FS4D varD) {

  gridH = Kokkos::create_mirror_view(gridD);
  varH = Kokkos::create_mirror_view(varD);
}

// struct inputConfig writeGrid(struct inputConfig cf, FSCAL *x, FSCAL *y,
// FSCAL *z,const char * fname){
void serialVTKWriter::writeGrid(struct inputConfig cf, const FS4D gridD,
                                const char *fname) {}
// struct inputConfig writeSPGrid(struct inputConfig cf, float *x, float *y,
// float *z, const char * fname){
void serialVTKWriter::writeSPGrid(struct inputConfig cf, const FS4D gridD,
                                  const char *fname) {}

void serialVTKWriter::writeRestart(struct inputConfig cf, std::unique_ptr<class rk_func>&f, int tdx,
                                   FSCAL time) {}

// void writeSolution(struct inputConfig cf, float *x, float *y, float *z, const
// FS4D deviceV, int tdx, FSCAL time){
void serialVTKWriter::writeSolution(struct inputConfig cf, std::unique_ptr<class rk_func>&mod,
                                    int tdx, FSCAL time) {

  Kokkos::deep_copy(varH, mod->var);
  Kokkos::deep_copy(gridH, mod->grid);

  stringstream ss;
  ss cv.pathName << "/" << "solution-" << setw(7) << setfill('0') << tdx << ".vtk";
  string fsname = ss.str();

  if (cf.rank == 0) {
    cout << c(0, YEL) << left << setw(22)
         << "    Writing Solution: " << c(0, NON) << c(0, CYA) << left
         << "'" + fsname + "'" << c(0, NON) << endl;
  }

  std::ofstream f;
  f.open(ss.str());
  f << "# vtk DataFile Version 4.2" << endl;
  f << cf.title.c_str() << endl;
  f << "BINARY\nDATASET STRUCTURED_GRID" << endl;
  // f << "ASCII\nDATASET STRUCTURED_GRID" << endl;
  // if (cf.ndim == 3){
  f << "DIMENSIONS " << cf.ni << " " << cf.nj << " " << cf.nk << endl;
  f << "POINTS " << cf.ni * cf.nj * cf.nk << " FSCAL" << endl;
  for (int k = 0; k < cf.nk; k++) {
    for (int j = 0; j < cf.nj; j++) {
      for (int i = 0; i < cf.ni; i++) {
        for (int v = 0; v < 3; v++) {
          char data[8], *pDouble = (char *)(FSCAL *)(&gridH(i, j, k, v));
          for (int b = 0; b < 8; ++b)
            data[b] = pDouble[7 - b];
          f.write(data, 8);
        }
        // f << gridH(i,j,k,0) << " " << gridH(i,j,k,1) << " " << gridH(i,j,k,2)
        // << endl;
      }
    }
  }
  //}else{
  //    f << "DIMENSIONS " << cf.ni << " " << cf.nj << endl;
  //    f << "POINTS " << cf.ni*cf.nj << " float" << endl;
  //    for (int j=0; j < cf.nj; j++) {
  //        for (int i=0; i < cf.ni; i++) {
  //            f << gridH(i,j,0,0) << " " << gridH(i,j,0,1) << endl;
  //        }
  //    }
  //}

  int ncells = cf.nci * cf.ncj * cf.nck;
  f << "CELL_DATA " << ncells << endl;
  f << "FIELD FieldData 1" << endl;
  f << "Variables " << cf.nv << " " << ncells << " FSCAL" << endl;
  for (int k = 0; k < cf.nck; k++) {
    for (int j = 0; j < cf.ncj; j++) {
      for (int i = 0; i < cf.nci; i++) {
        for (int v = 0; v < cf.nv; v++) {
          if (cf.ndim == 3) {
            char data[8], *pDouble = (char *)(FSCAL *)(&varH(
                              i + cf.ng, j + cf.ng, k + cf.ng, v));
            for (int b = 0; b < 8; ++b)
              data[b] = pDouble[7 - b];
            f.write(data, 8);
            // f.write((char*)&varH(i+cf.ng,j+cf.ng,k+cf.ng,v),sizeof(FSCAL));
            // f << setprecision(6) << scientific <<
            // varH(i+cf.ng,j+cf.ng,k+cf.ng,v) << " ";
          } else {
            char data[8],
                *pDouble =
                    (char *)(FSCAL *)(&varH(i + cf.ng, j + cf.ng, k, v));
            for (int b = 0; b < 8; ++b)
              data[b] = pDouble[7 - b];
            f.write(data, 8);
            // f.write((char*)&varH(i+cf.ng,j+cf.ng,k,v),sizeof(FSCAL));
            // f << setprecision(6) << scientific << varH(i+cf.ng,j+cf.ng,k,v)
            // << " ";
          }
        }
        // f << endl;
      }
    }
  }

  f.flush();
  f.close();

  // if (cf.numProcs == 1){
  //    std::ofstream myfile;
  //    myfile.open("output.txt");
  //    for (int i=cf.ng; i<cf.nci; ++i){
  //       int  idx = (cf.nci*cf.ncj)*0+cf.nci*3+i;
  //  //          myfile << x[idx] << ", " << hostV(i,3,0,0) << ", " <<
  //  hostV(i,3,0,1) << ", " << hostV(i,3,0,2) << ", " << hostV(i,3,0,3) <<
  //  std::endl;
  //        myfile << x[idx] << ", " << hostV(i,3,0,3) << std::endl;
  //    }
  //    myfile << std::endl;
  //    myfile.close();
  //}
}

void serialVTKWriter::readSolution(struct inputConfig cf, FS4D &gridD,
                                   FS4D &varD) {}
