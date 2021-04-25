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
#ifndef BLOCK_H
#define BLOCK_H

#include <string>
#include "rkfunction.hpp"
#include "mpi.h"
#include "kokkosTypes.hpp"
#include "hdf5.h"

template <typename T>
class blockWriter {
  public:
    blockWriter();
    blockWriter(struct inputConfig&,rk_func*,std::string,std::string,bool,size_t);
    blockWriter(struct inputConfig&,rk_func*,std::string,std::string,bool,size_t,
                                  std::vector<size_t>,std::vector<size_t>,std::vector<size_t>);
    void write(struct inputConfig cf, rk_func *f, int tdx, double time);
    size_t frq();

  private:
      
    std::vector<size_t> lStart;  // local starting index
    std::vector<size_t> lEnd;    // local ending cell index
    std::vector<size_t> lEndG;   // local ending grid index
    std::vector<size_t> lExt;    // local extent
    std::vector<size_t> lExtG;   // local grid extent
    std::vector<size_t> lOffset; // local offset
    std::vector<size_t> gStart;  // local starting index
    std::vector<size_t> gEnd;    // local ending index
    std::vector<size_t> gExt;    // global extent
    std::vector<size_t> gExtG;   // global grid extent
    std::vector<size_t> stride; // slice stride

    size_t freq;    // block write frequency

    bool writeVarx;
    
    size_t lElems;
    size_t lElemsG;

    MPI_Comm sliceComm;
    int myColor;

    std::vector<T> varData;
    std::vector<T> gridData;

    FS4DH gridH;
    FS4DH varH;
    FS4DH varxH;

    bool slicePresent;
    int offsetDelta;

    std::string path;
    std::string name;
    int pad;
    bool avg;

    void write_h5(hid_t, std::string, int, std::vector<size_t>, std::vector<size_t>, std::vector<size_t>, std::vector<T>&);

    void dataPack(int, int, const std::vector<size_t>&, const std::vector<size_t>&, const std::vector<size_t>&,
                            const std::vector<size_t>&, std::vector<T>&, const FS4DH&, const int, const bool);
    
    hid_t openHDF5ForWrite(MPI_Comm, MPI_Info, std::string);
    void close_h5(hid_t);
};

template <typename H>
hid_t getH5Type();

#endif
