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
#include "mpi.h"

template<typename T>
class h5Writer {

  public:
    h5Writer();
    void open(MPI_Comm comm, MPI_Info info, std::string fname);

    void close();

    void write(std::string, int, std::vector<size_t>, std::vector<size_t>, std::vector<size_t>, std::vector<T>&);

    void openGroup(std::string name);
    void closeGroup();

  private:
    std::string filename;
    MPI_Comm comm;
    hid_t file_id;
    hid_t group_id;
};
