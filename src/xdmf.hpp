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

#ifndef XDMF_H
#define XDMF_H

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <map>
#include <vector>
#include "kokkosTypes.hpp"

using namespace std;
void writeXMFDataItem(FILE*, string, int, vector<int> &);
void writeXMF(string, string, int, FSCAL, int, size_t*, vector<FSCAL>, vector<FSCAL>, int, bool, vector<string>, vector<string>, std::map<std::string,FS4D>&);

#endif
