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

#include "kokkosTypes.hpp"
#ifndef VTK_H
#define VTK_H

#include "input.hpp"
#include "rkfunction.hpp"
#include "writer.hpp"

class serialVTKWriter : public writer {

public:
  serialVTKWriter(struct inputConfig, FS4D gridD, FS4D varD);
  void writeGrid(struct inputConfig cf, const FS4D gridD, const char *fname);
  void writeSPGrid(struct inputConfig cf, const FS4D gridD, const char *fname);

  void writeSolution(struct inputConfig cf, std::unique_ptr<class rk_func>&f, int tdx, FSCAL time);
  void writeRestart(struct inputConfig cf, std::unique_ptr<class rk_func>&f, int tdx, FSCAL time);

  void readSolution(struct inputConfig cf, FS4D &deviceG, FS4D &deviceV);

private:
  FS4DH gridH;
  FS4DH varH;
};

#endif
