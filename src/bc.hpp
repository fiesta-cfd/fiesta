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
#ifndef BC_H
#define BC_H

#include "rkfunction.hpp"
#include "Kokkos_Core.hpp"
#include "kokkosTypes.hpp"
#include "input.hpp"
#include <string>
#ifdef HAVE_MPI
#include "mpi.hpp"
#endif

enum class BCType {outflow,reflective,noslip,hydrostatic};
/* void applyBCs(struct inputConfig cf, std::unique_ptr<class rk_func>&f); */
void applyBCs(struct inputConfig cf, class rk_func *f);
BCType parseBC(std::string name);

#endif
