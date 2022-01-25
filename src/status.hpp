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

#ifndef STATUS_H
#define STATUS_H

#include "Kokkos_Core.hpp"
#include "kokkosTypes.hpp"
#include "timer.hpp"
#include "rkfunction.hpp"
#include <iomanip>
#include <iostream>

void statusCheck(int cFlag, struct inputConfig cf, std::unique_ptr<class rk_func>&f, FSCAL time,
                                         Timer::fiestaTimer &wall, Timer::fiestaTimer &sim);

#endif
