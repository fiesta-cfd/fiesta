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

#ifndef RKFUNCTION_H
#define RKFUNCTION_H

#include "Kokkos_Core.hpp"
#include "kokkosTypes.hpp"
#include "input.hpp"
#include "timer.hpp"
#ifdef HAVE_MPI
#include "mpi.hpp"
#endif
#include <map>
#include <string>
#include <vector>
#include "block.hpp"

class rk_func {

public:
  rk_func(struct inputConfig &cf_);
  //rk_func(struct inputConfig &cf_, Kokkos::View<FSCAL *> &cd_);
  virtual ~rk_func() = default;

  virtual void compute() = 0;
  virtual void preStep() = 0;
  virtual void postStep() = 0;
  virtual void preSim() = 0;
  virtual void postSim() = 0;
  // virtual void compute(const FS4D & mvar, FS4D & mdvar) = 0;

  FS4D var;
  std::vector<std::string> varNames;
  FS4D varx;
  std::vector<std::string> varxNames;
  FS4D dvar;
  FS4D tmp1;
  FS4D grid;

  std::map<std::string, Timer::fiestaTimer> timers;
  policy_f ghostPol = policy_f({0, 0}, {1, 1});
  policy_f cellPol = policy_f({0, 0}, {1, 1});
  policy_f facePol = policy_f({0, 0}, {1, 1});

protected:
  struct inputConfig &cf;
  FS1D cd;
};

#endif
