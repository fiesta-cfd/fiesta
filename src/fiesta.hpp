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

#ifndef FIESTA2_HPP
#define FIESTA2_HPP

#include "rkfunction.hpp"
#include "writer.hpp"
#include "input.hpp"
#include <iostream>
#include <vector>
#include "block.hpp"
#include <memory>

#define FIESTA_RESTART_VERSION 2

namespace Fiesta {
    struct Simulation {
      std::unique_ptr<class rk_func> f;
      std::unique_ptr<blockWriter<FSCAL>> restartview;
      std::vector<blockWriter<float>> ioviews;
      fsconf cf;
    };

    struct inputConfig initialize(struct inputConfig&,int, char **);
    void initializeSimulation(struct inputConfig&, std::unique_ptr<class rk_func>&);
    void initializeSimulation(Simulation &sim);
    void reportTimers(struct inputConfig&, std::unique_ptr<class rk_func>&);
    void checkIO(Simulation &sim, size_t t);
    //void finalize(struct inputConfig &);
    void step(Simulation &sim, size_t t);
    void collectSignals(struct inputConfig &cf);
    void finalize();
    int fiestaTest(int,int);

};

#endif
