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

#ifndef TIMER_H
#define TIMER_H

#include "Kokkos_Core.hpp"
#include <string>

using namespace std;

class fiestaTimer {

public:
  fiestaTimer();
  fiestaTimer(int format);
  fiestaTimer(string description);
  string name;
  int format;
  void accumulate();
  double get();
  void clear();
  void reset();
  void start();
  void stop();
  string describe();
  string getf();
  string getf(int format);
  string checkf();
  string checkf(int format);
  double check();
  string formatTime(double tin);

private:
  double time;
  Kokkos::Timer *timer;
  string description;
};
#endif
