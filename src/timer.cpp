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

#include "timer.hpp"
#include "Kokkos_Core.hpp"
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

using namespace std;

fiestaTimer::fiestaTimer() {
  timer = new Kokkos::Timer();
  timer->reset();
  time = 0.0;
  description = "No Description";
}

fiestaTimer::fiestaTimer(string n_) : description(n_) {
  timer = new Kokkos::Timer();
  timer->reset();
  time = 0.0;
}

void fiestaTimer::accumulate() { time = time + timer->seconds(); }
double fiestaTimer::get() { return time; }

string fiestaTimer::getf() {
  return formatTime(time);
}

string fiestaTimer::getf(int format) {
  if (format == 0){
    return formatTime(time);
  }else{
    stringstream ss;
    ss.precision(format);
    ss.setf(ios::scientific);
    ss << time;
    return ss.str();
  }
}
void fiestaTimer::clear() { time = 0.0; }
void fiestaTimer::reset() { timer->reset(); }
void fiestaTimer::start() {
  timer->reset();
  time = 0.0;
}
void fiestaTimer::stop() { time = timer->seconds(); }
string fiestaTimer::checkf() { return formatTime(timer->seconds()); }
string fiestaTimer::checkf(int format) {
  if (format == 0){
    return formatTime(time);
  }else{
    stringstream ss;
    ss.precision(format);
    ss.setf(ios::scientific);
    ss << time;
    return ss.str();
  }
}
double fiestaTimer::check() { return timer->seconds(); }
string fiestaTimer::formatTime(double tin) {
  int d;
  int h;
  int m;
  double s;
  double temp;

  d = 0;
  h = 0;
  m = 0;
  s = 0.0;

  temp = tin;
  d = time / 864000;

  temp = temp - d * 864000;
  h = temp / 3600;

  temp = temp - h * 3600;
  m = temp / 60;

  temp = temp - m * 60;
  s = floor(temp);

  stringstream ss;

  ss.precision(0);
  ss.setf(ios::fixed);
  if (d > 0)
    ss << d << "d";
  if (h > 0)
    ss << h << "h";
  if (m > 0)
    ss << m << "m";
  ss << s << "s";

  return ss.str();
}
string fiestaTimer::describe() { return description; }
