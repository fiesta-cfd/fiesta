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
#include <string>
#include <map>

#ifndef DIAGNOSTICS_HPP
#define DIAGNOSTICS_HPP

class Diagnostics {
  public:
    Diagnostics(size_t ng_, size_t ni_, size_t nj_, size_t nk_, size_t nv_, int freq_);
    Diagnostics(Diagnostics&& d1);
    Diagnostics();
    Diagnostics& operator=(const Diagnostics& d1);
    Diagnostics& operator=(Diagnostics&& d1);
    void start(FS4D dvar);
    void check(std::string name, size_t t, FS4D dvar, std::map<std::string,FS4D> &dgvar);
    void stop(FS4D dvar);
  private:
    size_t ng,ni,nj,nk,nv,freq;
    FS4D diag;
};

#endif
