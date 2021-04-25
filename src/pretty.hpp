

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

#ifndef PRETTY_HPP
#define PRETTY_HPP

#include "unistd.h"
#include <string>

using std::string;

enum Colour { red, green, yellow, blue, magenta, cyan, reset, none};

struct ansiColors {
  ansiColors(const int c_):colorMode(c_){};
  string operator()(const Colour color){
    if ( (colorMode == 1) || ((colorMode == 2) && isatty(fileno(stdout))) ){
      switch (color) {
      case red:
        return "\033[0;31m";
      case green:
        return "\033[0;32m";
      case yellow:
        return "\033[0;33m";
      case blue:
        return "\033[0;34m";
      case magenta:
        return "\033[0;35m";
      case cyan:
        return "\033[0;36m";
      case reset:
        return "\033[0m";
      case none:
        return "\033[0m";
      default:
        return "\033[0m";
      }
    }else{
      return "";
    }
  }

  private:
  const int colorMode;
};

#endif
