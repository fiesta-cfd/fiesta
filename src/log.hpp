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

#ifndef LOG_HPP
#define LOG_HPP

#include "Kokkos_Core.hpp"
#include <iostream>
#include <string>
#include <iomanip>
#include <locale>
#include "unistd.h"
#include <sstream>
#include <fstream>
#include <ctime>
#include <chrono>
#include <algorithm>
#include "fmt/core.h"
#include "pretty.hpp"

using namespace std;
using fmt::format;

class Logger {
  public:
    Logger(int v_, int c_, int r_) : verbosity(v_), colorMode(c_), rank(r_) {
      timer = new Kokkos::Timer();
      print("Logger",format("Log started {}",getTime()),green);
      print("Logger",format("Log Level {}",verbosity),green);
      if (colorMode) print("Logger","Color logs enabled",green);
    }

    ~Logger() {
      print("Logger",format("Log ended {}",getTime()),green);
    }

    template<typename... Args>
    void debug(Args... args) {
      if (verbosity >= 5){
        stringstream message;
        (message << ... << args) << "";
        print("Debug",message.str(),magenta);
      }
    }
    
    template<typename... Args>
    void info(Args... args) {
      if (verbosity >= 4){
        stringstream message;
        (message << ... << args) << "";
        print("Info",message.str(),none);
      }
    }
    
    template<typename... Args>
    void message(Args... args) {
      if (verbosity >= 3){
        stringstream message;
        (message << ... << args) << "";
        print("Message",message.str(),blue);
      }
    }
    
    template<typename... Args>
    void warning(Args... args) {
      if (verbosity >= 2){
        stringstream message;
        (message << ... << args) << "";
        print("WARNING",message.str(),yellow);
      }
    }
    
    template<typename... Args>
    void error(Args... args) {
      if (verbosity >= 1){
        stringstream message;
        (message << ... << args) << "";
        print("ERROR",message.str(),red);
      }
    }
    
    string getTime(){
      auto currentTime = chrono::system_clock::now();
      time_t end_time = chrono::system_clock::to_time_t(currentTime);
      string message = ctime(&end_time);
      message.erase(remove(message.begin(), message.end(), '\n'), message.end());
      return message;
    }
    
  private:
    void print(std::string header, std::string message, Colour color){
      using namespace std;
      using namespace fmt;
      ansiColors c(colorMode);

      if (rank==0){
        cout << format("{}[{: >12.5f}] {: >7}: {}{}",c(color),timer->seconds(),header,message,c(reset)) << endl;
      }
    }
    Kokkos::Timer *timer;
    const int verbosity;
    const int rank;
    const int colorMode;
};

#endif
