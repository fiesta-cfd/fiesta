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

#ifndef LOG2_HPP
#define LOG2_HPP

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
#include "fmt/ranges.h"
#include "pretty.hpp"

using fmt::format;

namespace Fiesta {
  namespace Log {

    extern int verbosity;
    extern ansiColors *c;
    extern int rank;
    extern Kokkos::Timer *timer;
    
    template<typename... Args> inline
    void log(Colour color, std::string type, std::string logformat, Args&&... args){
      if (rank==0){
        std::string format = fmt::format("{}[{: >12.5f}] {: >7}: {}{}\n",(*c)(color),timer->seconds(),type,logformat,(*c)(reset));
        fmt::vprint(format,fmt::make_args_checked<Args...>(format,args...));
        std::cout << std::flush;
      }
    }

    template<typename... Args> inline
    void logRank(int ra, Colour color, std::string type, std::string logformat, Args&&... args){
      if (rank==ra){
        std::string format = fmt::format("{}[{: >12.5f}] {: >7} <{}>(: {}{}\n",(*c)(color),timer->seconds(),type,rank,logformat,(*c)(reset));
        fmt::vprint(format,fmt::make_args_checked<Args...>(format,args...));
        std::cout << std::flush;
      }
    }
    
    template<typename... Args> inline
    void logAll(Colour color, std::string type, std::string logformat, Args&&... args){
      std::string format = fmt::format("{}[{: >12.5f}] {: >7} <{}>: {}{}\n",(*c)(color),timer->seconds(),type,rank,logformat,(*c)(reset));
      fmt::vprint(format,fmt::make_args_checked<Args...>(format,args...));
    }
    
    inline
    void print(std::string header, std::string message, Colour color){
      if (rank==0){
        fmt::print("{}[{: >12.5f}] {: >7}: {}{}\n",(*c)(color),timer->seconds(),header,message,(*c)(reset));
      }
    }
    
    inline
    string getTime(){
      auto currentTime = std::chrono::system_clock::now();
      time_t end_time = std::chrono::system_clock::to_time_t(currentTime);
      string message = ctime(&end_time);
      message.erase(remove(message.begin(), message.end(), '\n'), message.end());
      return message;
    }
    
    inline
    void Logger(int v, int c_, int r){
      timer = new Kokkos::Timer();
      c = new ansiColors(c_);
      verbosity=v;
      rank=r;
    
      print("Logger",format("Log started {}",getTime()),green);
      print("Logger",format("Log Level {}",verbosity),green);
    }
    
    template<typename... Args> inline
    void debugAll(string format, Args&&... args) {
      if (verbosity>=5) logAll(magenta,"Debug",format,args...);
    }
    
    template<typename... Args> inline
    void infoAll(string format, Args&&... args) {
      if (verbosity>=4) logAll(none,"Info",format,args...);
    }
    
    template<typename... Args> inline
    void debug(string format, Args&&... args) {
      if (verbosity>=5) log(magenta,"Debug",format,args...);
    }

    template<typename... Args> inline
    void debugRank(int ra, string format, Args&&... args) {
      if (verbosity>=5) logRank(ra,magenta,"Debug",format,args...);
    }
    
    template<typename... Args> inline
    void debugWarning(string format, Args&&... args) {
      if (verbosity>=5) log(yellow,"Debug",format,args...);
    }

    template<typename... Args> inline
    void info(string format, Args&&... args) {
      if (verbosity>=4) log(none,"Info",format,args...);
    }

    template<typename... Args> inline
    void infoWarning(string format, Args&&... args) {
      if (verbosity>=4) log(yellow,"Info",format,args...);
    }
    
    template<typename... Args> inline
    void message(string format, Args&&... args) {
      if (verbosity>=3) log(blue,"Message",format,args...);
    }
    
    template<typename... Args> inline
    void warning(string format, Args&&... args) {
      if (verbosity>=2) log(yellow,"Warning",format,args...);
    }
    
    template<typename... Args> inline
    void error(string format, Args&&... args) {
      if (verbosity>=1) log(red,"Error",format,args...);
    }
  }
}
#endif
