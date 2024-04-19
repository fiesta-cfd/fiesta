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

#include <iostream>
#include <chrono>
#include <string>
#include "fmt/core.h"
#include "fmt/ranges.h"

using fmt::format;

namespace Log {

  inline int& verbosity(){
    static int verb = 5;
    return verb;
  }

  inline int& rank(){
    static int rk = 0;
    return rk;
  }

  inline bool& colorable(){
    static bool co = false;
    return co;
  }

  inline double timer(){
    static std::chrono::high_resolution_clock::time_point m_old = std::chrono::high_resolution_clock::now();
    std::chrono::high_resolution_clock::time_point m_new = std::chrono::high_resolution_clock::now();
    return std::chrono::duration_cast<std::chrono::duration<double>>(m_new - m_old).count();
  }

  enum Colors { red, green, yellow, blue, magenta, cyan, reset, none};
  inline std::string getColor(const Colors color){
    if (colorable()){
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

  template<typename... Args> inline
  void log(Colors color, std::string type, std::string logformat, Args&&... args){
    if (rank()==0){
      std::string format = fmt::format("{}[{: >12.5f}] {: >7}: {}{}\n",getColor(color),timer(),type,logformat,getColor(reset));
      fmt::vprint(format,fmt::make_args_checked<Args...>(format,args...));
      std::cout << std::flush;
    }
  }

  template<typename... Args> inline
  void logrank(int ra, Colors color, std::string type, std::string logformat, Args&&... args){
    if (rank()==ra){
      std::string format = fmt::format("{}[{: >12.5f}] {: >7} <{}>(: {}{}\n",getColor(color),timer(),type,rank(),logformat,getColor(reset));
      fmt::vprint(format,fmt::make_args_checked<Args...>(format,args...));
      std::cout << std::flush;
    }
  }
  
  template<typename... Args> inline
  void logAll(Colors color, std::string type, std::string logformat, Args&&... args){
    std::string format = fmt::format("{}[{: >12.5f}] {: >7} <{}>: {}{}\n",getColor(color),timer(),type,rank(),logformat,getColor(reset));
    fmt::vprint(format,fmt::make_args_checked<Args...>(format,args...));
  }
  
  inline
  void print(std::string header, std::string message, Colors color){
    if (rank()==0){
      fmt::print("{}[{: >12.5f}] {: >7}: {}{}\n",getColor(color),timer(),header,message,getColor(reset));
    }
  }
  
  inline
  std::string getTime(){
    auto currentTime = std::chrono::system_clock::now();
    time_t end_time = std::chrono::system_clock::to_time_t(currentTime);
    std::string message = ctime(&end_time);
    message.erase(remove(message.begin(), message.end(), '\n'), message.end());
    return message;
  }
  
  inline
  void Logger(int v, int c_, int r){
    colorable() = bool(c_);
    verbosity()=v;
    rank()=r;
  
    print("Logger",format("Log started {}",getTime()),green);
    print("Logger",format("Log Level {}",verbosity()),green);
  }
  
  template<typename... Args> inline
  void debugAll(std::string format, Args&&... args) {
    if (verbosity()>=5) logAll(magenta,"Debug",format,args...);
  }
  
  template<typename... Args> inline
  void infoAll(std::string format, Args&&... args) {
    if (verbosity()>=4) logAll(none,"Info",format,args...);
  }
  
  template<typename... Args> inline
  void debug(std::string format, Args&&... args) {
    if (verbosity()>=5) log(magenta,"Debug",format,args...);
  }

  template<typename... Args> inline
  void debugrank(int ra, std::string format, Args&&... args) {
    if (verbosity()>=5) logrank(ra,magenta,"Debug",format,args...);
  }
  
  template<typename... Args> inline
  void debugWarning(std::string format, Args&&... args) {
    if (verbosity()>=5) log(yellow,"Debug",format,args...);
  }

  template<typename... Args> inline
  void info(std::string format, Args&&... args) {
    if (verbosity()>=4) log(none,"Info",format,args...);
  }

  template<typename... Args> inline
  void infoWarning(std::string format, Args&&... args) {
    if (verbosity()>=4) log(yellow,"Info",format,args...);
  }
  
  template<typename... Args> inline
  void message(std::string format, Args&&... args) {
    if (verbosity()>=3) log(blue,"Message",format,args...);
  }
  
  template<typename... Args> inline
  void warning(std::string format, Args&&... args) {
    if (verbosity()>=2) log(yellow,"Warning",format,args...);
  }
  
  template<typename... Args> inline
  void error(std::string format, Args&&... args) {
    if (verbosity()>=1) log(red,"Error",format,args...);
  }
}
#endif
