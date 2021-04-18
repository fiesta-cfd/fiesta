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

using namespace std;
class Logger {
  public:
    Logger(int v_, int c_, int cl_, int r_, string f_) : verbosity(v_), colorFlag(c_), colorLogs(cl_), rank(r_), logFilename(f_) {
      timer = new Kokkos::Timer();
      logFile.open(logFilename,ios::trunc);
    
      stringstream message;
      message << "Log started " << getTime();
      print("Logger",message.str(),green);

      stringstream info;
      info << "Log level: " << verbosity;
      print("Logger",info.str(),green);

      if (colorFlag){
        stringstream cinfo;
        cinfo << "Color logs enabled";
        print("Logger",cinfo.str(),green);
      }

      if (colorLogs){
        stringstream clinfo;
        clinfo << "Color log file enabled";
        print("Logger",clinfo.str(),green);
      }
    }

    ~Logger() {
      stringstream message;
      message << "Log ended " << getTime();
      print("Logger",message.str(),green);

      logFile.close();
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
    enum logColor { red, green, yellow, blue, magenta, cyan, none };

    void print(std::string header, std::string message, logColor color){
      using namespace std;

      if (rank==0){
        stringstream logString;
        logString
            << "[" << setw(16) << setprecision(7) << fixed << right << timer->seconds() << "] "
            << setw(8) << right << header << ": "
            << message;

        // cout << setColor(color)
        //     << logString.str()
        //     << setColor(none)
        //     << endl;

        if (colorLogs){
          logFile << setColor(color)
                  << logString.str()
                  << setColor(none)
                  << endl;
        }else{
          logFile << logString.str() << endl;
        }
      }
    }

    string setColor(logColor color) {
      int use_color = 0;

      if (colorFlag == 1)
        use_color = 1;
      if ((colorFlag == 2) && isatty(fileno(stdout)))
        use_color = 1;

      if (use_color) {
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
        case none:
          return "\033[0m";
        default:
          return "\033[0m";
        }
      } else {
        return "";
      }
    }

    ofstream logFile;
    Kokkos::Timer *timer;
    const int colorLogs;
    const int verbosity;
    const int rank;
    const int colorFlag;
    const std::string logFilename;
};

#endif




