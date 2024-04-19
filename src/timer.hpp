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

#include <string>
#include <chrono>
#include <fmt/core.h>

namespace Timer{
  inline std::string format(double time){
    int d,h,m,s;
    double tmp;
    std::string ftime;

    tmp = time;
    d = time/86400;

    tmp = tmp - d*86400;
    h = tmp/3600;

    tmp = tmp - h*3600;
    m = tmp/60;

    s = tmp - m*60;

    ftime = fmt::format("{}s",s);
    if (m>0) ftime = fmt::format("{}m{}",m,ftime);
    if (h>0) ftime = fmt::format("{}h{}",h,ftime);
    if (d>0) ftime = fmt::format("{}d{}",d,ftime);

    return ftime;
  }

  class fiestaTimer {

  public:
    // Constructors
    fiestaTimer(){
      time = 0.0;
      m_old = std::chrono::high_resolution_clock::now();
      description = "No Description";
    }

    fiestaTimer(std::string desc_):description(desc_){
      time = 0.0;
      m_old = std::chrono::high_resolution_clock::now();
    }

    // Operators
    void reset(){
      m_old = std::chrono::high_resolution_clock::now();
    }

    void start(){
      m_old = std::chrono::high_resolution_clock::now();
      time = 0.0;
    }

    void accumulate(){
      m_new = std::chrono::high_resolution_clock::now();
      time += std::chrono::duration_cast<std::chrono::duration<double>>(m_new-m_old).count();
    }

    void stop(){
      m_new = std::chrono::high_resolution_clock::now();
      time += std::chrono::duration_cast<std::chrono::duration<double>>(m_new-m_old).count();
    }

    void clear(){
      time = 0.0;
    }

    // Retrievers
    double get(){
      return time;
    }

    std::string getf(){
      return format(time);
    }

    std::string checkf(){
      m_new = std::chrono::high_resolution_clock::now();
      return format(std::chrono::duration_cast<std::chrono::duration<double>>(m_new-m_old).count());
    }

    double check(){
      m_new = std::chrono::high_resolution_clock::now();
      return std::chrono::duration_cast<std::chrono::duration<double>>(m_new-m_old).count();
    }

    std::string describe(){
      return description;
    }

  private:
    double time;
    std::chrono::high_resolution_clock::time_point m_old;
    std::chrono::high_resolution_clock::time_point m_new;
    std::string description;
  };
}
#endif
