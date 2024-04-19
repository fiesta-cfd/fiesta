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

#ifndef SIGNAL_H
#define SIGNAL_H

#include <csignal>
#include "input.hpp"

class fiestaSignalHandler{
    static fiestaSignalHandler *instance;
    struct inputConfig &cf;

    fiestaSignalHandler(struct inputConfig &cf_):cf(cf_){};

  public:
    static fiestaSignalHandler *getInstance(struct inputConfig &cf){
      if (!instance)
        instance = new fiestaSignalHandler(cf);
      return instance;
    }

    void registerSignals(){
      // write restart
      signal(SIGUSR1,fiestaSignalHandler::sigusr1Handler);

      // write restart and exit
      signal(SIGURG,fiestaSignalHandler::sigurgHandler);
      signal(SIGUSR2,fiestaSignalHandler::sigurgHandler);

      //exit
      signal(SIGTERM,fiestaSignalHandler::sigtermHandler);
      signal(SIGINT,fiestaSignalHandler::sigtermHandler);
    }

    void sigurgFunction([[maybe_unused]] int signum){
      cf.exitFlag=1;
      cf.restartFlag=1;
    }
    void sigusr1Function([[maybe_unused]] int signum){
      cf.restartFlag=1;
    }
    void sigtermFunction([[maybe_unused]] int signum){
      cf.exitFlag=1;
    }

    static void sigurgHandler(int signum){
      instance->sigurgFunction(signum);
    }
    static void sigusr1Handler(int signum){
      instance->sigusr1Function(signum);
    }
    static void sigtermHandler(int signum){
      instance->sigtermFunction(signum);
    }
};

fiestaSignalHandler *fiestaSignalHandler::instance=0;
#endif
