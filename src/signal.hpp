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

    void sigurgFunction(int signum){
      cf.exitFlag=1;
      cf.restartFlag=1;
    }
    void sigusr1Function(int signum){
      cf.restartFlag=1;
    }
    void sigtermFunction(int signum){
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
