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
      signal(SIGINT,fiestaSignalHandler::sigintHandler);
      signal(SIGUSR1,fiestaSignalHandler::sigusr1Handler);
      signal(SIGTERM,fiestaSignalHandler::sigtermHandler);
    }

    void sigintFunction(int signum){
      //cf.log->warning("Recieved SIGINT:  Writing restart and exiting after timestep ",cf.t,".");
      cf.exitFlag=1;
      cf.restartFlag=1;
    }
    void sigusr1Function(int signum){
      //cf.log->message("Recieved SIGUSR1:  Writing restart after timestep ",cf.t,".");
      cf.restartFlag=1;
    }
    void sigtermFunction(int signum){
      //cf.log->error("Recieved SIGTERM:  Exiting after timestep ",cf.t,".");
      cf.exitFlag=1;
    }

    static void sigintHandler(int signum){
      instance->sigintFunction(signum);
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
