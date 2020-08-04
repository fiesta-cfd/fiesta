#ifndef TIMER_H
#define TIMER_H

#include "Kokkos_Core.hpp"
#include <string>

using namespace std;

class fiestaTimer {

public:
  fiestaTimer();
  fiestaTimer(string description);
  string name;
  void accumulate();
  double get();
  void clear();
  void reset();
  void start();
  void stop();
  string describe();
  string getf();
  string checkf();
  double check();
  string formatTime(double tin);

private:
  double time;
  Kokkos::Timer *timer;
  string description;
};
#endif
