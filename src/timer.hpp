#ifndef TIMER_H
#define TIMER_H

#include "Kokkos_Core.hpp"
#include <string>

using namespace std;

class fiestaTimer {

public:
  fiestaTimer();
  fiestaTimer(int format);
  fiestaTimer(string description);
  string name;
  int format;
  void accumulate();
  double get();
  void clear();
  void reset();
  void start();
  void stop();
  string describe();
  string getf();
  string getf(int format);
  string checkf();
  string checkf(int format);
  double check();
  string formatTime(double tin);

private:
  double time;
  Kokkos::Timer *timer;
  string description;
};
#endif
