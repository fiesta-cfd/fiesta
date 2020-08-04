#ifndef RKFUNCTION_H
#define RKFUNCTION_H

#include "Kokkos_Core.hpp"
#include "fiesta.hpp"
#include "input.hpp"
#include "particle.hpp"
#include "timer.hpp"
#include <map>
#include <string>

class rk_func {

public:
  rk_func(struct inputConfig &cf_, Kokkos::View<double *> &cd_);

  virtual void compute() = 0;
  virtual void preStep() = 0;
  virtual void postStep() = 0;
  virtual void preSim() = 0;
  virtual void postSim() = 0;
  // virtual void compute(const FS4D & mvar, FS4D & mdvar) = 0;

  FS4D var;
  FS4D dvar;
  FS4D tmp1;
  FS4D grid;
  FSP2D particles;   // Particle array
  FSP2DH particlesH; // Particle host array

  std::map<std::string, fiestaTimer> timers;
  policy_f ghostPol = policy_f({0, 0}, {1, 1});
  policy_f cellPol = policy_f({0, 0}, {1, 1});
  policy_f facePol = policy_f({0, 0}, {1, 1});

protected:
  struct inputConfig &cf;
  Kokkos::View<double *> mcd;
};

#endif
