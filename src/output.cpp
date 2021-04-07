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

#include "output.hpp"
#include "Kokkos_Core.hpp"
#include "input.hpp"
#include "unistd.h"
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <locale>

using namespace std;

// if color is enabled insert ascii color codes
string c(int cFlag, Color color) {
  int use_color = 0;
  if (cFlag == 1)
    use_color = 1;
  if ((cFlag == 2) && isatty(fileno(stdout)))
    use_color = 1;
  if (use_color) {
    switch (color) {
    case RED:
      return "\033[0;31m";
    case GRE:
      return "\033[0;32m";
    case YEL:
      return "\033[0;33m";
    case BLU:
      return "\033[0;34m";
    case MAG:
      return "\033[0;35m";
    case CYA:
      return "\033[0;36m";
    case NON:
      return "\033[0m";
    default:
      return "\033[0m";
    }
  } else {
    return "";
  }
}

void printSplash(int cFlag) {
  cout << " " << endl;
  cout << c(cFlag, RED) << R"(   _____ _           _         )" << c(cFlag, NON)
       << endl;
  cout << c(cFlag, RED) << R"(  |  ___(_) ___  ___| |_ __ _  )" << c(cFlag, NON)
       << endl;
  cout << c(cFlag, RED) << R"(  | |_  | |/ _ \/ __| __/ _` | )" << c(cFlag, NON)
       << endl;
  cout << c(cFlag, RED) << R"(  |  _| | |  __/\__ \ || (_| | )" << c(cFlag, NON)
       << endl;
  cout << c(cFlag, RED) << R"(  |_|   |_|\___||___/\__\__,_| )" << c(cFlag, NON)
       << endl;
  cout << c(cFlag, RED) << R"(                               )" << c(cFlag, NON)
       << endl;
  cout << " " << endl;
  cout << "Version:    " << c(cFlag, CYA) << "'" << FIESTA_VERSION << "'"
       << c(cFlag, NON) << endl;
  cout << "Build Type: " << c(cFlag, CYA) << "'" << FIESTA_OPTIONS << "'"
       << c(cFlag, NON) << endl;
  cout << "Build Time: " << c(cFlag, CYA) << "'" << FIESTA_BTIME << "'"
       << c(cFlag, NON) << endl
       << endl;
  cout << "-----------------------" << endl << endl;
}

void printConfig(struct inputConfig cf) {
if (cf.rank == 0){
  int cFlag = cf.colorFlag;

  cout.precision(2);

  // Kokkos::print_configuration(cout,1);

  cout << "Title: " << c(cFlag, BLU) << cf.title << c(cFlag, NON) << endl;
  cout << "Input File Name: " << c(cFlag, CYA) << "'" << cf.inputFname << "'"
       << c(cFlag, NON) << endl
       << endl;

  cout << "-----------------------" << endl << endl;

  // Restart Options
  cout << c(cFlag, GRE) << "Restart:" << c(cFlag, NON) << endl;
  if (cf.restart) {
    cout << "    Running from Restart File: " << c(cFlag, CYA) << "'"
         << cf.restartName << "'" << c(cFlag, NON) << endl;
    cout << "    Start Time:  " << scientific << c(cFlag, CYA) << cf.time
         << c(cFlag, NON) << endl;
    cout << "    Start Index: " << c(cFlag, CYA) << cf.tstart << c(cFlag, NON)
         << endl;
  } else
    cout << "    Not running from a restart" << endl;

  // Output Frequencies
  cout << c(cFlag, GRE) << "Output:" << c(cFlag, NON) << endl;
  cout   << "    Path: " << c(cFlag,BLU) << cf.pathName << c(cFlag,NON) << endl;
  if (cf.out_freq > 0)
    cout << setw(31) << left
         << "    Timestep output frequency: " << c(cFlag, CYA) << right
         << setw(0) << cf.out_freq << c(cFlag, NON) << endl;
  else
    cout << setw(20) << left << "    Timestep output " << c(cFlag, YEL)
         << "disabled" << c(cFlag, NON) << endl;
  if (cf.write_freq > 0)
    cout << setw(31) << left
         << "    Solution write frequency: " << c(cFlag, CYA) << right
         << setw(0) << cf.write_freq << c(cFlag, NON) << endl;
  else
    cout << setw(20) << left << "    CGNS writes " << c(cFlag, YEL)
         << "disabled" << c(cFlag, NON) << endl;
  if (cf.restart_freq > 0)
    cout << setw(31) << left << "    Restart write frequency: " << c(cFlag, CYA)
         << right << setw(0) << cf.restart_freq << c(cFlag, NON) << endl;
  else
    cout << setw(20) << left << "    Restart writes " << c(cFlag, YEL)
         << "disabled" << c(cFlag, NON) << endl;
  if (cf.stat_freq > 0)
    cout << setw(31) << left << "    Status check frequency: " << c(cFlag, CYA)
         << right << setw(0) << cf.stat_freq << c(cFlag, NON) << endl;
  else
    cout << setw(20) << left << "    Status checks " << c(cFlag, YEL)
         << "disabled" << c(cFlag, NON) << endl;

  // Discretization
  cout << c(cFlag, GRE) << "Discretization:" << c(cFlag, NON) << endl;
  cout << "    Running " << c(cFlag, CYA) << cf.numProcs << c(cFlag, NON)
       << " processes as (" << c(cFlag, CYA) << cf.xProcs << c(cFlag, NON)
       << "," << c(cFlag, CYA) << cf.yProcs << c(cFlag, NON) << ","
       << c(cFlag, CYA) << cf.zProcs << c(cFlag, NON) << ")" << endl;

  cout << "    tstart = " << c(cFlag, CYA) << cf.tstart << c(cFlag, NON)
       << ", nt = " << c(cFlag, CYA) << cf.nt << c(cFlag, NON)
       << ", tend = " << c(cFlag, CYA) << cf.tend << c(cFlag, NON)
       << ", dt = " << scientific << c(cFlag, CYA) << cf.dt << c(cFlag, NON)
       << endl;

  cout << "    Num Cells X = " << c(cFlag, CYA) << cf.glbl_nci << c(cFlag, NON)
       << ", dx = " << scientific << c(cFlag, CYA) << cf.dx << c(cFlag, NON)
       << endl;
  cout << "    Num Cells Y = " << c(cFlag, CYA) << cf.glbl_ncj << c(cFlag, NON)
       << ", dy = " << scientific << c(cFlag, CYA) << cf.dy << c(cFlag, NON)
       << endl;
  if (cf.ndim == 3)
    cout << "    Num Cells Z = " << c(cFlag, CYA) << cf.glbl_nck
         << c(cFlag, NON) << ", dz = " << scientific << c(cFlag, CYA) << cf.dz
         << c(cFlag, NON) << endl;
  cout << "    Num Cells Total: " << scientific << c(cFlag, CYA)
       << (double)cf.glbl_nck * cf.glbl_ncj * cf.glbl_nci << c(cFlag, NON)
       << endl;

  // Options
  cout << c(cFlag, GRE) << "Options:" << c(cFlag, NON) << endl;
  if (cf.scheme == 1)
    cout << "    Using 5th order weno scheme" << endl;
  if (cf.scheme == 2)
    cout << "    Using 4th order centered scheme" << endl;
  if (cf.scheme == 3)
    cout << "    Using Quick scheme" << endl;
  if (cf.visc)
    cout << setw(18) << left << "    Viscosity " << c(cFlag, YEL) << "enabled"
         << c(cFlag, NON) << endl;
  else
    cout << setw(18) << left << "    Viscosity " << c(cFlag, YEL) << "disabled"
         << c(cFlag, NON) << endl;
  // C-Equations
  if (cf.ceq)
    cout << setw(18) << left << "    C-Equation " << c(cFlag, YEL) << "enabled"
         << c(cFlag, NON) << endl;
  else
    cout << setw(18) << left << "    C-Equation " << c(cFlag, YEL) << "disabled"
         << c(cFlag, NON) << endl;
  // Noise Removal
  if (cf.noise == 1)
    cout << setw(18) << left << "    Noise Removal " << c(cFlag, YEL)
         << "enabled" << c(cFlag, NON) << endl;
  else
    cout << setw(18) << left << "    Noise Removal " << c(cFlag, YEL)
         << "disabled" << c(cFlag, NON) << endl;

  // Gas Properties
  cout << c(cFlag, GRE) << "Gas Properties:" << c(cFlag, NON) << endl;
  cout << "    Number of Species = " << c(cFlag, CYA) << cf.ns << c(cFlag, NON)
       << ":" << endl;
  for (int s = 0; s < cf.ns; ++s)
    cout << "    Species " << s + 1 << ", Gamma = " << setw(4)
         << setprecision(3) << fixed << c(cFlag, CYA) << cf.gamma[s]
         << c(cFlag, NON) << ", M = " << setw(6) << setprecision(4)
         << scientific << c(cFlag, CYA) << cf.M[s] << c(cFlag, NON)
         << ", mu = " << setw(6) << setprecision(4) << scientific
         << c(cFlag, CYA) << cf.mu[s] << c(cFlag, NON) << endl;
  // printf("    Species %d, Gamma = %4.2f, M =
  // %6.4f\n",s+1,cf.gamma[s],cf.M[s]);

  cout << endl << "-----------------------" << endl << endl << flush;
}
}
