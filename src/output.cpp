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
#include <string>
#include "pretty.hpp"
#include "fmt/core.h"

using std::cout;
using fmt::format;

void printSplash(int cm) {
  ansiColors k(cm);
  string splashFmt = format("{}{: >4}{{}}{}\n",k(red),"",k(reset));

  cout << format(splashFmt,R"(   _____ _           _         )");
  cout << format(splashFmt,R"(  |  ___(_) ___  ___| |_ __ _  )");
  cout << format(splashFmt,R"(  | |_  | |/ _ \/ __| __/ _` | )");
  cout << format(splashFmt,R"(  |  _| | |  __/\__ \ || (_| | )");
  cout << format(splashFmt,R"(  |_|   |_|\___||___/\__\__,_| )");
  cout << format(splashFmt,R"(                               )");

  cout << format("{: <12}{}'{}'{}\n","Version:",   k(blue),FIESTA_VERSION,k(reset));
  cout << format("{: <12}{}'{}'{}\n","Build Type:",k(blue),FIESTA_OPTIONS,k(reset));
  cout << format("{: <12}{}'{}'{}\n","Build Time:",k(blue),FIESTA_BTIME  ,k(reset));
  cout << std::endl;
}

void printConfig(struct inputConfig cf) {
  if (cf.rank == 0){
  
    ansiColors k(cf.colorFlag);
    string keyString   = format("{: >8}{{: <24}}{}'{{}}'{}\n","",k(blue),k(reset));
    string keyValue    = format("{: >8}{{: <24}}{}{{}}{}\n","",k(magenta),k(reset));
    string keyEnabled  = format("{: >8}{{: <24}}{}enabled{}\n","",k(green),k(reset));
    string keyDisabled = format("{: >8}{{: <24}}{}disabled{}\n","",k(yellow),k(reset));
    string keyYes      = format("{: >8}{{: <24}}{}yes{}\n","",k(green),k(reset));
    string keyNo       = format("{: >8}{{: <24}}{}no{}\n","",k(yellow),k(reset));
    string keyTupleInt = format("{: >8}{{: <24}}({}{{}}{},{}{{}}{},{}{{}}{})\n","",
                       k(magenta),k(reset),k(magenta),k(reset),k(magenta),k(reset));
    string keyTupleGen = format("{: >8}{{: <24}}({}{{:.3g}}{},{}{{:.3g}}{},{}{{:.3g}}{})\n","",
                       k(magenta),k(reset),k(magenta),k(reset),k(magenta),k(reset));
    string keyPair     = format("{: >8}{{: <24}}({}{{}}{},{}{{}}{})\n","",
                                  k(magenta),k(reset),k(magenta),k(reset));
  
    cout << format(keyString,"Title:",cf.title);
    cout << format(keyString,"Input file name:",cf.inputFname);
  
    if (cf.restart) {
      cout << format(keyString,"Running from:",cf.restartName);
      cout << format(keyValue,"Start Time:",cf.time);
      cout << format(keyValue,"Start Index",cf.tstart);
    }
  
    cout << format(keyString,"Solution Path:",cf.pathName);
    
    if (cf.out_freq > 0) cout << format(keyValue,"Progress Frequency:",cf.out_freq);
    else cout << format(keyDisabled,"Progress Indicator:");
  
    if (cf.write_freq > 0) cout << format(keyValue,"HDF5 Frequency:",cf.write_freq);
    else cout << format(keyDisabled,"HDF5:");
  
    if (cf.restart_freq > 0) cout << format(keyValue,"Restart frequency:",cf.restart_freq);
    else cout << format(keyDisabled,"Restart writes:");
  
    if (cf.stat_freq > 0) cout << format(keyValue,"Status frequency:",cf.stat_freq);
    else cout << format(keyDisabled,"Status reports:");

#ifdef HAVE_SINGLE
    cout << format(keyString,"Precision:","single");
#else
    cout << format(keyString,"Precision:","double");
#endif
  
    cout << format(keyValue,"Number of Processes:",cf.numProcs);
    cout << format(keyTupleInt,"MPI Discretization:",cf.xProcs,cf.yProcs,cf.zProcs);
  
    cout << format(keyValue,"tstart:",cf.tstart);
    cout << format(keyValue,"nt:",cf.nt);
    cout << format(keyValue,"tend:",cf.tend);
    cout << format(keyValue,"dt:",cf.dt);

    if (cf.ndim==3){
      cout << format(keyTupleInt,"Number of Cells:",cf.glbl_nci,cf.glbl_ncj,cf.glbl_nck);
      cout << format(keyTupleGen,"Cell Size:",cf.dx,cf.dy,cf.dz);
    }else{
      cout << format(keyPair,"Number of Cells:",cf.glbl_nci,cf.glbl_ncj);
      cout << format(keyPair,"Cell Size:",cf.dx,cf.dy);
    }
  
    if (cf.scheme == 1) cout << format(keyString,"Scheme:","weno5");
    if (cf.scheme == 2) cout << format(keyString,"Scheme:","centered4");
    if (cf.scheme == 3) cout << format(keyString,"Scheme:","quick");
  
    if (cf.visc) cout << format(keyEnabled,"Viscosity:");
    else cout << format(keyDisabled,"Viscosity:");
  
    if (cf.ceq) cout << format(keyEnabled,"C-Equations:");
    else cout << format(keyDisabled,"C-Equations:");
  
    if (cf.noise) cout << format(keyEnabled,"Noise Indicator:");
    else cout << format(keyDisabled,"Noise Indicator:");

    if (cf.buoyancy) cout << format(keyEnabled,"Buoyancy:");
    else cout << format(keyDisabled,"Buoyancy:");

    if (cf.chunkable) cout << format(keyYes,"HDF5 Chunks:");
    else cout << format(keyEnabled,"HDF5 Chunks:");

    if (cf.compressible) cout << format(keyYes,"HDF5 Compression:");
    else cout << format(keyDisabled,"HDF5 Compression:");
  
    cout << format(keyValue,"Number of species:",cf.ns);
    string val = format("{}{{}}{}",k(magenta),k(reset));
    for (int s = 0; s < cf.ns; ++s)
      cout << format("        {}: gamma={}, M={}, mu={}\n",cf.speciesName[s],
          format(val,cf.gamma[s]),
          format(val,cf.M[s]),
          format(val,cf.mu[s]));
  }
  cout << std::flush;
}
