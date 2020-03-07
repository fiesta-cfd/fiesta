#include "input.hpp"
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <locale>

using namespace std;

void printSplash(){
    cout << R"(--------------------------------------)" << endl;
    cout << R"(|     _____ _           _            |)" << endl;
    cout << R"(|    |  ___(_) ___  ___| |_ __ _     |)" << endl;
    cout << R"(|    | |_  | |/ _ \/ __| __/ _` |    |)" << endl;
    cout << R"(|    |  _| | |  __/\__ \ || (_| |    |)" << endl;
    cout << R"(|    |_|   |_|\___||___/\__\__,_|    |)" << endl;
    cout << R"(|                                    |)" << endl;
    cout << R"(--------------------------------------)" << endl;
}

void printConfig(struct inputConfig cf){

    cout.precision(2);

    std::cout << endl << "Input File Name: '" << cf.inputFname << "'" << endl << endl;

    // Restart Options
    cout << "Restart:" << endl;
    if (cf.restart){
        cout << "    Running from Restart File: '" << cf.sfName << "'" << endl;
        cout << "    Start Time is " << scientific << cf.time << endl;
        cout << "    Start Index is " << cf.tstart << endl;
    }else
        cout << "    Not running from a restart" << endl;

    // Output Frequencies
    cout << "Output:" << endl;
    if (cf.out_freq > 0)
        cout << "    Output frequency: " << cf.out_freq << endl;
    else
        cout << "    Output disabled" << endl;
    if (cf.write_freq > 0)
        cout << "    CGNS write frequency: " << cf.write_freq << endl;
    else
        cout << "    CGNS writes disabled" << endl;
    if (cf.restart_freq > 0)
        cout << "    Restart write frequency: " << cf.restart_freq << endl;
    else
        cout << "    Restart writes disabled" << endl;

    // Discretization
    cout << "Discretization:" << endl;
    cout << "    Running " << cf.numProcs << " processes as ("
         << cf.xProcs << "," << cf.yProcs << "," << cf.zProcs << ")" << endl;
    cout << "    tstart = "<< cf.tstart << ", nt = " << cf.nt << ", tend = "
         << cf.tend << ", dt = " << scientific << cf.dt << endl;
    cout << "    Num Cells X = " << cf.glbl_nci << ", dx = " << scientific << cf.dx << endl;
    cout << "    Num Cells Y = " << cf.glbl_ncj << ", dy = " << scientific << cf.dy << endl;
    if (cf.ndim == 3)
        cout << "    Num Cells Z = " << cf.glbl_nck << ", dz = " << scientific << cf.dz << endl;
    cout << "    Num Cells Total: " << scientific << (double)cf.glbl_nck*cf.glbl_ncj*cf.glbl_nci << endl;

    cout << "Options:" << endl;
    // C-Equations
    if (cf.ceq)
        cout << "    C-Equation enabled" << endl;
    else
        cout << "    C-Equation disabled" << endl;
    if (cf.visc)
        cout << "    Viscosity enabled" << endl;
    else
        cout << "    Viscosity disabled" << endl;

    // Gas Properties
    cout << "Gas Properties:" << endl;
    cout << "    Number of Species = " << cf.ns << ":" << endl;
    for (int s=0; s<cf.ns; ++s)
        cout << "    Species " << s+1 
             << ", Gamma = " << setw(4) << setprecision(2) << cf.gamma[s]
             << ", M = "     << setw(6) << setprecision(4) << cf.M[s]
             << endl;
        //printf("    Species %d, Gamma = %4.2f, M = %6.4f\n",s+1,cf.gamma[s],cf.M[s]);

    cout << endl << "-----------------------" << endl << flush;
}
