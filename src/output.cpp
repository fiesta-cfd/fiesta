#include "input.hpp"
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <locale>
#include "unistd.h"
#include "output.hpp"

using namespace std;

string c(Color color){
    if (isatty(STDOUT_FILENO)){
        switch (color){
            case RED :
                return "\033[0;31m";
            case GRE :
                return "\033[0;32m";
            case YEL :
                return "\033[0;33m";
            case BLU :
                return "\033[0;34m";
            case MAG :
                return "\033[0;35m";
            case CYA :
                return "\033[0;36m";
            case NON :
                return "\033[0m";
            default :
                return "\033[0m";
        }
    }else{
        return "";
    }
}

//#define c(CYA) "\033[0;36m"
//#define MAG "\033[0;35m"
//#define BLU "\033[0;34m"
//#define YEL "\033[0;33m"
//#define c(GRE) "\033[0;32m"
//#define RED "\033[0;31m"
//#define c(NON) "\033[0m"

void printSplash(){
    //cout << R"(--------------------------------------)" << endl;
    //cout << R"(|     _____ _           _            |)" << endl;
    //cout << R"(|    |  ___(_) ___  ___| |_ __ _     |)" << endl;
    //cout << R"(|    | |_  | |/ _ \/ __| __/ _` |    |)" << endl;
    //cout << R"(|    |  _| | |  __/\__ \ || (_| |    |)" << endl;
    //cout << R"(|    |_|   |_|\___||___/\__\__,_|    |)" << endl;
    //cout << R"(|                                    |)" << endl;
    //cout << R"(--------------------------------------)" << endl;
    cout << "----" << c(NON) << R"(------------------------------)" << c(NON) << "----" << endl;
    cout << "|   " << c(RED) << R"(  _____ _           _         )" << c(NON) << "   |" << endl;
    cout << "|   " << c(RED) << R"( |  ___(_) ___  ___| |_ __ _  )" << c(NON) << "   |" << endl;
    cout << "|   " << c(RED) << R"( | |_  | |/ _ \/ __| __/ _` | )" << c(NON) << "   |" << endl;
    cout << "|   " << c(RED) << R"( |  _| | |  __/\__ \ || (_| | )" << c(NON) << "   |" << endl;
    cout << "|   " << c(RED) << R"( |_|   |_|\___||___/\__\__,_| )" << c(NON) << "   |" << endl;
    cout << "|   " << c(RED) << R"(                              )" << c(NON) << "   |" << endl;
    cout << "----" << c(NON) << R"(------------------------------)" << c(NON) << "----" << endl;
}

void printConfig(struct inputConfig cf){

    cout.precision(2);

    std::cout << endl << "Input File Name: " << c(CYA) << "'" << cf.inputFname << "'" << c(NON) << "'" << endl << endl;

    // Restart Options
    cout << c(GRE) << "Restart:" << c(NON) << endl;
    if (cf.restart){
        cout << "    Running from Restart File: " << c(CYA) << "'" << cf.sfName << "'" << c(NON) << endl;
        cout << "    Start Time is " << scientific << cf.time << endl;
        cout << "    Start Index is " << cf.tstart << endl;
    }else
        cout << "    Not running from a restart" << endl;

    // Output Frequencies
    cout << c(GRE) << "Output:" << c(NON) << endl;
    if (cf.out_freq > 0)
        cout << setw(36) << left << "    Output frequency: " << c(CYA) << right << setw(5) << cf.out_freq << c(NON) << endl;
    else
        cout << setw(36) << left << "    Output " << YEL << "disabled" << c(NON) << endl;
    if (cf.write_freq > 0)
        cout << setw(36) << left << "    CGNS write frequency: "  << c(CYA) << right << setw(5) << cf.write_freq << c(NON) << endl;
    else
        cout << setw(36) << left << "    CGNS writes " << YEL << "disabled" << c(NON) << endl;
    if (cf.restart_freq > 0)
        cout << setw(36) << left << "    Restart write frequency: " << c(CYA) << right << setw(5) << cf.restart_freq << c(NON) << endl;
    else
        cout << setw(36) << left << "    Restart writes " << YEL << "disabled" << c(NON) << endl;

    // Discretization
    cout << c(GRE) << "Discretization:" << c(NON) << endl;
    cout << "    Running " << cf.numProcs << " processes as ("
         << cf.xProcs << "," << cf.yProcs << "," << cf.zProcs << ")" << endl;
    cout << "    tstart = "<< cf.tstart << ", nt = " << cf.nt << ", tend = "
         << cf.tend << ", dt = " << scientific << cf.dt << endl;
    cout << "    Num Cells X = " << cf.glbl_nci << ", dx = " << scientific << cf.dx << endl;
    cout << "    Num Cells Y = " << cf.glbl_ncj << ", dy = " << scientific << cf.dy << endl;
    if (cf.ndim == 3)
        cout << "    Num Cells Z = " << cf.glbl_nck << ", dz = " << scientific << cf.dz << endl;
    cout << "    Num Cells Total: " << scientific << (double)cf.glbl_nck*cf.glbl_ncj*cf.glbl_nci << endl;

    cout << c(GRE) << "Options:" << c(NON) << endl;
    if (cf.scheme == 1)
        cout << "    USing 5th order weno scheme" << endl;
    if (cf.scheme == 2)
        cout << "    Using 4th order centered scheme" << endl;
    if (cf.visc)
        cout << "    Viscosity enabled" << endl;
    else
        cout << "    Viscosity disabled" << endl;
    // C-Equations
    if (cf.ceq)
        cout << "    C-Equation enabled" << endl;
    else
        cout << "    C-Equation disabled" << endl;


    // Gas Properties
    cout << c(GRE) << "Gas Properties:" << c(NON) << endl;
    cout << "    Number of Species = " << cf.ns << ":" << endl;
    for (int s=0; s<cf.ns; ++s)
        cout << "    Species " << s+1 
             << ", Gamma = " << setw(4) << setprecision(3) << fixed      << c(CYA) << cf.gamma[s] << c(NON)
             << ", M = "     << setw(6) << setprecision(4) << scientific << c(CYA) << cf.M[s] << c(NON)
             << ", mu = "    << setw(6) << setprecision(4) << scientific << c(CYA) << cf.mu[s] << c(NON)
             << endl;
        //printf("    Species %d, Gamma = %4.2f, M = %6.4f\n",s+1,cf.gamma[s],cf.M[s]);

    cout << endl << "-----------------------" << endl << flush;
}
