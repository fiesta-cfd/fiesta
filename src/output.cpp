#include "input.hpp"
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <locale>
#include "unistd.h"
#include "output.hpp"

using namespace std;

string c(Color color){
    //if (isatty(fileno(stdout))){
    ////if (isatty(STDOUT_FILENO)){
    //    switch (color){
    //        case RED :
    //            return "\033[0;31m";
    //        case GRE :
    //            return "\033[0;32m";
    //        case YEL :
    //            return "\033[0;33m";
    //        case BLU :
    //            return "\033[0;34m";
    //        case MAG :
    //            return "\033[0;35m";
    //        case CYA :
    //            return "\033[0;36m";
    //        case NON :
    //            return "\033[0m";
    //        default :
    //            return "\033[0m";
    //    }
    //}else{
        return "";
    //}
}

//#define c(CYA) "\033[0;36m"
//#define MAG "\033[0;35m"
//#define BLU "\033[0;34m"
//#define YEL "\033[0;33m"
//#define c(GRE) "\033[0;32m"
//#define RED "\033[0;31m"
//#define c(NON) "\033[0m"

void printSplash(){
    //cout << "----" << c(NON) << R"(------------------------------)" << c(NON) << "----" << endl;
    //cout << "|   " << c(RED) << R"(  _____ _           _         )" << c(NON) << "   |" << endl;
    //cout << "|   " << c(RED) << R"( |  ___(_) ___  ___| |_ __ _  )" << c(NON) << "   |" << endl;
    //cout << "|   " << c(RED) << R"( | |_  | |/ _ \/ __| __/ _` | )" << c(NON) << "   |" << endl;
    //cout << "|   " << c(RED) << R"( |  _| | |  __/\__ \ || (_| | )" << c(NON) << "   |" << endl;
    //cout << "|   " << c(RED) << R"( |_|   |_|\___||___/\__\__,_| )" << c(NON) << "   |" << endl;
    //cout << "|   " << c(RED) << R"(                              )" << c(NON) << "   |" << endl;
    //cout << "----" << c(NON) << R"(------------------------------)" << c(NON) << "----" << endl;

    cout << " " << endl;
    cout << c(RED) << R"(   _____ _           _         )" << c(NON) << endl;
    cout << c(RED) << R"(  |  ___(_) ___  ___| |_ __ _  )" << c(NON) << endl;
    cout << c(RED) << R"(  | |_  | |/ _ \/ __| __/ _` | )" << c(NON) << endl;
    cout << c(RED) << R"(  |  _| | |  __/\__ \ || (_| | )" << c(NON) << endl;
    cout << c(RED) << R"(  |_|   |_|\___||___/\__\__,_| )" << c(NON) << endl;
    cout << c(RED) << R"(                               )" << c(NON) << endl;
    cout << " " << endl;
    cout << "Version: " << c(CYA) << "'" << FIESTA_VERSION << "'" << c(NON) << endl;
    cout << "Type:    " << c(CYA) << "'" << FIESTA_OPTIONS << "'" << c(NON) << endl << endl;
    cout << "-----------------------" << endl << endl;
}

void printConfig(struct inputConfig cf){

    cout.precision(2);

    cout << "Title: " << c(BLU) << cf.title << c(NON) << endl;
    cout << "Input File Name: " << c(CYA) << "'" << cf.inputFname << "'" << c(NON) << endl << endl;

    cout << "-----------------------" << endl << endl;

    // Restart Options
    cout << c(GRE) << "Restart:" << c(NON) << endl;
    if (cf.restart){
        cout << "    Running from Restart File: " << c(CYA) << "'" << cf.sfName << "'" << c(NON) << endl;
        cout << "    Start Time:  " << scientific << c(CYA) << cf.time << c(NON) << endl;
        cout << "    Start Index: " << c(CYA) << cf.tstart << c(NON) << endl;
    }else
        cout << "    Not running from a restart" << endl;

    // Output Frequencies
    cout << c(GRE) << "Output:" << c(NON) << endl;
    if (cf.out_freq > 0)
        cout << setw(30) << left << "    Progress Output frequency: " << c(CYA) << right << setw(0) << cf.out_freq << c(NON) << endl;
    else
        cout << setw(19) << left << "    Progress Output " << c(YEL) << "disabled" << c(NON) << endl;
    if (cf.write_freq > 0)
        cout << setw(30) << left << "    Solution write frequency: "  << c(CYA) << right << setw(0) << cf.write_freq << c(NON) << endl;
    else
        cout << setw(19) << left << "    CGNS writes " << c(YEL) << "disabled" << c(NON) << endl;
    if (cf.restart_freq > 0)
        cout << setw(30) << left << "    Restart write frequency: " << c(CYA) << right << setw(0) << cf.restart_freq << c(NON) << endl;
    else
        cout << setw(19) << left << "    Restart writes " << c(YEL) << "disabled" << c(NON) << endl;

    // Discretization
    cout << c(GRE) << "Discretization:" << c(NON) << endl;
    cout << "    Running " << c(CYA) << cf.numProcs << c(NON) << " processes as ("
         << c(CYA) << cf.xProcs << c(NON) << ","
         << c(CYA) << cf.yProcs << c(NON) << ","
         << c(CYA) << cf.zProcs << c(NON) << ")" << endl;

    cout << "    tstart = "<< c(CYA) << cf.tstart << c(NON) << ", nt = " << c(CYA) << cf.nt << c(NON) << ", tend = "
         << c(CYA) << cf.tend << c(NON) << ", dt = " << scientific << c(CYA) << cf.dt << c(NON) << endl;

    cout << "    Num Cells X = " << c(CYA) << cf.glbl_nci << c(NON) << ", dx = " << scientific << c(CYA) << cf.dx <<  c(NON) << endl;
    cout << "    Num Cells Y = " << c(CYA) << cf.glbl_ncj << c(NON) << ", dy = " << scientific << c(CYA) << cf.dy <<  c(NON) << endl;
    if (cf.ndim == 3)
        cout << "    Num Cells Z = " << c(CYA) << cf.glbl_nck << c(NON) << ", dz = " << scientific << c(CYA) << cf.dz << c(NON)<< endl;
    cout << "    Num Cells Total: " << scientific << c(CYA) << (double)cf.glbl_nck*cf.glbl_ncj*cf.glbl_nci << c(NON) << endl;

    cout << c(GRE) << "Options:" << c(NON) << endl;
    if (cf.scheme == 1)
        cout << "    Using 5th order weno scheme" << endl;
    if (cf.scheme == 2)
        cout << "    Using 4th order centered scheme" << endl;
    if (cf.visc)
        cout << setw(15) << left << "    Viscosity " << c(YEL) << "enabled" <<c(NON) <<  endl;
    else
        cout << setw(15) << left << "    Viscosity " << c(YEL) << "disabled" << c(NON) << endl;
    // C-Equations
    if (cf.ceq)
        cout << setw(15) << left << "    C-Equation " << c(YEL) << "enabled" << c(NON) << endl;
    else
        cout << setw(15) << left << "    C-Equation " << c(YEL) << "disabled" << c(NON) << endl;


    // Gas Properties
    cout << c(GRE) << "Gas Properties:" << c(NON) << endl;
    cout << "    Number of Species = " << c(CYA) << cf.ns << c(NON) << ":" << endl;
    for (int s=0; s<cf.ns; ++s)
        cout << "    Species " << s+1 
             << ", Gamma = " << setw(4) << setprecision(3) << fixed      << c(CYA) << cf.gamma[s] << c(NON)
             << ", M = "     << setw(6) << setprecision(4) << scientific << c(CYA) << cf.M[s] << c(NON)
             << ", mu = "    << setw(6) << setprecision(4) << scientific << c(CYA) << cf.mu[s] << c(NON)
             << endl;
        //printf("    Species %d, Gamma = %4.2f, M = %6.4f\n",s+1,cf.gamma[s],cf.M[s]);

    cout << endl << "-----------------------" << endl << endl << flush;
}
