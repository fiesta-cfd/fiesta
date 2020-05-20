#include "fiesta.hpp"
#include "input.hpp"
#ifndef NOMPI
#include "mpi.hpp"
#include "cgns.hpp"
#include <mpi.h>
#endif
#include "bc.hpp"
#include "Kokkos_Core.hpp"
#include "debug.hpp"
#include "cart3d.hpp"
#include "cart2d.hpp"
#include "rkfunction.hpp"
#include <iostream>
#include <cstdio>
#include <ctime>
#include "output.hpp"
#include <iomanip>
#include "timer.hpp"
#include "status.hpp"
#include <set>

using namespace std;

void fnExit1(void){
    Kokkos::finalize();
}

int main(int argc, char* argv[]){

    // Create variables for command line information
    int cFlag,vFlag;
    string fName;

    // Get the command line options including input file name if supplied.
    getCommandlineOptions(argc, argv, vFlag,cFlag,fName);

    // If the --version command line option is present, then just print version
    // and compilation information and exit.
    if (vFlag){
        cout << "Fiesta" << endl;
        cout << "Version:    '" << FIESTA_VERSION << "'" << endl;
        cout << "Build Type: '" << FIESTA_OPTIONS << "'" << endl;
        cout << "Build Time: '" << FIESTA_BTIME   << "'" << endl;
        exit(EXIT_SUCCESS);
    }
    
    // Create and start some timers
    fiestaTimer totalTimer;     // Total program timer
    fiestaTimer initTimer;      // Timer for program initialization
    fiestaTimer solWriteTimer;  // Timer for solution file write times
    fiestaTimer simTimer;       // Timer for simulation time steppin
    fiestaTimer loadTimer;      // Timer for I.C. generation or restart read
    fiestaTimer resWriteTimer;  // Timer for restart file writes
    fiestaTimer gridTimer;      // Timer for grid generation
    fiestaTimer writeTimer;     // Timer for 

    // Initialize command MPI and get rank.
    int temp_rank;
    temp_rank = 0;
#ifndef NOMPI
    MPI_Init(NULL,NULL);
    MPI_Comm_rank(MPI_COMM_WORLD,&temp_rank);
#endif

    // Print spash screen with logo, version number and compilation info
    // to signify process startup.
    if (temp_rank == 0)
        printSplash(cFlag);

    // Initialize kokkos and set kokkos finalize as exit function.
    Kokkos::initialize(argc, argv);
    atexit(fnExit1);

    // Execute input file and generate simulation configuration
    struct inputConfig cf;
    cf = executeConfiguration(fName);

#ifndef NOMPI
    // Perform MPI decomposition and assign cells to ranks
    cf = mpi_init(cf);

    // create mpi buffers
    mpiBuffers m(cf);
#endif

    // Print out configuration
    if (cf.rank == 0)
        printConfig(cf,cFlag);

    // Create and copy minimal configuration array for data needed
    // withing Kokkos kernels.
    Kokkos::View<double*> cd("deviceCF",5+cf.ns*3);
    typename Kokkos::View<double*>::HostMirror hostcd = Kokkos::create_mirror_view(cd);
    Kokkos::deep_copy(hostcd, cd);
    hostcd(0) = cf.ns;  // number of gas species
    hostcd(1) = cf.dx;  // cell size
    hostcd(2) = cf.dy;  // cell size
    hostcd(3) = cf.dz;  // cell size
    hostcd(4) = cf.nv;  // number of flow variables
    
    // include gas properties for each gas species
    int sdx = 5;
    for (int s=0; s<cf.ns; ++s){
        hostcd(sdx) = cf.gamma[s];      // ratio of specific heats
        hostcd(sdx+1) = cf.R/cf.M[s];   // species gas comstant
        hostcd(sdx+2) = cf.mu[s];       // kinematic viscosity
        sdx += 3;
    }
    Kokkos::deep_copy(cd,hostcd);       // copy congifuration array to device

    /*** Choose Scheme ***/
    rk_func *f;
    if (cf.ndim == 3){
        f = new hydroc3d_func(cf,cd);
    }else{
        f = new hydro2dvisc_func(cf,cd);
    }

    // Create writer object
#ifndef NOMPI
    cgnsWriter w(cf,f->grid,f->var);
#endif

    // If not restarting, generate initial conditions and grid
    if (cf.restart == 0){
        if (cf.rank == 0) cout << c(cFlag,GRE) << "Generating Initial Conditions:" << c(cFlag,NON) << endl;
        loadTimer.start();
        loadInitialConditions(cf,f->var);
        loadTimer.stop();
        if (cf.rank == 0) cout << "    Generated in: " << c(cFlag,CYA) << loadTimer.getf() << c(cFlag,NON) << endl << endl;

        if (cf.rank == 0) cout << c(cFlag,GRE) << "Generating Grid:" << c(cFlag,NON) << endl;
        gridTimer.start();
        loadGrid(cf,f->grid);
        gridTimer.stop();
        if (cf.rank == 0) cout << "    Generated in: " << c(cFlag,CYA) << gridTimer.getf() << c(cFlag,NON) << endl << endl;
    }

    // execution policy for Runge Kutta update kernels, should be moved to header
    typedef Kokkos::MDRangePolicy<Kokkos::Rank<4>> policy_1;
    
    // Variables to track time step information
    double time = cf.time;
    int tstart = cf.tstart;
    

    // Read Restart or Write initial conditions
#ifndef NOMPI
    if (cf.restart == 1){
        if (cf.rank == 0) cout << c(cFlag,GRE) << "Loading Restart File:" << c(cFlag,NON) << endl;
        loadTimer.reset();
        w.readSolution(cf,f->grid,f->var);
        loadTimer.stop();
        cout << "    Loaded in: " << setprecision(2) << c(cFlag,CYA) << loadTimer.get() << "s" << c(cFlag,NON) << endl;
    }else{
        if (cf.rank == 0)
            if (cf.write_freq >0 || cf.restart_freq>0)
                cout << c(cFlag,GRE) << "Writing Initial Conditions:" << c(cFlag,NON) << endl;
        writeTimer.start();
        if (cf.write_freq >0){
            solWriteTimer.reset();
            w.writeSolution(cf,f->grid,f->var,0,0.00);
            solWriteTimer.accumulate();
        }
        if (cf.restart_freq >0){
            resWriteTimer.reset();
            w.writeRestart(cf,f->grid,f->var,0,0.00);
            resWriteTimer.accumulate();
        }
        writeTimer.stop();
        if (cf.rank == 0)
            if (cf.write_freq >0 || cf.restart_freq>0)
                cout << "    Wrote in: " << c(cFlag,CYA) << writeTimer.getf() << c(cFlag,NON) << endl;
    }

#endif

    // stop initialization timer
    initTimer.stop();

    // notify simulation start
    if (cf.rank == 0){
        cout << endl << "-----------------------" << endl << endl;
        cout << c(cFlag,GRE) << "Starting Simulation:" << c(cFlag,NON) << endl;
    }

    // reset simulation timer
    simTimer.reset();

    // pre simulation hook
    f->preSim();

  // // // // // // // //  \\ \\ \\ \\ \\ \\ \\ \\
 // // // // // // MAIN TIME LOOP \\ \\ \\ \\ \\ \\
// // // // // // // // //\\ \\ \\ \\ \\ \\ \\ \\ \\

    for (int t=tstart; t<cf.tend; ++t){
        time = time + cf.dt;

        //pre timestep hook
        f->preStep();

        // apply boundary conditions
#ifndef NOMPI
        applyBCs(cf,f->var,m);
#else
        applyBCs(cf,f->var);
#endif

        // First Stage Compute
        f->compute();
        
        // assign temporary variables
        FS4D mytmp = f->tmp1;
        FS4D myvar = f->var;
        FS4D mydvar = f->dvar;

        // First stage update
        Kokkos::parallel_for("Loop1", policy_1({0,0,0,0},{cf.ngi, cf.ngj, cf.ngk, cf.nvt}),
               KOKKOS_LAMBDA  (const int i, const int j, const int k, const int v) {
            mytmp(i,j,k,v) = myvar(i,j,k,v);
            myvar(i,j,k,v) = myvar(i,j,k,v) + 0.5*cf.dt*mydvar(i,j,k,v);
        });

        // apply boundary conditions
#ifndef NOMPI
        applyBCs(cf,f->var,m);
#else
        applyBCs(cf,f->var);
#endif

        // Seconf stage compute
        f->compute();

        //assign temporary variables
        mytmp = f->tmp1;
        myvar = f->var;
        mydvar = f->dvar;

        // Second stage update
        Kokkos::parallel_for("Loop2", policy_1({0,0,0,0},{cf.ngi, cf.ngj, cf.ngk, cf.nvt}),
               KOKKOS_LAMBDA  (const int i, const int j, const int k, const int v) {
            myvar(i,j,k,v) = myvar(i,j,k,v) + cf.dt*mydvar(i,j,k,v);
        });

        //post timestep hook
        f->postStep();
        
        // Print time step info if necessary
        if (cf.rank==0){
            if (cf.out_freq > 0)
                if ((t+1) % cf.out_freq == 0)
                    cout << c(cFlag,YEL) << left << setw(15) << "    Iteration:" << c(cFlag,NON) 
                         << c(cFlag,CYA) << right << setw(0) << t+1 << c(cFlag,NON) << "/" << c(cFlag,CYA) << left << setw(0) << cf.tend << c(cFlag,NON) << ", "
                         << c(cFlag,CYA) << right << setw(0) << setprecision(3) << scientific << time << "s" << c(cFlag,NON) << endl;
        }
#ifndef NOMPI
        // Write solution file if necessary
        if (cf.write_freq > 0){
            if ((t+1) % cf.write_freq == 0){
                f->timers["solWrite"].reset();
                w.writeSolution(cf,f->grid,f->var,t+1,time);
                Kokkos::fence();
                f->timers["solWrite"].accumulate();
            }
        }
        // Write restart file if necessary
        if (cf.restart_freq > 0){
            if ((t+1) % cf.restart_freq == 0){
                f->timers["resWrite"].reset();
                w.writeRestart(cf,f->grid,f->var,t+1,time);
                Kokkos::fence();
                f->timers["resWrite"].accumulate();
            }
        }
#endif
        // Print status chevk if necessary
        if (cf.stat_freq > 0){
            if ((t+1) % cf.stat_freq == 0){
                f->timers["statCheck"].reset();
                statusCheck(cFlag,cf,f->var,t+1,time,totalTimer);
                f->timers["statCheck"].accumulate();
            }
        }
    } // End main time loop

    // post simulation hook
    f->postSim();

    //stop simulation timer
    simTimer.stop();

    // notify simulation complete
    if (cf.rank == 0)
        cout << c(cFlag,GRE) << "Simulation Complete!" << c(cFlag,NON) << endl;


    // Sort computer timers
    typedef std::function<bool(std::pair<std::string, fiestaTimer>, std::pair<std::string, fiestaTimer>)> Comparator;
    Comparator compFunctor =
            [](std::pair<std::string, fiestaTimer> elem1 ,std::pair<std::string, fiestaTimer> elem2)
            {
                return elem1.second.get() > elem2.second.get();
            };
    std::set<std::pair<std::string, fiestaTimer>, Comparator> stmr(
            f->timers.begin(), f->timers.end(), compFunctor);
 
    // stop total program timer
    totalTimer.stop();

    // prin timer values
    if (cf.rank == 0){
        cout << endl << "-----------------------" << endl << endl;
        cout.precision(2);
        cout << c(cFlag,GRE) << left  << setw(36) << "Total Time:"     << c(cFlag,NON) 
             << c(cFlag,CYA) << right << setw(13) << totalTimer.getf() << c(cFlag,NON) << endl << endl;

        cout << c(cFlag,GRE) << left  << setw(36) << "  Setup Time:"    << c(cFlag,NON)
             << c(cFlag,CYA) << right << setw(13) << initTimer.getf() << c(cFlag,NON) << endl;
        if(cf.restart == 1) 
            cout << c(cFlag,NON) << left  << setw(36) << "    Restart Read:" << c(cFlag,NON)
                 << c(cFlag,CYA) << right << setw(13) << loadTimer.getf()                    << c(cFlag,NON) << endl << endl;
        else{
            cout << c(cFlag,NON) << left  << setw(36) << "    Initial Condition Generation:" << c(cFlag,NON)
                 << c(cFlag,CYA) << right << setw(13) << loadTimer.getf()                    << c(cFlag,NON) << endl;
            cout << c(cFlag,NON) << left  << setw(36) << "    Grid Generation:" << c(cFlag,NON)
                 << c(cFlag,CYA) << right << setw(13) << gridTimer.getf()                    << c(cFlag,NON) << endl;
            cout << c(cFlag,NON) << left  << setw(36) << "    Initial Consition WriteTime:" << c(cFlag,NON)
                 << c(cFlag,CYA) << right << setw(13) << writeTimer.getf()                    << c(cFlag,NON) << endl << endl;
        }

        cout << c(cFlag,GRE) << left  << setw(36) << "  Simulation Time:" << c(cFlag,NON)
             << c(cFlag,CYA) << right << setw(13) << simTimer.getf() << c(cFlag,NON) << endl;
        for (auto tmr : stmr){
            cout << c(cFlag,NON) << left  << setw(36) << "    "+tmr.second.describe()+":" << c(cFlag,NON)
                 << c(cFlag,CYA) << right << setw(13) << tmr.second.getf()                << c(cFlag,NON) << endl;
        }
        cout << " " << endl;
    }

    
    // finalize mpi
#ifndef NOMPI
    MPI_Finalize();
#endif

    //exit success
    return 0;
}
