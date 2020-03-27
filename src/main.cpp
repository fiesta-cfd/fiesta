#include "fiesta.hpp"
#include "input.hpp"
#include "mpi.hpp"
#include "cgns.hpp"
#include "bc.hpp"
#include "Kokkos_Core.hpp"
#include <mpi.h>
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
    // INITIALIZE
    MPI_Init(NULL,NULL);

    int temp_rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&temp_rank);
    if (temp_rank == 0)
        printSplash();

    Kokkos::initialize(argc, argv);
    
    fiestaTimer totalTimer;
    fiestaTimer initTimer;
    fiestaTimer solWriteTimer;
    fiestaTimer simTimer;
    fiestaTimer loadTimer;
    fiestaTimer resWriteTimer;
    fiestaTimer gridTimer;
    fiestaTimer writeTimer;


    atexit(fnExit1);

    struct inputConfig cf;

    // CONFIGURE

    cf = executeConfiguration(argc,argv);

    cf = mpi_init(cf);
    MPI_Barrier(cf.comm);

    //int cv = 0;
    //if (cf.ceq == 1)
    //    cv = 5;
        
    Kokkos::View<double*> cd("deviceCF",5+cf.ns*3);
    typename Kokkos::View<double*>::HostMirror hostcd = Kokkos::create_mirror_view(cd);
    Kokkos::deep_copy(hostcd, cd);
    hostcd(0) = cf.ns;
    hostcd(1) = cf.dx;
    hostcd(2) = cf.dy;
    hostcd(3) = cf.dz;
    hostcd(4) = cf.nv;
    
    int sdx = 5;
    for (int s=0; s<cf.ns; ++s){
        hostcd(sdx) = cf.gamma[s];
        hostcd(sdx+1) = cf.R/cf.M[s];
        hostcd(sdx+2) = cf.mu[s];
        sdx += 3;
    }
    Kokkos::deep_copy(cd,hostcd);

    MPI_Barrier(cf.comm);
    /*** Output runtime information ***/
    if (cf.rank == 0)
        printConfig(cf);

    /*** Choose Scheme ***/
    rk_func *f;
    if (cf.ndim == 3){
        f = new hydroc3d_func(cf,cd);
    }else{
        f = new hydro2dvisc_func(cf,cd);
    }

    cgnsWriter w(cf,f->grid,f->var);

    MPI_Barrier(cf.comm);
    if (cf.restart == 0){
        if (cf.rank == 0) cout << c(GRE) << "Generating Initial Conditions:" << c(NON) << endl;
        loadTimer.start();
        loadInitialConditions(cf,f->var);
        loadTimer.stop();
        if (cf.rank == 0) cout << "    Generated in: " << c(CYA) << loadTimer.getf() << c(NON) << endl << endl;

        if (cf.rank == 0) cout << c(GRE) << "Generating Grid:" << c(NON) << endl;
        gridTimer.start();
        loadGrid(cf,f->grid);
        gridTimer.stop();
        if (cf.rank == 0) cout << "    Generated in: " << c(CYA) << gridTimer.getf() << c(NON) << endl << endl;
    }

    typedef Kokkos::MDRangePolicy<Kokkos::Rank<4>> policy_1;
    
    double time = cf.time;
    int tstart = cf.tstart;
    
    MPI_Barrier(cf.comm);

    /*** Read Restart or Write initial conditions ***/
    if (cf.restart == 1){
        if (cf.rank == 0) cout << c(GRE) << "Loading Restart File:" << c(NON) << endl;
        loadTimer.reset();
        w.readSolution(cf,f->grid,f->var);
        loadTimer.stop();
        cout << "    Loaded in: " << setprecision(2) << c(CYA) << loadTimer.get() << "s" << c(NON) << endl;
    }else{
        if (cf.rank == 0)
            if (cf.write_freq >0 || cf.restart_freq>0)
                cout << c(GRE) << "Writing Initial Conditions:" << c(NON) << endl;
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
                cout << "    Wrote in: " << c(CYA) << writeTimer.getf() << c(NON) << endl;
    }


    // create mpi buffers
    mpiBuffers m(cf);

    if (cf.rank == 0){
        cout << endl << "-----------------------" << endl << endl;
        cout << c(GRE) << "Starting Simulation:" << c(NON) << endl;
    }

    MPI_Barrier(cf.comm);

    initTimer.stop();
    simTimer.reset();

    for (int t=tstart; t<cf.tend; ++t){
        time = time + cf.dt;

        /****** Low Storage Runge-Kutta 2nd order ******/
        //K1 = f(myV)
        applyBCs(cf,f->var,m);
        f->compute();
        //f->compute(myV,K1);
        
        //tmp = myV + k1/2
        FS4D mytmp = f->tmp1;
        FS4D myvar = f->var;
        FS4D mydvar = f->dvar;
        Kokkos::parallel_for("Loop1", policy_1({0,0,0,0},{cf.ngi, cf.ngj, cf.ngk, cf.nvt}),
               KOKKOS_LAMBDA  (const int i, const int j, const int k, const int v) {
            mytmp(i,j,k,v) = myvar(i,j,k,v);
            myvar(i,j,k,v) = myvar(i,j,k,v) + 0.5*cf.dt*mydvar(i,j,k,v);
        });

        applyBCs(cf,f->var,m);
        f->compute();

        mytmp = f->tmp1;
        myvar = f->var;
        mydvar = f->dvar;
        Kokkos::parallel_for("Loop2", policy_1({0,0,0,0},{cf.ngi, cf.ngj, cf.ngk, cf.nvt}),
               KOKKOS_LAMBDA  (const int i, const int j, const int k, const int v) {
            myvar(i,j,k,v) = myvar(i,j,k,v) + cf.dt*mydvar(i,j,k,v);
        });
        
        /****** Output Control ******/
        if (cf.rank==0){
            if (cf.out_freq > 0)
                if ((t+1) % cf.out_freq == 0)
                    cout << c(YEL) << left << setw(15) << "    Iteration:" << c(NON) 
                         << c(CYA) << right << setw(0) << t+1 << c(NON) << "/" << c(CYA) << left << setw(0) << cf.tend << c(NON) << ", "
                         << c(CYA) << right << setw(0) << setprecision(3) << scientific << time << "s" << c(NON) << endl;
                    //printf("    Iteration: %d/%d, Sim Time: %.2e\n",t+1,cf.tend,time);
        }
        if (cf.write_freq > 0){
            if ((t+1) % cf.write_freq == 0){
                f->timers["solWrite"].reset();
                w.writeSolution(cf,f->grid,f->var,t+1,time);
                f->timers["solWrite"].accumulate();
            }
        }
        if (cf.restart_freq > 0){
            if ((t+1) % cf.restart_freq == 0){
                f->timers["resWrite"].reset();
                w.writeRestart(cf,f->grid,f->var,t+1,time);
                f->timers["resWrite"].accumulate();
            }
        }
        if (cf.stat_freq > 0){
            if ((t+1) % cf.stat_freq == 0){
                f->timers["statCheck"].reset();
                statusCheck(cf,f->var,t+1,time,totalTimer);
                f->timers["statCheck"].accumulate();
            }
        }
    }
    simTimer.stop();
    if (cf.rank == 0)
        cout << c(GRE) << "Simulation Complete!" << c(NON) << endl;

    MPI_Barrier(cf.comm);

    typedef std::function<bool(std::pair<std::string, fiestaTimer>, std::pair<std::string, fiestaTimer>)> Comparator;
 
    Comparator compFunctor =
            [](std::pair<std::string, fiestaTimer> elem1 ,std::pair<std::string, fiestaTimer> elem2)
            {
                return elem1.second.get() > elem2.second.get();
            };
 
    std::set<std::pair<std::string, fiestaTimer>, Comparator> stmr(
            f->timers.begin(), f->timers.end(), compFunctor);
 
    totalTimer.stop();
    if (cf.rank == 0){
        cout << endl << "-----------------------" << endl << endl;
        cout.precision(2);
        cout << c(GRE) << left  << setw(36) << "Total Time:"     << c(NON) 
             << c(CYA) << right << setw(13) << totalTimer.getf() << c(NON) << endl << endl;

        cout << c(GRE) << left  << setw(36) << "  Setup Time:"    << c(NON)
             << c(CYA) << right << setw(13) << initTimer.getf() << c(NON) << endl;
        if(cf.restart == 1) 
            cout << c(NON) << left  << setw(36) << "    Restart Read:" << c(NON)
                 << c(CYA) << right << setw(13) << loadTimer.getf()                    << c(NON) << endl << endl;
        else{
            cout << c(NON) << left  << setw(36) << "    Initial Condition Generation:" << c(NON)
                 << c(CYA) << right << setw(13) << loadTimer.getf()                    << c(NON) << endl;
            cout << c(NON) << left  << setw(36) << "    Grid Generation:" << c(NON)
                 << c(CYA) << right << setw(13) << gridTimer.getf()                    << c(NON) << endl;
            cout << c(NON) << left  << setw(36) << "    Initial Consition WriteTime:" << c(NON)
                 << c(CYA) << right << setw(13) << writeTimer.getf()                    << c(NON) << endl << endl;
        }

        //cout << c(GRE) << left  << setw(36) << "Write Times:"  << c(NON)
        //     << c(CYA) << right << setw(8)  << simTimer.getf() << c(NON) << endl;
        //cout << c(NON) << left  << setw(36) << "    Solution write Time:" << c(NON)
        //     << c(CYA) << right << setw(8)  << solWriteTimer.getf()   << c(NON) << endl;
        //cout << c(NON) << left  << setw(36) << "    Restart write Time:" << c(NON)
        //     << c(CYA) << right << setw(8)  << resWriteTimer.getf()      << c(NON) << endl << endl;

        cout << c(GRE) << left  << setw(36) << "  Simulation Time:" << c(NON)
             << c(CYA) << right << setw(13) << simTimer.getf() << c(NON) << endl;
        for (auto tmr : stmr){
            cout << c(NON) << left  << setw(36) << "    "+tmr.second.describe()+":" << c(NON)
                 << c(CYA) << right << setw(13) << tmr.second.getf()                << c(NON) << endl;
        }
        cout << " " << endl;
    }

    
    MPI_Finalize();
    //Kokkos::finalize();
    return 0;
}
