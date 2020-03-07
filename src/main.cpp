#include "fiesta.hpp"
#include "input.hpp"
#include "mpi.hpp"
#include "cgns.hpp"
#include "bc.hpp"
#include "Kokkos_Core.hpp"
#include <mpi.h>
#include "debug.hpp"
//#include "hydroc3d.hpp"
#include "hydro2d.hpp"
#include "hydro2dvisc.hpp"
#include "rkfunction.hpp"
#include <iostream>
#include <cstdio>
#include <ctime>
#include "output.hpp"
#include <iomanip>
#include "timer.hpp"

using namespace std;


void fnExit1(void){
    Kokkos::finalize();
}

int main(int argc, char* argv[]){
    MPI_Init(NULL,NULL);

    int temp_rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&temp_rank);
    if (temp_rank == 0)
        printSplash();

    Kokkos::initialize(argc, argv);
    
    fiestaTimer totalTimer("Total Runtime");
    fiestaTimer initTimer("Initialization Timer");
    fiestaTimer solWriteTimer("Solution Write Timer");
    fiestaTimer simTimer("Simulation Timer");
    fiestaTimer loadTimer("Initial Condition Generate Timer");
    fiestaTimer readTimer("Restart Read Timer");
    fiestaTimer resWriteTimer("Restart Write Timer");


    atexit(fnExit1);

    int idx;

    struct inputConfig cf;

    cf = executeConfiguration(argc,argv);

    cf = mpi_init(cf);
    MPI_Barrier(cf.comm);

    int cv = 0;
    if (cf.ceq == 1)
        cv = 5;
        
    Kokkos::View<double*> cd("deviceCF",5+cf.ns*2);
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
        sdx += 2;
    }
    Kokkos::deep_copy(cd,hostcd);

    MPI_Barrier(cf.comm);
    /*** Output runtime information ***/
    if (cf.rank == 0)
        printConfig(cf);

    /*** Choose Scheme ***/
    rk_func *f;
    //if (cf.ndim == 3){
    //    f = new hydroc3d_func(cf,cd);
    //}else{
        if (cf.visc == 1)
            f = new hydro2dvisc_func(cf,cd);
        else
            f = new hydro2d_func(cf,cd);
    //}
    
    MPI_Barrier(cf.comm);
    if (cf.restart == 0){
        if (cf.rank == 0) printf("\nGenerating Initial Conditions...\n");
        loadTimer.reset();
        loadInitialConditions(cf,f->var);
        loadTimer.stop();
    }

    if (cf.rank == 0){
        printf("    Generating Initial Conditions in %.1fs\n",loadTimer.get());
    }

    /* allocate grid coordinate and flow variables */
    double *x = (double*)malloc(cf.ni*cf.nj*cf.nk*sizeof(double));
    double *y = (double*)malloc(cf.ni*cf.nj*cf.nk*sizeof(double));
    double *z = (double*)malloc(cf.ni*cf.nj*cf.nk*sizeof(double));

    float *xSP = (float*)malloc(cf.ni*cf.nj*cf.nk*sizeof(float));
    float *ySP = (float*)malloc(cf.ni*cf.nj*cf.nk*sizeof(float));
    float *zSP = (float*)malloc(cf.ni*cf.nj*cf.nk*sizeof(float));
    

    /* calculate rank local grid coordinates */
    for (int k=0; k<cf.nk; ++k){
        for (int j=0; j<cf.nj; ++j){
            for (int i=0; i<cf.ni; ++i){
                idx = (cf.ni*cf.nj)*k+cf.ni*j+i;
                x[idx] = (cf.iStart + i)*cf.dx;
                y[idx] = (cf.jStart + j)*cf.dy;
                z[idx] = (cf.kStart + k)*cf.dz;

                xSP[idx] = (float)x[idx];
                ySP[idx] = (float)y[idx];
                zSP[idx] = (float)z[idx];
            }
        }
    }
    
    typedef Kokkos::MDRangePolicy<Kokkos::Rank<4>> policy_1;
    
    double time = cf.time;
    int tstart = cf.tstart;
    
    MPI_Barrier(cf.comm);

    /*** Read Restart or Write initial conditions ***/
    if (cf.restart == 1){
        if (cf.rank == 0) printf("\nLoading Restart File...\n");
        readTimer.reset();
        readSolution(cf,f->var);
        readTimer.stop();
        cout << "    Loaded restart file in " << setprecision(2) << readTimer.get() << "s" << endl;
    }else{
        if (cf.rank == 0)
            if (cf.write_freq >0 || cf.restart_freq>0)
                printf("\nWriting Initial Conditions...\n");
        if (cf.write_freq >0){
            solWriteTimer.reset();
            writeSolution(cf,xSP,ySP,zSP,f->var,0,0.00);
            solWriteTimer.accumulate();
        }
        if (cf.restart_freq >0){
            resWriteTimer.reset();
            writeRestart(cf,x,y,z,f->var,0,0.00);
            resWriteTimer.accumulate();
        }
        if (cf.rank == 0)
            if (cf.write_freq >0 || cf.restart_freq>0)
                cout << "    Wrote initial conditions in "
                     << solWriteTimer.get() + resWriteTimer.get()
                     << "s" << endl;
    }


    // create mpi buffers
    mpiBuffers m(cf);

    if (cf.rank == 0){
        cout << endl << "-----------------------" << endl << endl;
        cout << "Starting Simulation..." << endl;
    }
    MPI_Barrier(cf.comm);
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
        Kokkos::parallel_for("Loop1", policy_1({0,0,0,0},{cf.ngi, cf.ngj, cf.ngk, cf.nv+cv}),
               KOKKOS_LAMBDA  (const int i, const int j, const int k, const int v) {
            mytmp(i,j,k,v) = myvar(i,j,k,v);
            myvar(i,j,k,v) = myvar(i,j,k,v) + 0.5*cf.dt*mydvar(i,j,k,v);
        });

        applyBCs(cf,f->var,m);
        f->compute();

        mytmp = f->tmp1;
        myvar = f->var;
        mydvar = f->dvar;
        Kokkos::parallel_for("Loop2", policy_1({0,0,0,0},{cf.ngi, cf.ngj, cf.ngk, cf.nv+cv}),
               KOKKOS_LAMBDA  (const int i, const int j, const int k, const int v) {
            myvar(i,j,k,v) = myvar(i,j,k,v) + cf.dt*mydvar(i,j,k,v);
        });
        
        /****** Output Control ******/
        if (cf.rank==0){
            if (cf.out_freq > 0)
                if ((t+1) % cf.out_freq == 0)
                    printf("    Iteration: %d/%d, Sim Time: %.2e\n",t+1,cf.tend,time);
        }
        if (cf.write_freq > 0){
            if ((t+1) % cf.write_freq == 0){
                solWriteTimer.reset();
                writeSolution(cf,xSP,ySP,zSP,f->var,t+1,time);
                solWriteTimer.accumulate();
            }
        }
        if (cf.restart_freq > 0){
            if ((t+1) % cf.restart_freq == 0){
                resWriteTimer.reset();
                writeRestart(cf,x,y,z,f->var,t+1,time);
                resWriteTimer.accumulate();
            }
        }
    }
    simTimer.stop();
    if (cf.rank == 0)
        cout << "Simulation Complete" << endl;

    MPI_Barrier(cf.comm);

    totalTimer.stop();
    if (cf.rank == 0){
        cout << endl << "-----------------------" << endl << endl;
        cout.precision(2);
        cout << "Total Time: " << totalTimer.getf() << endl;
        cout << "Setup Time: " << initTimer.get() << endl;
        cout << "Sim Time: " << simTimer.get() << endl;
        cout << "    Solution write Time: " << solWriteTimer.get() << endl;
        cout << "    Restart write Time: " << resWriteTimer.get() << endl;
    }

    
    MPI_Finalize();
    //Kokkos::finalize();
    return 0;
}
