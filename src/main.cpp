#include "input.hpp"
#include "mpi.hpp"
#include "cgns.hpp"
#include "bc.hpp"
#include "Kokkos_Core.hpp"
#include <mpi.h>
#include "debug.hpp"
#include "hydroc3d.hpp"
#include "hydro2d.hpp"
#include "rkfunction.hpp"
#include <iostream>
#include <cstdio>
#include <ctime>

void fnExit1(void){
    Kokkos::finalize();
}

int main(int argc, char* argv[]){
    MPI_Init(NULL,NULL);

    int temp_rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&temp_rank);
    if (temp_rank == 0){
        printf("--------------------------------------\n");
        printf("|     _____ _           _            |\n");
        printf("|    |  ___(_) ___  ___| |_ __ _     |\n");
        printf("|    | |_  | |/ _ \\/ __| __/ _` |    |\n");
        printf("|    |  _| | |  __/\\__ \\ || (_| |    |\n");
        printf("|    |_|   |_|\\___||___/\\__\\__,_|    |\n");
        printf("|                                    |\n");
        printf("--------------------------------------\n");
    }
    //if (temp_rank == 0) printf("\n----------  Fiesta ----------\n");
    if (temp_rank == 0) printf("\nInitializing...\n");

    Kokkos::initialize(argc, argv);

    atexit(fnExit1);

    int idx;

    struct inputConfig cf;

    cf = executeConfiguration(argv[1]);

    cf = mpi_init(cf);
    MPI_Barrier(cf.comm);

    int cv = 0;
    if (cf.ceq == 1)
        cv = 5;

    Kokkos::View<double****> myV("myV",cf.ngi,cf.ngj,cf.ngk,cf.nv+cv);
    Kokkos::View<double****> tmp("RK_tmp",cf.ngi,cf.ngj,cf.ngk,cf.nv+cv);
    Kokkos::View<double****> K1("RK_K1",cf.ngi,cf.ngj,cf.ngk,cf.nv+cv);
    Kokkos::View<double****> K2("RK_K2",cf.ngi,cf.ngj,cf.ngk,cf.nv+cv);
        
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
    if (cf.rank == 0){
        printf("%s\n",cf.inputFname);
        if (cf.restart)
            printf("Running from Restart File: %d | %s\n",cf.restart,cf.sfName);
        printf("-----------------------\n");
        printf("Running %d processes as (%d,%d,%d)\n",cf.numProcs,cf.xProcs,cf.yProcs,cf.zProcs);
        printf("nt = %d, dt = %.2e\n",cf.nt,cf.dt);
        printf("Num Cells X = %d, dx = %.2e\n",cf.glbl_nci,cf.dx);
        printf("Num Cells Y = %d, dy = %.2e\n",cf.glbl_ncj,cf.dy);
        if (cf.ndim == 3)
            printf("Num Cells Z = %d, dz = %.2e\n",cf.glbl_nck,cf.dz);
        if (cf.ceq)
            printf("C-Equation Enables");
        else
            printf("C-Equation Disabled\n");
        printf("Number of Species = %d:\n",cf.ns);
        for (int s=0; s<cf.ns; ++s)
            printf("    Species %d, Gamma = %4.2f, M = %6.4f\n",s+1,cf.gamma[s],cf.M[s]);
        printf("-----------------------\n");
    }
    
    MPI_Barrier(cf.comm);
    if (cf.rank == 0) printf("\nLoading Initial Conditions...\n");
    double total_time,mean_time;
    std::clock_t start;
    start = std::clock();
    if (cf.restart == 0)
        loadInitialConditions(cf,myV);
    total_time = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

    if (cf.rank == 0){
        printf("Loaded Initial Conditions in %.1fs\n",total_time);
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
        if (cf.rank == 0) printf("\nReading Restart File...\n");
        readSolution(cf,myV);
    }else{
        if (cf.rank == 0)
            if (cf.write_freq >0 || cf.restart_freq>0)
                printf("\nWriting Initial Conditions...\n");
        if (cf.write_freq >0)
            writeSolution(cf,xSP,ySP,zSP,myV,0,0.00);
        if (cf.restart_freq >0)
            writeRestart(cf,x,y,z,myV,0,0.00);
    }

    /*** Choose Scheme ***/
    rk_func *f;
    if (cf.ndim == 3){
        f = new hydroc3d_func(cf,cd);
    }else{
        f = new hydro2d_func(cf,cd);
    }

    // create mpi buffers
    mpiBuffers m(cf);

    if (cf.rank == 0) printf("\nStarting Simulation...\n");
    MPI_Barrier(cf.comm);
    start = std::clock();
    MPI_Barrier(cf.comm);

    for (int t=tstart; t<cf.nt; ++t){
        time = time + cf.dt;

        /****** Low Storage Runge-Kutta 2nd order ******/
        //K1 = f(myV)
        applyBCs(cf,myV,m);
        f->compute(myV,K1);
        
        //tmp = myV + k1/2
        Kokkos::parallel_for("Loop1", policy_1({0,0,0,0},{cf.ngi, cf.ngj, cf.ngk, cf.nv+cv}),
               KOKKOS_LAMBDA  (const int i, const int j, const int k, const int v) {
            tmp(i,j,k,v) = myV(i,j,k,v) + cf.dt*K1(i,j,k,v)/2;
        });

        //K2 = f(tmp)
        applyBCs(cf,tmp,m);
        f->compute(tmp,K2);

        //myV = myV + K2
        Kokkos::parallel_for("Loop2", policy_1({0,0,0,0},{cf.ngi, cf.ngj, cf.ngk, cf.nv+cv}),
               KOKKOS_LAMBDA  (const int i, const int j, const int k, const int v) {
            myV(i,j,k,v) = myV(i,j,k,v) + cf.dt*K2(i,j,k,v);
        });
        
        
        /****** Output Control ******/
        if (cf.rank==0){
            if (cf.out_freq > 0)
                if ((t+1) % cf.out_freq == 0)
                    printf("Iteration: %d/%d, Sim Time: %.2e\n",t+1,cf.nt,time);
        }
        if (cf.write_freq > 0)
            if ((t+1) % cf.write_freq == 0)
                writeSolution(cf,xSP,ySP,zSP,myV,t+1,time);
        if (cf.restart_freq > 0)
            if ((t+1) % cf.restart_freq == 0)
                writeRestart(cf,x,y,z,myV,t+1,time);
    }

    MPI_Barrier(cf.comm);
    total_time = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    mean_time = total_time/cf.nt;

    if (cf.rank == 0){
        printf("\nSimulation Complete\n");
        printf("Total Time: %.2fs\nMean Time Step: %.2es\n",total_time,mean_time);
    }

    //if (cf.rank==0) printf("\n");
    
    MPI_Finalize();
    //Kokkos::finalize();
    return 0;
}
