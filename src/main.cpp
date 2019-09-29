#include "input.hpp"
#include "mpi_init.hpp"
#include "cgns.hpp"
#include "bc.hpp"
#include "Kokkos_Core.hpp"
#include <mpi.h>
#include "lsdebug.hpp"
#include "weno_function.hpp"

void fnExit1(void){
    Kokkos::finalize();
}

int main(int argc, char* argv[]){
    MPI_Init(NULL,NULL);
    Kokkos::initialize(argc, argv);

    atexit(fnExit1);

    int idx;

    struct inputConfig cf;

    cf = executeConfiguration(argv[1]);

    cf = mpi_init(cf);

    Kokkos::View<double****> myV("testView",cf.ngi,cf.ngj,cf.ngk,cf.nv+4);
    Kokkos::View<double****> tmp("testView",cf.ngi,cf.ngj,cf.ngk,cf.nv+4);
    Kokkos::View<double****> K1("testView",cf.ngi,cf.ngj,cf.ngk,cf.nv+4);
    Kokkos::View<double****> K2("testView",cf.ngi,cf.ngj,cf.ngk,cf.nv+4);
    Kokkos::View<double*> cd("deviceCF",4+cf.ns*2);
    typename Kokkos::View<double*>::HostMirror hostcd = Kokkos::create_mirror_view(cd);
    Kokkos::deep_copy(hostcd, cd);
    hostcd(0) = cf.ns;
    hostcd(1) = cf.dx;
    hostcd(2) = cf.dy;
    hostcd(3) = cf.dz;
    
    int sdx = 4;
    for (int s=0; s<cf.ns; ++s){
        hostcd(sdx) = cf.gamma[s];
        hostcd(sdx+1) = cf.R/cf.M[s];
        sdx += 2;
    }
    Kokkos::deep_copy(cd,hostcd);

    if (cf.restart == 0)
        loadInitialConditions(cf,myV);
    
    if (cf.rank == 0){
        printf("loboSHOK - %s\n",cf.inputFname);
        printf("Restart File: %d | %s\n",cf.restart,cf.sfName);
        printf("-----------------------\n");
        printf("Running %d processes as (%d,%d,%d)\n",cf.numProcs,cf.xProcs,cf.yProcs,cf.zProcs);
        printf("nt = %d, dt = %f\n",cf.nt,cf.dt);
        printf("glbl_ni = %d, dx = %f\n",cf.glbl_ni,cf.dx);
        printf("glbl_nj = %d, dy = %f\n",cf.glbl_nj,cf.dy);
        printf("glbl_nk = %d, dz = %f\n",cf.glbl_nk,cf.dz);
        printf("Number of Species = %d:\n",cf.ns);
        for (int s=0; s<cf.ns; ++s)
            printf("    Species %d, Gamma = %4.2f, M = %6.4f\n",s+1,cf.gamma[s],cf.M[s]);
        printf("-----------------------\n");
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
    
    //typedef typename Kokkos::View<double****> ViewType;
    typedef Kokkos::MDRangePolicy<Kokkos::Rank<4>> policy_1;
    //Kokkos::MDRangePolicy<Kokkos::Rank<3>> policy_1({0,0,0},{cf.nci, cf.ncj, cf.nck});
    
    double time = cf.time;
    int tstart = cf.tstart;
    
    MPI_Barrier(cf.comm);

    if (cf.restart == 1){
        readSolution(cf,myV);
    }else{
        writeSolution(cf,xSP,ySP,zSP,myV,0,0.00);
        writeRestart(cf,x,y,z,myV,0,0.00);
    }

    for (int t=tstart; t<cf.nt; ++t){
        time = time + cf.dt;

//        applyBCs(cf,myV);
//        rhs_func f1(cf.dt,cf,myV,K1, cd);
//        f1();
//
//        //myV = myV + K2
//        Kokkos::parallel_for("Loop2", policy_1({cf.ng,cf.ng,cf.ng,0},{cf.ngi-cf.ng,cf.ngj-cf.ng,cf.ngk-cf.ng,cf.nv}),
//               KOKKOS_LAMBDA __device__ (const int i, const int j, const int k, const int v) {
//            myV(i,j,k,v) = myV(i,j,k,v) + cf.dt*K1(i,j,k,v);
//        });

        //tmp = myV
        Kokkos::deep_copy(tmp,myV);

        //K1 = dt*f(tmp) dt is member data to f()
        applyBCs(cf,tmp);
        weno_func f1(cf,tmp,K1, cd);
        f1();
        
        //tmp = myV + k1/2
        Kokkos::parallel_for("Loop1", policy_1({0,0,0,0},{cf.ngi, cf.ngj, cf.ngk, cf.nv+4}),
               KOKKOS_LAMBDA __device__ (const int i, const int j, const int k, const int v) {
            tmp(i,j,k,v) = myV(i,j,k,v) + cf.dt*K1(i,j,k,v)/2;
        });

        //K2 = dt*f(tmp) dt is member data to f()
        applyBCs(cf,tmp);
        weno_func f2(cf,tmp,K2, cd);
        f2();

        //myV = myV + K2
        Kokkos::parallel_for("Loop2", policy_1({0,0,0,0},{cf.ngi, cf.ngj, cf.ngk, cf.nv+4}),
               KOKKOS_LAMBDA __device__ (const int i, const int j, const int k, const int v) {
            myV(i,j,k,v) = myV(i,j,k,v) + cf.dt*K2(i,j,k,v);
        });
        
        if (cf.rank==0){
            if ((t+1) % cf.out_freq == 0)
                printf("%d/%d, %f\n",t+1,cf.nt,time);
        }
        if ((t+1) % cf.write_freq == 0)
            writeSolution(cf,xSP,ySP,zSP,myV,t+1,time);
        if ((t+1) % cf.restart_freq == 0)
            writeRestart(cf,x,y,z,myV,t+1,time);
    }

    //if (cf.rank==0) printf("\n");
    //MPI_Barrier(cf.comm);
    //MYDBGR
    
    MPI_Finalize();
    //Kokkos::finalize();
    return 0;
}
