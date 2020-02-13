#include "input.hpp"
#include "mpi.hpp"
#include "cgns.hpp"
#include "bc.hpp"
#include "Kokkos_Core.hpp"
#include <mpi.h>
#include "lsdebug.hpp"
#include "wenoc3d.hpp"
#include "weno2d.hpp"
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
    if (temp_rank == 0) printf("\n----------  Fiesta ----------\n");
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

    int ncells = cf.ngi * cf.ngj * cf.ngk; // missnamed, not global, just overal size
    //int nxborder;
    //if (cf.ndim == 2)
    //    nborder = cf.ngi * cf.ngj - ((cf.ngi - 2) * (cf.ngj - 2));
    //else
    //    nborder = cf.ngi * cf.ngj * cf.ngk - ((cf.ngi - 2) * (cf.ngj - 2) * (cf.ngk - 2));
    //int nxface = cf.ngi - 1;
    //int nyface = cf.ngj - 1;
    //if (cf.ngk > 1)
    //    int nzface = cf.ngk - 1;
    //else
    //    int nzface = 1;

    //if (cf.rank == 0) 
    //    printf("%d %d %d %d\n", cf.ngi, cf.ngj, cf.ngk, ncells);

    Kokkos::View<double**> myV("myV",cf.ngi*cf.ngj*cf.ngk,cf.nv+cv);
    Kokkos::View<double**> tmp("RK_tmp",cf.ngi*cf.ngj*cf.ngk,cf.nv+cv);
    Kokkos::View<double**> K1("RK_K1",cf.ngi*cf.ngj*cf.ngk,cf.nv+cv);
    Kokkos::View<double**> K2("RK_K2",cf.ngi*cf.ngj*cf.ngk,cf.nv+cv);
        
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

    //printf("%d %d %d\n", cf.ngi, cf.ngj, cf.ngk);

    // Quick preliminary mesh mapping
    Kokkos::View<int*> nlft("nlft",cf.ngi * cf.ngj * cf.ngk);
    Kokkos::View<int*> nrht("nrht",cf.ngi * cf.ngj * cf.ngk);
    Kokkos::View<int*> nbot("nbot",cf.ngi * cf.ngj * cf.ngk);
    Kokkos::View<int*> ntop("ntop",cf.ngi * cf.ngj * cf.ngk);
    Kokkos::View<int*> nbak("nbak",cf.ngi * cf.ngj * cf.ngk);
    Kokkos::View<int*> nfrt("nfrt",cf.ngi * cf.ngj * cf.ngk);
    Kokkos::View<int*>::HostMirror h_nlft = Kokkos::create_mirror_view(nlft);
    Kokkos::View<int*>::HostMirror h_nrht = Kokkos::create_mirror_view(nrht);
    Kokkos::View<int*>::HostMirror h_nbot = Kokkos::create_mirror_view(nbot);
    Kokkos::View<int*>::HostMirror h_ntop = Kokkos::create_mirror_view(ntop);
    Kokkos::View<int*>::HostMirror h_nbak = Kokkos::create_mirror_view(nbak);
    Kokkos::View<int*>::HostMirror h_nfrt = Kokkos::create_mirror_view(nfrt);

    //Kokkos::View<int*> realcell("realcell",cf.ngi * cf.ngj * cf.ngk);
    //Kokkos::View<int*>::HostMirror h_realcell = Kokkos::create_mirror_view(realcell);
    //Kokkos::View<int*> facecell("facecell",cf.ngi * cf.ngj * cf.ngk);
    //Kokkos::View<int*>::HostMirror h_facecell = Kokkos::create_mirror_view(facecell);
    //Kokkos::View<int*>::HostMirror h_bordcell = Kokkos::create_mirror_view(bordcell);
    Kokkos::View<int*> celltype("celltype",cf.ngi * cf.ngj * cf.ngk);
    //Kokkos::View<int*> bordcell("bordcell",cf.ngi * cf.ngj * cf.ngk);
    //Kokkos::View<int*> corntype("corntype",cf.ngi * cf.ngj * cf.ngk);
    Kokkos::View<int*>::HostMirror h_celltype = Kokkos::create_mirror_view(celltype);
    //Kokkos::View<int*>::HostMirror h_bordcell = Kokkos::create_mirror_view(bordcell);
    //Kokkos::View<int*>::HostMirror h_corncell = Kokkos::create_mirror_view(corncell);

    /*
    Kokkos::View<int*> xface2cell_lower("xface2cell_lower",nxface);
    Kokkos::View<int*> xface2cell_upper("xface2cell_upper",nxface);
    Kokkos::View<int*> yface2cell_lower("yface2cell_lower",nyface);
    Kokkos::View<int*> yface2cell_upper("yface2cell_upper",nyface);
    Kokkos::View<int*> zface2cell_lower("zface2cell_lower",nzface);
    Kokkos::View<int*> zface2cell_upper("zface2cell_upper",nzface);
    Kokkos::View<int*>::HostMirror h_xface2cell_lower = Kokkos::create_mirror_view(xface2cell_lower);
    Kokkos::View<int*>::HostMirror h_xface2cell_upper = Kokkos::create_mirror_view(xface2cell_upper);
    Kokkos::View<int*>::HostMirror h_yface2cell_lower = Kokkos::create_mirror_view(yface2cell_lower);
    Kokkos::View<int*>::HostMirror h_yface2cell_upper = Kokkos::create_mirror_view(yface2cell_upper);
    Kokkos::View<int*>::HostMirror h_zface2cell_lower = Kokkos::create_mirror_view(zface2cell_lower);
    Kokkos::View<int*>::HostMirror h_zface2cell_upper = Kokkos::create_mirror_view(zface2cell_upper);
    */
    
    for (int k = 0; k < cf.ngk; k++) {
        for (int j = 0; j < cf.ngj; j++) {
            for (int i = 0; i < cf.ngi; i++) {
                int loc_idx = k * cf.ngj * cf.ngi + j * cf.ngi + i;
                // for i direction
                if (i > 0)
                    h_nlft(loc_idx) = k * cf.ngj * cf.ngi + j * cf.ngi + i - 1;
                else
                    h_nlft(loc_idx) = loc_idx;
                if (i < cf.ngi-1)
                    h_nrht(loc_idx) = k * cf.ngj * cf.ngi + j * cf.ngi + i + 1;
                else
                    h_nrht(loc_idx) = loc_idx;
                // for j direction
                if (j > 0)
                    h_nbot(loc_idx) = k * cf.ngj * cf.ngi + (j - 1) * cf.ngi + i;
                else
                    h_nbot(loc_idx) = loc_idx;
                if (j < cf.ngj-1)
                    h_ntop(loc_idx) = k * cf.ngj * cf.ngi + (j + 1) * cf.ngi + i;
                else
                    h_ntop(loc_idx) = loc_idx;
                // for k direction
                if (k > 0)
                    h_nbak(loc_idx) = (k - 1) * cf.ngj * cf.ngi + j * cf.ngi + i;
                else
                    h_nbak(loc_idx) = loc_idx;
                if (k < cf.ngk-1)
                    h_nfrt(loc_idx) = (k + 1) * cf.ngj * cf.ngi + j * cf.ngi + i;
                else
                    h_nfrt(loc_idx) = loc_idx;
                // for boundary/halo cell
                if (i < cf.ng || i >= cf.ngi - cf.ng || j < cf.ng || j >= cf.ngj - cf.ng) { // || k < cf.ng || k >= cf.ngk - cf.ng)
                    h_celltype(loc_idx) = BOUNDARY_CELL;
                }
                // for border cells (cells touching ghost/halo cells)
                else if (i == cf.ng || i == cf.ngi - cf.ng - 1 || j == cf.ng || j == cf.ngj - cf.ng - 1) { // || k == cf.ng || k == cf.ngk - cf.ng - 1) {
                    h_celltype(loc_idx) = BORDER_CELL;
                }
                // for real cells
                else {
                    h_celltype(loc_idx) = REAL_CELL;
                }
            }
        }
    }
            //printf("%d %d %d %d %d %d %d\n", cf.rank, cf.xMinus, cf.xPlus, cf.yMinus, cf.yPlus, cf.zMinus, cf.zPlus);

    //for (int n = 0; n < ncells; n++) {
    //    if (h_celltype(n) == BORDER_CELL && h_celltype(h_ntop(n)) == BOUNDARY_CELL) printf("%d\n", n);
    //}


    Kokkos::deep_copy(nlft, h_nlft);
    Kokkos::deep_copy(nrht, h_nrht);
    Kokkos::deep_copy(nbot, h_nbot);
    Kokkos::deep_copy(ntop, h_ntop);
    Kokkos::deep_copy(nbak, h_nbak);
    Kokkos::deep_copy(nfrt, h_nfrt);
    //Kokkos::deep_copy(realcell, h_realcell);
    //Kokkos::deep_copy(facecell, h_facecell);
    Kokkos::deep_copy(celltype, h_celltype);
    //Kokkos::deep_copy(bordcell, h_bordcell);

    //printf("Here1\n");

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

    //printf("Here2\n");
    
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
    
    //printf("Here3\n");
    

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
        /*
        if (cf.write_freq >0)
            writeSolution(cf,xSP,ySP,zSP,myV,0,0.00);
        if (cf.restart_freq >0)
            writeRestart(cf,x,y,z,myV,0,0.00);
        */
    }

    //printf("Here4\n");
    
    /*** Choose Scheme ***/
    rk_func *f;
    if (cf.ndim == 3){
        //f = new wenoc3d_func(cf,cd);
    }else{
        f = new weno2d_func(cf,cd);
    }

    //printf("Here5\n");
    
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
        applyBCs(cf,myV,m,nlft,nrht,nbot,ntop,nbak,nfrt,celltype,ncells);
        //if (t == 0) printf("Here6\n");
    
        f->compute(myV,nlft,nrht,nbot,ntop,nbak,nfrt,celltype,K1,ncells);
        //if (t == 0) printf("Here7\n");
    
        
        //tmp = myV + k1/2
        Kokkos::parallel_for("Loop1", policy_1({0,0,0,0},{cf.ngi, cf.ngj, cf.ngk, cf.nv+cv}),
               KOKKOS_LAMBDA  (const int i, const int j, const int k, const int v) {
            tmp(i+j*cf.ngi+k*cf.ngj*cf.ngi,v) = myV(i+j*cf.ngi+k*cf.ngj*cf.ngi,v) + cf.dt*K1(i+j*cf.ngi+k*cf.ngj*cf.ngi,v)/2;
        });

        //K2 = f(tmp)
        applyBCs(cf,tmp,m,nlft,nrht,nbot,ntop,nbak,nfrt,celltype,ncells);
        f->compute(tmp,nlft,nrht,nbot,ntop,nbak,nfrt,celltype,K2,ncells);

        //myV = myV + K2
        Kokkos::parallel_for("Loop2", policy_1({0,0,0,0},{cf.ngi, cf.ngj, cf.ngk, cf.nv+cv}),
               KOKKOS_LAMBDA  (const int i, const int j, const int k, const int v) {
            myV(i+j*cf.ngi+k*cf.ngj*cf.ngi,v) = myV(i+j*cf.ngi+k*cf.ngj*cf.ngi,v) + cf.dt*K2(i+j*cf.ngi+k*cf.ngj*cf.ngi,v);
        });
        
        //if (t == 0) printf("Here8\n");
    
        
        /****** Output Control ******/
        if (cf.rank==0){
            if (cf.out_freq > 0)
                if ((t+1) % cf.out_freq == 0)
                    printf("Iteration: %d/%d, Sim Time: %.2e\n",t+1,cf.nt,time);
        }
        /*
        if (cf.write_freq > 0)
            if ((t+1) % cf.write_freq == 0)
                writeSolution(cf,xSP,ySP,zSP,myV,t+1,time);
        if (cf.restart_freq > 0)
            if ((t+1) % cf.restart_freq == 0)
                writeRestart(cf,x,y,z,myV,t+1,time);
        */
    }

    MPI_Barrier(cf.comm);
    total_time = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    mean_time = total_time/cf.nt;

    
    Kokkos::View<double**>::HostMirror h_myV = Kokkos::create_mirror_view(myV);
    Kokkos::deep_copy(h_myV, myV);

    if (cf.rank == 0) {
        for (int var = 0; var < cf.nv + cv; var++) {
            for (int kk = 0; kk < cf.ngk; kk++) {
                for (int jj = cf.ng; jj < cf.ngj - cf.ng; jj++) {
                    for (int ii = cf.ng; ii < cf.ngi - cf.ng; ii++) {
                        int loc_idx = ii + jj * cf.ngi + kk * cf.ngi * cf.ngj;
                        printf("%f\n", h_myV(loc_idx,var));
                    }
                }
            }
        }
    }
        

    if (cf.rank == 0){
        printf("\nSimulation Complete\n");
        printf("Total Time: %.2fs\nMean Time Step: %.2es\n",total_time,mean_time);
    }

    //if (cf.rank==0) printf("\n");
    
    MPI_Finalize();
    //Kokkos::finalize();
    return 0;
}
