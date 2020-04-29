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
#include "l7/l7.h"
#include "mesh/mesh.h"
#include "mesh/partition.h"
#ifdef HAVE_GRAPHICS
#include <unistd.h>
#include "graphics/display.h"
#endif
#include "MallocPlus/MallocPlus.h"

#define HAVE_MPI
#define DYNAMIC_REZONE

enum partition_method initial_order,  //  Initial order of mesh.
                      cycle_reorder;  //  Order of mesh every cycle.
bool localStencil;

static Mesh *mesh;           //  Object containing mesh information
MallocPlus state_memory;     // MallocPlus state array
double *X = NULL;
double *Y = NULL;
double *E = NULL;
double *D = NULL;

int do_quo_setup, lttrace_on;

void fnExit1(void){
    Kokkos::finalize();
}

void memory_reset_ptrs(void){
   X = (double *)state_memory.get_memory_ptr("X");
   Y = (double *)state_memory.get_memory_ptr("Y");
   E = (double *)state_memory.get_memory_ptr("E");
   D = (double *)state_memory.get_memory_ptr("D");
}

void state_reorder(vector<int> iorder) {
    X = state_memory.memory_reorder(X, &iorder[0]);
    Y = state_memory.memory_reorder(Y, &iorder[0]);
    E = state_memory.memory_reorder(E, &iorder[0]);
    D = state_memory.memory_reorder(D, &iorder[0]);
}

int main(int argc, char* argv[]){

    // Needed for code to compile correctly on the Mac
    int mype=0;
    int numpe=0;


    //MPI_Init(NULL,NULL);

    L7_Init(&mype, &numpe, &argc, argv, do_quo_setup, lttrace_on);

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

    int nx = cf.ngi - cf.ng*2;
    int ny = cf.ngj - cf.ng*2;
    int levmx = 2;
    mesh = new Mesh(nx, ny, levmx, 2, 1.0, 1.0, 1, 1, 0);
    mesh->init(nx, ny, 0);
    size_t &ncells_loc   = mesh->ncells;
    size_t new_ncells   = 0;
    printf("ncells %d\n", ncells_loc);

    mesh->calc_celltype(ncells_loc);
    mesh->calc_neighbors_local();
    for (int check = 0; check < mesh->ncells; check++) {
        //printf("%d) %d\n", check, mesh->celltype[check]);
    }


    Kokkos::View<double**> myV("myV",ncells_loc,cf.nv+cv);
    Kokkos::View<double**> tmp("RK_tmp",ncells_loc,cf.nv+cv);
    Kokkos::View<double**> K1("RK_K1",ncells_loc,cf.nv+cv);
    Kokkos::View<double**> K2("RK_K2",ncells_loc,cf.nv+cv);
    Kokkos::View<int*> nlft("nlft",ncells_loc);
    Kokkos::View<int*> nrht("nrht",ncells_loc);
    Kokkos::View<int*> nbot("nbot",ncells_loc);
    Kokkos::View<int*> ntop("ntop",ncells_loc);
    Kokkos::View<int*> nbak("nbak",ncells_loc);
    Kokkos::View<int*> nfrt("nfrt",ncells_loc);
    Kokkos::View<int*> celltype("celltype",ncells_loc);
    Kokkos::View<int*>::HostMirror h_nlft = Kokkos::create_mirror_view(nlft);
    Kokkos::View<int*>::HostMirror h_nrht = Kokkos::create_mirror_view(nrht);
    Kokkos::View<int*>::HostMirror h_nbot = Kokkos::create_mirror_view(nbot);
    Kokkos::View<int*>::HostMirror h_ntop = Kokkos::create_mirror_view(ntop);
    Kokkos::View<int*>::HostMirror h_nbak = Kokkos::create_mirror_view(nbak);
    Kokkos::View<int*>::HostMirror h_nfrt = Kokkos::create_mirror_view(nfrt);
    Kokkos::View<int*>::HostMirror h_celltype = Kokkos::create_mirror_view(celltype);
        
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

    //Kokkos::View<int*> realcell("realcell",cf.ngi * cf.ngj * cf.ngk);
    //Kokkos::View<int*>::HostMirror h_realcell = Kokkos::create_mirror_view(realcell);
    //Kokkos::View<int*> facecell("facecell",cf.ngi * cf.ngj * cf.ngk);
    //Kokkos::View<int*>::HostMirror h_facecell = Kokkos::create_mirror_view(facecell);
    //Kokkos::View<int*>::HostMirror h_bordcell = Kokkos::create_mirror_view(bordcell);
    //Kokkos::View<int*> celltype("celltype",ncells);
    //Kokkos::View<int*> bordcell("bordcell",cf.ngi * cf.ngj * cf.ngk);
    //Kokkos::View<int*> corntype("corntype",cf.ngi * cf.ngj * cf.ngk);
    //Kokkos::View<int*>::HostMirror h_celltype = Kokkos::create_mirror_view(celltype);
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
    
    /*
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
                    h_celltype(loc_idx) = BOUNDARY_CELL_;
                }
                // for border cells (cells touching ghost/halo cells)
                else if (i == cf.ng || i == cf.ngi - cf.ng - 1 || j == cf.ng || j == cf.ngj - cf.ng - 1) { // || k == cf.ng || k == cf.ngk - cf.ng - 1) {
                    h_celltype(loc_idx) = BORDER_CELL_;
                }
                // for real cells
                else {
                    h_celltype(loc_idx) = REAL_CELL_;
                }
            }
        }
    }
    */
            //printf("%d %d %d %d %d %d %d\n", cf.rank, cf.xMinus, cf.xPlus, cf.yMinus, cf.yPlus, cf.zMinus, cf.zPlus);

    //for (int n = 0; n < ncells; n++) {
    //    if (h_celltype(n) == BORDER_CELL_ && h_celltype(h_ntop(n)) == BOUNDARY_CELL_) printf("%d\n", n);
    //}


    //Kokkos::deep_copy(nlft, h_nlft);
    //Kokkos::deep_copy(nrht, h_nrht);
    //Kokkos::deep_copy(nbot, h_nbot);
    //Kokkos::deep_copy(ntop, h_ntop);
    //Kokkos::deep_copy(nbak, h_nbak);
    //Kokkos::deep_copy(nfrt, h_nfrt);
    //Kokkos::deep_copy(realcell, h_realcell);
    //Kokkos::deep_copy(facecell, h_facecell);
    //Kokkos::deep_copy(celltype, h_celltype);
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
    double *x = (double*)malloc(ncells_loc*sizeof(double));
    double *y = (double*)malloc(ncells_loc*sizeof(double));
    double *z = (double*)malloc(ncells_loc*sizeof(double));

    float *xSP = (float*)malloc(ncells_loc*sizeof(float));
    float *ySP = (float*)malloc(ncells_loc*sizeof(float));
    float *zSP = (float*)malloc(ncells_loc*sizeof(float));

    double *g_x, *g_y, *g_dx, *g_dy;
    float xwinmin, xwinmax, ywinmin, ywinmax;
    
    //printf("Here3\n");
    

    /* calculate rank local grid coordinates */
    //for (int k=0; k<cf.nk; ++k){
    //    for (int j=0; j<cf.nj; ++j){
    //        for (int i=0; i<cf.ni; ++i){
    //            idx = (cf.ni*cf.nj)*k+cf.ni*j+i;
    /*
            for (int idx = 0; idx < ncells_loc; idx++) {
                x[idx] = mesh->x[idx];
                y[idx] = mesh->x[idx];
                z[idx] = 1.0;

                xSP[idx] = (float)x[idx];
                ySP[idx] = (float)y[idx];
                zSP[idx] = (float)z[idx];
            }
    */
    //        }
    //    }
    //}
    
    double time = cf.time;
    int tstart = cf.tstart;
    
    MPI_Barrier(cf.comm);

    /*** Read Restart or Write initial conditions ***/
    if (cf.restart == 1){
        if (cf.rank == 0) printf("\nReading Restart File...\n");
        //readSolution(cf,myV);
    }else{
        if (cf.rank == 0)
            if (cf.write_freq >0 || cf.restart_freq>0)
                printf("\nWriting Initial Conditions...\n");
        //if (cf.write_freq >0)
        //    writeSolution(cf,xSP,ySP,zSP,myV,0,0.00,mesh);
        //if (cf.restart_freq >0)
            //writeRestart(cf,x,y,z,myV,0,0.00);
    }

    //printf("Here4\n");


    //double *X = NULL;
    //double *Y = NULL;
    //double *E = NULL;
    //double *D = NULL;
    int flags = 0;
    flags = (RESTART_DATA | REZONE_DATA | LOAD_BALANCE_MEMORY);
    // x momentum
    X = (double *)state_memory.memory_malloc(ncells, sizeof(double), "X", flags);
    // y momentum
    Y = (double *)state_memory.memory_malloc(ncells, sizeof(double), "Y", flags);
    // total energy
    E = (double *)state_memory.memory_malloc(ncells, sizeof(double), "E", flags);
    // density
    D = (double *)state_memory.memory_malloc(ncells, sizeof(double), "D", flags);

    // temporary copy over of state variables, should be done with pointers.
    Kokkos::View<double**>::HostMirror h_myV = Kokkos::create_mirror_view(myV);
    Kokkos::deep_copy(h_myV, myV);
        
    int conv_idx = 0;
    for (int k = 0; k < 1; k++) {
        for (int j = cf.ng; j < cf.ngj - cf.ng; j++) {
            for (int i = cf.ng; i < cf.ngi - cf.ng; i++) {
                int loc_idx = k * cf.ngj * cf.ngi + j * cf.ngi + i;
                while(mesh->celltype[conv_idx] != REAL_CELL) {
                    conv_idx++;
                }
                X[conv_idx] = h_myV(loc_idx, 0);
                Y[conv_idx] = h_myV(loc_idx, 1);
                E[conv_idx] = h_myV(loc_idx, 2);
                D[conv_idx] = h_myV(loc_idx, 3);
                //printf("%d) X %f Y %f E %f D %f\n", conv_idx, X[conv_idx], Y[conv_idx], E[conv_idx], D[conv_idx]);
                conv_idx++;
            }
        }
    }
    for (int ic = 0; ic < ncells_loc; ic++) {
        if (mesh->celltype[ic] == LEFT_BOUNDARY) {
            X[ic] = X[mesh->nrht[ic]];
            Y[ic] = Y[mesh->nrht[ic]];
            E[ic] = E[mesh->nrht[ic]];
            D[ic] = D[mesh->nrht[ic]];
        }
        else if (mesh->celltype[ic] == RIGHT_BOUNDARY) {
            X[ic] = X[mesh->nlft[ic]];
            Y[ic] = Y[mesh->nlft[ic]];
            E[ic] = E[mesh->nlft[ic]];
            D[ic] = D[mesh->nlft[ic]];
        }
        else if (mesh->celltype[ic] == BOTTOM_BOUNDARY) {
            X[ic] = X[mesh->ntop[ic]];
            Y[ic] = Y[mesh->ntop[ic]];
            E[ic] = E[mesh->ntop[ic]];
            D[ic] = D[mesh->ntop[ic]];
        }
        else if (mesh->celltype[ic] == TOP_BOUNDARY) {
            X[ic] = X[mesh->nbot[ic]];
            Y[ic] = Y[mesh->nbot[ic]];
            E[ic] = E[mesh->nbot[ic]];
            D[ic] = D[mesh->nbot[ic]];
        }
        //printf("%d) X %f Y %f E %f D %f\n", ic, h_myV(ic,0), h_myV(ic,1), h_myV(ic,2), h_myV(ic,3));
    }
    //Kokkos::deep_copy(myV,h_myV);

#ifndef DYNAMIC_REZONE
    vector<char_t>      mpot;

    mpot.resize(mesh->ncells_ghost);

    mesh->set_bounds(ncells_loc);
    int icount, jcount;
    new_ncells = mesh->state_calc_refine_potential(mpot, icount, jcount, D);
    printf("New cells %d\n\n\n\n\n\n", new_ncells);

    //  Resize the mesh, inserting cells where refinement is necessary.
    mesh->rezone_all(icount, jcount, mpot, 1, state_memory);
    memory_reset_ptrs();

    vector<char_t>().swap(mpot);

    mesh->ncells = new_ncells;
    ncells_loc   = new_ncells;

    mesh->do_load_balance_local(ncells_loc, NULL, state_memory);
    memory_reset_ptrs();
    mesh->set_bounds(ncells_loc);

    mesh->proc.resize(ncells_loc);
    if (icount)
    {  
        vector<int> index(ncells_loc);
        mesh->partition_cells(numpe, index, cycle_reorder);
        state_reorder(index);
        memory_reset_ptrs();
    }

    mesh->calc_neighbors_local();

    mesh->partition_measure();

    mesh->set_bounds(ncells_loc);

    Kokkos::resize(myV, ncells_loc, cf.nv+cv);
    Kokkos::resize(h_myV, ncells_loc, cf.nv+cv);
    Kokkos::resize(tmp, ncells_loc, cf.nv+cv);
    Kokkos::resize(K1, ncells_loc, cf.nv+cv);
    Kokkos::resize(K2, ncells_loc, cf.nv+cv);
    Kokkos::resize(nlft, ncells_loc);
    Kokkos::resize(nrht, ncells_loc);
    Kokkos::resize(nbot, ncells_loc);
    Kokkos::resize(ntop, ncells_loc);
    Kokkos::resize(nbak, ncells_loc);
    Kokkos::resize(nfrt, ncells_loc);
    Kokkos::resize(celltype, ncells_loc);
    Kokkos::resize(h_nlft, ncells_loc);
    Kokkos::resize(h_nrht, ncells_loc);
    Kokkos::resize(h_nbot, ncells_loc);
    Kokkos::resize(h_ntop, ncells_loc);
    Kokkos::resize(h_nbak, ncells_loc);
    Kokkos::resize(h_nfrt, ncells_loc);
    Kokkos::resize(h_celltype, ncells_loc);

    for (int ic = 0; ic < ncells_loc; ic++) {
        if (mesh->celltype[ic] == LEFT_BOUNDARY) {
            X[ic] = -X[mesh->nrht[ic]];
            Y[ic] = -Y[mesh->nrht[ic]];
            E[ic] = -E[mesh->nrht[ic]];
            D[ic] = -D[mesh->nrht[ic]];
        }
        else if (mesh->celltype[ic] == RIGHT_BOUNDARY) {
            X[ic] = -X[mesh->nlft[ic]];
            Y[ic] = -Y[mesh->nlft[ic]];
            E[ic] = -E[mesh->nlft[ic]];
            D[ic] = -D[mesh->nlft[ic]];
        }
        else if (mesh->celltype[ic] == BOTTOM_BOUNDARY) {
            X[ic] = -X[mesh->ntop[ic]];
            Y[ic] = -Y[mesh->ntop[ic]];
            E[ic] = -E[mesh->ntop[ic]];
            D[ic] = -D[mesh->ntop[ic]];
        }
        else if (mesh->celltype[ic] == TOP_BOUNDARY) {
            X[ic] = -X[mesh->nbot[ic]];
            Y[ic] = -Y[mesh->nbot[ic]];
            E[ic] = -E[mesh->nbot[ic]];
            D[ic] = -D[mesh->nbot[ic]];
        }
        h_nlft(ic) = mesh->nlft[ic];
        h_nrht(ic) = mesh->nrht[ic];
        h_nbot(ic) = mesh->nbot[ic];
        h_ntop(ic) = mesh->ntop[ic];
        h_celltype(ic) = mesh->celltype[ic];
        h_myV(ic,0) = X[ic];
        h_myV(ic,1) = Y[ic];
        h_myV(ic,2) = E[ic];
        h_myV(ic,3) = D[ic];
    }

    Kokkos::deep_copy(nlft, h_nlft);
    Kokkos::deep_copy(nrht, h_nrht);
    Kokkos::deep_copy(nbot, h_nbot);
    Kokkos::deep_copy(ntop, h_ntop);
    //Kokkos::deep_copy(nbak, h_nbak);
    //Kokkos::deep_copy(nfrt, h_nfrt);
    Kokkos::deep_copy(celltype, h_celltype);
    Kokkos::deep_copy(myV, h_myV);



#endif


#ifdef HAVE_GRAPHICS
    mesh->calc_spatial_coordinates(0);
    g_x = (double *) malloc(ncells_loc*sizeof(double));
    g_y = (double *) malloc(ncells_loc*sizeof(double));
    g_dx = (double *) malloc(ncells_loc*sizeof(double));
    g_dy = (double *) malloc(ncells_loc*sizeof(double));
    for (int ic = 0; ic < ncells_loc; ic++) {
        g_x[ic] = mesh->x[ic]; 
        g_y[ic] = mesh->y[ic]; 
        g_dx[ic] = mesh->dx[ic]; 
        g_dy[ic] = mesh->dy[ic]; 
        //printf("%d %lf %lf %d %d %d\n", ic, mesh->x[ic], mesh->y[ic], mesh->level[ic], mesh->i[ic], mesh->j[ic]);
    }
    printf("\n");
    xwinmin = mesh->xmin-2.0;
    xwinmax = (float)(mesh->xmax+2.0);
    ywinmin = mesh->ymin-12.0+2.0;
    ywinmax = (float)(mesh->ymax+2.0);

    set_display_mysize(ncells_loc);
    set_display_cell_coordinates_double(g_x, g_dx, g_y, g_dy);
    set_display_cell_data_double(E);

    set_display_window(xwinmin,xwinmax,ywinmin,ywinmax);
    set_display_viewmode(1);
    set_display_outline(1);
    init_display(&argc, argv, "Fiesta AMR");
    draw_scene();
    sleep(10);
#endif




    //for (int ic = 0; ic < ncells_loc; ic++) {
        //printf("%d) lvl %d celltype %d X %f Y %f E %f D %f\n", ic, mesh->level[ic], mesh->celltype[ic], h_myV(ic,0), h_myV(ic,1), h_myV(ic,2), h_myV(ic,3));
        //printf("%d) lvl %d celltype %d i %d j %d\n", ic, mesh->level[ic], mesh->celltype[ic], mesh->i[ic], mesh->j[ic]);
    //}





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

    typedef Kokkos::MDRangePolicy<Kokkos::Rank<4>> policy_1;
    typedef Kokkos::MDRangePolicy<Kokkos::Rank<2>> policy_f1;
    //policy_f1 unstruct = policy_f1({0,0},{ncells_loc,cf.nv+cv});
    
//#ifdef TEST_DONE
    for (int t=tstart; t<cf.nt; ++t){
        time = time + cf.dt;

#ifdef DYNAMIC_REZONE
    vector<char_t>      mpot;

    mpot.resize(mesh->ncells_ghost);

    mesh->set_bounds(ncells_loc);
    int icount, jcount;
    new_ncells = mesh->state_calc_refine_potential(mpot, icount, jcount, D);
    printf("%d\n\n", new_ncells);

    //  Resize the mesh, inserting cells where refinement is necessary.
    mesh->rezone_all(icount, jcount, mpot, 1, state_memory);
    memory_reset_ptrs();

    vector<char_t>().swap(mpot);

    mesh->ncells = new_ncells;
    ncells_loc   = new_ncells;

    mesh->do_load_balance_local(ncells_loc, NULL, state_memory);
    memory_reset_ptrs();
    mesh->set_bounds(ncells_loc);

    mesh->proc.resize(ncells_loc);
    if (icount)
    {  
        vector<int> index(ncells_loc);
        mesh->partition_cells(numpe, index, cycle_reorder);
        state_reorder(index);
        memory_reset_ptrs();
    }

    mesh->calc_neighbors_local();

    mesh->partition_measure();

    mesh->set_bounds(ncells_loc);

    //printf("%d %d %d %d\n", mesh->lev_ibegin[2], mesh->lev_iend[2], mesh->lev_jbegin[2], mesh->lev_jend[2]);
    for (int ic = 0; ic < ncells_loc; ic++) {
        //printf("%d) %d %d %d %d\n", ic, mesh->celltype[ic], mesh->level[ic], mesh->i[ic], mesh->j[ic]);
    }

    Kokkos::resize(myV, ncells_loc, cf.nv+cv);
    Kokkos::resize(h_myV, ncells_loc, cf.nv+cv);
    Kokkos::resize(tmp, ncells_loc, cf.nv+cv);
    Kokkos::resize(K1, ncells_loc, cf.nv+cv);
    Kokkos::resize(K2, ncells_loc, cf.nv+cv);
    Kokkos::resize(nlft, ncells_loc);
    Kokkos::resize(nrht, ncells_loc);
    Kokkos::resize(nbot, ncells_loc);
    Kokkos::resize(ntop, ncells_loc);
    Kokkos::resize(nbak, ncells_loc);
    Kokkos::resize(nfrt, ncells_loc);
    Kokkos::resize(celltype, ncells_loc);
    Kokkos::resize(h_nlft, ncells_loc);
    Kokkos::resize(h_nrht, ncells_loc);
    Kokkos::resize(h_nbot, ncells_loc);
    Kokkos::resize(h_ntop, ncells_loc);
    Kokkos::resize(h_nbak, ncells_loc);
    Kokkos::resize(h_nfrt, ncells_loc);
    Kokkos::resize(h_celltype, ncells_loc);

    for (int ic = 0; ic < ncells_loc; ic++) {
        /*
        if (mesh->celltype[ic] == LEFT_BOUNDARY) {
            X[ic] = -X[mesh->nrht[ic]];
            Y[ic] = -Y[mesh->nrht[ic]];
            E[ic] = -E[mesh->nrht[ic]];
            D[ic] = -D[mesh->nrht[ic]];
        }
        else if (mesh->celltype[ic] == RIGHT_BOUNDARY) {
            X[ic] = -X[mesh->nlft[ic]];
            Y[ic] = -Y[mesh->nlft[ic]];
            E[ic] = -E[mesh->nlft[ic]];
            D[ic] = -D[mesh->nlft[ic]];
        }
        else if (mesh->celltype[ic] == BOTTOM_BOUNDARY) {
            X[ic] = -X[mesh->ntop[ic]];
            Y[ic] = -Y[mesh->ntop[ic]];
            E[ic] = -E[mesh->ntop[ic]];
            D[ic] = -D[mesh->ntop[ic]];
        }
        else if (mesh->celltype[ic] == TOP_BOUNDARY) {
            X[ic] = -X[mesh->nbot[ic]];
            Y[ic] = -Y[mesh->nbot[ic]];
            E[ic] = -E[mesh->nbot[ic]];
            D[ic] = -D[mesh->nbot[ic]];
        }
        */
        h_nlft(ic) = mesh->nlft[ic];
        h_nrht(ic) = mesh->nrht[ic];
        h_nbot(ic) = mesh->nbot[ic];
        h_ntop(ic) = mesh->ntop[ic];
        h_celltype(ic) = mesh->celltype[ic];
        h_myV(ic,0) = X[ic];
        h_myV(ic,1) = Y[ic];
        h_myV(ic,2) = E[ic];
        h_myV(ic,3) = D[ic];
        //printf("%d) X %f Y %f E %f D %f\n", ic, h_myV(ic,0), h_myV(ic,1), h_myV(ic,2), h_myV(ic,3));
    }

    Kokkos::deep_copy(nlft, h_nlft);
    Kokkos::deep_copy(nrht, h_nrht);
    Kokkos::deep_copy(nbot, h_nbot);
    Kokkos::deep_copy(ntop, h_ntop);
    //Kokkos::deep_copy(nbak, h_nbak);
    //Kokkos::deep_copy(nfrt, h_nfrt);
    Kokkos::deep_copy(celltype, h_celltype);
    Kokkos::deep_copy(myV, h_myV);

#endif




        policy_f1 unstruct = policy_f1({0,0},{ncells_loc,cf.nv+cv});

        /****** Low Storage Runge-Kutta 2nd order ******/
        //K1 = f(myV)
        //applyBCs(cf,myV,m,nlft,nrht,nbot,ntop,nbak,nfrt,celltype,ncells_loc);
        Kokkos::parallel_for("temp_bc", unstruct,
               KOKKOS_LAMBDA  (const int ic, const int v) {
            if (celltype(ic) == LEFT_BOUNDARY)
                myV(ic,v) = -myV(nrht(ic),v);
            else if (celltype(ic) == RIGHT_BOUNDARY)
                myV(ic,v) = -myV(nlft(ic),v);
            else if (celltype(ic) == BOTTOM_BOUNDARY)
                myV(ic,v) = -myV(ntop(ic),v);
            else if (celltype(ic) == TOP_BOUNDARY)
                myV(ic,v) = -myV(nbot(ic),v);
        });
        //if (t == 0) printf("Here6\n");
    
        f->compute(myV,nlft,nrht,nbot,ntop,nbak,nfrt,celltype,K1,ncells_loc);
        //if (t == 0) printf("Here7\n");
    
        
        //tmp = myV + k1/2
        //Kokkos::parallel_for("Loop1", policy_1({0,0,0,0},{cf.ngi, cf.ngj, cf.ngk, cf.nv+cv}),
        //       KOKKOS_LAMBDA  (const int i, const int j, const int k, const int v) {
        //    tmp(i+j*cf.ngi+k*cf.ngj*cf.ngi,v) = myV(i+j*cf.ngi+k*cf.ngj*cf.ngi,v) + cf.dt*K1(i+j*cf.ngi+k*cf.ngj*cf.ngi,v)/2;
        Kokkos::parallel_for("Loop1", unstruct,
               KOKKOS_LAMBDA  (const int ic, const int v) {
            tmp(ic,v) = myV(ic,v) + cf.dt*K1(ic,v)/2;
        });

        //printf("\n\n*******************************************************************\n\n");

        //K2 = f(tmp)
        //applyBCs(cf,tmp,m,nlft,nrht,nbot,ntop,nbak,nfrt,celltype,ncells_loc);
        Kokkos::parallel_for("temp_bc", unstruct,
               KOKKOS_LAMBDA  (const int ic, const int v) {
            if (celltype(ic) == LEFT_BOUNDARY)
                tmp(ic,v) = -tmp(nrht(ic),v);
            else if (celltype(ic) == RIGHT_BOUNDARY)
                tmp(ic,v) = -tmp(nlft(ic),v);
            else if (celltype(ic) == BOTTOM_BOUNDARY)
                tmp(ic,v) = -tmp(ntop(ic),v);
            else if (celltype(ic) == TOP_BOUNDARY)
                tmp(ic,v) = -tmp(nbot(ic),v);
        });
        f->compute(tmp,nlft,nrht,nbot,ntop,nbak,nfrt,celltype,K2,ncells_loc);

        //myV = myV + K2
        //Kokkos::parallel_for("Loop2", policy_1({0,0,0,0},{cf.ngi, cf.ngj, cf.ngk, cf.nv+cv}),
        //       KOKKOS_LAMBDA  (const int i, const int j, const int k, const int v) {
        //    myV(i+j*cf.ngi+k*cf.ngj*cf.ngi,v) = myV(i+j*cf.ngi+k*cf.ngj*cf.ngi,v) + cf.dt*K2(i+j*cf.ngi+k*cf.ngj*cf.ngi,v);
        Kokkos::parallel_for("Loop2", unstruct,
               KOKKOS_LAMBDA  (const int ic, const int v) {
            myV(ic,v) = myV(ic,v) + cf.dt*K2(ic,v);
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

#ifdef DYNAMIC_REZONE
    Kokkos::deep_copy(h_myV, myV);
    for (int ic = 0; ic < ncells_loc; ic++) {
        X[ic] = h_myV(ic, 0);
        Y[ic] = h_myV(ic, 1);
        E[ic] = h_myV(ic, 2);
        D[ic] = h_myV(ic, 3);
    }
#endif




    // Real-time graphics
    ///*
#ifdef HAVE_GRAPHICS
    mesh->calc_spatial_coordinates(0);
    free(g_x);
    free(g_dx);
    free(g_y);
    free(g_dy);
    g_x = (double *) malloc(ncells_loc*sizeof(double));
    g_y = (double *) malloc(ncells_loc*sizeof(double));
    g_dx = (double *) malloc(ncells_loc*sizeof(double));
    g_dy = (double *) malloc(ncells_loc*sizeof(double));
    for (int ic = 0; ic < ncells_loc; ic++) {
        g_x[ic] = mesh->x[ic]; 
        g_y[ic] = mesh->y[ic]; 
        g_dx[ic] = mesh->dx[ic]; 
        g_dy[ic] = mesh->dy[ic]; 
        //printf("%d %lf %lf %d %d %d\n", ic, mesh->x[ic], mesh->y[ic], mesh->level[ic], mesh->i[ic], mesh->j[ic]);
    }
    xwinmin = mesh->xmin-2.0;
    xwinmax = (float)(mesh->xmax+2.0);
    ywinmin = mesh->ymin-12.0+2.0;
    ywinmax = (float)(mesh->ymax+2.0);

    set_display_mysize(ncells_loc);
    set_display_cell_coordinates_double(g_x, g_dx, g_y, g_dy);
    set_display_cell_data_double(E);

    set_display_window(xwinmin,xwinmax,ywinmin,ywinmax);
    provide_sim_progress(time, t);
    draw_scene(); 
    sleep(5);
#endif
    //*/
   
    }
//#endif

    MPI_Barrier(cf.comm);
    total_time = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    mean_time = total_time/cf.nt;

    
    //Kokkos::View<double**>::HostMirror h_myV = Kokkos::create_mirror_view(myV);


    //Kokkos::deep_copy(h_myV, myV);
    int real_cell_cnt = 0;
    /*
    if (cf.rank == 0) {
        //    for (int kk = 0; kk < cf.ngk; kk++) {
        //        for (int jj = cf.ng; jj < cf.ngj - cf.ng; jj++) {
        //            for (int ii = cf.ng; ii < cf.ngi - cf.ng; ii++) {
        //                int loc_idx = ii + jj * cf.ngi + kk * cf.ngi * cf.ngj;
        for (int ic = 0; ic < ncells_loc; ic++) {           
            if (h_celltype(ic) != REAL_CELL) continue;
            real_cell_cnt ++;
            printf("%d) %d ", ic, h_celltype(ic));
            for (int var = 0; var < cf.nv + cv; var++) {
                printf("%f ",  h_myV(ic,var));
            }
            printf("\n");
        }
        printf("%d\n", real_cell_cnt);
        //            }
        //        }
        //    }
    }
    */
    
    mesh->terminate();


    if (cf.rank == 0){
        printf("\nSimulation Complete\n");
        printf("Total Time: %.2fs\nMean Time Step: %.2es\n",total_time,mean_time);
    }

    //if (cf.rank==0) printf("\n");
    
    MPI_Finalize();
    //Kokkos::finalize();
    return 0;
}
