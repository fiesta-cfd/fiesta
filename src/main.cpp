#include "input.hpp"
#include "mpi_init.hpp"
#include "cgns.hpp"
#include "bc.hpp"
#include "Kokkos_Core.hpp"
#include <mpi.h>
#include "lsdebug.hpp"

//template < class ViewType >
struct rhs_func {

    Kokkos::View<double****> var;
    Kokkos::View<double****> dvar;
    Kokkos::View<double*> cd;

    typedef Kokkos::MDRangePolicy<Kokkos::Rank<3>> policy_f;
    struct inputConfig cf;
    double dt;

    rhs_func(double dt_, struct inputConfig &cf_, const Kokkos::View<double****> & u_, Kokkos::View<double****> & k_, Kokkos::View<double*> & cd_)
        : dt(dt_), cf(cf_) , var(u_), dvar(k_) , cd(cd_){}

    void operator()() const {
        Kokkos::View<double****> my_dvar = dvar;
        Kokkos::View<double****> my_var = var;
        Kokkos::View<double*> my_cd = cd;

        Kokkos::parallel_for("Loopf", policy_f({0,0,0},{cf.ngi, cf.ngj, cf.ngk}),
        KOKKOS_LAMBDA __device__ (const int i, const int j, const int k) {
            
            int ns = (int)my_cd(0);
            double rho = 0;

            for (int s=0; s<my_cd(0); ++s)
                rho = rho + my_var(i,j,k,4+s);

            //rim3 = 
            //uim3 = u(i-3,j,k,0);
            
            //for (int v=0; v<6; ++v)
            //    my_dvar(i,j,k,v) = my_var(i,j,k,v)/my_cd(1);
            my_dvar(i,j,k,0) = my_var(i,j,k,0)/my_cd(1);
            
        });

    }
};

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

    Kokkos::View<double****> myV("testView",cf.ngi,cf.ngj,cf.ngk,cf.nv);
    Kokkos::View<double****> tmp("testView",cf.ngi,cf.ngj,cf.ngk,cf.nv);
    Kokkos::View<double****> K1("testView",cf.ngi,cf.ngj,cf.ngk,cf.nv);
    Kokkos::View<double****> K2("testView",cf.ngi,cf.ngj,cf.ngk,cf.nv);
    Kokkos::View<double*> cd("deviceCF",4+cf.ns*2);
    typename Kokkos::View<double*>::HostMirror hostcd = Kokkos::create_mirror_view(cd);
    Kokkos::deep_copy(hostcd, cd);
    hostcd(0) = cf.nv;
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

    loadInitialConditions(cf,myV);
    
    applyBCs(cf,myV);

    if (cf.rank == 0){
        printf("loboSHOK - %s\n",cf.inputFname);
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
    

    /* calculate rank local grid coordinates */
    for (int k=0; k<cf.nk; ++k){
        for (int j=0; j<cf.nj; ++j){
            for (int i=0; i<cf.ni; ++i){
                idx = (cf.ni*cf.nj)*k+cf.ni*j+i;
                x[idx] = (cf.iStart + i)*cf.dx;
                y[idx] = (cf.jStart + j)*cf.dy;
                z[idx] = (cf.kStart + k)*cf.dz;
            }
        }
    }
    
    //typedef typename Kokkos::View<double****> ViewType;
    typedef Kokkos::MDRangePolicy<Kokkos::Rank<4>> policy_1;
    //Kokkos::MDRangePolicy<Kokkos::Rank<3>> policy_1({0,0,0},{cf.nci, cf.ncj, cf.nck});
    
    double time = 0.0;
    
    MPI_Barrier(cf.comm);
    
    writeSolution(cf,x,y,z,myV,0,0.00);

    for (int t=0; t<cf.nt; ++t){
        time = time + cf.dt;

        //tmp = myV
        Kokkos::deep_copy(tmp,myV);

        //K1 = dt*f(tmp) dt is member data to f()
        applyBCs(cf,tmp);
        rhs_func f1(cf.dt,cf,tmp,K1, cd);
        f1();
        
        //tmp = myV + k1/2
        Kokkos::parallel_for("Loop1", policy_1({0,0,0,0},{cf.ngi, cf.ngj, cf.ngk, cf.nv}),
               KOKKOS_LAMBDA __device__ (const int i, const int j, const int k, const int v) {
            tmp(i,j,k,v) = myV(i,j,k,v) + K1(i,j,k,v)/2;
        });

        //K2 = dt*f(tmp) dt is member data to f()
        applyBCs(cf,tmp);
        rhs_func f2(cf.dt,cf,tmp,K2, cd);
        f2();

        //myV = myV + K2
        Kokkos::parallel_for("Loop2", policy_1({0,0,0,0},{cf.ngi, cf.ngj, cf.ngk, cf.nv}),
               KOKKOS_LAMBDA __device__ (const int i, const int j, const int k, const int v) {
            myV(i,j,k,v) = myV(i,j,k,v) + K2(i,j,k,v);
        });
        
        if (t % cf.out_freq == 0)
            writeSolution(cf,x,y,z,myV,t+1,time);
        if (cf.rank==0) printf("%d/%d, %f\n",t+1,cf.nt,time);
    }

    //if (cf.rank==0) printf("\n");
    //MPI_Barrier(cf.comm);
    //MYDBGR
    
    MPI_Finalize();
    //Kokkos::finalize();
    return 0;
}
