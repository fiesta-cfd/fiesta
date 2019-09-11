#include "input.hpp"
#include "mpi_init.hpp"
#include "cgns.hpp"
#include "Kokkos_Core.hpp"
#include <mpi.h>
#include "lsdebug.hpp"

//template < class ViewType >
struct rhs_func {

    Kokkos::View<double****> u;
    Kokkos::View<double****> k;

    typedef Kokkos::MDRangePolicy<Kokkos::Rank<4>> policy_f;
    struct inputConfig cf;
    double dt;

    rhs_func(double dt_, struct inputConfig &cf_, const Kokkos::View<double****> & u_, const Kokkos::View<double****> & k_)
            : dt(dt_), cf(cf_) , u(u_), k(k_) {}

    void operator()() const {
        Kokkos::View<double****> myK("testK",cf.nci,cf.ncj,cf.nck,cf.nv);
        Kokkos::View<double****> myU("testK",cf.nci,cf.ncj,cf.nck,cf.nv);
        myK = k;
        myU = u;
        Kokkos::parallel_for("Loopf", policy_f({0,0,0,0},{cf.nci, cf.ncj, cf.nck, cf.nv}),
        KOKKOS_LAMBDA __device__ (const int i, const int j, const int k, const int v) {
            myK(i,j,k,v) = 1; //u(i,j,k,v)*dt;
        });
    }
};

//void f(double dt, struct inputConfig &cf, const Kokkos::View<double****> & u, const Kokkos::View<double****> & k) {
//    typedef Kokkos::MDRangePolicy<Kokkos::Rank<4>> policy_f;
//    Kokkos::parallel_for("Loopf", policy_f({0,0,0,0},{cf.nci, cf.ncj, cf.nck, cf.nv}),
//    KOKKOS_LAMBDA __device__ (const int i, const int j, const int k, const int v) {
//        k(i,j,k,v) = 0; //u(i,j,k,v)*dt;
//    });
//}

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

    Kokkos::View<double****> myV("testView",cf.nci,cf.ncj,cf.nck,cf.nv);
    Kokkos::View<double****> tmp("testView",cf.nci,cf.ncj,cf.nck,cf.nv);
    Kokkos::View<double****> K1("testView",cf.nci,cf.ncj,cf.nck,cf.nv);
    Kokkos::View<double****> K2("testView",cf.nci,cf.ncj,cf.nck,cf.nv);
    
    loadInitialConditions(cf,myV);

    if (cf.rank == 0){
        printf("loboSHOK - %s\n",cf.inputFname);
        printf("-----------------------\n");
        printf("Running %d processes as (%d.%d,%d)\n",cf.numProcs,cf.xProcs,cf.yProcs,cf.zProcs);
        printf("nt = %d, dt = %f\n",cf.nt,cf.dt);
        printf("glbl_ni = %d, dx = %f\n",cf.glbl_ni,cf.dx);
        printf("glbl_nj = %d, dy = %f\n",cf.glbl_nj,cf.dy);
        printf("glbl_nk = %d, dz = %f\n",cf.glbl_nk,cf.dz);
        printf("Gamma = %4.2f\n",cf.gamma);
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
    writeSolution(cf,x,y,z,myV,0,0.00);
    for (int t=0; t<cf.nt; ++t){
        time = time + cf.dt;

        //tmp = myV
        Kokkos::deep_copy(tmp,myV);
        //K1 = dt*f(tmp) dt is member data to f()
        rhs_func f1(cf.dt,cf,tmp,K1);
        f1();
        //f(cf.dt,cf,tmp,K1);
        //tmp = myV + k1/2
        Kokkos::parallel_for("Loop1", policy_1({0,0,0,0},{cf.nci, cf.ncj, cf.nck, cf.nv}),
               KOKKOS_LAMBDA __device__ (const int i, const int j, const int k, const int v) {
            tmp(i,j,k,v) = myV(i,j,k,v) + K1(i,j,k,v)/2;
        });
        //K2 = dt*f(tmp) dt is member data to f()
        rhs_func f2(cf.dt,cf,tmp,K2);
        f2();
        //f(cf.dt,cf,tmp,K2);
        //myV = myV + K2
        Kokkos::parallel_for("Loop2", policy_1({0,0,0,0},{cf.nci, cf.ncj, cf.nck, cf.nv}),
               KOKKOS_LAMBDA __device__ (const int i, const int j, const int k, const int v) {
            myV(i,j,k,v) = myV(i,j,k,v) + K2(i,j,k,v);
        });
        
        if (t % cf.out_freq == 0)
            writeSolution(cf,x,y,z,myV,t+1,time);
    }

    MPI_Finalize();
    //Kokkos::finalize();
    return 0;
}
