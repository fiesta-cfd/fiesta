#include "input.hpp"
#include "mpi_init.hpp"
#include "cgns.hpp"
#include "bc.hpp"
#include "Kokkos_Core.hpp"
#include <mpi.h>
#include "lsdebug.hpp"

//template < class ViewType >
struct rhs_func {

    Kokkos::View<double****> mvar;
    Kokkos::View<double****> mdvar;
    Kokkos::View<double*> mcd;

    typedef Kokkos::MDRangePolicy<Kokkos::Rank<3>> policy_f;
    struct inputConfig cf;
    double dt;

    rhs_func(double dt_, struct inputConfig &cf_, const Kokkos::View<double****> & u_, Kokkos::View<double****> & k_, Kokkos::View<double*> & cd_)
        : dt(dt_), cf(cf_) , mvar(u_), mdvar(k_) , mcd(cd_){}

    void operator()() const {
        Kokkos::View<double****> dvar = mdvar;
        Kokkos::View<double****> var = mvar;
        Kokkos::View<double*> cd = mcd;

        Kokkos::parallel_for("Loopf", policy_f({0,0,0},{cf.ngi, cf.ngj, cf.ngk}),
        KOKKOS_LAMBDA __device__ (const int i, const int j, const int k) {
            
            int ns = (int)cd(0);

            double Cp = 0;
            double Cv = 0;
            double gamma,gammas,Rs;

            double ur,ul,vr,vl,wr,wl;
            double dxp,dyp,dzp;

            double rho[13] = {0,0,0,0,0,0,0,0,0,0,0,0,0};
            double uvel[13] = {0,0,0,0,0,0,0,0,0,0,0,0,0};
            double vvel[13] = {0,0,0,0,0,0,0,0,0,0,0,0,0};
            double wvel[13] = {0,0,0,0,0,0,0,0,0,0,0,0,0};
            double pres[13] = {0,0,0,0,0,0,0,0,0,0,0,0,0};

            for (int s=0; s<ns; ++s){
                rho[0]  = rho[0]  + var(i  ,j  ,k  ,4+s);
                rho[1]  = rho[1]  + var(i-1,j  ,k  ,4+s);
                rho[2]  = rho[2]  + var(i+1,j  ,k  ,4+s);
                rho[3]  = rho[3]  + var(i  ,j-1,k  ,4+s);
                rho[4]  = rho[4]  + var(i  ,j+1,k  ,4+s);
                rho[5]  = rho[5]  + var(i  ,j  ,k-1,4+s);
                rho[6]  = rho[6]  + var(i  ,j  ,k+1,4+s);
                rho[7]  = rho[7]  + var(i-2,j  ,k  ,4+s);
                rho[8]  = rho[8]  + var(i+2,j  ,k  ,4+s);
                rho[9]  = rho[9]  + var(i  ,j-2,k  ,4+s);
                rho[10] = rho[10] + var(i  ,j+2,k  ,4+s);
                rho[11] = rho[11] + var(i  ,j  ,k-2,4+s);
                rho[12] = rho[12] + var(i  ,j  ,k+2,4+s);
            }

            uvel[0]  = var(i  ,j  ,k  ,0)/rho[0];
            uvel[1]  = var(i-1,j  ,k  ,0)/rho[1];
            uvel[2]  = var(i+1,j  ,k  ,0)/rho[2];
            uvel[3]  = var(i  ,j-1,k  ,0)/rho[3];
            uvel[4]  = var(i  ,j+1,k  ,0)/rho[4];
            uvel[5]  = var(i  ,j  ,k-1,0)/rho[5];
            uvel[6]  = var(i  ,j  ,k+1,0)/rho[6];
            uvel[7]  = var(i-2,j  ,k  ,0)/rho[7];
            uvel[8]  = var(i+2,j  ,k  ,0)/rho[8];
            uvel[9]  = var(i  ,j-2,k  ,0)/rho[9];
            uvel[10] = var(i  ,j+2,k  ,0)/rho[10];
            uvel[11] = var(i  ,j  ,k-2,0)/rho[11];
            uvel[12] = var(i  ,j  ,k+2,0)/rho[12];

            vvel[0]  = var(i  ,j  ,k  ,1)/rho[0];
            vvel[1]  = var(i-1,j  ,k  ,1)/rho[1];
            vvel[2]  = var(i+1,j  ,k  ,1)/rho[2];
            vvel[3]  = var(i  ,j-1,k  ,1)/rho[3];
            vvel[4]  = var(i  ,j+1,k  ,1)/rho[4];
            vvel[5]  = var(i  ,j  ,k-1,1)/rho[5];
            vvel[6]  = var(i  ,j  ,k+1,1)/rho[6];
            vvel[7]  = var(i-2,j  ,k  ,1)/rho[7];
            vvel[8]  = var(i+2,j  ,k  ,1)/rho[8];
            vvel[9]  = var(i  ,j-2,k  ,1)/rho[9];
            vvel[10] = var(i  ,j+2,k  ,1)/rho[10];
            vvel[11] = var(i  ,j  ,k-2,1)/rho[11];
            vvel[12] = var(i  ,j  ,k+2,1)/rho[12];

            wvel[0]  = var(i  ,j  ,k  ,2)/rho[0];
            wvel[1]  = var(i-1,j  ,k  ,2)/rho[1];
            wvel[2]  = var(i+1,j  ,k  ,2)/rho[2];
            wvel[3]  = var(i  ,j-1,k  ,2)/rho[3];
            wvel[4]  = var(i  ,j+1,k  ,2)/rho[4];
            wvel[5]  = var(i  ,j  ,k-1,2)/rho[5];
            wvel[6]  = var(i  ,j  ,k+1,2)/rho[6];
            wvel[7]  = var(i-2,j  ,k  ,2)/rho[7];
            wvel[8]  = var(i+2,j  ,k  ,2)/rho[8];
            wvel[9]  = var(i  ,j-2,k  ,2)/rho[9];
            wvel[10] = var(i  ,j+2,k  ,2)/rho[10];
            wvel[11] = var(i  ,j  ,k-2,2)/rho[11];
            wvel[12] = var(i  ,j  ,k+2,2)/rho[12];
            
            for (int s=0; s<ns; ++s){
                gammas = cd(4+2*s);
                Rs = cd(4+2*s+1);
                Cp = Cp + (gammas*Rs)/(gammas-1);
                Cv = Cv + Rs/(gammas-1);
            }

            gamma = Cp/Cv;

            pres[0]  = (gamma-1)*(var(i  ,j  ,k  ,3)-0.5*rho[0]  *(uvel[0] *uvel[0]  + vvel[0] *vvel[0]  + wvel[0] *wvel[0] ));
            pres[1]  = (gamma-1)*(var(i-1,j  ,k  ,3)-0.5*rho[1]  *(uvel[1] *uvel[1]  + vvel[1] *vvel[1]  + wvel[1] *wvel[1] ));
            pres[2]  = (gamma-1)*(var(i+1,j  ,k  ,3)-0.5*rho[2]  *(uvel[2] *uvel[2]  + vvel[2] *vvel[2]  + wvel[2] *wvel[2] ));
            pres[3]  = (gamma-1)*(var(i  ,j-1,k  ,3)-0.5*rho[3]  *(uvel[3] *uvel[3]  + vvel[3] *vvel[3]  + wvel[3] *wvel[3] ));
            pres[4]  = (gamma-1)*(var(i  ,j+1,k  ,3)-0.5*rho[4]  *(uvel[4] *uvel[4]  + vvel[4] *vvel[4]  + wvel[4] *wvel[4] ));
            pres[5]  = (gamma-1)*(var(i  ,j  ,k-1,3)-0.5*rho[5]  *(uvel[5] *uvel[5]  + vvel[5] *vvel[5]  + wvel[5] *wvel[5] ));
            pres[6]  = (gamma-1)*(var(i  ,j  ,k+1,3)-0.5*rho[6]  *(uvel[6] *uvel[6]  + vvel[6] *vvel[6]  + wvel[6] *wvel[6] ));
            pres[7]  = (gamma-1)*(var(i-2,j  ,k  ,3)-0.5*rho[7]  *(uvel[7] *uvel[7]  + vvel[7] *vvel[7]  + wvel[7] *wvel[7] ));
            pres[8]  = (gamma-1)*(var(i+2,j  ,k  ,3)-0.5*rho[8]  *(uvel[8] *uvel[8]  + vvel[8] *vvel[8]  + wvel[8] *wvel[8] ));
            pres[9]  = (gamma-1)*(var(i  ,j-2,k  ,3)-0.5*rho[9]  *(uvel[9] *uvel[9]  + vvel[9] *vvel[9]  + wvel[9] *wvel[9] ));
            pres[10] = (gamma-1)*(var(i  ,j+2,k  ,3)-0.5*rho[10] *(uvel[10]*uvel[10] + vvel[10]*vvel[10] + wvel[10]*wvel[10]));
            pres[11] = (gamma-1)*(var(i  ,j  ,k-2,3)-0.5*rho[11] *(uvel[11]*uvel[11] + vvel[11]*vvel[11] + wvel[11]*wvel[11]));
            pres[12] = (gamma-1)*(var(i  ,j  ,k+2,3)-0.5*rho[12] *(uvel[12]*uvel[12] + vvel[12]*vvel[12] + wvel[12]*wvel[12]));

            ur = (-uvel[2]  + 7.0*uvel[0] + 7.0*uvel[1] - uvel[7])/12.0;
            ul = (-uvel[8]  + 7.0*uvel[2] + 7.0*uvel[0] - uvel[1])/12.0;
            
            vr = (-uvel[4]  + 7.0*vvel[0] + 7.0*vvel[3] - vvel[9])/12.0;
            vl = (-uvel[10] + 7.0*vvel[4] + 7.0*vvel[0] - vvel[3])/12.0;
            
            wr = (-uvel[6]  + 7.0*wvel[0] + 7.0*wvel[5] - wvel[11])/12.0;
            wl = (-uvel[12] + 7.0*wvel[6] + 7.0*wvel[0] - wvel[5] )/12.0;

            dxp = (pres[7]  - 8.0*pres[1] + 8.0*pres[2] - pres[8] )/(12.0*cd(1));
            dyp = (pres[9]  - 8.0*pres[3] + 8.0*pres[4] - pres[10])/(12.0*cd(2));
            dzp = (pres[11] - 8.0*pres[5] + 8.0*pres[6] - pres[12])/(12.0*cd(3));

            //do weno here

            dvar(i,j,k,0) = var(i,j,k,0)/cd(1);
            dvar(i,j,k,1) = var(i,j,k,1)/cd(1);
            dvar(i,j,k,2) = var(i,j,k,2)/cd(1);
            dvar(i,j,k,3) = var(i,j,k,3)/cd(1);
            for (int s=0; s<ns; ++s)
                dvar(i,j,k,4+s) = var(i,j,k,4+s)/cd(1);
            
            
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
