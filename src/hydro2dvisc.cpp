#include "input.hpp"
#include "mpi.hpp"
#include "cgns.hpp"
#include "bc.hpp"
#include "Kokkos_Core.hpp"
#include <mpi.h>
#include "debug.hpp"
#include "hydro2dvisc.hpp"
#include "flux.hpp"

struct calculateRhoAndPressure2dv {
    FS4D var;
    FS2D p;
    FS2D T;
    FS2D rho;
    Kokkos::View<double*> cd;

    calculateRhoAndPressure2dv (FS4D var_, FS2D p_, FS2D rho_, FS2D T_, Kokkos::View<double*> cd_)
         : var(var_), p(p_), rho(rho_), T(T_), cd(cd_) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const int i, const int j) const {

        int ns = (int)cd(0);
        int nv = (int)cd(4);
        double gamma, gammas, Rs;
        double Cp = 0;
        double Cv = 0;

        rho(i,j) = 0.0;

        // Total Density for this cell
        for (int s=0; s<ns; ++s){
            rho(i,j) = rho(i,j) + var(i,j,0,3+s);
        }

        // Calculate mixture ratio of specific heats
        for (int s=0; s<ns; ++s){
            gammas = cd(5+2*s);
            Rs = cd(5+2*s+1);

            // accumulate mixture heat capacity by mass fraction weights
            Cp = Cp + (var(i,j,0,3+s)/rho(i,j))*( gammas*Rs/(gammas-1) );
            Cv = Cv + (var(i,j,0,3+s)/rho(i,j))*( Rs/(gammas-1) );
        }
        gamma = Cp/Cv;

        // calculate pressure assuming perfect gas
        p(i,j) = (gamma-1)*( var(i,j,0,2) - (0.5/rho(i,j))
                  *(var(i,j,0,0)*var(i,j,0,0) + var(i,j,0,1)*var(i,j,0,1)) );

        T(i,j) = p(i,j)/( (Cp-Cv)*rho(i,j) );

        //if (i==2004)
        //    printf("%f, %f, %f, %f\n",T(2001,j),T(2002,j),T(2003,j),T(2004,j));
    }
};


struct advect2dv {
    
    FS4D dvar;
    FS2D fluxx;
    FS2D fluxy;
    Kokkos::View<double*> cd;
    int v;

    advect2dv (FS4D dvar_, FS2D fx_, FS2D fy_, Kokkos::View<double*> cd_, int v_)
        : dvar(dvar_), fluxx(fx_), fluxy(fy_), cd(cd_), v(v_) {}
    
    KOKKOS_INLINE_FUNCTION
    void operator()(const int i, const int j) const {

        double a;

        a = -( (fluxx(i,j) - fluxx(i-1,j))
                          +(fluxy(i,j) - fluxy(i,j-1)) );

        //if (i==2002 && j==3)
        //    printf("adv %d: %f\n",v,a);
        dvar(i,j,0,v) = a;
    }
};

struct applyPressure2dv {
    
    FS4D dvar;
    FS2D p;
    Kokkos::View<double*> cd;

    applyPressure2dv (FS4D dvar_, FS2D p_, Kokkos::View<double*> cd_)
        : dvar(dvar_), p(p_), cd(cd_) {}
    
    KOKKOS_INLINE_FUNCTION
    void operator()(const int i, const int j) const {
        double dx = cd(1);
        double dy = cd(2);
        // calculate pressure gradient across cell in each direction using 4th order
        // central difference
        double dxp = ( p(i-2,j) - 8.0*p(i-1,j) + 8.0*p(i+1,j) - p(i+2,j) )/(12.0*dx);
        double dyp = ( p(i,j-2) - 8.0*p(i,j-1) + 8.0*p(i,j+1) - p(i,j+2) )/(12.0*dy);

        //apply pressure gradient term to right hand side of Euler equation dV/dt = ...
        double p1,p2;
        p1 = dvar(i,j,0,0) - dxp;
        p2 = dvar(i,j,0,1) - dyp;

        //if (i==1999 && j==3)
            //printf("pres %f %f\n",dxp,dyp);

        dvar(i,j,0,0) = p1;
        dvar(i,j,0,1) = p2;
    }
};

struct calculateStressTensor2dv {
    
    FS4D var;
    FS2D rho;
    FS4D stressx;
    FS4D stressy;
    FS2D T;
    Kokkos::View<double*> cd;

    calculateStressTensor2dv (FS4D var_, FS2D rho_, FS2D T_, FS4D strx_, FS4D stry_, Kokkos::View<double*> cd_)
        : var(var_), rho(rho_), T(T_), stressx(strx_), stressy(stry_), cd(cd_) {}
    
    KOKKOS_INLINE_FUNCTION
    void operator()(const int i, const int j) const {
        int ns = (int)cd(0);
        double dx = cd(1);
        double dy = cd(2);
        double mu1 = 2.928e-5;
        double mu2 = 1.610e-5;
        double dudx,dvdy,dudy,dvdx;

        double muij = (var(i,j,0,3)*mu1 + var(i,j,0,4)*mu2)/rho(i,j);
        double muip = (var(i+1,j,0,3)*mu1 + var(i+1,j,0,4)*mu2)/rho(i+1,j);
        double mujp = (var(i,j+1,0,3)*mu1 + var(i,j+1,0,4)*mu2)/rho(i,j+1);

        double mur = (muip + muij)/2;
        double mut = (mujp + muij)/2;

        //double rhox = ( -     rho(i+2,j) + 7.0*rho(i+1,j)
        //         + 7.0*rho(i  ,j) -     rho(i-1,j) )/12.0;

        //double rhoy = ( -     rho(i,j+2) + 7.0*rho(i,j+1)
        //         + 7.0*rho(i  ,j) -     rho(i,j-1) )/12.0;

        //xface
        dudx = (var(i+1,j,0,0)/rho(i+1,j) - var(i,j,0,0)/rho(i,j))/dx;

        dvdy = ( (var(i+1,j+1,0,1)/rho(i+1,j+1) + var(i,j+1,0,1)/rho(i,j+1))
                -(var(i+1,j-1,0,1)/rho(i+1,j-1) + var(i,j-1,0,1)/rho(i,j-1)) )/(4.0*dy);

        dvdx = (var(i+1,j,0,1)/rho(i+1,j) - var(i,j,0,1)/rho(i,j))/dx;

        dudy = ( (var(i+1,j+1,0,0)/rho(i+1,j+1) + var(i,j+1,0,0)/rho(i,j+1))
                -(var(i+1,j-1,0,0)/rho(i+1,j-1) + var(i,j-1,0,0)/rho(i,j-1)) )/(4.0*dy);
        
        //stressx(i,j,0,0) = (4.0/3.0)*mu*dudx;
        //stressx(i,j,0,0) = mu*dudx;
        //stressx(i,j,0,0) = (2.0/3.0)*mu*(2.0*dudx-dvdy);
        //stressx(i,j,0,1) = 0.0;
        //stressx(i,j,1,0) = 0.0;
        //stressx(i,j,1,1) = 0.0;
        stressx(i,j,1,1) = (2.0/3.0)*mur*(2.0*dvdy-dudx);
        stressx(i,j,0,1) = mur*(dudy+dvdx);
        stressx(i,j,1,0) =stressx(i,j,0,1);

        //if (stressx(i,j,0,0) > 0 || stressx(i,j,0,0) < 0)
        //    printf("%f\n",stressx(i,j,0,0));
        //if (dudx >0 || dudx < 0)
        //    printf("dudx: %.25f\n",dudx);
        //if (i==2002 && j==3)
        //    printf("FACE: %f = %f * %f * %f\n",stressx(2002,3,0,0),dudx,rhox,mu);

        //yface
        dudx = ( (var(i+1,j+1,0,0)/rho(i+1,j+1) + var(i+1,j,0,0)/rho(i+1,j))
                -(var(i-1,j+1,0,0)/rho(i-1,j+1) + var(i-1,j,0,0)/rho(i-1,j)) )/(4.0*dx);

        dvdy = (var(i,j+1,0,1)/rho(i,j+1) - var(i,j,0,1)/rho(i,j))/dy;

        dvdx = ( (var(i+1,j+1,0,1)/rho(i+1,j+1) + var(i+1,j,0,1)/rho(i+1,j))
                -(var(i-1,j+1,0,1)/rho(i-1,j+1) + var(i-1,j,0,1)/rho(i-1,j)) )/(4.0*dx);

        dudy = (var(i,j+1,0,0)/rho(i,j+1) - var(i,j,0,0)/rho(i,j))/dy;

        //stressy(i,j,0,0) = (2.0/3.0)*rhoy*mu*(2.0*dudx);
        //stressy(i,j,0,0) = 0.0;
        //stressy(i,j,0,0) = (2.0/3.0)*mu*(2.0*dudx-dvdy);
        //stressy(i,j,0,1) = 0.0;
        //stressy(i,j,1,0) = 0.0;
        //stressy(i,j,1,1) = 0.0;
        stressy(i,j,1,1) = (2.0/3.0)*mut*(2.0*dvdy-dudx);
        stressy(i,j,0,1) = mut*(dudy+dvdx);
        stressy(i,j,1,0) = stressy(i,j,0,1);
    }
};
struct calculateHeatFlux2dv {
    
    FS4D var;
    FS2D rho;
    FS2D qx;
    FS2D qy;
    FS2D T;
    Kokkos::View<double*> cd;

    double k = 0.018;

    calculateHeatFlux2dv (FS4D var_, FS2D rho_, FS2D T_, FS2D qx_, FS2D qy_, Kokkos::View<double*> cd_)
        : var(var_), rho(rho_), T(T_), qx(qx_), qy(qy_), cd(cd_) {}
    
    KOKKOS_INLINE_FUNCTION
    void operator()(const int i, const int j) const {
        double dx = cd(1);
        double dy = cd(2);

        qx(i,j) = -k*(T(i+1,j)-T(i,j))/dx;
        qy(i,j) = -k*(T(i,j+1)-T(i,j))/dy;

        //if (i==2002 && j==3)
        //    printf("%f\n",qx(2002,3));

    }
};

struct applyViscousTerm2dv {
    
    FS4D dvar;
    FS4D var;
    FS2D rho;
    FS2D qx;
    FS2D qy;
    FS4D stressx;
    FS4D stressy;
    Kokkos::View<double*> cd;

    applyViscousTerm2dv (FS4D dvar_, FS4D var_, FS2D rho_, FS4D strx_, FS4D stry_, FS2D qx_, FS2D qy_, Kokkos::View<double*> cd_)
        : dvar(dvar_), var(var_), rho(rho_), stressx(strx_), stressy(stry_), qx(qx_), qy(qy_), cd(cd_) {}
    
    KOKKOS_INLINE_FUNCTION
    void operator()(const int i, const int j) const {
        double dx = cd(1);
        double dy = cd(2);
        double a,b,c1,c2;

        double ur = (var(i+1,j  ,0,0)/rho(i+1,j  ) + var(i  ,j  ,0,0)/rho(i  ,j  ))/2.0;
        double ul = (var(i  ,j  ,0,0)/rho(i  ,j  ) + var(i-1,j  ,0,0)/rho(i-1,j  ))/2.0;
        double vr = (var(i+1,j  ,0,1)/rho(i+1,j  ) + var(i  ,j  ,0,1)/rho(i  ,j  ))/2.0;
        double vl = (var(i  ,j  ,0,1)/rho(i  ,j  ) + var(i-1,j  ,0,1)/rho(i-1,j  ))/2.0;

        double ut = (var(i  ,j+1,0,0)/rho(i  ,j+1) + var(i  ,j  ,0,0)/rho(i  ,j  ))/2.0;
        double ub = (var(i  ,j  ,0,0)/rho(i  ,j  ) + var(i  ,j-1,0,0)/rho(i  ,j-1))/2.0;
        double vt = (var(i  ,j+1,0,1)/rho(i  ,j+1) + var(i  ,j  ,0,1)/rho(i  ,j  ))/2.0;
        double vb = (var(i  ,j  ,0,1)/rho(i  ,j  ) + var(i  ,j-1,0,1)/rho(i  ,j-1))/2.0;

        a = (stressx(i,j,0,0)-stressx(i-1,j,0,0))/dx
            +(stressx(i,j,0,1)-stressx(i-1,j,0,1))/dy;

        b = (stressx(i,j,1,0)-stressx(i,j-1,1,0))/dx
            +(stressx(i,j,1,1)-stressx(i,j-1,1,1))/dy;

        c1 =  (ur*stressx(i,j,0,0)-ul*stressx(i-1,j,0,0))/dx
             +(vl*stressy(i,j,0,1)-vl*stressy(i,j-1,0,1))/dx
             +(ut*stressy(i,j,1,0)-ub*stressy(i,j-1,1,0))/dy
             +(vt*stressy(i,j,1,1)-vb*stressy(i,j-1,1,1))/dy;

        c2 =  (qx(i,j)-qx(i-1,j))/dx
             +(qy(i,j)-qy(i,j-1))/dy;

        //if (i==2002 && j==3)
        //    printf("div %f: %f, %f\n",a1,stressx(i-1,j,0,0),stressx(i,j,0,0));
        dvar(i,j,0,0) = dvar(i,j,0,0) + a;
        dvar(i,j,0,1) = dvar(i,j,0,1) + b;
        dvar(i,j,0,2) = dvar(i,j,0,2) + c1 - c2;

        //if (i == 2000)
        //    printf("%f, %f, %f\n",a1,a2,a3);
    }
};

hydro2dvisc_func::hydro2dvisc_func(struct inputConfig &cf_, Kokkos::View<double*> & cd_):rk_func(cf_,cd_){
    
    var     = Kokkos::View<double****,FS_LAYOUT>("var",cf.ngi,cf.ngj,cf.ngk,cf.nv); // Primary Variable Array
    tmp1    = Kokkos::View<double****,FS_LAYOUT>("tmp1",cf.ngi,cf.ngj,cf.ngk,cf.nv); // Temporary Variable Arrayr1
    tmp2    = Kokkos::View<double****,FS_LAYOUT>("tmp2",cf.ngi,cf.ngj,cf.ngk,cf.nv); // Temporary Variable Array2
    dvar    = Kokkos::View<double****,FS_LAYOUT>("dvar",cf.ngi,cf.ngj,cf.ngk,cf.nv); // RHS Output
    p       = Kokkos::View<double**,FS_LAYOUT>("p",cf.ngi,cf.ngj);          // Pressure
    T       = Kokkos::View<double**,FS_LAYOUT>("T",cf.ngi,cf.ngj);          // Pressure
    rho     = Kokkos::View<double**,FS_LAYOUT>("rho",cf.ngi,cf.ngj);      // Total Density
    qx      = Kokkos::View<double**,FS_LAYOUT>("qx",cf.ngi,cf.ngj);  // Weno Fluxes in X direction
    qy      = Kokkos::View<double**,FS_LAYOUT>("qy",cf.ngi,cf.ngj);  // Weno Fluxes in X direction
    fluxx   = Kokkos::View<double**,FS_LAYOUT>("fluxx",cf.ngi,cf.ngj);  // Weno Fluxes in X direction
    fluxy   = Kokkos::View<double**,FS_LAYOUT>("fluxy",cf.ngi,cf.ngj);  // Weno Fluxes in Y direction
    stressx = Kokkos::View<double****,FS_LAYOUT>("stressx",cf.ngi,cf.ngj,2,2);  // stress tensor on x faces
    stressy = Kokkos::View<double****,FS_LAYOUT>("stressy",cf.ngi,cf.ngj,2,2);  // stress tensor on y faces
    cd = mcd;

    timers["flux"] = fiestaTimer("Flux Calculation");
    timers["pressgrad"] = fiestaTimer("Pressure Gradient Calculation");
    timers["calcSecond"] = fiestaTimer("Secondary Variable Calculation");
    timers["stress"] = fiestaTimer("Stress Tensor Computation");
    timers["qflux"] = fiestaTimer("Heat Flux Calculation");
    timers["visc"] = fiestaTimer("Viscous Term Calculation");

};

void hydro2dvisc_func::compute(){

    policy_f ghost_pol = policy_f({0,0},{cf.ngi, cf.ngj});
    policy_f cell_pol  = policy_f({cf.ng,cf.ng},{cf.ngi-cf.ng, cf.ngj-cf.ng});
    policy_f face_pol  = policy_f({cf.ng-1,cf.ng-1},{cf.ngi-cf.ng, cf.ngj-cf.ng});


    // Calcualte Total Density and Pressure Fields
    timers["calcSecond"].reset();
    Kokkos::parallel_for( ghost_pol, calculateRhoAndPressure2dv(var,p,rho,T,cd) );
    timers["calcSecond"].accumulate();

    // Calculate and apply weno fluxes for each variable
    for (int v=0; v<cf.nv; ++v){
        timers["flux"].reset();
        Kokkos::parallel_for( face_pol, computeFluxCentered2D(var,p,rho,fluxx,fluxy,cd,v) );
        Kokkos::parallel_for( cell_pol, advect2dv(dvar,fluxx,fluxy,cd,v) );
        timers["flux"].accumulate();
    }

    // Apply Pressure Gradient Term
    timers["pressgrad"].reset();
    Kokkos::parallel_for( cell_pol, applyPressure2dv(dvar,p,cd) );
    timers["pressgrad"].accumulate();

    timers["stress"].reset();
    Kokkos::parallel_for( face_pol, calculateStressTensor2dv(var,rho,T,stressx,stressy,cd) );
    timers["stress"].accumulate();

    timers["qflux"].reset();
    Kokkos::parallel_for( face_pol, calculateHeatFlux2dv(var,rho,T,qx,qy,cd) );
    timers["qflux"].accumulate();

    timers["visc"].reset();
    Kokkos::parallel_for( cell_pol, applyViscousTerm2dv(dvar,var,rho,stressx,stressy,qx,qy,cd) );
    timers["visc"].accumulate();
}
