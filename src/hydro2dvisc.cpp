#include "input.hpp"
#include "mpi.hpp"
#include "cgns.hpp"
#include "bc.hpp"
#include "Kokkos_Core.hpp"
#include <mpi.h>
#include "debug.hpp"
#include "hydro2dvisc.hpp"

struct calculateRhoAndPressure2dv {
    typedef typename Kokkos::View<double****> V4D;
    typedef typename Kokkos::View<double**> V2D;
    V4D var;
    V2D p;
    V2D rho;
    Kokkos::View<double*> cd;

    calculateRhoAndPressure2dv (V4D var_, V2D p_, V2D rho_, Kokkos::View<double*> cd_)
         : var(var_), p(p_), rho(rho_), cd(cd_) {}

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
    }
};

struct advect2dv {
    
    typedef typename Kokkos::View<double****> V4D;
    typedef typename Kokkos::View<double**> V2D;
    V4D dvar;
    V4D var;
    V2D p;
    V2D rho;
    Kokkos::View<double*> cd;
    int v;

    advect2dv (V4D dvar_, V4D var_, V2D p_, V2D rho_, Kokkos::View<double*> cd_, int v_)
                         : dvar(dvar_), var(var_), p(p_), rho(rho_), cd(cd_), v(v_) {}
    
    KOKKOS_INLINE_FUNCTION
    void operator()(const int i, const int j) const {

        double dx = cd(1);
        double dy = cd(2);
        double ur,vr,w,f1,f2,f3,f4;

        ur = var(i,j,0,0)/rho(i,j);
        vr = var(i,j,0,1)/rho(i,j);

        //for each direction
        for (int idx=0; idx<2; ++idx){
            //get stencil data.  the flux for the energy equation includes pressure so only add 
            //pressure for the energy variable (index 2 for 2d problem)
            if (idx == 0){
                f1 = var(i-2,j,0,v) + (v==2)*p(i-2,j);
                f2 = var(i-1,j,0,v) + (v==2)*p(i-1,j);
                f3 = var(i+1,j,0,v) + (v==2)*p(i+1,j);
                f4 = var(i+2,j,0,v) + (v==2)*p(i+2,j);
            } 
            if (idx == 1) {
                f1 = var(i,j-2,0,v) + (v==2)*p(i,j-2);
                f2 = var(i,j-1,0,v) + (v==2)*p(i,j-1);
                f3 = var(i,j+1,0,v) + (v==2)*p(i,j+1);
                f4 = var(i,j+2,0,v) + (v==2)*p(i,j+2);
            }

            w = (f1 - 8.0*f2 + 8.0*f3 - f4)/12.0;
            //calculate advection flux
            if (idx == 0){
                dvar(i,j,0,v) = ur*w/dx;
            }
            if (idx == 1) {
                dvar(i,j,0,v) = vr*w/dy;
            }

        }
    }
};

struct applyPressure2dv {
    
    typedef typename Kokkos::View<double****> V4D;
    typedef typename Kokkos::View<double**> V2D;
    V4D dvar;
    V2D p;
    Kokkos::View<double*> cd;

    applyPressure2dv (V4D dvar_, V2D p_, Kokkos::View<double*> cd_)
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
        dvar(i,j,0,0) = dvar(i,j,0,0) - dxp;
        dvar(i,j,0,1) = dvar(i,j,0,1) - dyp;
    }
};

struct calculateStressTensor2dv {
    
    typedef typename Kokkos::View<double****> V4D;
    typedef typename Kokkos::View<double**> V2D;
    V4D var;
    V2D rho;
    V4D stressx;
    V4D stressy;
    Kokkos::View<double*> cd;

    calculateStressTensor2dv (V4D var_, V2D rho_, V4D strx_, V4D stry_, Kokkos::View<double*> cd_)
        : var(var_), rho(rho_), stressx(strx_), stressy(stry_), cd(cd_) {}
    
    KOKKOS_INLINE_FUNCTION
    void operator()(const int i, const int j) const {
        double dx = cd(1);
        double dy = cd(2);
        double mu = 2.295e-5;
        double dudx,dvdy;

        //xface
        dudx = (var(i+1,j,0,0)/rho(i,j) - var(i,j,0,0)/rho(i,j))/dx;
        dvdy = ( (var(i+1,j+1,0,1)/rho(i,j) + var(i,j+1,0,1)/rho(i,j))
                -(var(i+1,j-1,0,1)/rho(i,j) + var(i,j-1,0,1)/rho(i,j)) )/(4*dy);
        
        stressx(i,j,0,0) = (2/3)*mu*(2*dudx-dvdy);
        stressx(i,j,1,1) = (2/3)*mu*(2*dvdy-dudx);
        stressx(i,j,0,1) = 0;
        stressx(i,j,1,0) = 0;

        //yface
        dudx = ( (var(i+1,j+1,0,0)/rho(i,j) + var(i+1,j,0,0)/rho(i,j))
                -(var(i-1,j+1,0,0)/rho(i,j) + var(i-1,j,0,0)/rho(i,j)) )/(4*dx);
        dvdy = (var(i,j+1,0,1)/rho(i,j) - var(i,j,0,1)/rho(i,j))/dy;

        stressy(i,j,0,0) = (2/3)*mu*(2*dudx-dvdy);
        stressy(i,j,1,1) = (2/3)*mu*(2*dvdy-dudx);
        stressy(i,j,0,1) = 0;
        stressy(i,j,1,0) = 0;
    }
};

struct applyViscousTerm2dv {
    
    typedef typename Kokkos::View<double****> V4D;
    typedef typename Kokkos::View<double**> V2D;
    V4D dvar;
    V4D var;
    V2D rho;
    V4D stressx;
    V4D stressy;
    Kokkos::View<double*> cd;

    applyViscousTerm2dv (V4D dvar_, V4D var_, V2D rho_, V4D strx_, V4D stry_, Kokkos::View<double*> cd_)
        : dvar(dvar_), var(var_), rho(rho_), stressx(strx_), stressy(stry_), cd(cd_) {}
    
    KOKKOS_INLINE_FUNCTION
    void operator()(const int i, const int j) const {
        double dx = cd(1);
        double dy = cd(2);

        double ur = (var(i+1,j,0,0)/rho(i,j) - var(i,j,0,0)/rho(i,j))/2;
        double ul = (var(i,j,0,0)/rho(i,j) - var(i-1,j,0,0)/rho(i,j))/2;
        double vt = (var(i,j+1,0,1)/rho(i,j) - var(i,j,0,1)/rho(i,j))/2;
        double vb = (var(i,j,0,1)/rho(i,j) - var(i,j-1,0,1)/rho(i,j))/2;

        dvar(i,j,0,0) += (stressx(i,j,0,0)+stressx(i,j,0,0))/dx;
        dvar(i,j,0,1) += (stressx(i,j,1,1)+stressx(i,j,1,1))/dy;
        dvar(i,j,0,2) += (ur*stressx(i,j,0,0)-ul*stressx(i-1,j,0,0))/dx
                        +(vt*stressy(i,j,0,1)-vb*stressy(i-1,j,0,0))/dy;
    }
};


hydro2dvisc_func::hydro2dvisc_func(struct inputConfig &cf_, Kokkos::View<double*> & cd_):rk_func(cf_,cd_){};

void hydro2dvisc_func::compute(const Kokkos::View<double****> & mvar, Kokkos::View<double****> & mdvar){

    // Typename acronyms for 2D and 4D variables
    typedef typename Kokkos::View<double****> V4D;
    typedef typename Kokkos::View<double**> V2D;

    // Copy input and output views
    V4D var = mvar;
    V4D dvar = mdvar;

    // Copy Configuration Data
    Kokkos::View<double*> cd = mcd;

    /*** Temporary Views ***/
    V2D p("p",cf.ngi,cf.ngj);          // Pressure
    V2D rho("rho",cf.ngi,cf.ngj);      // Total Density
    V2D wenox("wenox",cf.ngi,cf.ngj);  // Weno Fluxes in X direction
    V2D wenoy("wenoy",cf.ngi,cf.ngj);  // Weno Fluxes in Y direction
    V4D stressx("stressx",cf.ngi,cf.ngj,2,2);  // stress tensor on x faces
    V4D stressy("stressy",cf.ngi,cf.ngj,2,2);  // stress tensor on y faces

    /*** Range Policies ***/

    // Physical and Ghost cells
    policy_f ghost_pol = policy_f({0,0},{cf.ngi, cf.ngj});
    // Physical Cells only
    policy_f cell_pol  = policy_f({cf.ng,cf.ng},{cf.ngi-cf.ng, cf.ngj-cf.ng});
    // Cell Faces
    policy_f weno_pol  = policy_f({cf.ng-1,cf.ng-1},{cf.ngi-cf.ng, cf.ngj-cf.ng});


    /**** WENO ****/

    // Calcualte Total Density and Pressure Fields
    Kokkos::parallel_for( ghost_pol, calculateRhoAndPressure2dv(var,p,rho,cd) );

    // Calculate and apply weno fluxes for each variable
    for (int v=0; v<cf.nv; ++v){
        Kokkos::parallel_for( cell_pol, advect2dv(dvar,var,p,rho,cd,v) );
    }

    // Apply Pressure Gradient Term
    Kokkos::parallel_for( cell_pol, applyPressure2dv(dvar,p,cd) );

    Kokkos::parallel_for( weno_pol, calculateStressTensor2dv(var,rho,stressx,stressy,cd) );

    Kokkos::parallel_for( cell_pol, applyViscousTerm2dv(dvar,var,rho,stressx,stressy,cd) );
}
