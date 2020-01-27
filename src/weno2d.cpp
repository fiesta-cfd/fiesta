#include "input.hpp"
#include "mpi.hpp"
#include "cgns.hpp"
#include "bc.hpp"
#include "Kokkos_Core.hpp"
#include <mpi.h>
#include "lsdebug.hpp"
#include "weno2d.hpp"

struct calculateRhoAndPressure2d {
    typedef typename Kokkos::View<double****> V4D;
    typedef typename Kokkos::View<double**> V2D;
    V4D var;
    V2D p;
    V2D rho;
    Kokkos::View<double*> cd;

    calculateRhoAndPressure2d (V4D var_, V2D p_, V2D rho_, Kokkos::View<double*> cd_)
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

struct calculateWenoFluxes2d {
    
    typedef typename Kokkos::View<double****> V4D;
    typedef typename Kokkos::View<double**> V2D;
    V4D var;
    V2D p;
    V2D rho;
    V2D wenox;
    V2D wenoy;
    Kokkos::View<double*> cd;
    int v;
    double eps = 0.000001;

    calculateWenoFluxes2d (V4D var_, V2D p_, V2D rho_,
                         V2D wenox_, V2D wenoy_, Kokkos::View<double*> cd_, int v_)
                         : var(var_), p(p_), rho(rho_),
                           wenox(wenox_), wenoy(wenoy_), cd(cd_), v(v_) {}
    
    KOKKOS_INLINE_FUNCTION
    void operator()(const int i, const int j) const {

        int ns = (int)cd(0);
        double ur,vr,w,b1,b2,b3,w1,w2,w3,p1,p2,p3,f1,f2,f3,f4,f5;
        double dx = cd(1);
        double dy = cd(2);

        //calculate cell face velocities (in positive direction) with 4th order interpolation
        //velocity is momentum divided by total density
        ur = ( -     var(i+2,j,0,0)/rho(i+2,j) + 7.0*var(i+1,j,0,0)/rho(i+1,j)
               + 7.0*var(i  ,j,0,0)/rho(i  ,j) -     var(i-1,j,0,0)/rho(i-1,j) )/12.0;

        vr = ( -     var(i,j+2,0,1)/rho(i,j+2) + 7.0*var(i,j+1,0,1)/rho(i,j+1)
               + 7.0*var(i,j  ,0,1)/rho(i,j  ) -     var(i,j-1,0,1)/rho(i,j-1) )/12.0;

        //for each direction
        for (int idx=0; idx<2; ++idx){
            //get stencil data.  the flux for the energy equation includes pressure so only add 
            //pressure for the energy variable (index 2 for 2d problem)
            if (idx == 0){
                if (ur < 0.0){
                    f1 = var(i+3,j,0,v) + (v==2)*p(i+3,j);
                    f2 = var(i+2,j,0,v) + (v==2)*p(i+2,j);
                    f3 = var(i+1,j,0,v) + (v==2)*p(i+1,j);
                    f4 = var(i  ,j,0,v) + (v==2)*p(i  ,j);
                    f5 = var(i-1,j,0,v) + (v==2)*p(i-1,j);
                }else{
                    f1 = var(i-2,j,0,v) + (v==2)*p(i-2,j);
                    f2 = var(i-1,j,0,v) + (v==2)*p(i-1,j);
                    f3 = var(i  ,j,0,v) + (v==2)*p(i  ,j);
                    f4 = var(i+1,j,0,v) + (v==2)*p(i+1,j);
                    f5 = var(i+2,j,0,v) + (v==2)*p(i+2,j);
                }
            } 
            if (idx == 1) {
                if (vr < 0.0){
                    f1 = var(i,j+3,0,v) + (v==2)*p(i,j+3);
                    f2 = var(i,j+2,0,v) + (v==2)*p(i,j+2);
                    f3 = var(i,j+1,0,v) + (v==2)*p(i,j+1);
                    f4 = var(i,j  ,0,v) + (v==2)*p(i,j  );
                    f5 = var(i,j-1,0,v) + (v==2)*p(i,j-1);
                }else{
                    f1 = var(i,j-2,0,v) + (v==2)*p(i,j-2);
                    f2 = var(i,j-1,0,v) + (v==2)*p(i,j-1);
                    f3 = var(i,j  ,0,v) + (v==2)*p(i,j  );
                    f4 = var(i,j+1,0,v) + (v==2)*p(i,j+1);
                    f5 = var(i,j+2,0,v) + (v==2)*p(i,j+2);
                }
            }

            // calculate weights and other weno stuff
            b1 = (13/12)*pow((f1-2.0*f2+f3),2.0) + (0.25)*pow((f1-4.0*f2+3.0*f3),2.0);
            b2 = (13/12)*pow((f2-2.0*f3+f4),2.0) + (0.25)*pow((f2-f4),2.0);
            b3 = (13/12)*pow((f3-2.0*f4+f5),2.0) + (0.25)*pow((3.0*f3-4.0*f4+f5),2.0);

            w1 = (0.1)/pow((eps+b1),2.0);
            w2 = (0.6)/pow((eps+b2),2.0);
            w3 = (0.3)/pow((eps+b3),2.0);

            p1 = ( 1.0/3.0)*f1 + (-7.0/6.0)*f2 + (11.0/6.0)*f3;
            p2 = (-1.0/6.0)*f2 + ( 5.0/6.0)*f3 + ( 1.0/3.0)*f4;
            p3 = ( 1.0/3.0)*f3 + ( 5.0/6.0)*f4 + (-1.0/6.0)*f5;

            w = (w1*p1+w2*p2+w3*p3)/(w1+w2+w3);

            //calculate weno flux
            if (idx == 0){
                wenox(i,j) = ur*w/dx;
            }
            if (idx == 1) {
                wenoy(i,j) = vr*w/dy;
            }
        }
    }
};

struct applyWenoFluxes2d {
    
    typedef typename Kokkos::View<double****> V4D;
    typedef typename Kokkos::View<double**> V2D;
    V4D dvar;
    V2D wenox;
    V2D wenoy;
    int v;
    

    applyWenoFluxes2d (V4D dvar_, V2D wenox_, V2D wenoy_, int v_)
        : dvar(dvar_), wenox(wenox_), wenoy(wenoy_), v(v_) {}
    
    KOKKOS_INLINE_FUNCTION
    void operator()(const int i, const int j) const {

        //apply weno fluxes to right hand side of Euler equation dV/dt = ...
        dvar(i,j,0,v) = -( (wenox(i,j) - wenox(i-1,j))
                          +(wenoy(i,j) - wenoy(i,j-1)) );
    }
};

struct applyPressure2d {
    
    typedef typename Kokkos::View<double****> V4D;
    typedef typename Kokkos::View<double**> V2D;
    V4D dvar;
    V2D p;
    Kokkos::View<double*> cd;

    applyPressure2d (V4D dvar_, V2D p_, Kokkos::View<double*> cd_)
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

weno2d_func::weno2d_func(struct inputConfig &cf_, Kokkos::View<double*> & cd_):rk_func(cf_,cd_){};

void weno2d_func::compute(const Kokkos::View<double****> & mvar, Kokkos::View<double****> & mdvar){

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

    /*** Range Policies ***/

    // Physical and Ghost cells
    policy_f ghost_pol = policy_f({0,0},{cf.ngi, cf.ngj});
    // Physical Cells only
    policy_f cell_pol  = policy_f({cf.ng,cf.ng},{cf.ngi-cf.ng, cf.ngj-cf.ng});
    // Cell Faces
    policy_f weno_pol  = policy_f({cf.ng-1,cf.ng-1},{cf.ngi-cf.ng, cf.ngj-cf.ng});


    /**** WENO ****/

    // Calcualte Total Density and Pressure Fields
    Kokkos::parallel_for( ghost_pol, calculateRhoAndPressure2d(var,p,rho,cd) );

    // Calculate and apply weno fluxes for each variable
    for (int v=0; v<cf.nv; ++v){
        Kokkos::parallel_for( weno_pol, calculateWenoFluxes2d(var,p,rho,wenox,wenoy,cd,v) );

        Kokkos::parallel_for( cell_pol, applyWenoFluxes2d(dvar,wenox,wenoy,v) );
    }

    // Apply Pressure Gradient Term
    Kokkos::parallel_for( cell_pol, applyPressure2d(dvar,p,cd) );
}