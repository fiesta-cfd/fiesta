#include "fiesta.hpp"
#include "input.hpp"
#include "mpi.hpp"
#include "cgns.hpp"
#include "bc.hpp"
#include "Kokkos_Core.hpp"
#include <mpi.h>
#include "debug.hpp"
#include "hydro2d.hpp"
#include <map>
#include <string>
#include "timer.hpp"
#include "flux.hpp"

struct calculateRhoAndPressure2d {
    FS4D var;
    FS2D p;
    FS2D rho;
    Kokkos::View<double*> cd;

    calculateRhoAndPressure2d (FS4D var_, FS2D p_, FS2D rho_, Kokkos::View<double*> cd_)
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

struct applyWenoFluxes2d {
    
    FS4D dvar;
    FS2D wenox;
    FS2D wenoy;
    int v;
    

    applyWenoFluxes2d (FS4D dvar_, FS2D wenox_, FS2D wenoy_, int v_)
        : dvar(dvar_), wenox(wenox_), wenoy(wenoy_), v(v_) {}
    
    KOKKOS_INLINE_FUNCTION
    void operator()(const int i, const int j) const {

        //apply weno fluxes to right hand side of Euler equation dV/dt = ...
        dvar(i,j,0,v) = -( (wenox(i,j) - wenox(i-1,j))
                          +(wenoy(i,j) - wenoy(i,j-1)) );
    }
};

struct applyPressure2d {
    
    FS4D dvar;
    FS2D p;
    Kokkos::View<double*> cd;

    applyPressure2d (FS4D dvar_, FS2D p_, Kokkos::View<double*> cd_)
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

hydro2d_func::hydro2d_func(struct inputConfig &cf_, Kokkos::View<double*> & cd_):rk_func(cf_,cd_){

    using namespace std;
    var     = Kokkos::View<double****,FS_LAYOUT>("var",cf.ngi,cf.ngj,cf.ngk,cf.nv); // Primary Variable Array
    tmp1    = Kokkos::View<double****,FS_LAYOUT>("tmp1",cf.ngi,cf.ngj,cf.ngk,cf.nv); // Temporary Variable Arrayr1
    tmp2    = Kokkos::View<double****,FS_LAYOUT>("tmp2",cf.ngi,cf.ngj,cf.ngk,cf.nv); // Temporary Variable Array2
    dvar    = Kokkos::View<double****,FS_LAYOUT>("dvar",cf.ngi,cf.ngj,cf.ngk,cf.nv); // RHS Output
    p       = Kokkos::View<double**,FS_LAYOUT>("p",cf.ngi,cf.ngj);             // Pressure
    rho     = Kokkos::View<double**,FS_LAYOUT>("rho",cf.ngi,cf.ngj);           // Total Density
    wenox   = Kokkos::View<double**,FS_LAYOUT>("wenox",cf.ngi,cf.ngj);           // Total Density
    wenoy   = Kokkos::View<double**,FS_LAYOUT>("wenoy",cf.ngi,cf.ngj);           // Total Density
    cd = mcd;

    //timers.insert( make_pair("weno",fiestaTimer()) );
    //timers.emplace("weno",fiestaTimer());
    timers["flux"] = fiestaTimer("Flux Calculation");
    timers["pressgrad"] = fiestaTimer("Pressure Gradient Calculation");
    timers["calcRhoP"] = fiestaTimer("Secondary Variable Calculation");
    //timers.insert( make_pair<string, class fiestaTimer>("weno",fiestaTimer()) );

};

//void hydro2d_func::compute(const FS4D & mvar, FS4D & mdvar){
void hydro2d_func::compute(){
    policy_f ghost_pol = policy_f({0,0},{cf.ngi, cf.ngj});
    policy_f cell_pol  = policy_f({cf.ng,cf.ng},{cf.ngi-cf.ng, cf.ngj-cf.ng});
    policy_f face_pol  = policy_f({cf.ng-1,cf.ng-1},{cf.ngi-cf.ng, cf.ngj-cf.ng});

    /**** WENO ****/

    // Calcualte Total Density and Pressure Fields
    timers["calcRhoP"].reset();
    Kokkos::parallel_for( ghost_pol, calculateRhoAndPressure2d(var,p,rho,cd) );
    timers["calcRhoP"].accumulate();

    // Calculate and apply weno fluxes for each variable
    for (int v=0; v<cf.nv; ++v){
        timers["flux"].reset();

        Kokkos::parallel_for( face_pol, calculateFluxWeno2D(var,p,rho,wenox,wenoy,cd,v) );
        Kokkos::parallel_for( cell_pol, applyWenoFluxes2d(dvar,wenox,wenoy,v) );

        timers["flux"].accumulate();
    }

    // Apply Pressure Gradient Term
    timers["pressgrad"].reset();
    Kokkos::parallel_for( cell_pol, applyPressure2d(dvar,p,cd) );
    timers["pressgrad"].accumulate();
}
