#include "input.hpp"
#include "mpi.hpp"
#include "cgns.hpp"
#include "bc.hpp"
#include "Kokkos_Core.hpp"
#include <mpi.h>
#include "lsdebug.hpp"
#include "weno2d.hpp"

struct calculateRhoAndPressure2d {
    typedef typename Kokkos::View<double**> V4D;
    typedef typename Kokkos::View<double*> V2D;
    typedef typename Kokkos::View<int*> I2D;
    V4D var;
    V2D p;
    V2D rho;
    I2D nlft;
    I2D nrht;
    I2D nbot;
    I2D ntop;
    Kokkos::View<double*> cd;
    //int ngi, ngj;

    calculateRhoAndPressure2d (V4D var_, I2D nlft_, I2D nrht_, I2D nbot_, I2D ntop_, V2D p_, V2D rho_, Kokkos::View<double*> cd_)
         : var(var_), nlft(nlft_), nrht(nrht_), nbot(nbot_), ntop(ntop_), p(p_), rho(rho_), cd(cd_) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const int idx) const {

        int ns = (int)cd(0);
        int nv = (int)cd(4);
        double gamma, gammas, Rs;
        double Cp = 0;
        double Cv = 0;
        //int idx = i + j*ngi;

        rho(idx) = 0.0;

        // Total Density for this cell
        for (int s=0; s<ns; ++s){
            rho(idx) = rho(idx) + var(idx,3+s);
        }

        // Calculate mixture ratio of specific heats
        for (int s=0; s<ns; ++s){
            gammas = cd(5+2*s);
            Rs = cd(5+2*s+1);

            // accumulate mixture heat capacity by mass fraction weights
            Cp = Cp + (var(idx,3+s)/rho(idx))*( gammas*Rs/(gammas-1) );
            Cv = Cv + (var(idx,3+s)/rho(idx))*( Rs/(gammas-1) );
        }
        gamma = Cp/Cv;

        // calculate pressure assuming perfect gas
        p(idx) = (gamma-1)*( var(idx,2) - (0.5/rho(idx))
                  *(var(idx,0)*var(idx,0) + var(idx,1)*var(idx,1)) );
    }
};

struct calculateWenoFluxes2d {
    
    typedef typename Kokkos::View<double**> V4D;
    typedef typename Kokkos::View<double*> V2D;
    typedef typename Kokkos::View<int*> I2D;
    V4D var;
    V2D p;
    V2D rho;
    V2D wenox;
    V2D wenoy;
    I2D nlft;
    I2D nrht;
    I2D nbot;
    I2D ntop;
    I2D cell_type;
    Kokkos::View<double*> cd;
    int v;
    double eps = 0.000001;
    int ngi, ngj;

    calculateWenoFluxes2d (V4D var_, I2D nlft_, I2D nrht_, I2D nbot_, I2D ntop_, I2D cell_type_, V2D p_, V2D rho_,
                         V2D wenox_, V2D wenoy_, Kokkos::View<double*> cd_, int v_, int ngi_, int ngj_)
                         : var(var_), nlft(nlft_), nrht(nrht_), nbot(nbot_), ntop(ntop_), cell_type(cell_type_), p(p_), rho(rho_),
                           wenox(wenox_), wenoy(wenoy_), cd(cd_), v(v_), ngi(ngi_), ngj(ngj_) {}
    
    KOKKOS_INLINE_FUNCTION
    void operator()(const int index) const {
      if (cell_type(index) == REAL_CELL || cell_type(nrht(index)) == BORDER_CELL || cell_type(ntop(index)) == BORDER_CELL  || cell_type(index) == BORDER_CELL) {
        int ns = (int)cd(0);
        double ur,vr,w,b1,b2,b3,w1,w2,w3,p1,p2,p3,f1,f2,f3,f4,f5;
        double dx = cd(1);
        double dy = cd(2);
        //int index = i + j*ngi;

        //calculate cell face velocities (in positive direction) with 4th order interpolation
        //velocity is momentum divided by total density
        ur = ( -     var(nrht(nrht(index))+0,0)/rho(nrht(nrht(index))) + 7.0*var(nrht(index),0)/rho(nrht(index))
               + 7.0*var(index,0)/rho(index)                           -     var(nlft(index),0)/rho(nlft(index)) )/12.0;

        vr = ( -     var(ntop(ntop(index)),1)/rho(ntop(ntop(index))) + 7.0*var(ntop(index),1)/rho(ntop(index))
               + 7.0*var(index,1)/rho(index) -     var(nbot(index),1)/rho(nbot(index)) )/12.0;

        //for each direction
        for (int idx=0; idx<2; ++idx){
            //get stencil data.  the flux for the energy equation includes pressure so only add 
            //pressure for the energy variable (index 2 for 2d problem)
            if (idx == 0){
                if (ur < 0.0){
                    f1 = var(nrht(nrht(nrht(index)))+0,v) + (v==2)*p(nrht(nrht(nrht(index))));
                    f2 = var(nrht(nrht(index)),v) + (v==2)*p(nrht(nrht(index)));
                    f3 = var(nrht(index),v) + (v==2)*p(nrht(index));
                    f4 = var(index,v) + (v==2)*p(index);
                    f5 = var(nlft(index),v) + (v==2)*p(nlft(index));
                }else{
                    f1 = var(nlft(nlft(index)),v) + (v==2)*p(nlft(nlft(index)));
                    f2 = var(nlft(index),v) + (v==2)*p(nlft(index));
                    f3 = var(index,v) + (v==2)*p(index);
                    f4 = var(nrht(index),v) + (v==2)*p(nrht(index));
                    f5 = var(nrht(nrht(index)),v) + (v==2)*p(nrht(nrht(index)));
                }
            } 
            if (idx == 1) {
                if (vr < 0.0){
                    f1 = var(ntop(ntop(ntop(index))),v) + (v==2)*p(ntop(ntop(ntop(index))));
                    f2 = var(ntop(ntop(index)),v) + (v==2)*p(ntop(ntop(index)));
                    f3 = var(ntop(index),v) + (v==2)*p(ntop(index));
                    f4 = var(index,v) + (v==2)*p(index);
                    f5 = var(nbot(index),v) + (v==2)*p(nbot(index));
                }else{
                    f1 = var(nbot(nbot(index)),v) + (v==2)*p(nbot(nbot(index)));
                    f2 = var(nbot(index),v) + (v==2)*p(nbot(index));
                    f3 = var(index,v) + (v==2)*p(index);
                    f4 = var(ntop(index),v) + (v==2)*p(ntop(index));
                    f5 = var(ntop(ntop(index)),v) + (v==2)*p(ntop(ntop(index)));
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
                wenox(index) = ur*w/dx;
            }
            if (idx == 1) {
                wenoy(index) = vr*w/dy;
            }
        }
      }
    }
};

struct applyWenoFluxes2d {
    
    typedef typename Kokkos::View<double**> V4D;
    typedef typename Kokkos::View<double*> V2D;
    typedef typename Kokkos::View<int*> I2D;
    V4D dvar;
    V2D wenox;
    V2D wenoy;
    I2D nlft;
    I2D nrht;
    I2D nbot;
    I2D ntop;
    I2D cell_type;
    int v;
    //int ngi, ngj;
    

    applyWenoFluxes2d (V4D dvar_, I2D nlft_, I2D nrht_, I2D nbot_, I2D ntop_, I2D cell_type_, V2D wenox_, V2D wenoy_, int v_)
        : dvar(dvar_), nlft(nlft_), nrht(nrht_), nbot(nbot_), ntop(ntop_), cell_type(cell_type_), wenox(wenox_), wenoy(wenoy_), v(v_) {}
    
    KOKKOS_INLINE_FUNCTION
    void operator()(const int idx) const {

        //int idx = i + j*26;
        //apply weno fluxes to right hand side of Euler equation dV/dt = ...
        // if real
        if (cell_type(idx) != BOUNDARY_CELL) {
            dvar(idx,v) = -( (wenox(idx) - wenox(nlft(idx)))
                              +(wenoy(idx) - wenoy(nbot(idx))) );
        }
    }
};

struct applyPressure2d {
    
    typedef typename Kokkos::View<double**> V4D;
    typedef typename Kokkos::View<double*> V2D;
    typedef typename Kokkos::View<int*> I2D;
    V4D dvar;
    V2D p;
    I2D nlft;
    I2D nrht;
    I2D nbot;
    I2D ntop;
    I2D cell_type;
    Kokkos::View<double*> cd;
    //int ngi, ngj;

    applyPressure2d (V4D dvar_, I2D nlft_, I2D nrht_, I2D nbot_, I2D ntop_, I2D cell_type_, V2D p_, Kokkos::View<double*> cd_)
        : dvar(dvar_), nlft(nlft_), nrht(nrht_), nbot(nbot_), ntop(ntop_), cell_type(cell_type_), p(p_), cd(cd_) {}
    
    KOKKOS_INLINE_FUNCTION
    void operator()(const int idx) const {
        //int idx = i + j*26;
        // if real
        if (cell_type(idx) != BOUNDARY_CELL) {
            double dx = cd(1);
            double dy = cd(2);
            // calculate pressure gradient across cell in each direction using 4th order
            // central difference
            double dxp = ( p(nlft(nlft(idx))) - 8.0*p(nlft(idx)) + 8.0*p(nrht(idx)) - p(nrht(nrht(idx))) )/(12.0*dx);
            double dyp = ( p(nbot(nbot(idx))) - 8.0*p(nbot(idx)) + 8.0*p(ntop(idx)) - p(ntop(ntop(idx))) )/(12.0*dy);

            //apply pressure gradient term to right hand side of Euler equation dV/dt = ...
            dvar(idx,0) = dvar(idx,0) - dxp;
            dvar(idx,1) = dvar(idx,1) - dyp;
        }
    }
};

weno2d_func::weno2d_func(struct inputConfig &cf_, Kokkos::View<double*> & cd_):rk_func(cf_,cd_){};

void weno2d_func::compute(const Kokkos::View<double**> & mvar, Kokkos::View<int*> nlft, Kokkos::View<int*> nrht, Kokkos::View<int*> nbot, Kokkos::View<int*> ntop, Kokkos::View<int*> nbak, Kokkos::View<int*> nfrt, Kokkos::View<int*> cell_type, Kokkos::View<double**> & mdvar, int ncells){

    // Typename acronyms for 2D and 4D variables
    typedef typename Kokkos::View<double**> V4D;
    typedef typename Kokkos::View<double*> V2D;
    typedef typename Kokkos::View<int*> I2D;

    // Copy input and output views
    V4D var = mvar;
    V4D dvar = mdvar;
    I2D lft = nlft;
    I2D rht = nrht;
    I2D bot = nbot;
    I2D top = ntop;

    // Copy Configuration Data
    Kokkos::View<double*> cd = mcd;

    /*** Temporary Views ***/
    V2D p("p",cf.ngi*cf.ngj);          // Pressure
    V2D rho("rho",cf.ngi*cf.ngj);      // Total Density
    V2D wenox("wenox",cf.ngi*cf.ngj);  // Weno Fluxes in X direction
    V2D wenoy("wenoy",cf.ngi*cf.ngj);  // Weno Fluxes in Y direction

    /*** Range Policies ***/

    // Physical and Ghost cells
    policy_f1 ghost_pol = policy_f1({0},{ncells});
    // Physical Cells only
    policy_f cell_pol  = policy_f({cf.ng,cf.ng},{cf.ngi-cf.ng, cf.ngj-cf.ng});
    // Cell Faces
    policy_f weno_pol  = policy_f({cf.ng-1,cf.ng-1},{cf.ngi-cf.ng, cf.ngj-cf.ng});


    /**** WENO ****/

    // Calcualte Total Density and Pressure Fields
    Kokkos::parallel_for( ghost_pol, calculateRhoAndPressure2d(var,lft,rht,bot,top,p,rho,cd) );

    // Calculate and apply weno fluxes for each variable
    for (int v=0; v<cf.nv; ++v){
        Kokkos::parallel_for( ghost_pol, calculateWenoFluxes2d(var,lft,rht,bot,top,cell_type,p,rho,wenox,wenoy,cd,v,cf.ngi,cf.ngj) );

        Kokkos::parallel_for( ghost_pol, applyWenoFluxes2d(dvar,lft,rht,bot,top,cell_type,wenox,wenoy,v) );
    }

    // Apply Pressure Gradient Term
    Kokkos::parallel_for( ghost_pol, applyPressure2d(dvar,lft,rht,bot,top,cell_type,p,cd) );
}
