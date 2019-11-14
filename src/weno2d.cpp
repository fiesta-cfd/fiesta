#include "input.hpp"
#include "mpi.hpp"
#include "cgns.hpp"
#include "bc.hpp"
#include "Kokkos_Core.hpp"
#include <mpi.h>
#include "lsdebug.hpp"
#include "weno2d.hpp"

struct calculateRhoAndPressure2D {
    typedef typename Kokkos::View<double****> V4D;
    typedef typename Kokkos::View<double**> V2D;
    V4D var;
    V2D p;
    V2D rho;
    Kokkos::View<double*> cd;

    calculateRhoAndPressure2D (V4D var_, V2D p_, V2D rho_, Kokkos::View<double*> cd_)
         : var(var_), p(p_), rho(rho_), cd(cd_) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const int i, const int j) const {

        int ns = (int)cd(0);
        int nv = (int)cd(4);
        double gamma, gammas, Rs;
        double Cp = 0;
        double Cv = 0;

        rho(i,j) = 0.0;

        for (int s=0; s<ns; ++s){
            rho(i,j) = rho(i,j) + var(i,j,0,3+s);
        }

        for (int s=0; s<ns; ++s){
            gammas = cd(5+2*s);
            Rs = cd(5+2*s+1);

            Cp = Cp + (var(i,j,0,3+s)/rho(i,j))*( gammas*Rs/(gammas-1) );
            Cv = Cv + (var(i,j,0,3+s)/rho(i,j))*( Rs/(gammas-1) );
        }

        gamma = Cp/Cv;

        p(i,j) = (gamma-1)*( var(i,j,0,2) - (0.5/rho(i,j))
                  *(var(i,j,0,0)*var(i,j,0,0) + var(i,j,0,1)*var(i,j,0,1)) );
    }
};

struct calculateWenoFluxes2D {
    
    typedef typename Kokkos::View<double****> V4D;
    typedef typename Kokkos::View<double**> V2D;
    V4D var;
    V2D p;
    V2D rho;
    V2D wenox;
    V2D wenoy;
    Kokkos::View<double*> cd;
    double eps = 0.000001;

    calculateWenoFluxes2D (V4D var_, V2D p_, V2D rho_,
                         V2D wenox_, V2D wenoy_, Kokkos::View<double*> cd_)
                         : var(var_), p(p_), rho(rho_),
                           wenox(wenox_), wenoy(wenoy_), cd(cd_) {}
    
    KOKKOS_INLINE_FUNCTION
    void operator()(const int i, const int j, const int v) const {

        int ns = (int)cd(0);
        double ur,vr,w,b1,b2,b3,w1,w2,w3,p1,p2,p3,f1,f2,f3,f4,f5;

        ur = ( -     var(i+2,j,0,0)/rho(i+2,j) + 7.0*var(i+1,j,0,0)/rho(i+1,j)
               + 7.0*var(i  ,j,0,0)/rho(i  ,j) -     var(i-1,j,0,0)/rho(i-1,j) )/12.0;

        vr = ( -     var(i,j+2,0,1)/rho(i,j+2) + 7.0*var(i,j+1,0,1)/rho(i,j+1)
               + 7.0*var(i,j  ,0,1)/rho(i,j  ) -     var(i,j-1,0,1)/rho(i,j-1) )/12.0;

        for (int idx=0; idx<2; ++idx){
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

            if (idx == 0){
                wenox(i,j) = ur*w/cd(1);
            }
            if (idx == 1) {
                wenoy(i,j) = vr*w/cd(2);
            }
        }
    }
};

struct applyWenoFluxes2D {
    
    typedef typename Kokkos::View<double****> V4D;
    typedef typename Kokkos::View<double**> V2D;
    V4D dvar;
    V2D wenox;
    V2D wenoy;
    

    applyWenoFluxes2D (V4D dvar_, V2D wenox_, V2D wenoy_)
        : dvar(dvar_), wenox(wenox_), wenoy(wenoy_) {}
    
    KOKKOS_INLINE_FUNCTION
    void operator()(const int i, const int j, const int v) const {

        dvar(i,j,0,v) = -( (wenox(i,j) - wenox(i-1,j))
                          +(wenoy(i,j) - wenoy(i,j-1)) );
    }
};

struct applyPressure2D {
    
    typedef typename Kokkos::View<double****> V4D;
    typedef typename Kokkos::View<double**> V2D;
    V4D dvar;
    V2D p;
    Kokkos::View<double*> cd;

    applyPressure2D (V4D dvar_, V2D p_, Kokkos::View<double*> cd_)
        : dvar(dvar_), p(p_), cd(cd_) {}
    
    KOKKOS_INLINE_FUNCTION
    void operator()(const int i, const int j) const {
        double dxp = ( p(i-2,j) - 8.0*p(i-1,j) + 8.0*p(i+1,j) - p(i+2,j) )/(12.0*cd(1));
        double dyp = ( p(i,j-2) - 8.0*p(i,j-1) + 8.0*p(i,j+1) - p(i,j+2) )/(12.0*cd(2));

        dvar(i,j,0,0) = dvar(i,j,0,0) - dxp;
        dvar(i,j,0,1) = dvar(i,j,0,1) - dyp;
    }
};

weno2d_func::weno2d_func(struct inputConfig &cf_, const Kokkos::View<double****> & u_,
                     Kokkos::View<double****> & k_, Kokkos::View<double*> & cd_)
                     : cf(cf_) , mvar(u_), mdvar(k_) , mcd(cd_) {};

void weno2d_func::operator()() {

    // Typename acronyms for 3D and 4D variables
    typedef typename Kokkos::View<double****> V4D;
    typedef typename Kokkos::View<double**> V2D;

    // copy input views
    V4D dvar = mdvar;
    V4D var = mvar;

    // create configuration data view
    Kokkos::View<double*> cd = mcd;

    // create temprary views
    V2D p("p",cf.ngi,cf.ngj);
    V2D rho("rho",cf.ngi,cf.ngj);
    V2D wenox("wenox",cf.ngi,cf.ngj);
    V2D wenoy("wenoy",cf.ngi,cf.ngj);

    // create range policies
    policy_f2 ghost_pol2 = policy_f2({0,0},{cf.ngi, cf.ngj});
    policy_f2 cell_pol2  = policy_f2({cf.ng,cf.ng},{cf.ngi-cf.ng, cf.ngj-cf.ng});
    policy_f2 weno_pol2  = policy_f2({cf.ng-1,cf.ng-1},{cf.ngi-cf.ng, cf.ngj-cf.ng});
    policy_fv2 cell_pol_v2 = policy_fv2({cf.ng,cf.ng},{cf.ngi-cf.ng, cf.ngj-cf.ng,cf.nv});
    policy_fv2 weno_pol_v2 = policy_fv2({cf.ng-1,cf.ng-1},{cf.ngi-cf.ng, cf.ngj-cf.ng,cf.nv});


    /**** WENO ****/
    Kokkos::parallel_for( ghost_pol2, calculateRhoAndPressure2D(var,p,rho,cd) );

    Kokkos::parallel_for( weno_pol_v2, calculateWenoFluxes2D(var,p,rho,wenox,wenoy,cd) );

    Kokkos::parallel_for( cell_pol_v2, applyWenoFluxes2D(dvar,wenox,wenoy) );

    Kokkos::parallel_for( cell_pol2, applyPressure2D(dvar,p,cd) );
}
