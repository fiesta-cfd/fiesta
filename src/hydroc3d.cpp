#include "fiesta.hpp"
#include "input.hpp"
#include "mpi.hpp"
#include "cgns.hpp"
#include "bc.hpp"
#include "Kokkos_Core.hpp"
#include <mpi.h>
#include "debug.hpp"
#include "hydroc3d.hpp"

struct calculateRhoAndPressure {
    FS4D var;
    FS3D p;
    FS3D rho;
    Kokkos::View<double*> cd;

    calculateRhoAndPressure (FS4D var_, FS3D p_, FS3D rho_, Kokkos::View<double*> cd_)
         : var(var_), p(p_), rho(rho_), cd(cd_) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const int i, const int j, const int k) const {

        int ns = (int)cd(0);
        int nv = (int)cd(4);
        double gamma, gammas, Rs;
        double Cp = 0;
        double Cv = 0;

        rho(i,j,k) = 0.0;

        for (int s=0; s<ns; ++s){
            rho(i,j,k) = rho(i,j,k) + var(i,j,k,4+s);
        }

        for (int s=0; s<ns; ++s){
            gammas = cd(5+2*s);
            Rs = cd(5+2*s+1);

            Cp = Cp + (var(i,j,k,4+s)/rho(i,j,k))*( gammas*Rs/(gammas-1) );
            Cv = Cv + (var(i,j,k,4+s)/rho(i,j,k))*( Rs/(gammas-1) );
        }

        gamma = Cp/Cv;

        p(i,j,k) = (gamma-1)*( var(i,j,k,3) - (0.5/rho(i,j,k))
                  *(var(i,j,k,0)*var(i,j,k,0) + var(i,j,k,1)*var(i,j,k,1) + var(i,j,k,2)*var(i,j,k,2)) );
    }
};

struct calculateWenoFluxes {
    
    FS4D var;
    FS3D p;
    FS3D rho;
    FS3D wenox;
    FS3D wenoy;
    FS3D wenoz;
    Kokkos::View<double*> cd;
    int v;
    double eps = 0.000001;

    calculateWenoFluxes (FS4D var_, FS3D p_, FS3D rho_,
                         FS3D wenox_, FS3D wenoy_, FS3D wenoz_, Kokkos::View<double*> cd_, int v_)
                         : var(var_), p(p_), rho(rho_),
                           wenox(wenox_), wenoy(wenoy_), wenoz(wenoz_), cd(cd_), v(v_) {}
    
    KOKKOS_INLINE_FUNCTION
    void operator()(const int i, const int j, const int k) const {

        int ns = (int)cd(0);
        double ur,vr,wr,w,b1,b2,b3,w1,w2,w3,p1,p2,p3,f1,f2,f3,f4,f5;

        ur = ( -     var(i+2,j,k,0)/rho(i+2,j,k) + 7.0*var(i+1,j,k,0)/rho(i+1,j,k)
               + 7.0*var(i  ,j,k,0)/rho(i  ,j,k) -     var(i-1,j,k,0)/rho(i-1,j,k) )/12.0;

        vr = ( -     var(i,j+2,k,1)/rho(i,j+2,k) + 7.0*var(i,j+1,k,1)/rho(i,j+1,k)
               + 7.0*var(i,j  ,k,1)/rho(i,j  ,k) -     var(i,j-1,k,1)/rho(i,j-1,k) )/12.0;

        wr = ( -     var(i,j,k+2,2)/rho(i,j,k+2) + 7.0*var(i,j,k+1,2)/rho(i,j,k+1)
               + 7.0*var(i,j,k  ,2)/rho(i,j,k  ) -     var(i,j,k-1,2)/rho(i,j,k-1) )/12.0;
        
        for (int idx=0; idx<3; ++idx){
            if (idx == 0){
                if (ur < 0.0){
                    f1 = var(i+3,j,k,v) + (v==3)*p(i+3,j,k);
                    f2 = var(i+2,j,k,v) + (v==3)*p(i+2,j,k);
                    f3 = var(i+1,j,k,v) + (v==3)*p(i+1,j,k);
                    f4 = var(i  ,j,k,v) + (v==3)*p(i  ,j,k);
                    f5 = var(i-1,j,k,v) + (v==3)*p(i-1,j,k);
                }else{                                
                    f1 = var(i-2,j,k,v) + (v==3)*p(i-2,j,k);
                    f2 = var(i-1,j,k,v) + (v==3)*p(i-1,j,k);
                    f3 = var(i  ,j,k,v) + (v==3)*p(i  ,j,k);
                    f4 = var(i+1,j,k,v) + (v==3)*p(i+1,j,k);
                    f5 = var(i+2,j,k,v) + (v==3)*p(i+2,j,k);
                }
            } 
            if (idx == 1) {
                if (vr < 0.0){
                    f1 = var(i,j+3,k,v) + (v==3)*p(i,j+3,k);
                    f2 = var(i,j+2,k,v) + (v==3)*p(i,j+2,k);
                    f3 = var(i,j+1,k,v) + (v==3)*p(i,j+1,k);
                    f4 = var(i,j  ,k,v) + (v==3)*p(i,j  ,k);
                    f5 = var(i,j-1,k,v) + (v==3)*p(i,j-1,k);
                }else{                                  
                    f1 = var(i,j-2,k,v) + (v==3)*p(i,j-2,k);
                    f2 = var(i,j-1,k,v) + (v==3)*p(i,j-1,k);
                    f3 = var(i,j  ,k,v) + (v==3)*p(i,j  ,k);
                    f4 = var(i,j+1,k,v) + (v==3)*p(i,j+1,k);
                    f5 = var(i,j+2,k,v) + (v==3)*p(i,j+2,k);
                }
            }
            if (idx == 2){
                if (wr < 0.0){
                    f1 = var(i,j,k+3,v) + (v==3)*p(i,j,k+3);
                    f2 = var(i,j,k+2,v) + (v==3)*p(i,j,k+2);
                    f3 = var(i,j,k+1,v) + (v==3)*p(i,j,k+1);
                    f4 = var(i,j,k  ,v) + (v==3)*p(i,j,k  );
                    f5 = var(i,j,k-1,v) + (v==3)*p(i,j,k-1);
                }else{                                    
                    f1 = var(i,j,k-2,v) + (v==3)*p(i,j,k-2);
                    f2 = var(i,j,k-1,v) + (v==3)*p(i,j,k-1);
                    f3 = var(i,j,k  ,v) + (v==3)*p(i,j,k  );
                    f4 = var(i,j,k+1,v) + (v==3)*p(i,j,k+1);
                    f5 = var(i,j,k+2,v) + (v==3)*p(i,j,k+2);
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
                wenox(i,j,k) = ur*w/cd(1);
            }
            if (idx == 1) {
                wenoy(i,j,k) = vr*w/cd(2);
            }
            if (idx == 2) {
                wenoz(i,j,k) = wr*w/cd(3);
            }
        }
    }
};

struct applyWenoFluxes {
    
    FS4D dvar;
    FS3D wenox;
    FS3D wenoy;
    FS3D wenoz;
    int v;

    applyWenoFluxes (FS4D dvar_, FS3D wenox_, FS3D wenoy_, FS3D wenoz_, int v_)
        : dvar(dvar_), wenox(wenox_), wenoy(wenoy_), wenoz(wenoz_), v(v_) {}
    
    KOKKOS_INLINE_FUNCTION
    void operator()(const int i, const int j, const int k) const {

        dvar(i,j,k,v) = -( (wenox(i,j,k) - wenox(i-1,j,k))
                          +(wenoy(i,j,k) - wenoy(i,j-1,k))
                          +(wenoz(i,j,k) - wenoz(i,j,k-1)) );
    }
};

struct applyPressure {
    
    FS4D dvar;
    FS3D p;
    Kokkos::View<double*> cd;

    applyPressure (FS4D dvar_, FS3D p_, Kokkos::View<double*> cd_)
        : dvar(dvar_), p(p_), cd(cd_) {}
    
    KOKKOS_INLINE_FUNCTION
    void operator()(const int i, const int j, const int k) const {
        double dxp = ( p(i-2,j,k) - 8.0*p(i-1,j,k) + 8.0*p(i+1,j,k) - p(i+2,j,k) )/(12.0*cd(1));
        double dyp = ( p(i,j-2,k) - 8.0*p(i,j-1,k) + 8.0*p(i,j+1,k) - p(i,j+2,k) )/(12.0*cd(2));
        double dzp = ( p(i,j,k-2) - 8.0*p(i,j,k-1) + 8.0*p(i,j,k+1) - p(i,j,k+2) )/(12.0*cd(3));

        dvar(i,j,k,0) = dvar(i,j,k,0) - dxp;
        dvar(i,j,k,1) = dvar(i,j,k,1) - dyp;
        dvar(i,j,k,2) = dvar(i,j,k,2) - dzp;
    }
};

struct maxWaveSpeed {
    
    FS4D var;
    FS3D p;
    FS3D rho;
    Kokkos::View<double*> cd;

    maxWaveSpeed (FS4D var_, FS3D p_, FS3D rho_, Kokkos::View<double*> cd_)
        : var(var_), p(p_), rho(rho_), cd(cd_) {}
    
    KOKKOS_INLINE_FUNCTION
    void operator()(const int i, const int j, const int k, double& lmax) const {
        int ns = (int)cd(0);
        double gamma, gammas, Rs;
        double Cp = 0;
        double Cv = 0;

        double a, s;

        for (int s=0; s<ns; ++s){
            gammas = cd(5+2*s);
            Rs = cd(5+2*s+1);

            Cp = Cp + (var(i,j,k,4+s)/rho(i,j,k))*( gammas*Rs/(gammas-1) );
            Cv = Cv + (var(i,j,k,4+s)/rho(i,j,k))*( Rs/(gammas-1) );
        }

        gamma = Cp/Cv;

        a = sqrt(gamma*p(i,j,k)/rho(i,j,k));
        s = a + sqrt(var(i,j,k,0)*var(i,j,k,0) + var(i,j,k,1)*var(i,j,k,1) + var(i,j,k,2)*var(i,j,k,2))/rho(i,j,k);

        if (s > lmax)
            lmax = s;
    }
};

struct calculateRhoGrad {

    FS4D var;
    FS3D rho;
    FS4D gradRho;
    Kokkos::View<double*> cd;

    calculateRhoGrad (FS4D var_, FS3D rho_, FS4D gradRho_, Kokkos::View<double*> cd_)
        : var(var_), rho(rho_), gradRho(gradRho_), cd(cd_) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const int i, const int j, const int k) const {

        int indicator = 0;

        double u1 = var(i-2,j,k,0)/rho(i-2,j,k);
        double u2 = var(i-1,j,k,0)/rho(i-1,j,k);
        double u3 = var(i+1,j,k,0)/rho(i+1,j,k);
        double u4 = var(i+2,j,k,0)/rho(i+2,j,k);

        double v1 = var(i,j-2,k,0)/rho(i,j-2,k);
        double v2 = var(i,j-1,k,0)/rho(i,j-1,k);
        double v3 = var(i,j+1,k,0)/rho(i,j+1,k);
        double v4 = var(i,j+2,k,0)/rho(i,j+2,k);

        double w1 = var(i,j,k-2,0)/rho(i,j,k-2);
        double w2 = var(i,j,k-1,0)/rho(i,j,k-1);
        double w3 = var(i,j,k+1,0)/rho(i,j,k+1);
        double w4 = var(i,j,k+2,0)/rho(i,j,k+2);

        double ex1 = var(i-2,j,k,3)/rho(i,j,k-2)-0.5*(1/(rho(i-2,j,k)*rho(i-2,j,k)))*
            (var(i-2,j,k,0)*var(i-2,j,k,0) + var(i-2,j,k,1)*var(i-2,j,k,1) + var(i-2,j,k,2)*var(i-2,j,k,2));
        double ex2 = var(i-1,j,k,3)/rho(i,j,k-1)-0.5*(1/(rho(i-1,j,k)*rho(i-1,j,k)))*
            (var(i-1,j,k,0)*var(i-1,j,k,0) + var(i-1,j,k,1)*var(i-1,j,k,1) + var(i-1,j,k,2)*var(i-1,j,k,2));
        double ex3 = var(i+1,j,k,3)/rho(i,j,k+1)-0.5*(1/(rho(i+1,j,k)*rho(i+1,j,k)))*
            (var(i+1,j,k,0)*var(i+1,j,k,0) + var(i+1,j,k,1)*var(i+1,j,k,1) + var(i+1,j,k,2)*var(i+1,j,k,2));
        double ex4 = var(i+2,j,k,3)/rho(i,j,k+2)-0.5*(1/(rho(i+2,j,k)*rho(i+2,j,k)))*
            (var(i+2,j,k,0)*var(i+2,j,k,0) + var(i+2,j,k,1)*var(i+2,j,k,1) + var(i+2,j,k,2)*var(i+2,j,k,2));

        double ey1 = var(i,j-2,k,3)/rho(i,j,k-2)-0.5*(1/(rho(i,j-2,k)*rho(i,j-2,k)))*
            (var(i,j-2,k,0)*var(i,j-2,k,0) + var(i,j-2,k,1)*var(i,j-2,k,1) + var(i,j-2,k,2)*var(i,j-2,k,2));
        double ey2 = var(i,j-1,k,3)/rho(i,j,k-1)-0.5*(1/(rho(i,j-1,k)*rho(i,j-1,k)))*
            (var(i,j-1,k,0)*var(i,j-1,k,0) + var(i,j-1,k,1)*var(i,j-1,k,1) + var(i,j-1,k,2)*var(i,j-1,k,2));
        double ey3 = var(i,j+1,k,3)/rho(i,j,k+1)-0.5*(1/(rho(i,j+1,k)*rho(i,j+1,k)))*
            (var(i,j+1,k,0)*var(i,j+1,k,0) + var(i,j+1,k,1)*var(i,j+1,k,1) + var(i,j+1,k,2)*var(i,j+1,k,2));
        double ey4 = var(i,j+2,k,3)/rho(i,j,k+2)-0.5*(1/(rho(i,j+2,k)*rho(i,j+2,k)))*
            (var(i,j+2,k,0)*var(i,j+2,k,0) + var(i,j+2,k,1)*var(i,j+2,k,1) + var(i,j+2,k,2)*var(i,j+2,k,2));

        double ez1 = var(i,j,k-2,3)/rho(i,j,k-2)-0.5*(1/(rho(i,j,k-2)*rho(i,j,k-2)))*
            (var(i,j,k-2,0)*var(i,j,k-2,0) + var(i,j,k-2,1)*var(i,j,k-2,1) + var(i,j,k-2,2)*var(i,j,k-2,2));
        double ez2 = var(i,j,k-1,3)/rho(i,j,k-1)-0.5*(1/(rho(i,j,k-1)*rho(i,j,k-1)))*
            (var(i,j,k-1,0)*var(i,j,k-1,0) + var(i,j,k-1,1)*var(i,j,k-1,1) + var(i,j,k-1,2)*var(i,j,k-1,2));
        double ez3 = var(i,j,k+1,3)/rho(i,j,k+1)-0.5*(1/(rho(i,j,k+1)*rho(i,j,k+1)))*
            (var(i,j,k+1,0)*var(i,j,k+1,0) + var(i,j,k+1,1)*var(i,j,k+1,1) + var(i,j,k+1,2)*var(i,j,k+1,2));
        double ez4 = var(i,j,k+2,3)/rho(i,j,k+2)-0.5*(1/(rho(i,j,k+2)*rho(i,j,k+2)))*
            (var(i,j,k+2,0)*var(i,j,k+2,0) + var(i,j,k+2,1)*var(i,j,k+2,1) + var(i,j,k+2,2)*var(i,j,k+2,2));


        double dxr = (rho(i-2,j,k) - 8.0*rho(i-1,j,k) + 8.0*rho(i+1,j,k) - rho(i+2,j,k))/(12.0*cd(1));
        double dyr = (rho(i,j-2,k) - 8.0*rho(i,j-1,k) + 8.0*rho(i,j+1,k) - rho(i,j+2,k))/(12.0*cd(2));
        double dzr = (rho(i,j,k-2) - 8.0*rho(i,j,k-1) + 8.0*rho(i,j,k+1) - rho(i,j,k+2))/(12.0*cd(3));

        double dxe = (ex1 - 8.0*ex2 + 8.0*ex3 - ex4)/(12.0*cd(1));
        double dye = (ey1 - 8.0*ey2 + 8.0*ey3 - ey4)/(12.0*cd(2));
        double dze = (ez1 - 8.0*ez2 + 8.0*ez3 - ez4)/(12.0*cd(3));

        double dxu = (u1 - 8.0*u2 + 8.0*u3 - u4)/(12.0*cd(1));
        double dyv = (v1 - 8.0*v2 + 8.0*v3 - v4)/(12.0*cd(2));
        double dzw = (w1 - 8.0*w2 + 8.0*w3 - w4)/(12.0*cd(3));

        double n1 = dxr;
        double n2 = dyr;
        double n3 = dzr;

        double rgrad = sqrt(dxr*dxr+dyr*dyr+dzr*dzr);
        double divu = dxu+dyv+dzw;

        double dnednr = (n1*dxe+n2*dye+n3*dze)*(n1*dxr+n2*dyr+n3*dzr);

        //compression switch
        if (dnednr <= 0)
            indicator = 1;
        else
            indicator = 0;
        
        // detect shock front (equation 5a)
        if (divu <= 0)
            gradRho(i,j,k,0) = (1-indicator)*rgrad;
            //gradRho(i,j,k,0) = (1-indicator)*divu*rgrad;
        else
            gradRho(i,j,k,0) = 0;

        // detect contact surface
        gradRho(i,j,k,1) = indicator*rgrad;

        // gradient components
        gradRho(i,j,k,2) = dxr;
        gradRho(i,j,k,3) = dyr;
        gradRho(i,j,k,4) = dzr;
    }
};

struct maxGradFunctor {

    FS4D gradRho;
    int n;

    maxGradFunctor(FS4D gradRho_, int n_) : gradRho(gradRho_), n(n_) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const int i, const int j, const int k, double& lmax) const {
        if (abs(gradRho(i,j,k,n)) > lmax)
            lmax = abs(gradRho(i,j,k,n));
    }
};

struct updateCeq {

    FS4D dvar;
    FS4D var;
    FS4D gradRho;
    double maxS,kap,eps;
    Kokkos::View<double*> cd;

    updateCeq (FS4D dvar_, FS4D var_, FS4D gradRho_, double maxS_, Kokkos::View<double*> cd_,
                      double kap_, double eps_)
        : dvar(dvar_), var(var_), gradRho(gradRho_), maxS(maxS_), cd(cd_), kap(kap_), eps(eps_) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const int i, const int j, const int k) const {
        double dx [3];
        dx[0] = cd(1);
        dx[1] = cd(2);
        dx[2] = cd(3);

        // calculate cequation variable indices based on number of species (c variables come after species densities)
        int nv = (int)cd(0) + 4;
        int nc  = nv;

        double dc_left [3][3];
        double dc_right[3][3];
        double lap;

        // average cell size
        double dxmag = pow(dx[0]*dx[1]*dx[2],1.0/3.0);
        //double dxmag = sqrt(dx*dx + dy*dy + dz*dz);
          
        for (int n=0; n<5; ++n){
            // left face
            dc_left[0][0]  = (var(i  ,j  ,k  ,nc+n) - var(i-1,j  ,k  ,nc+n))/dx[0];
            dc_left[0][1]  = (var(i-1,j+1,k  ,nc+n) + var(i  ,j+1,k  ,nc+n)
                             +var(i-1,j-1,k  ,nc+n) + var(i  ,j-1,k  ,nc+n))/(4*dx[1]);
            dc_left[0][2]  = (var(i-1,j  ,k+1,nc+n) + var(i  ,j  ,k+1,nc+n)
                             +var(i-1,j  ,k-1,nc+n) + var(i  ,j  ,k-1,nc+n))/(4*dx[2]);

            // bottom face
            dc_left[1][0]  = (var(i-1,j  ,k  ,nc+n) + var(i+1,j  ,k  ,nc+n)
                             +var(i-1,j-1,k  ,nc+n) + var(i+1,j-1,k  ,nc+n))/(4*dx[0]);
            dc_left[1][1]  = (var(i  ,j  ,k  ,nc+n) - var(i  ,j-1,k  ,nc+n))/dx[1];
            dc_left[1][2]  = (var(i  ,j  ,k-1,nc+n) + var(i  ,j  ,k+1,nc+n)
                             +var(i  ,j-1,k-1,nc+n) + var(i  ,j-1,k+1,nc+n))/(4*dx[2]);

            // back (hind) face 
            dc_left[2][0]  = (var(i-1,j  ,k  ,nc+n) + var(i+1,j  ,k  ,nc+n)
                             +var(i-1,j  ,k-1,nc+n) + var(i+1,j  ,k-1,nc+n))/(4*dx[0]);
            dc_left[2][1]  = (var(i  ,j-1,k  ,nc+n) + var(i  ,j+1,k  ,nc+n)
                             +var(i  ,j-1,k-1,nc+n) + var(i  ,j+1,k-1,nc+n))/(4*dx[1]);
            dc_left[2][2]  = (var(i  ,j  ,k  ,nc+n) - var(i  ,j  ,k-1,nc+n))/dx[2];

            // right face
            dc_right[0][0] = (var(i+1,j  ,k  ,nc+n) - var(i  ,j  ,k  ,nc+n))/dx[0];
            dc_right[0][1] = (var(i  ,j+1,k  ,nc+n) + var(i+1,j+1,k  ,nc+n)
                             +var(i  ,j-1,k  ,nc+n) + var(i+1,j-1,k  ,nc+n))/(4*dx[1]);
            dc_right[0][2] = (var(i  ,j  ,k+1,nc+n) + var(i+1,j  ,k+1,nc+n)
                             +var(i  ,j  ,k-1,nc+n) + var(i+1,j  ,k-1,nc+n))/(4*dx[2]);

            //top face
            dc_right[1][0] = (var(i-1,j+1,k  ,nc+n) + var(i+1,j+1,k  ,nc+n)
                             +var(i-1,j  ,k  ,nc+n) + var(i+1,j  ,k  ,nc+n))/(4*dx[0]);
            dc_right[1][1] = (var(i  ,j+1,k  ,nc+n) - var(i  ,j  ,k  ,nc+n))/dx[1];
            dc_right[1][2] = (var(i  ,j+1,k-1,nc+n) + var(i  ,j+1,k+1,nc+n)
                             +var(i  ,j  ,k-1,nc+n) + var(i  ,j  ,k+1,nc+n))/(4*dx[2]);

            //front face
            dc_right[2][0] = (var(i-1,j  ,k+1,nc+n) + var(i+1,j  ,k+1,nc+n)
                             +var(i-1,j  ,k  ,nc+n) + var(i+1,j  ,k  ,nc+n))/(4*dx[0]);
            dc_right[2][1] = (var(i  ,j-1,k+1,nc+n) + var(i  ,j+1,k+1,nc+n)
                             +var(i  ,j-1,k  ,nc+n) + var(i  ,j+1,k  ,nc+n))/(4*dx[1]);
            dc_right[2][2] = (var(i  ,j  ,k+1,nc+n) - var(i  ,j  ,k  ,nc+n))/dx[2];


            // update ceq right hand side
            lap = 0;
            for (int d=0; d<3; ++d){
                for (int f=0; f<3; ++f){
                    lap = lap + (dc_right[d][f] - dc_left[d][f])/dx[d];
                }
            }
            dvar(i,j,k,nc+n) = maxS/(eps*dxmag)*(gradRho(i,j,k,n) - var(i,j,k,nc+n)) + kap*maxS*dxmag*lap;
        }
    }
};

struct calculateCeqFlux {

    FS4D var;
    FS3D rho;
    FS6D mFlux;  //(m,n,i,j,k,dir)
    FS4D cFlux; //(i,j,k,dir)
    Kokkos::View<double*> cd;

    calculateCeqFlux (FS4D var_, FS3D rho_, FS6D mFlux_, FS4D cFlux_, Kokkos::View<double*> cd_)
        : var(var_), rho(rho_), mFlux(mFlux_), cFlux(cFlux_), cd(cd_) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const int i, const int j, const int k) const {

        int nv = (int)cd(0) + 4;
        int nc  = nv;
        int nch = nv+1;
        int nc1  = nv+2;

        int ip=0;
        int jp=0;
        int kp=0;

        double c_left;
        double c_right;
        double ch_left;
        double ch_right;

        double cn_left  [3];
        double cn_right [3];
        double cmag_left  = 0;
        double cmag_right = 0;

        double m_left;
        double m_right;

        double rho_left;
        double rho_right;

        //for each direction (i=0, j=1, k=2)
        for (int dir=0; dir <3; ++dir){
            if (dir==0) ip=1;
            if (dir==1) jp=1;
            if (dir==2) kp=1;


            // left is current cell, right is cell in positive direction

            // get right and left components of isotropic C
            c_left  = var(i   ,j   ,k   ,nc);
            c_right = var(i+ip,j+jp,k+kp,nc);

            // get right and left components of anisotropic C
            ch_left  = var(i   ,j   ,k   ,nch);
            ch_right = var(i+ip,j+jp,k+kp,nch);

            rho_left  = rho(i   ,j   ,k   );
            rho_right = rho(i+ip,j+jp,k+kp);

            //get right and left values of each directional C
            for (int idx=0; idx<3; ++idx){
                cn_left[idx]  = var(i   ,j   ,k   ,nc1+idx);
                cn_right[idx] = var(i+ip,j+jp,k+kp,nc1+idx);
            }

            // calculate magnitude of directional C
            for (int idx=0; idx<3; ++idx) cmag_left += cn_left[idx]*cn_left[idx];
            //cmag_left  = sqrt(cmag_left);
            for (int idx=0; idx<3; ++idx) cmag_right += cn_right[idx]*cn_right[idx];
            //cmag_right = sqrt(cmag_right);

            if (cmag_left<=0.000001 || cmag_right<=0.000001){
                for (int m=0; m<3; ++m){
                    for (int n=0; n<3; ++n){
                        mFlux(m,n,i,j,k,dir) = 0.0;
                        cFlux(i,j,k,dir) = 0.0;
                    }
                }
            }else{
                //tensor components
                for (int m=0; m<3; ++m){
                    for (int n=0; n<3; ++n){

                        //dirac delta
                        double d = 0;
                        if (m==n) d=1;
                        
                        // calculate right and left tensor components
                        m_left  = d - (cn_left[m]*cn_left[n]/cmag_left);
                        m_right = d - (cn_right[m]*cn_right[n]/cmag_right);

                        //include isotropic c
                        m_left  = m_left*ch_left;
                        m_right = m_right*ch_right;

                        //include density
                        m_left  = m_left *rho_left;
                        m_right = m_right*rho_right;

                        //find flux
                        mFlux(m,n,i,j,k,dir) = (m_right - m_left)/2.0;
                    }
                }

                // calcualte isotropic C flux
                cFlux(i,j,k,dir) = (c_right*rho_right - c_left*rho_left)/2.0;

                ip=0;
                jp=0;
                kp=0;
                
            }
        }
    }
};

struct applyCeq {

    FS4D dvar;
    FS4D var;
    FS3D rho;
    FS6D mFlux;  //(m,n,i,j,k,dir)
    FS4D cFlux; //(i,j,k,dir)
    double alpha,beta,betae;
    Kokkos::View<double*> cd;

    applyCeq (FS4D dvar_, FS4D var_, FS3D rho_, FS6D mFlux_, FS4D cFlux_, Kokkos::View<double*> cd_,
                      double alpha_, double beta_, double betae_)
        : dvar(dvar_), var(var_), rho(rho_), mFlux(mFlux_), cFlux(cFlux_),
          cd(cd_), alpha(alpha_), beta(beta_), betae(betae_) {} 

    KOKKOS_INLINE_FUNCTION
    void operator()(const int i, const int j, const int k) const {

        int nv = (int)cd(0) + 4;
        double diffu;

        double du_right [3][3]; // face, direction
        double du_left [3][3];

        double dx [3];

        dx[0] = cd(1);
        dx[1] = cd(2);
        dx[2] = cd(3);

        double an; //anisotropic part
        double is; //isotropic part

        int ip,jp,kp;

        // for each velocity component and energy
        for (int n=0; n<4; ++n){
            // left face
            du_left[0][0] = (var(i  ,j  ,k  ,n)/rho(i  ,j  ,k  ) - var(i-1,j  ,k  ,n)/rho(i-1,j  ,k  ))/cd(1);
            du_left[0][1] = (var(i-1,j+1,k  ,n)/rho(i-1,j+1,k  ) + var(i  ,j+1,k  ,n)/rho(i  ,j+1,k  )
                            +var(i-1,j-1,k  ,n)/rho(i-1,j-1,k  ) + var(i  ,j-1,k  ,n)/rho(i  ,j-1,k  ))/(4*cd(2));
            du_left[0][2] = (var(i-1,j  ,k+1,n)/rho(i-1,j  ,k+1) + var(i  ,j  ,k+1,n)/rho(i  ,j  ,k+1)
                            +var(i-1,j  ,k-1,n)/rho(i-1,j  ,k-1) + var(i  ,j  ,k-1,n)/rho(i  ,j  ,k-1))/(4*cd(3));

            // bottom face
            du_left[1][0] = (var(i-1,j  ,k  ,n)/rho(i-1,j  ,k  ) + var(i+1,j  ,k  ,n)/rho(i+1,j  ,k  )
                            +var(i-1,j-1,k  ,n)/rho(i-1,j-1,k  ) + var(i+1,j-1,k  ,n)/rho(i+1,j-1,k  ))/(4*cd(1));
            du_left[1][1] = (var(i  ,j  ,k  ,n)/rho(i  ,j  ,k  ) - var(i  ,j-1,k  ,n)/rho(i  ,j-1,k  ))/cd(2);
            du_left[1][2] = (var(i  ,j  ,k-1,n)/rho(i  ,j  ,k-1) + var(i  ,j  ,k+1,n)/rho(i  ,j  ,k+1)
                            +var(i  ,j-1,k-1,n)/rho(i  ,j-1,k-1) + var(i  ,j-1,k+1,n)/rho(i  ,j-1,k+1))/(4*cd(3));

            // back (hind) face 
            du_left[2][0] = (var(i-1,j  ,k  ,n)/rho(i-1,j  ,k  ) + var(i+1,j  ,k  ,n)/rho(i+1,j  ,k  )
                            +var(i-1,j  ,k-1,n)/rho(i-1,j  ,k-1) + var(i+1,j  ,k-1,n)/rho(i+1,j  ,k-1))/(4*cd(1));
            du_left[2][1] = (var(i  ,j-1,k  ,n)/rho(i  ,j-1,k  ) + var(i  ,j+1,k  ,n)/rho(i  ,j+1,k  )
                            +var(i  ,j-1,k-1,n)/rho(i  ,j-1,k-1) + var(i  ,j+1,k-1,n)/rho(i  ,j+1,k-1))/(4*cd(2));
            du_left[2][2] = (var(i  ,j  ,k  ,n)/rho(i  ,j  ,k  ) - var(i  ,j  ,k-1,n)/rho(i  ,j  ,k-1))/cd(3);

            // right face
            du_right[0][0] = (var(i+1,j  ,k  ,n)/rho(i+1,j  ,k  ) - var(i  ,j  ,k  ,n)/rho(i  ,j  ,k  ))/cd(1);
            du_right[0][1] = (var(i  ,j+1,k  ,n)/rho(i  ,j+1,k  ) + var(i+1,j+1,k  ,n)/rho(i+1,j+1,k  )
                             +var(i  ,j-1,k  ,n)/rho(i  ,j-1,k  ) + var(i+1,j-1,k  ,n)/rho(i+1,j-1,k  ))/(4*cd(2));
            du_right[0][2] = (var(i  ,j  ,k+1,n)/rho(i  ,j  ,k+1) + var(i+1,j  ,k+1,n)/rho(i+1,j  ,k+1)
                             +var(i  ,j  ,k-1,n)/rho(i  ,j  ,k-1) + var(i+1,j  ,k-1,n)/rho(i+1,j  ,k-1))/(4*cd(3));

            //top face
            du_right[1][0]= (var(i-1,j+1,k  ,n)/rho(i-1,j+1,k  ) + var(i+1,j+1,k  ,n)/rho(i+1,j+1,k  )
                            +var(i-1,j  ,k  ,n)/rho(i-1,j  ,k  ) + var(i+1,j  ,k  ,n)/rho(i+1,j  ,k  ))/(4*cd(1));
            du_right[1][1]= (var(i  ,j+1,k  ,n)/rho(i  ,j+1,k  ) - var(i  ,j  ,k  ,n)/rho(i  ,j  ,k  ))/cd(2);
            du_right[1][2]= (var(i  ,j+1,k-1,n)/rho(i  ,j+1,k-1) + var(i  ,j+1,k+1,n)/rho(i  ,j+1,k+1)
                            +var(i  ,j  ,k-1,n)/rho(i  ,j  ,k-1) + var(i  ,j  ,k+1,n)/rho(i  ,j  ,k+1))/(4*cd(3));

            //front face
            du_right[2][0] = (var(i-1,j  ,k+1,n)/rho(i-1,j  ,k+1) + var(i+1,j  ,k+1,n)/rho(i+1,j  ,k+1)
                             +var(i-1,j  ,k  ,n)/rho(i-1,j  ,k  ) + var(i+1,j  ,k  ,n)/rho(i+1,j  ,k  ))/(4*cd(1));
            du_right[2][1] = (var(i  ,j-1,k+1,n)/rho(i  ,j-1,k+1) + var(i  ,j+1,k+1,n)/rho(i  ,j+1,k+1)
                             +var(i  ,j-1,k  ,n)/rho(i  ,j-1,k  ) + var(i  ,j+1,k  ,n)/rho(i  ,j+1,k  ))/(4*cd(2));
            du_right[2][2] = (var(i  ,j  ,k+1,n)/rho(i  ,j  ,k+1) - var(i  ,j  ,k  ,n)/rho(i  ,j  ,k  ))/cd(3);

            an = 0;
            is = 0;
            ip = 0;
            jp = 0;
            kp = 0;

            if (n < 3){
                diffu = 0;
                for (int d=0; d<3; ++d){ // each direction for divergence
                    if (d==0) ip = 1;
                    if (d==1) jp = 1;
                    if (d==2) kp = 1;

                    for (int f=0; f<3; ++f){ //each direction for gradient
                        an = an + ( mFlux(d,f,i,j,k,d)*du_right[d][f] - mFlux(d,f,i-ip,j-jp,k-kp,d)*du_left[d][f] )/dx[d];
                        is = is + ( cFlux(i,j,k,d)*du_right[d][f] - cFlux(i-ip,j-jp,k-kp,d)*du_left[d][f] )/dx[d];
                    }

                    ip = 0;
                    jp = 0;
                    kp = 0;
                    }

                diffu = alpha*an + beta*is;

                dvar(i,j,k,n) = dvar(i,j,k,n) + diffu;
            }else{
                diffu = 0;
                for (int d=0; d<3; ++d){ // each direction for divergence
                    if (d==0) ip = 1;
                    if (d==1) jp = 1;
                    if (d==2) kp = 1;

                    for (int f=0; f<3; ++f){ // each direction for gradient 
                        diffu = diffu + ( cFlux(i,j,k,d)*du_right[d][f] - cFlux(i-ip,j-jp,k-kp,d)*du_left[d][f] )/dx[d];
                    }

                    ip = 0;
                    jp = 0;
                    kp = 0;
                }
                dvar(i,j,k,n) = dvar(i,j,k,n) + betae*diffu;
            }
        }
    }
};

hydroc3d_func::hydroc3d_func(struct inputConfig &cf_, Kokkos::View<double*> & cd_):rk_func(cf_,cd_) {};

void hydroc3d_func::compute(const FS4D & mvar, FS4D & mdvar) {

    // Typename acronyms for 3D and 4D variables

    // copy input views
    FS4D dvar = mdvar;
    FS4D var = mvar;

    // create configuration data view
    Kokkos::View<double*> cd = mcd;

    // create temprary views
    FS3D p("p",cf.ngi,cf.ngj,cf.ngk);
    FS3D rho("rho",cf.ngi,cf.ngj,cf.ngk);
    FS3D wenox("wenox",cf.ngi,cf.ngj,cf.ngk);
    FS3D wenoy("wenoy",cf.ngi,cf.ngj,cf.ngk);
    FS3D wenoz("wenoz",cf.ngi,cf.ngj,cf.ngk);
    FS4D gradRho("gradRho",cf.ngi,cf.ngj,cf.ngk,5);
    FS6D mFlux("mFlux",3,3,cf.ngi,cf.ngj,cf.ngk,3);
    FS4D cFlux("cFlux",cf.ngi,cf.ngj,cf.ngk,3);

    // create range policies
    policy_f ghost_pol = policy_f({0,0,0},{cf.ngi, cf.ngj, cf.ngk});
    policy_f cell_pol  = policy_f({cf.ng,cf.ng,cf.ng},{cf.ngi-cf.ng, cf.ngj-cf.ng, cf.ngk-cf.ng});
    policy_f weno_pol  = policy_f({cf.ng-1,cf.ng-1,cf.ng-1},{cf.ngi-cf.ng, cf.ngj-cf.ng, cf.ngk-cf.ng});

    double myMaxS,maxS;

    double myMaxC, maxC;
    double myMaxCh, maxCh;
    double myMaxC1, maxC1;
    double myMaxC2, maxC2;
    double myMaxC3, maxC3;
    double mu,alpha,beta,betae;

    /**** WENO ****/
    Kokkos::parallel_for( ghost_pol, calculateRhoAndPressure(var,p,rho,cd) );

    for (int v=0; v<cf.nv; ++v){
        Kokkos::parallel_for( weno_pol, calculateWenoFluxes(var,p,rho,wenox,wenoy,wenoz,cd,v) );

        Kokkos::parallel_for( cell_pol, applyWenoFluxes(dvar,wenox,wenoy,wenoz,v) );
    }

    Kokkos::parallel_for( cell_pol, applyPressure(dvar,p,cd) );

    if (cf.ceq != 0){
        /**** CEQ ****/
        //find mac wavespeed
        Kokkos::parallel_reduce(cell_pol,maxWaveSpeed(var,p,rho,cd), Kokkos::Max<double>(myMaxS));
        MPI_Allreduce(&myMaxS,&maxS,1,MPI_DOUBLE,MPI_MAX,cf.comm);

        //calculate density and energy gradients and find ceq source terms
        Kokkos::parallel_for(cell_pol,calculateRhoGrad(var,rho,gradRho,cd));

        // update c equations
        Kokkos::parallel_for(cell_pol,
            updateCeq(dvar,var,gradRho,maxS,cd,cf.kap,cf.eps));

        /**** Apply CEQ ****/
        Kokkos::parallel_reduce(cell_pol,maxGradFunctor(var,cf.nv+0), Kokkos::Max<double>(myMaxC));
        MPI_Allreduce(&myMaxC,&maxC,1,MPI_DOUBLE,MPI_MAX,cf.comm);

        Kokkos::parallel_reduce(cell_pol,maxGradFunctor(var,cf.nv+1), Kokkos::Max<double>(myMaxCh));
        MPI_Allreduce(&myMaxCh,&maxCh,1,MPI_DOUBLE,MPI_MAX,cf.comm);

        Kokkos::parallel_reduce(cell_pol,maxGradFunctor(var,cf.nv+2), Kokkos::Max<double>(myMaxC1));
        MPI_Allreduce(&myMaxC1,&maxC1,1,MPI_DOUBLE,MPI_MAX,cf.comm);

        Kokkos::parallel_reduce(cell_pol,maxGradFunctor(var,cf.nv+3), Kokkos::Max<double>(myMaxC2));
        MPI_Allreduce(&myMaxC2,&maxC2,1,MPI_DOUBLE,MPI_MAX,cf.comm);

        Kokkos::parallel_reduce(cell_pol,maxGradFunctor(var,cf.nv+4), Kokkos::Max<double>(myMaxC3));
        MPI_Allreduce(&myMaxC3,&maxC3,1,MPI_DOUBLE,MPI_MAX,cf.comm);

        mu = maxC1;
        if (maxC2 > mu)
            mu = maxC2;
        if (maxC3 > mu)
            mu = maxC3;

        alpha = 0.0;
        beta = 0.0;
        betae = 0.0;

        double dxmag = sqrt(cf.dx*cf.dx+cf.dy*cf.dy+cf.dz*cf.dz);
        if (mu>0 && maxCh>0)
            alpha = (dxmag/(mu*mu*maxCh))*cf.alpha;
        if (maxC>0){
            beta = (dxmag/maxC)*cf.beta;
            betae = (dxmag/maxC)*cf.betae;
        }

        //if (cf.rank == 0)
        //    printf("%.4f,%.4f,%.4f,%.4f,%.4f\n",maxC,maxCh,maxC1,maxC2,maxC3);
        //    //printf("alpha = %f, beta = %f, betae = %f\n",alpha,beta,betae);
        
        Kokkos::parallel_for(weno_pol,
            calculateCeqFlux(var,rho,mFlux,cFlux,cd));

        Kokkos::parallel_for( cell_pol, applyCeq(dvar,var,rho,mFlux,cFlux,cd,alpha,beta,betae) );
    }
}
