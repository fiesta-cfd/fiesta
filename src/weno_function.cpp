#include "input.hpp"
#include "mpi_init.hpp"
#include "cgns.hpp"
#include "bc.hpp"
#include "Kokkos_Core.hpp"
#include <mpi.h>
#include "lsdebug.hpp"
#include "weno_function.hpp"

struct calculateRhoAndPressure {
    typedef typename Kokkos::View<double****> V4D;
    typedef typename Kokkos::View<double***> V3D;
    V4D var;
    V3D p;
    V3D rho;
    Kokkos::View<double*> cd;

    calculateRhoAndPressure (V4D var_, V3D p_, V3D rho_, Kokkos::View<double*> cd_)
         : var(var_), p(p_), rho(rho_), cd(cd_) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const int i, const int j, const int k) const {

        int ns = (int)cd(0);
        double gamma, gammas, Rs;
        double Cp = 0;
        double Cv = 0;

        rho(i,j,k) = 0.0;

        for (int s=0; s<ns; ++s){
            rho(i,j,k) = rho(i,j,k) + var(i,j,k,4+s);
        }

        for (int s=0; s<ns; ++s){
            gammas = cd(4+2*s);
            Rs = cd(4+2*s+1);

            Cp = Cp + (var(i,j,k,4+s)/rho(i,j,k))*( gammas*Rs/(gammas-1) );
            Cv = Cv + (var(i,j,k,4+s)/rho(i,j,k))*( Rs/(gammas-1) );
        }

        gamma = Cp/Cv;

        p(i,j,k) = (gamma-1)*( var(i,j,k,3) - (0.5/rho(i,j,k))
                  *(var(i,j,k,0)*var(i,j,k,0) + var(i,j,k,1)*var(i,j,k,1) + var(i,j,k,2)*var(i,j,k,2)) );
    }
};

struct calculateWenoFluxes {
    
    typedef typename Kokkos::View<double****> V4D;
    typedef typename Kokkos::View<double***> V3D;
    V4D var;
    V3D p;
    V3D rho;
    V3D wenox;
    V3D wenoy;
    V3D wenoz;
    Kokkos::View<double*> cd;
    int v;
    double eps = 0.000001;

    calculateWenoFluxes (V4D var_, V3D p_, V3D rho_,
                         V3D wenox_, V3D wenoy_, V3D wenoz_, Kokkos::View<double*> cd_, int v_)
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
    
    typedef typename Kokkos::View<double****> V4D;
    typedef typename Kokkos::View<double***> V3D;
    V4D dvar;
    V3D wenox;
    V3D wenoy;
    V3D wenoz;
    int v;

    applyWenoFluxes (V4D dvar_, V3D wenox_, V3D wenoy_, V3D wenoz_, int v_)
        : dvar(dvar_), wenox(wenox_), wenoy(wenoy_), wenoz(wenoz_), v(v_) {}
    
    KOKKOS_INLINE_FUNCTION
    void operator()(const int i, const int j, const int k) const {

        dvar(i,j,k,v) = -( (wenox(i,j,k) - wenox(i-1,j,k))
                          +(wenoy(i,j,k) - wenoy(i,j-1,k))
                          +(wenoz(i,j,k) - wenoz(i,j,k-1)) );
    }
};

struct applyPressure {
    
    typedef typename Kokkos::View<double****> V4D;
    typedef typename Kokkos::View<double***> V3D;
    V4D dvar;
    V3D p;
    Kokkos::View<double*> cd;

    applyPressure (V4D dvar_, V3D p_, Kokkos::View<double*> cd_)
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
    
    typedef typename Kokkos::View<double****> V4D;
    typedef typename Kokkos::View<double***> V3D;
    V4D var;
    V3D p;
    V3D rho;
    Kokkos::View<double*> cd;

    maxWaveSpeed (V4D var_, V3D p_, V3D rho_, Kokkos::View<double*> cd_)
        : var(var_), p(p_), rho(rho_), cd(cd_) {}
    
    KOKKOS_INLINE_FUNCTION
    void operator()(const int i, const int j, const int k, double& lmax) const {
        int ns = (int)cd(0);
        double gamma, gammas, Rs;
        double Cp = 0;
        double Cv = 0;

        double a, s;

        for (int s=0; s<ns; ++s){
            gammas = cd(4+2*s);
            Rs = cd(4+2*s+1);

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

    typedef typename Kokkos::View<double****> V4D;
    typedef typename Kokkos::View<double***> V3D;
    V4D var;
    V3D rho;
    V4D gradRho;
    Kokkos::View<double*> cd;

    calculateRhoGrad (V4D var_, V3D rho_, V4D gradRho_, Kokkos::View<double*> cd_)
        : var(var_), rho(rho_), gradRho(gradRho_), cd(cd_) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const int i, const int j, const int k) const {

        double dxr = (rho(i-2,j,k) - 8.0*rho(i-1,j,k) + 8.0*rho(i+1,j,k) - rho(i+2,j,k))/(12.0*cd(1));
        double dyr = (rho(i,j-2,k) - 8.0*rho(i,j-1,k) + 8.0*rho(i,j+1,k) - rho(i,j+2,k))/(12.0*cd(2));
        double dzr = (rho(i,j,k-2) - 8.0*rho(i,j,k-1) + 8.0*rho(i,j,k+1) - rho(i,j,k+2))/(12.0*cd(3));

        gradRho(i,j,k,0) = dxr*dxr+dyr*dyr+dzr*dzr;
        gradRho(i,j,k,1) = dyr-dzr;
        gradRho(i,j,k,2) = dzr-dxr;
        gradRho(i,j,k,3) = dxr-dyr;

    }
};

struct maxGradFunctor {

    typedef typename Kokkos::View<double****> V4D;
    V4D gradRho;
    int n;

    maxGradFunctor(V4D gradRho_, int n_) : gradRho(gradRho_), n(n_) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const int i, const int j, const int k, double& lmax) const {
        if (abs(gradRho(i,j,k,n)) > lmax)
            lmax = abs(gradRho(i,j,k,n));
    }
};
struct calculateCeqFlux {

    typedef typename Kokkos::View<double****> V4D;
    typedef typename Kokkos::View<double***> V3D;
    V4D dvar;
    V4D var;
    V4D gradRho;
    double maxS,maxGradRho,maxTau1,maxTau2,maxTau3,kap,eps;
    Kokkos::View<double*> cd;

    calculateCeqFlux (V4D dvar_, V4D var_, V4D gradRho_, double maxS_, double maxGradRho_, double maxTau1_,
                      double maxTau2_, double maxTau3_, Kokkos::View<double*> cd_, double kap_, double eps_)
        : dvar(dvar_), var(var_), gradRho(gradRho_), maxS(maxS_), maxGradRho(maxGradRho_), maxTau1(maxTau1_),
          maxTau2(maxTau2_),  maxTau3(maxTau3_), cd(cd_), kap(kap_), eps(eps_) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const int i, const int j, const int k) const {
        int nv = (int)cd(0) + 4;
        int ng = nv;
        int nt1 = nv+1;
        int nt2 = nv+2;
        int nt3 = nv+3;

        double G  = 0.0;
        double T1 = 0.0;
        double T2 = 0.0;
        double T3 = 0.0;

        double dcx,dc1x,dc2x,dc3x;
        double dcy,dc1y,dc2y,dc3y;
        double dcz,dc1z,dc2z,dc3z;

        double dxmag = sqrt(cd(1)*cd(1) + cd(2)*cd(2) + cd(3)*cd(3));

        dcx  = (var(i+1,j,k,ng ) - 2*var(i,j,k,ng ) + var(i-1,j,k,ng ))/cd(1);
        dc1x = (var(i+1,j,k,nt1) - 2*var(i,j,k,nt1) + var(i-1,j,k,nt1))/cd(1);
        dc2x = (var(i+1,j,k,nt2) - 2*var(i,j,k,nt2) + var(i-1,j,k,nt2))/cd(1);
        dc3x = (var(i+1,j,k,nt3) - 2*var(i,j,k,nt3) + var(i-1,j,k,nt3))/cd(1);

        dcy  = (var(i,j+1,k,ng ) - 2*var(i,j,k,ng ) + var(i,j-1,k,ng ))/cd(2);
        dc1y = (var(i,j+1,k,nt1) - 2*var(i,j,k,nt1) + var(i,j-1,k,nt1))/cd(2);
        dc2y = (var(i,j+1,k,nt2) - 2*var(i,j,k,nt2) + var(i,j-1,k,nt2))/cd(2);
        dc3y = (var(i,j+1,k,nt3) - 2*var(i,j,k,nt3) + var(i,j-1,k,nt3))/cd(2);

        dcz  = (var(i,j,k+1,ng ) - 2*var(i,j,k,ng ) + var(i,j,k-1,ng ))/cd(3);
        dc1z = (var(i,j,k+1,nt1) - 2*var(i,j,k,nt1) + var(i,j,k-1,nt1))/cd(3);
        dc2z = (var(i,j,k+1,nt2) - 2*var(i,j,k,nt2) + var(i,j,k-1,nt2))/cd(3);
        dc3z = (var(i,j,k+1,nt3) - 2*var(i,j,k,nt3) + var(i,j,k-1,nt3))/cd(3);

        if (maxGradRho > 0.0)
            G = gradRho(i,j,k,0)/maxGradRho;
        if (maxTau1 > 0.0)
            T1 = gradRho(i,j,k,1)/maxTau1;
        if (maxTau2 > 0.0)
            T2 = gradRho(i,j,k,2)/maxTau2;
        if (maxTau3 > 0.0)
            T3 = gradRho(i,j,k,3)/maxTau3;

        dvar(i,j,k,ng ) = maxS/(eps*dxmag)*(G  - var(i,j,k,ng )) + kap*maxS*dxmag*(dcx +dcy +dcz);
        dvar(i,j,k,nt1) = maxS/(eps*dxmag)*(T1 - var(i,j,k,nt1)) + kap*maxS*dxmag*(dc1x+dc1y+dc1z);
        dvar(i,j,k,nt2) = maxS/(eps*dxmag)*(T2 - var(i,j,k,nt2)) + kap*maxS*dxmag*(dc2x+dc2y+dc2z);
        dvar(i,j,k,nt3) = maxS/(eps*dxmag)*(T3 - var(i,j,k,nt3)) + kap*maxS*dxmag*(dc3x+dc3y+dc3z);
    }
};

struct applyCeq {

    typedef typename Kokkos::View<double****> V4D;
    typedef typename Kokkos::View<double***> V3D;
    V4D dvar;
    V4D var;
    V3D rho;
    double maxC,mu,beta;
    Kokkos::View<double*> cd;

    applyCeq (V4D dvar_, V4D var_, V3D rho_, Kokkos::View<double*> cd_,
                      double maxC_, double mu_, double beta_)
        : dvar(dvar_), var(var_), rho(rho_), cd(cd_), maxC(maxC_), mu(mu_), beta(beta_) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const int i, const int j, const int k) const {

        int nv = (int)cd(0) + 4;
        double diffu;
        double dxmag = sqrt(cd(1)*cd(1) + cd(2)*cd(2) + cd(3)*cd(3));

        double mbeta = 0;
        if (mu>1 && maxC>1)
             mbeta = beta*(dxmag*dxmag/(mu*mu*maxC));

        double rl,rb,rh,rr,rt,rf;
        double cl,cb,ch,cr,ct,cf;
        double c1l,c1b,c1h,c1r,c1t,c1f;
        double c2l,c2b,c2h,c2r,c2t,c2f;
        double c3l,c3b,c3h,c3r,c3t,c3f;
        double duxl,duxb,duxh,duxr,duxt,duxf;
        double duyl,duyb,duyh,duyr,duyt,duyf;
        double duzl,duzb,duzh,duzr,duzt,duzf;

        rl = (rho(i  ,j  ,k  ) - rho(i-1,j  ,k  ))/2.0;
        rb = (rho(i  ,j  ,k  ) - rho(i  ,j-1,k  ))/2.0;
        rh = (rho(i  ,j  ,k  ) - rho(i  ,j  ,k-1))/2.0;

        rr = (rho(i+1,j  ,k  ) - rho(i  ,j  ,k  ))/2.0;
        rt = (rho(i  ,j+1,k  ) - rho(i  ,j  ,k  ))/2.0;
        rf = (rho(i  ,j  ,k+1) - rho(i  ,j  ,k  ))/2.0;

        cl = (var(i  ,j  ,k  ,nv  ) - var(i-1,j  ,k  ,nv  ))/2.0;
        cb = (var(i  ,j  ,k  ,nv  ) - var(i  ,j-1,k  ,nv  ))/2.0;
        ch = (var(i  ,j  ,k  ,nv  ) - var(i  ,j  ,k-1,nv  ))/2.0;

        cr = (var(i+1,j  ,k  ,nv  ) - var(i  ,j  ,k  ,nv  ))/2.0;
        ct = (var(i  ,j+1,k  ,nv  ) - var(i  ,j  ,k  ,nv  ))/2.0;
        cf = (var(i  ,j  ,k+1,nv  ) - var(i  ,j  ,k  ,nv  ))/2.0;

        c1l = (var(i  ,j  ,k  ,nv+1) - var(i-1,j  ,k  ,nv+1))/2.0;
        c1b = (var(i  ,j  ,k  ,nv+1) - var(i  ,j-1,k  ,nv+1))/2.0;
        c1h = (var(i  ,j  ,k  ,nv+1) - var(i  ,j  ,k-1,nv+1))/2.0;

        c1r = (var(i+1,j  ,k  ,nv+1) - var(i  ,j  ,k  ,nv+1))/2.0;
        c1t = (var(i  ,j+1,k  ,nv+1) - var(i  ,j  ,k  ,nv+1))/2.0;
        c1f = (var(i  ,j  ,k+1,nv+1) - var(i  ,j  ,k  ,nv+1))/2.0;

        c2l = (var(i  ,j  ,k  ,nv+2) - var(i-1,j  ,k  ,nv+2))/2.0;
        c2b = (var(i  ,j  ,k  ,nv+2) - var(i  ,j-1,k  ,nv+2))/2.0;
        c2h = (var(i  ,j  ,k  ,nv+2) - var(i  ,j  ,k-1,nv+2))/2.0;

        c2r = (var(i+1,j  ,k  ,nv+2) - var(i  ,j  ,k  ,nv+2))/2.0;
        c2t = (var(i  ,j+1,k  ,nv+2) - var(i  ,j  ,k  ,nv+2))/2.0;
        c2f = (var(i  ,j  ,k+1,nv+2) - var(i  ,j  ,k  ,nv+2))/2.0;

        c3l = (var(i  ,j  ,k  ,nv+3) - var(i-1,j  ,k  ,nv+3))/2.0;
        c3b = (var(i  ,j  ,k  ,nv+3) - var(i  ,j-1,k  ,nv+3))/2.0;
        c3h = (var(i  ,j  ,k  ,nv+3) - var(i  ,j  ,k-1,nv+3))/2.0;

        c3r = (var(i+1,j  ,k  ,nv+3) - var(i  ,j  ,k  ,nv+3))/2.0;
        c3t = (var(i  ,j+1,k  ,nv+3) - var(i  ,j  ,k  ,nv+3))/2.0;
        c3f = (var(i  ,j  ,k+1,nv+3) - var(i  ,j  ,k  ,nv+3))/2.0;


        for (int n=0; n<3; ++n){
            duxl = (var(i  ,j  ,k  ,n)/rho(i  ,j  ,k  ) - var(i-1,j  ,k  ,n)/rho(i-1,j  ,k  ))/cd(1);
            duyl = (var(i-1,j+1,k  ,n)/rho(i-1,j+1,k  ) + var(i  ,j+1,k  ,n)/rho(i  ,j+1,k  )
                   +var(i-1,j-1,k  ,n)/rho(i-1,j-1,k  ) + var(i  ,j-1,k  ,n)/rho(i  ,j-1,k  ))/(4*cd(2));
            duzl = (var(i-1,j  ,k+1,n)/rho(i-1,j  ,k+1) + var(i  ,j  ,k+1,n)/rho(i  ,j  ,k+1)
                   +var(i-1,j  ,k-1,n)/rho(i-1,j  ,k-1) + var(i  ,j  ,k-1,n)/rho(i  ,j  ,k-1))/(4*cd(3));

            duxb = (var(i-1,j  ,k  ,n)/rho(i-1,j  ,k  ) + var(i+1,j  ,k  ,n)/rho(i+1,j  ,k  )
                   +var(i-1,j-1,k  ,n)/rho(i-1,j-1,k  ) + var(i+1,j-1,k  ,n)/rho(i+1,j-1,k  ))/(4*cd(1));
            duyb = (var(i  ,j  ,k  ,n)/rho(i  ,j  ,k  ) - var(i  ,j-1,k  ,n)/rho(i  ,j-1,k  ))/cd(2);
            duzb = (var(i  ,j  ,k-1,n)/rho(i  ,j  ,k-1) + var(i  ,j  ,k+1,n)/rho(i  ,j  ,k+1)
                   +var(i  ,j-1,k-1,n)/rho(i  ,j-1,k-1) + var(i  ,j-1,k+1,n)/rho(i  ,j-1,k+1))/(4*cd(3));

            duxh = (var(i-1,j  ,k  ,n)/rho(i-1,j  ,k  ) + var(i+1,j  ,k  ,n)/rho(i+1,j  ,k  )
                   +var(i-1,j  ,k-1,n)/rho(i-1,j  ,k-1) + var(i+1,j  ,k-1,n)/rho(i+1,j  ,k-1))/(4*cd(1));
            duyh = (var(i  ,j-1,k  ,n)/rho(i  ,j-1,k  ) + var(i  ,j+1,k  ,n)/rho(i  ,j+1,k  )
                   +var(i  ,j-1,k-1,n)/rho(i  ,j-1,k-1) + var(i  ,j+1,k-1,n)/rho(i  ,j+1,k-1))/(4*cd(2));
            duzh = (var(i  ,j  ,k  ,n)/rho(i  ,j  ,k  ) - var(i  ,j  ,k-1,n)/rho(i  ,j  ,k-1))/cd(3);

            duxr = (var(i+1,j  ,k  ,n)/rho(i+1,j  ,k  ) - var(i  ,j  ,k  ,n)/rho(i  ,j  ,k  ))/cd(1);
            duyr = (var(i  ,j+1,k  ,n)/rho(i  ,j+1,k  ) + var(i+1,j+1,k  ,n)/rho(i+1,j+1,k  )
                   +var(i  ,j-1,k  ,n)/rho(i  ,j-1,k  ) + var(i+1,j-1,k  ,n)/rho(i+1,j-1,k  ))/(4*cd(2));
            duzr = (var(i  ,j  ,k+1,n)/rho(i  ,j  ,k+1) + var(i+1,j  ,k+1,n)/rho(i+1,j  ,k+1)
                   +var(i  ,j  ,k-1,n)/rho(i  ,j  ,k-1) + var(i+1,j  ,k-1,n)/rho(i+1,j  ,k-1))/(4*cd(3));

            duxt = (var(i-1,j+1,k  ,n)/rho(i-1,j+1,k  ) + var(i+1,j+1,k  ,n)/rho(i+1,j+1,k  )
                   +var(i-1,j  ,k  ,n)/rho(i-1,j  ,k  ) + var(i+1,j  ,k  ,n)/rho(i+1,j  ,k  ))/(4*cd(1));
            duyt = (var(i  ,j+1,k  ,n)/rho(i  ,j+1,k  ) - var(i  ,j  ,k  ,n)/rho(i  ,j  ,k  ))/cd(2);
            duzt = (var(i  ,j+1,k-1,n)/rho(i  ,j+1,k-1) + var(i  ,j+1,k+1,n)/rho(i  ,j+1,k+1)
                   +var(i  ,j  ,k-1,n)/rho(i  ,j  ,k-1) + var(i  ,j  ,k+1,n)/rho(i  ,j  ,k+1))/(4*cd(3));

            duxf = (var(i-1,j  ,k+1,n)/rho(i-1,j  ,k+1) + var(i+1,j  ,k  ,n)/rho(i+1,j  ,k+1)
                   +var(i-1,j  ,k  ,n)/rho(i-1,j  ,k  ) + var(i+1,j  ,k-1,n)/rho(i+1,j  ,k  ))/(4*cd(1));
            duyf = (var(i  ,j-1,k+1,n)/rho(i  ,j-1,k+1) + var(i  ,j+1,k  ,n)/rho(i  ,j+1,k+1)
                   +var(i  ,j-1,k  ,n)/rho(i  ,j-1,k  ) + var(i  ,j+1,k-1,n)/rho(i  ,j+1,k  ))/(4*cd(2));
            duzf = (var(i  ,j  ,k+1,n)/rho(i  ,j  ,k+1) - var(i  ,j  ,k-1,n)/rho(i  ,j  ,k  ))/cd(3);

            diffu = mbeta*(-(rl*cl*c1l*(c1l*duxl+c2l*duyl+c3l*duzl) - rr*cr*c1r*(c1r*duxr+c2r*duyr+c3r*duzr))/cd(1)
                           -(rb*cb*c1b*(c1b*duxb+c2b*duyb+c3b*duzb) - rt*ct*c1t*(c1t*duxt+c2t*duyt+c3t*duzt))/cd(2)
                           -(rh*ch*c1h*(c1h*duxh+c2h*duyh+c3h*duzh) - rf*cf*c1f*(c1f*duxf+c2f*duyf+c3f*duzf))/cd(3));

            dvar(i,j,k,n) = dvar(i,j,k,n) + diffu;
            //dvar(i,j,k,n) = dvar(i,j,k,n) + 0;
        }
    }
};

weno_func::weno_func(struct inputConfig &cf_, const Kokkos::View<double****> & u_,
                     Kokkos::View<double****> & k_, Kokkos::View<double*> & cd_)
                     : cf(cf_) , mvar(u_), mdvar(k_) , mcd(cd_) {};

void weno_func::operator()() {

    // Typename acronyms for 3D and 4D variables
    typedef typename Kokkos::View<double****> V4D;
    typedef typename Kokkos::View<double***> V3D;

    // copy input views
    V4D dvar = mdvar;
    V4D var = mvar;

    // create configuration data view
    Kokkos::View<double*> cd = mcd;

    // create temprary views
    V3D p("p",cf.ngi,cf.ngj,cf.ngk);
    V3D rho("rho",cf.ngi,cf.ngj,cf.ngk);
    V3D wenox("wenox",cf.ngi,cf.ngj,cf.ngk);
    V3D wenoy("wenoy",cf.ngi,cf.ngj,cf.ngk);
    V3D wenoz("wenoz",cf.ngi,cf.ngj,cf.ngk);
    V4D gradRho("gradRho",cf.ngi,cf.ngj,cf.ngk,4);

    // create range policies
    policy_f ghost_pol = policy_f({0,0,0},{cf.ngi, cf.ngj, cf.ngk});
    policy_f cell_pol  = policy_f({cf.ng,cf.ng,cf.ng},{cf.ngi-cf.ng, cf.ngj-cf.ng, cf.ngk-cf.ng});
    policy_f weno_pol  = policy_f({cf.ng-1,cf.ng-1,cf.ng-1},{cf.ngi-cf.ng, cf.ngj-cf.ng, cf.ngk-cf.ng});

    double myMaxS,maxS;
    double myMaxRhoGrad,maxRhoGrad;
    double myMaxTau1,maxTau1;
    double myMaxTau2,maxTau2;
    double myMaxTau3,maxTau3;

    double myMaxC, maxC;
    double myMaxC1, maxC1;
    double myMaxC2, maxC2;
    double myMaxC3, maxC3;
    double mu;

    // WENO
    Kokkos::parallel_for( ghost_pol, calculateRhoAndPressure(var,p,rho,cd) );

    for (int v=0; v<cf.nv; ++v){
        Kokkos::parallel_for( weno_pol, calculateWenoFluxes(var,p,rho,wenox,wenoy,wenoz,cd,v) );

        Kokkos::parallel_for( cell_pol, applyWenoFluxes(dvar,wenox,wenoy,wenoz,v) );
    }

    Kokkos::parallel_for( cell_pol, applyPressure(dvar,p,cd) );

    if (cf.ceq != 0){
        // CEQ
        Kokkos::parallel_reduce(cell_pol,maxWaveSpeed(var,p,rho,cd), Kokkos::Max<double>(myMaxS));
        MPI_Allreduce(&myMaxS,&maxS,1,MPI_DOUBLE,MPI_MAX,cf.comm);

        Kokkos::parallel_for(cell_pol,calculateRhoGrad(var,rho,gradRho,cd));

        Kokkos::parallel_reduce(ghost_pol,maxGradFunctor(gradRho,0), Kokkos::Max<double>(myMaxRhoGrad));
        MPI_Allreduce(&myMaxRhoGrad,&maxRhoGrad,1,MPI_DOUBLE,MPI_MAX,cf.comm);

        Kokkos::parallel_reduce(ghost_pol,maxGradFunctor(gradRho,1), Kokkos::Max<double>(myMaxTau1));
        MPI_Allreduce(&myMaxTau1,&maxTau1,1,MPI_DOUBLE,MPI_MAX,cf.comm);

        Kokkos::parallel_reduce(ghost_pol,maxGradFunctor(gradRho,2), Kokkos::Max<double>(myMaxTau2));
        MPI_Allreduce(&myMaxTau2,&maxTau2,1,MPI_DOUBLE,MPI_MAX,cf.comm);

        Kokkos::parallel_reduce(ghost_pol,maxGradFunctor(gradRho,3), Kokkos::Max<double>(myMaxTau3));
        MPI_Allreduce(&myMaxTau3,&maxTau3,1,MPI_DOUBLE,MPI_MAX,cf.comm);

        Kokkos::parallel_for(cell_pol,
            calculateCeqFlux(dvar,var,gradRho,maxS,maxRhoGrad,maxTau1,maxTau2,maxTau3,cd,cf.kap,cf.eps));

        // Apply CEQ
        Kokkos::parallel_reduce(ghost_pol,maxGradFunctor(var,cf.nv+0), Kokkos::Max<double>(myMaxC));
        MPI_Allreduce(&myMaxC,&maxC,1,MPI_DOUBLE,MPI_MAX,cf.comm);

        Kokkos::parallel_reduce(ghost_pol,maxGradFunctor(var,cf.nv+1), Kokkos::Max<double>(myMaxC1));
        MPI_Allreduce(&myMaxC1,&maxC1,1,MPI_DOUBLE,MPI_MAX,cf.comm);

        Kokkos::parallel_reduce(ghost_pol,maxGradFunctor(var,cf.nv+2), Kokkos::Max<double>(myMaxC2));
        MPI_Allreduce(&myMaxC2,&maxC2,1,MPI_DOUBLE,MPI_MAX,cf.comm);

        Kokkos::parallel_reduce(ghost_pol,maxGradFunctor(var,cf.nv+3), Kokkos::Max<double>(myMaxC3));
        MPI_Allreduce(&myMaxC3,&maxC3,1,MPI_DOUBLE,MPI_MAX,cf.comm);

        mu = maxC1;
        if (maxC2 > mu)
            mu = maxC2;
        if (maxC3 > mu)
            mu = maxC3;

        Kokkos::parallel_for( cell_pol, applyCeq(dvar,var,rho,cd,maxC,mu,cf.beta) );
    }
}