#include "input.hpp"
#ifndef NOMPI
#include "mpi.hpp"
#include "cgns.hpp"
#include <mpi.h>
#endif
#include "bc.hpp"
#include "Kokkos_Core.hpp"
#include "debug.hpp"
#include "cart2d.hpp"
#include "flux.hpp"
#include "ceq.hpp"
#include <iomanip>
#include <fstream>

struct calculateGravity {
    FS4D dvar;
    FS4D var;
    FS2D rho;
    double g, gx, gy;

    calculateGravity (FS4D dvar_, FS4D var_, FS2D rho_, double g_, double gx_, double gy_)
         : dvar(dvar_), var(var_), rho(rho_), g(g_), gx(gx_), gy(gy_) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const int i, const int j) const {

        double u = var(i,j,0,0)/rho(i,j);
        double v = var(i,j,0,1)/rho(i,j);
        double rhop = rho(i,j) - 1.0;
        double eps = 1e-6;

        if (rhop >= eps || rhop <= -eps){
            double f = -g*rhop;;

            //dvar(i,j,0,0) += f*gx;
            //dvar(i,j,0,1) += f*gy;
            //dvar(i,j,0,2) += u*f*gx + v*f*gy;

            dvar(i,j,0,1) += f;
            dvar(i,j,0,2) += v*f;

            //printf("%f, %f\n",f,v);
        }
    }
};

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
            gammas = cd(6+3*s);
            Rs = cd(6+3*s+1);

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
        //double mu1 = 2.928e-5;
        //double mu2 = 1.610e-5;
        double dudx,dvdy,dudy,dvdx;

        double muij = 0.0;
        double muip = 0.0;
        double mujp = 0.0;

        for (int s=0; s<ns; ++s){
            muij += var(i, j, 0,3+s)*cd(6+3*s+2)/rho(i,j);
            muip += var(i+1,j,0,3+s)*cd(6+3*s+2)/rho(i+1,j);
            mujp += var(i,j+1,0,3+s)*cd(6+3*s+2)/rho(i,j+1);
        }

        //double muij = (var(i,j,0,3)*mu1 + var(i,j,0,4)*mu2)/rho(i,j);
        //double muip = (var(i+1,j,0,3)*mu1 + var(i+1,j,0,4)*mu2)/rho(i+1,j);
        //double mujp = (var(i,j+1,0,3)*mu1 + var(i,j+1,0,4)*mu2)/rho(i,j+1);

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
        
        //stressx(i,j,0,0) = (4.0/3.0)*mur*dudx;
        //stressx(i,j,0,0) = mu*dudx;
        stressx(i,j,0,0) = (2.0/3.0)*mur*(2.0*dudx-dvdy);
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

        //stressy(i,j,0,0) = (2.0/3.0)*rhoy*mut*(2.0*dudx);
        //stressy(i,j,0,0) = 0.0;
        stressy(i,j,0,0) = (2.0/3.0)*mut*(2.0*dudx-dvdy);
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
    
    grid    = Kokkos::View<double****,FS_LAYOUT>("coords", cf.ni, cf.nj, cf.nk, 3);
    var     = Kokkos::View<double****,FS_LAYOUT>("var",    cf.ngi,cf.ngj,cf.ngk,cf.nvt); // Primary Variable Array
    tmp1    = Kokkos::View<double****,FS_LAYOUT>("tmp1",   cf.ngi,cf.ngj,cf.ngk,cf.nvt); // Temporary Variable Arrayr1
    //tmp2    = Kokkos::View<double****,FS_LAYOUT>("tmp2",   cf.ngi,cf.ngj,cf.ngk,cf.nvt); // Temporary Variable Array2
    dvar    = Kokkos::View<double****,FS_LAYOUT>("dvar",   cf.ngi,cf.ngj,cf.ngk,cf.nvt); // RHS Output
    p       = Kokkos::View<double**  ,FS_LAYOUT>("p",      cf.ngi,cf.ngj);                 // Pressure
    T       = Kokkos::View<double**  ,FS_LAYOUT>("T",      cf.ngi,cf.ngj);                 // Temperature
    rho     = Kokkos::View<double**  ,FS_LAYOUT>("rho",    cf.ngi,cf.ngj);                 // Total Density
    fluxx   = Kokkos::View<double**  ,FS_LAYOUT>("fluxx",  cf.ngi,cf.ngj);                 // Advective Fluxes in X direction
    fluxy   = Kokkos::View<double**  ,FS_LAYOUT>("fluxy",  cf.ngi,cf.ngj);                 // Advective Fluxes in Y direction
    if (cf.visc == 1){
        qx      = Kokkos::View<double**  ,FS_LAYOUT>("qx",     cf.ngi,cf.ngj);                 // Heat Fluxes in X direction
        qy      = Kokkos::View<double**  ,FS_LAYOUT>("qy",     cf.ngi,cf.ngj);                 // Heat Fluxes in X direction
        stressx = Kokkos::View<double****,FS_LAYOUT>("stressx",cf.ngi,cf.ngj,2,2);             // stress tensor on x faces
        stressy = Kokkos::View<double****,FS_LAYOUT>("stressy",cf.ngi,cf.ngj,2,2);             // stress tensor on y faces
    }
    if (cf.ceq == 1){
        gradRho = Kokkos::View<double*** ,FS_LAYOUT>("gradRho",cf.ngi,cf.ngj,4);               // stress tensor on y faces
    }
    if (cf.noise == 1){
        noise   = Kokkos::View<int   **  ,FS_LAYOUT>("noise"  ,cf.ngi,cf.ngj);
    }
#ifdef NOMPI
    if (cf.particle == 1){
        particles = Kokkos::View<particleStruct*>("particles",cf.p_np);
    }
#endif
    cd = mcd;

    timers["flux"]       = fiestaTimer("Flux Calculation");
    timers["pressgrad"]  = fiestaTimer("Pressure Gradient Calculation");
    timers["calcSecond"] = fiestaTimer("Secondary Variable Calculation");
    timers["solWrite"] = fiestaTimer("Solution Write Time");
    timers["resWrite"] = fiestaTimer("Restart Write Time");
    timers["statCheck"] = fiestaTimer("Status Check");
    timers["rk"] = fiestaTimer("Runge Stage Update");
    if (cf.gravity == 1){
        timers["gravity"]     = fiestaTimer("Gravity Term");
    }
    if (cf.visc == 1){
        timers["stress"]     = fiestaTimer("Stress Tensor Computation");
        timers["qflux"]      = fiestaTimer("Heat Flux Calculation");
        timers["visc"]       = fiestaTimer("Viscous Term Calculation");
    }
    if (cf.ceq == 1){
        timers["ceq"]       = fiestaTimer("C-Equation");
    }
    if (cf.noise == 1){
        timers["noise"]         = fiestaTimer("Noise Removal");
    }
    if (cf.particle == 1){
        timers["padvect"] = fiestaTimer("Particle Advection");
        timers["pwrite"] = fiestaTimer("Particle Write");
        timers["psetup"] = fiestaTimer("Particle Setup");
    }
};

void hydro2dvisc_func::preStep(){

}

void hydro2dvisc_func::postStep(){

    if (cf.noise == 1){
        int M = 0;
        int N = 0;
        double coff;

        if ((cf.nci-1) % 2 == 0)
            M = (cf.nci-1)/2;
        else
            M = cf.nci/2;

        if ((cf.ncj-1) % 2 == 0)
            N = (cf.ncj-1)/2;
        else
            N = cf.ncj/2;

        //double etat = 5.0e-3;
        //double dh = 5.0e-4;
        //double doff = 0.1;

        policy_f noise_pol = policy_f({0,0},{M,N});
        policy_f cell_pol  = policy_f({cf.ng,cf.ng},{cf.ngi-cf.ng, cf.ngj-cf.ng});

        if (cf.ceq == 1){
            double maxCh;
            Kokkos::parallel_reduce(cell_pol,maxCvar2D(var,1,cd), Kokkos::Max<double>(maxCh));
#ifndef NOMPI
            MPI_Allreduce(&maxCh,&maxCh,1,MPI_DOUBLE,MPI_MAX,cf.comm);
#endif
            coff = cf.n_coff*maxCh;

        }else{
            coff = 0.0;
        }


        timers["noise"].reset();
        for (int v=0; v<2; ++v){
        //int v = 2;
            Kokkos::parallel_for( noise_pol, detectNoise2D(var,noise,cf.n_dh,coff,cd,v) );
            for (int tau=0; tau<cf.n_nt; ++tau){
                Kokkos::parallel_for( cell_pol,  removeNoise2D(dvar,var,noise,cf.n_eta,cd,v) );
                Kokkos::parallel_for( cell_pol,  updateNoise2D(dvar,var,v) );
            }
        }
        Kokkos::fence();
        timers["noise"].accumulate();
    } //end noise

#ifdef NOMPI
    if (cf.particle == 1){

        Kokkos::View<particleStruct*>::HostMirror particlesH = Kokkos::create_mirror_view(particles);
        Kokkos::deep_copy(particlesH, particles);

        if (cf.write_freq > 0){
            if ((cf.t) % cf.write_freq == 0){
                timers["pwrite"].reset();
                stringstream ss;
                ss << "particle-" << setw(7) << setfill('0') << cf.t << ".vtk";
                ofstream f;
                //f.open("particle.vtk");
                f.open(ss.str());
                f << "# vtk DataFile Version 4.2" << endl;
                f << "Test Particles" << endl;
                f << "ASCII" << endl;
                f << "DATASET POLYDATA" << endl;
                f << "POINTS " << cf.p_np << " float" << endl;
                for (int p=0; p<cf.p_np; ++p){
                    f << particlesH(p).x << " " << particlesH(p).y << " " << "0.0" << endl;
                }
                f << "VERTICES " << cf.p_np << " " << cf.p_np*2 <<endl;
                for (int p=0; p<cf.p_np; ++p){
                    f << "1 " << p << endl;
                }
                f << "POINT_DATA " << cf.p_np << endl;
                f << "SCALARS state float" << endl;
                f << "LOOKUP_TABLE default" << endl;
                for (int p=0; p<cf.p_np; ++p){
                    f << particlesH(p).state << endl;
                }
                f << "SCALARS ci float" << endl;
                f << "LOOKUP_TABLE default" << endl;
                for (int p=0; p<cf.p_np; ++p){
                    f << particlesH(p).ci << endl;
                }
                f << "SCALARS cj float" << endl;
                f << "LOOKUP_TABLE default" << endl;
                for (int p=0; p<cf.p_np; ++p){
                    f << particlesH(p).cj << endl;
                }
                f.flush();
                f.close();
                Kokkos::fence();
                timers["pwrite"].accumulate();
                //for (int p=0; p<cf.p_np; ++p){
                //    cout << p << ":(" << particlesH(p).ci << ", " << particlesH(p).cj << ") ";
                //}
                //cout << endl;
            }
        } // end particle write

        // calculate density and othe secondary variables
        policy_f ghost_pol = policy_f({0,0},{cf.ngi, cf.ngj});

        // Calcualte Total Density and Pressure Fields
        timers["calcSecond"].reset();
        Kokkos::parallel_for( ghost_pol, calculateRhoAndPressure2dv(var,p,rho,T,cd) );
        Kokkos::fence();
        timers["calcSecond"].accumulate();

        // advect particles
        timers["padvect"].reset();
        Kokkos::parallel_for( cf.p_np, advectParticles2D(var,rho,grid,particles,cf.dt,cf.nci,cf.ncj,cf.ng) );

//        // find new cells
//        policy_f grid_pol  = policy_f({0,0},{cf.nci, cf.ncj});
//        for (int p=0; p<cf.p_np; ++p){
//            Kokkos::parallel_for( grid_pol, findInitialCell2D(grid,particles,p,cf.ng) );
//        }
        Kokkos::fence();
        timers["padvect"].accumulate();
    }
#endif

}
void hydro2dvisc_func::preSim(){

    if (cf.particle == 1){
        timers["psetup"].reset();
        Kokkos::View<particleStruct*>::HostMirror particlesH = Kokkos::create_mirror_view(particles);

        double ymax = 6.0;
        double xmax = 2.0;
        double dpx = xmax/((double)(cf.p_np-1));
        double dpy = 9.0/((double)(cf.p_np-1));
        for (int p=0; p<cf.p_np; ++p){
            particlesH(p).state = 1;
            particlesH(p).x = 1.0;
            particlesH(p).y = 0.5+p*dpy;
            //cout << p << " " << particlesH(p).state;
            //cout      << " " << particlesH(p).x;
            //cout      << " " << particlesH(p).y;
            //cout << endl;
        }

        // find initial cell id
        policy_f grid_pol  = policy_f({0,0},{cf.nci,cf.ncj});
        for (int p=0; p<cf.p_np; ++p){
            Kokkos::parallel_for( grid_pol, findInitialCell2D(grid,particles,p,cf.ng) );
        }
        Kokkos::fence();
        timers["psetup"].accumulate();

        Kokkos::deep_copy(particles, particlesH);

        if (cf.write_freq > 0){
            stringstream ss;
            ss << "particle-" << setw(7) << setfill('0') << cf.t << ".vtk";
            ofstream f;
            //f.open("particle.vtk");
            f.open(ss.str());
            f << "# vtk DataFile Version 4.2" << endl;
            f << "Test Particles" << endl;
            f << "ASCII" << endl;
            f << "DATASET POLYDATA" << endl;
            f << "POINTS " << cf.p_np << " float" << endl;
            for (int p=0; p<cf.p_np; ++p){
                f << particlesH(p).x << " " << particlesH(p).y << " " << "0.0" << endl;
            }
            f << "VERTICES " << cf.p_np << " " << cf.p_np*2 <<endl;
            for (int p=0; p<cf.p_np; ++p){
                f << "1 " << p << endl;
            }
            f << "POINT_DATA " << cf.p_np << endl;
            f << "SCALARS state float" << endl;
            f << "LOOKUP_TABLE default" << endl;
            for (int p=0; p<cf.p_np; ++p){
                f << particlesH(p).state << endl;
            }
            f << "SCALARS ci float" << endl;
            f << "LOOKUP_TABLE default" << endl;
            for (int p=0; p<cf.p_np; ++p){
                f << particlesH(p).ci << endl;
            }
            f << "SCALARS cj float" << endl;
            f << "LOOKUP_TABLE default" << endl;
            for (int p=0; p<cf.p_np; ++p){
                f << particlesH(p).cj << endl;
            }
        } // end initial write
    }
        

}
void hydro2dvisc_func::postSim(){}


void hydro2dvisc_func::compute(){

    policy_f ghost_pol = policy_f({0,0},{cf.ngi, cf.ngj});
    policy_f cell_pol  = policy_f({cf.ng,cf.ng},{cf.ngi-cf.ng, cf.ngj-cf.ng});
    policy_f face_pol  = policy_f({cf.ng-1,cf.ng-1},{cf.ngi-cf.ng, cf.ngj-cf.ng});


    // Calcualte Total Density and Pressure Fields
    timers["calcSecond"].reset();
    Kokkos::parallel_for( ghost_pol, calculateRhoAndPressure2dv(var,p,rho,T,cd) );
    Kokkos::fence();
    timers["calcSecond"].accumulate();

    // Calculate and apply weno fluxes for each variable
    for (int v=0; v<cf.nv; ++v){
        timers["flux"].reset();
        if (cf.scheme == 3){
            Kokkos::parallel_for( face_pol, computeFluxQuick2D(var,p,rho,fluxx,fluxy,cd,v) );
        }else if (cf.scheme == 2){
            Kokkos::parallel_for( face_pol, computeFluxCentered2D(var,p,rho,fluxx,fluxy,cd,v) );
        }else{
            Kokkos::parallel_for( face_pol, computeFluxWeno2D(var,p,rho,fluxx,fluxy,cd,v) );
        }
        Kokkos::fence();
        Kokkos::parallel_for( cell_pol, advect2dv(dvar,fluxx,fluxy,cd,v) );
        Kokkos::fence();
        timers["flux"].accumulate();
    }

    // Apply Pressure Gradient Term
    timers["pressgrad"].reset();
    Kokkos::parallel_for( cell_pol, applyPressure2dv(dvar,p,cd) );
    Kokkos::fence();
    timers["pressgrad"].accumulate();

    if (cf.gravity == 1){
        timers["gravity"].reset();
        Kokkos::parallel_for( cell_pol, calculateGravity(dvar,var,rho,cf.g_accel,cf.g_vec[0],cf.g_vec[1]) );
        Kokkos::fence();
        timers["gravity"].accumulate();
    }

    if (cf.visc == 1){
        timers["stress"].reset();
        Kokkos::parallel_for( face_pol, calculateStressTensor2dv(var,rho,T,stressx,stressy,cd) );
        Kokkos::fence();
        timers["stress"].accumulate();

        timers["qflux"].reset();
        Kokkos::parallel_for( face_pol, calculateHeatFlux2dv(var,rho,T,qx,qy,cd) );
        Kokkos::fence();
        timers["qflux"].accumulate();

        timers["visc"].reset();
        Kokkos::parallel_for( cell_pol, applyViscousTerm2dv(dvar,var,rho,stressx,stressy,qx,qy,cd) );
        Kokkos::fence();
        timers["visc"].accumulate();
    }

    if (cf.ceq == 1){
        double maxS;
        double maxG;
        double maxGh;
        double maxTau1;
        double maxTau2;
        double mu;
        double maxC;
        double maxCh;
        double maxT1;
        double maxT2;

        timers["ceq"].reset();
        Kokkos::parallel_reduce(cell_pol,maxWaveSpeed2D(var,p,rho,cd), Kokkos::Max<double>(maxS));
        Kokkos::parallel_reduce(cell_pol,maxGradRho2D(gradRho,1), Kokkos::Max<double>(maxG));
        Kokkos::parallel_reduce(cell_pol,maxGradRho2D(gradRho,1), Kokkos::Max<double>(maxGh));
        Kokkos::parallel_reduce(cell_pol,maxGradRho2D(gradRho,1), Kokkos::Max<double>(maxTau1));
        Kokkos::parallel_reduce(cell_pol,maxGradRho2D(gradRho,1), Kokkos::Max<double>(maxTau2));
        Kokkos::parallel_reduce(cell_pol,maxCvar2D(var,0,cd), Kokkos::Max<double>(maxC));
        Kokkos::parallel_reduce(cell_pol,maxCvar2D(var,1,cd), Kokkos::Max<double>(maxCh));
        Kokkos::parallel_reduce(cell_pol,maxCvar2D(var,2,cd), Kokkos::Max<double>(maxT1));
        Kokkos::parallel_reduce(cell_pol,maxCvar2D(var,3,cd), Kokkos::Max<double>(maxT2));

#ifndef NOMPI
        MPI_Allreduce(&maxS,&maxS,1,MPI_DOUBLE,MPI_MAX,cf.comm);
        MPI_Allreduce(&maxG,&maxG,1,MPI_DOUBLE,MPI_MAX,cf.comm);
        MPI_Allreduce(&maxGh,&maxGh,1,MPI_DOUBLE,MPI_MAX,cf.comm);
        MPI_Allreduce(&maxTau1,&maxTau1,1,MPI_DOUBLE,MPI_MAX,cf.comm);
        MPI_Allreduce(&maxTau1,&maxTau2,1,MPI_DOUBLE,MPI_MAX,cf.comm);
        MPI_Allreduce(&maxC,&maxC,1,MPI_DOUBLE,MPI_MAX,cf.comm);
        MPI_Allreduce(&maxCh,&maxCh,1,MPI_DOUBLE,MPI_MAX,cf.comm);
        MPI_Allreduce(&maxT1,&maxT1,1,MPI_DOUBLE,MPI_MAX,cf.comm);
        MPI_Allreduce(&maxT2,&maxT2,1,MPI_DOUBLE,MPI_MAX,cf.comm);
#endif


        if (maxG < 1e-6) maxG = 1e-6;
        if (maxGh < 1e-6) maxGh = 1e-6;
        if (maxTau1 < 1e-6) maxTau1 = 1e-6;
        if (maxTau2 < 1e-6) maxTau2 = 1e-6;
        if (maxC < 1e-6) maxC = 1e-6;
        if (maxCh < 1e-6) maxCh = 1e-6;
        if (maxT1 < 1e-6) maxT1 = 1e-6;
        if (maxT2 < 1e-6) maxT2 = 1e-6;

        mu = maxT1;
        if (maxT2 > mu) mu = maxT2;

        //printf("### DEBUG ### %f : %f : %f : %f : %f\n",maxC,maxCh,maxTau1,maxTau2,mu);

        Kokkos::parallel_for(cell_pol, calculateRhoGrad2D(var,rho,gradRho,cd));

        Kokkos::parallel_for(cell_pol, updateCeq2D(dvar,var,gradRho,maxS,cd,cf.kap,cf.eps,maxG,maxGh,maxTau1,maxTau2));

        Kokkos::parallel_for(cell_pol, applyCeq2D(dvar,var,rho,cf.beta,cf.betae,cf.alpha,maxC,maxCh,mu,cd));

        Kokkos::fence();
        timers["ceq"].accumulate();
    }
}
