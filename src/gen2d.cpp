#include "input.hpp"
#ifndef NOMPI
#include "mpi.hpp"
#include "cgns.hpp"
#include <mpi.h>
#endif
#include "bc.hpp"
#include "Kokkos_Core.hpp"
#include "debug.hpp"
#include "gen2d.hpp"
#include "metric.hpp"
#include "flux.hpp"
#include "ceq.hpp"
#include <iomanip>
#include <fstream>

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

struct computeTransformedVelocity2D {
    FS4D var;
    FS4D metrics;
    FS2D rho;
    FS2D J;
    FS3D vel;

    computeTransformedVelocity2D (FS4D var_, FS4D m_, FS2D j_, FS2D r_, FS3D v_)
         : var(var_), metrics(m_), J(j_), rho(r_), vel(v_) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const int i, const int j) const {

        vel(i,j,0) = ( metrics(i,j,0,0)*var(i,j,0,0)/rho(i,j) + metrics(i,j,0,1)*var(i,j,0,1)/rho(i,j) );///J(i,j);
        vel(i,j,1) = ( metrics(i,j,1,0)*var(i,j,0,0)/rho(i,j) + metrics(i,j,1,1)*var(i,j,0,1)/rho(i,j) );///J(i,j);
        //printf("%d, %d, %f, %f, %f, %f\n",i,j, metrics(i,j,0,0), metrics(i,j,0,1), metrics(i,j,1,0), metrics(i,j,1,1));
        //printf("%d, %d, %f, %f\n",i,j,vel(i,j,0),vel(i,j,1));
    }
};


struct advectGen2D {
    
    FS4D dvar;
    FS2D J;
    FS2D fluxx;
    FS2D fluxy;
    Kokkos::View<double*> cd;
    int v;

    advectGen2D (FS4D dvar_, FS2D j_, FS2D fx_, FS2D fy_, Kokkos::View<double*> cd_, int v_)
        : dvar(dvar_), J(j_), fluxx(fx_), fluxy(fy_), cd(cd_), v(v_) {}
    
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

struct applyGenPressure2D {
    
    FS4D dvar;
    FS2D p,J;
    FS4D m;
    FS1D cd;

    applyGenPressure2D (FS4D dvar_, FS4D m_, FS2D j_, FS2D p_, FS1D cd_)
        : dvar(dvar_), m(m_), J(j_), p(p_), cd(cd_) {}
    
    KOKKOS_INLINE_FUNCTION
    void operator()(const int i, const int j) const {
        
        double dxipxix;
        double detpetx;
        double dxipxiy;
        double detpety;

        //dxipxix = ( p(i-2,j)*m(i-2,j,0,0)/J(i-2,j) - 8.0*p(i-1,j)*m(i-1,j,0,0)/J(i-1,j) + 8.0*p(i+1,j)*m(i+1,j,0,0)/J(i+1,j) - p(i+2,j)*m(i+2,j,0,0)/J(i+2,j) )/(12.0);
        //detpetx = ( p(i,j-2)*m(i,j-2,1,0)/J(i,j-2) - 8.0*p(i,j-1)*m(i,j-1,1,0)/J(i,j-1) + 8.0*p(i,j+1)*m(i,j+1,1,0)/J(i,j+1) - p(i,j+2)*m(i,j+2,1,0)/J(i,j+2) )/(12.0);
        //dxipxiy = ( p(i-2,j)*m(i-2,j,0,1)/J(i-2,j) - 8.0*p(i-1,j)*m(i-1,j,0,1)/J(i-1,j) + 8.0*p(i+1,j)*m(i+1,j,0,1)/J(i+1,j) - p(i+2,j)*m(i+2,j,0,1)/J(i+2,j) )/(12.0);
        //detpety = ( p(i,j-2)*m(i,j-2,1,1)/J(i,j-2) - 8.0*p(i,j-1)*m(i,j-1,1,1)/J(i,j-1) + 8.0*p(i,j+1)*m(i,j+1,1,1)/J(i,j+1) - p(i,j+2)*m(i,j+2,1,1)/J(i,j+2) )/(12.0);

        dxipxix = ( p(i-2,j)*m(i-2,j,0,0) - 8.0*p(i-1,j)*m(i-1,j,0,0) + 8.0*p(i+1,j)*m(i+1,j,0,0) - p(i+2,j)*m(i+2,j,0,0) )/(12.0);
        detpetx = ( p(i,j-2)*m(i,j-2,1,0) - 8.0*p(i,j-1)*m(i,j-1,1,0) + 8.0*p(i,j+1)*m(i,j+1,1,0) - p(i,j+2)*m(i,j+2,1,0) )/(12.0);
        dxipxiy = ( p(i-2,j)*m(i-2,j,0,1) - 8.0*p(i-1,j)*m(i-1,j,0,1) + 8.0*p(i+1,j)*m(i+1,j,0,1) - p(i+2,j)*m(i+2,j,0,1) )/(12.0);
        detpety = ( p(i,j-2)*m(i,j-2,1,1) - 8.0*p(i,j-1)*m(i,j-1,1,1) + 8.0*p(i,j+1)*m(i,j+1,1,1) - p(i,j+2)*m(i,j+2,1,1) )/(12.0);

        //dvar(i,j,0,0) = dvar(i,j,0,0) - (dxipxix + detpetx);
        //dvar(i,j,0,1) = dvar(i,j,0,1) - (dxipxiy + detpety);
        dvar(i,j,0,0) = ( dvar(i,j,0,0) - (dxipxix + detpetx) );
        dvar(i,j,0,1) = ( dvar(i,j,0,1) - (dxipxiy + detpety) );
        //printf("%d, %d, %f, %f\n",i,j,dvar(i,j,0,0),dvar(i,j,0,1));
    }
};

gen2d_func::gen2d_func(struct inputConfig &cf_, Kokkos::View<double*> & cd_):rk_func(cf_,cd_){
    
    grid    = Kokkos::View<double****,FS_LAYOUT>("coords", cf.ni, cf.nj, cf.nk, 3);
    metrics = Kokkos::View<double****,FS_LAYOUT>("metrics", cf.ngi, cf.ngj, 2, 2);
    J       = Kokkos::View<double**  ,FS_LAYOUT>("jacobian", cf.ngi, cf.ngj);
    var     = Kokkos::View<double****,FS_LAYOUT>("var",    cf.ngi,cf.ngj,cf.ngk,cf.nvt); // Primary Variable Array
    tmp1    = Kokkos::View<double****,FS_LAYOUT>("tmp1",   cf.ngi,cf.ngj,cf.ngk,cf.nvt); // Temporary Variable Arrayr1
    dvar    = Kokkos::View<double****,FS_LAYOUT>("dvar",   cf.ngi,cf.ngj,cf.ngk,cf.nvt); // RHS Output
    tvel    = Kokkos::View<double*** ,FS_LAYOUT>("vel",    cf.ngi,cf.ngj,2); // RHS Output
    p       = Kokkos::View<double**  ,FS_LAYOUT>("p",      cf.ngi,cf.ngj);                 // Pressure
    T       = Kokkos::View<double**  ,FS_LAYOUT>("T",      cf.ngi,cf.ngj);                 // Temperature
    rho     = Kokkos::View<double**  ,FS_LAYOUT>("rho",    cf.ngi,cf.ngj);                 // Total Density
    fluxx   = Kokkos::View<double**  ,FS_LAYOUT>("fluxx",  cf.ngi,cf.ngj);                 // Advective Fluxes in X direction
    fluxy   = Kokkos::View<double**  ,FS_LAYOUT>("fluxy",  cf.ngi,cf.ngj);                 // Advective Fluxes in Y direction
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
    timers["calcMatrics"] = fiestaTimer("Metric Calculations");
    if (cf.noise == 1){
        timers["noise"]         = fiestaTimer("Noise Removal");
    }
    if (cf.particle == 1){
        timers["padvect"] = fiestaTimer("Particle Advection");
        timers["pwrite"] = fiestaTimer("Particle Write");
        timers["psetup"] = fiestaTimer("Particle Setup");
    }
};

void gen2d_func::preStep(){

}

void gen2d_func::postStep(){

    //policy_f3 var_pol = policy_f3({cf.ng,cf.ng,0},{cf.ngi-cf.ng,cf.ngj-cf.ng,cf.nv});
    //Kokkos::parallel_for("postJac",var_pol,KOKKOS_LAMBDA (const int i, const int j, const int v){
    //    var(i,j,0,v) = J(i,j)*var(i,j,0,v);
//  //      printf("%d, %d, %d, %f, %f\n",i,j,v,var(i,j,0,v),J(i,j));
    //});

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

        // execution policy for all cells including ghost cells
        policy_f ghost_pol = policy_f({0,0},{cf.ngi, cf.ngj});

        // Calcualte Total Density and Pressure Fields
        timers["calcSecond"].reset();
        Kokkos::parallel_for( ghost_pol, calculateRhoAndPressure2dv(var,p,rho,T,cd) );
        Kokkos::fence();
        timers["calcSecond"].accumulate();

        // advect particles
        timers["padvect"].reset();
        Kokkos::parallel_for( cf.p_np, advectParticles2D(var,rho,grid,particles,cf.dt,cf.ng) );
        Kokkos::fence();
        timers["padvect"].accumulate();
    }
#endif

}
void gen2d_func::preSim(){

    // execution policy for all cells including ghost cells
    policy_f cell_pol = policy_f({cf.ng,cf.ng},{cf.ngi-cf.ng, cf.ngj-cf.ng});

    // compute metrics
    timers["calcMetrics"].reset();
    Kokkos::parallel_for( cell_pol, computeMetrics2D(metrics,J,grid) );
    Kokkos::fence();
    timers["calcMetrics"].accumulate();

    Kokkos::parallel_for("metricbcx",policy_f({0,0},{cf.ng,cf.ngi}),
        KOKKOS_LAMBDA (const int g, const int i){
            metrics(i,cf.ng-g-1,0,0) = metrics(i,cf.ng+g,0,0);
            metrics(i,cf.ng-g-1,0,1) = metrics(i,cf.ng+g,0,1);
            metrics(i,cf.ng-g-1,1,0) = metrics(i,cf.ng+g,1,0);
            metrics(i,cf.ng-g-1,1,1) = metrics(i,cf.ng+g,1,1);
            J(i,cf.ng-g-1) = J(i,cf.ng+g);
            metrics(i,cf.ngi-1-g,0,0) = metrics(i,cf.nci+g,0,0);
            metrics(i,cf.ngi-1-g,0,1) = metrics(i,cf.nci+g,0,1);
            metrics(i,cf.ngi-1-g,1,0) = metrics(i,cf.nci+g,1,0);
            metrics(i,cf.ngi-1-g,1,1) = metrics(i,cf.nci+g,1,1);
            J(i,cf.ngi-1-g) = J(i,cf.nci+g);
    });
    Kokkos::parallel_for("metricbcy",policy_f({0,0},{cf.ng,cf.ngj}),
        KOKKOS_LAMBDA (const int g, const int j){
            metrics(cf.ng-g-1,j,0,0) = metrics(cf.ng+g,j,0,0);
            metrics(cf.ng-g-1,j,0,1) = metrics(cf.ng+g,j,0,1);
            metrics(cf.ng-g-1,j,1,0) = metrics(cf.ng+g,j,1,0);
            metrics(cf.ng-g-1,j,1,1) = metrics(cf.ng+g,j,1,1);
            J(cf.ng-g-1,j) = J(cf.ng+g,j);
            metrics(cf.ngj-1-g,j,0,0) = metrics(cf.ncj+g,j,0,0);
            metrics(cf.ngj-1-g,j,0,1) = metrics(cf.ncj+g,j,0,1);
            metrics(cf.ngj-1-g,j,1,0) = metrics(cf.ncj+g,j,1,0);
            metrics(cf.ngj-1-g,j,1,1) = metrics(cf.ncj+g,j,1,1);
            J(cf.ngj-1-g,j) = J(cf.ncj+g,j);
    });

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
void gen2d_func::postSim(){}


void gen2d_func::compute(){

    policy_f ghost_pol = policy_f({0,0},{cf.ngi, cf.ngj});
    policy_f cell_pol  = policy_f({cf.ng,cf.ng},{cf.ngi-cf.ng, cf.ngj-cf.ng});
    policy_f face_pol  = policy_f({cf.ng-1,cf.ng-1},{cf.ngi-cf.ng, cf.ngj-cf.ng});


    // Calcualte Total Density and Pressure Fields
    timers["calcSecond"].reset();
    Kokkos::parallel_for( ghost_pol, calculateRhoAndPressure2dv(var,p,rho,T,cd) );
    Kokkos::parallel_for( ghost_pol, computeTransformedVelocity2D(var,metrics,J,rho,tvel) );
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
            Kokkos::parallel_for( face_pol, computeFluxWeno2D(var,p,rho,tvel,fluxx,fluxy,cd,v) );
        }
        Kokkos::fence();
        Kokkos::parallel_for( cell_pol, advectGen2D(dvar,J,fluxx,fluxy,cd,v) );
        Kokkos::fence();
        timers["flux"].accumulate();
    }

    // Apply Pressure Gradient Term
    timers["pressgrad"].reset();
    Kokkos::parallel_for( cell_pol, applyGenPressure2D(dvar,metrics,J,p,cd) );
    Kokkos::fence();
    timers["pressgrad"].accumulate();
}
