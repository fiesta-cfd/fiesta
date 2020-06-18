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
    FS3D vel;

    computeTransformedVelocity2D (FS4D var_, FS4D m_, FS2D r_, FS3D v_)
         : var(var_), metrics(m_), rho(r_), vel(v_) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const int i, const int j) const {

        vel(i,j,0) = ( metrics(i,j,0,0)*var(i,j,0,0)/rho(i,j) + metrics(i,j,0,1)*var(i,j,0,1)/rho(i,j) );
        vel(i,j,1) = ( metrics(i,j,1,0)*var(i,j,0,0)/rho(i,j) + metrics(i,j,1,1)*var(i,j,0,1)/rho(i,j) );
    }
};


struct advectGen2D {
    
    FS4D dvar;
    FS2D fluxx;
    FS2D fluxy;
    FS1D cd;
    int v;

    advectGen2D (FS4D dvar_, FS2D fx_, FS2D fy_, FS1D cd_, int v_)
        : dvar(dvar_), fluxx(fx_), fluxy(fy_), cd(cd_), v(v_) {}
    
    KOKKOS_INLINE_FUNCTION
    void operator()(const int i, const int j) const {

        double a;

        a = -( (fluxx(i,j) - fluxx(i-1,j))
                          +(fluxy(i,j) - fluxy(i,j-1)) );

        dvar(i,j,0,v) = a;
    }
};

struct applyGenPressure2D {
    
    FS4D dvar;
    FS2D p;
    FS4D m;
    FS1D cd;

    applyGenPressure2D (FS4D dvar_, FS4D m_, FS2D p_, FS1D cd_)
        : dvar(dvar_), m(m_), p(p_), cd(cd_) {}
    
    KOKKOS_INLINE_FUNCTION
    void operator()(const int i, const int j) const {
        
        double dxipxix;
        double detpetx;
        double dxipxiy;
        double detpety;

        dxipxix = ( p(i-2,j)*m(i-2,j,0,0) - 8.0*p(i-1,j)*m(i-1,j,0,0) + 8.0*p(i+1,j)*m(i+1,j,0,0) - p(i+2,j)*m(i+2,j,0,0) )/(12.0);
        detpetx = ( p(i,j-2)*m(i,j-2,1,0) - 8.0*p(i,j-1)*m(i,j-1,1,0) + 8.0*p(i,j+1)*m(i,j+1,1,0) - p(i,j+2)*m(i,j+2,1,0) )/(12.0);
        dxipxiy = ( p(i-2,j)*m(i-2,j,0,1) - 8.0*p(i-1,j)*m(i-1,j,0,1) + 8.0*p(i+1,j)*m(i+1,j,0,1) - p(i+2,j)*m(i+2,j,0,1) )/(12.0);
        detpety = ( p(i,j-2)*m(i,j-2,1,1) - 8.0*p(i,j-1)*m(i,j-1,1,1) + 8.0*p(i,j+1)*m(i,j+1,1,1) - p(i,j+2)*m(i,j+2,1,1) )/(12.0);

        dvar(i,j,0,0) = ( dvar(i,j,0,0) - (dxipxix + detpetx) );
        dvar(i,j,0,1) = ( dvar(i,j,0,1) - (dxipxiy + detpety) );
    }
};

gen2d_func::gen2d_func(struct inputConfig &cf_, Kokkos::View<double*> & cd_):rk_func(cf_,cd_){
    
    grid    = Kokkos::View<double****,FS_LAYOUT>("coords", cf.ni, cf.nj, cf.nk, 3);
    metrics = Kokkos::View<double****,FS_LAYOUT>("metrics", cf.ngi, cf.ngj, 2, 2);
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
    Kokkos::parallel_for( cell_pol, computeMetrics2D(metrics,grid) );

#ifndef NOMPI
    //mpi exchange of metrics
    ls = Kokkos::View<double****,FS_LAYOUT>("leftSend",cf.ng,cf.ngj,2,2);
    lr = Kokkos::View<double****,FS_LAYOUT>("leftRecv",cf.ng,cf.ngj,2,2);
    rs = Kokkos::View<double****,FS_LAYOUT>("rightSend",cf.ng,cf.ngj,2,2);
    rr = Kokkos::View<double****,FS_LAYOUT>("rightRecv",cf.ng,cf.ngj,2,2);
    bs = Kokkos::View<double****,FS_LAYOUT>("bottomSend",cf.ngi,cf.ng,2,2);
    br = Kokkos::View<double****,FS_LAYOUT>("bottomRecv",cf.ngi,cf.ng,2,2);
    ts = Kokkos::View<double****,FS_LAYOUT>("topSend",cf.ngi,cf.ng,2,2);
    tr = Kokkos::View<double****,FS_LAYOUT>("topRecv",cf.ngi,cf.ng,2,2);
    
    lsH = Kokkos::create_mirror_view(ls);
    lrH = Kokkos::create_mirror_view(lr);
    rsH = Kokkos::create_mirror_view(rs);
    rrH = Kokkos::create_mirror_view(rr);
    bsH = Kokkos::create_mirror_view(bs);
    brH = Kokkos::create_mirror_view(br);
    tsH = Kokkos::create_mirror_view(ts);
    trH = Kokkos::create_mirror_view(tr);

    policy_f4 xPol = policy_f4({0,0,0,0},{cf.ng,cf.ngj,2,2});
    policy_f4 yPol = policy_f4({0,0,0,0},{cf.ngi,cf.ng,2,2});

    Kokkos::parallel_for( xPol, KOKKOS_LAMBDA (const int i, const int j, const int k, const int l){
        ls(i,j,k,l) = metrics(cf.ng+i,j,k,l);
        rs(i,j,k,l) = metrics(cf.nci+i,j,k,l);
    });
    Kokkos::parallel_for( yPol, KOKKOS_LAMBDA (const int i, const int j, const int k, const int l){
        bs(i,j,k,l) = metrics(i,cf.ng+j,k,l);
        ts(i,j,k,l) = metrics(i,cf.ncj+j,k,l);
    });
    Kokkos::deep_copy(lsH,ls);
    Kokkos::deep_copy(rsH,rs);
    Kokkos::deep_copy(bsH,bs);
    Kokkos::deep_copy(tsH,ts);

    MPI_Request reqs[8];

    MPI_Isend(lsH.data(), cf.ng*cf.ngj*4, MPI_DOUBLE, cf.xMinus,           0, cf.comm, &reqs[0]);
    MPI_Irecv(lrH.data(), cf.ng*cf.ngj*4, MPI_DOUBLE, cf.xMinus, MPI_ANY_TAG, cf.comm, &reqs[1]);

    MPI_Isend(rsH.data(), cf.ng*cf.ngj*4, MPI_DOUBLE, cf.xPlus,           0, cf.comm, &reqs[2]);
    MPI_Irecv(rrH.data(), cf.ng*cf.ngj*4, MPI_DOUBLE, cf.xPlus, MPI_ANY_TAG, cf.comm, &reqs[3]);

    MPI_Isend(bsH.data(), cf.ngi*cf.ng*4, MPI_DOUBLE, cf.yMinus,           0, cf.comm, &reqs[4]);
    MPI_Irecv(brH.data(), cf.ngi*cf.ng*4, MPI_DOUBLE, cf.yMinus, MPI_ANY_TAG, cf.comm, &reqs[5]);

    MPI_Isend(tsH.data(), cf.ngi*cf.ng*4, MPI_DOUBLE, cf.yPlus,           0, cf.comm, &reqs[6]);
    MPI_Irecv(trH.data(), cf.ngi*cf.ng*4, MPI_DOUBLE, cf.yPlus, MPI_ANY_TAG, cf.comm, &reqs[7]);

    MPI_Waitall(8, reqs, MPI_STATUS_IGNORE);

    Kokkos::deep_copy(lr,lrH);
    Kokkos::deep_copy(rr,rrH);
    Kokkos::deep_copy(br,brH);
    Kokkos::deep_copy(tr,trH);

    Kokkos::parallel_for( xPol, KOKKOS_LAMBDA (const int i, const int j, const int k, const int l){
        metrics(i,j,k,l) = lr(i,j,k,l);
        metrics(cf.ngi-cf.ng+i,j,k,l) = rr(i,j,k,l);
    });
    Kokkos::parallel_for( yPol, KOKKOS_LAMBDA (const int i, const int j, const int k, const int l){
        metrics(i,j,k,l) = br(i,j,k,l);
        metrics(i,cf.ngj-cf.ng+j,k,l) = tr(i,j,k,l);
    });
#endif

    if (cf.xMinus < 0){
        Kokkos::parallel_for("metricbcl",policy_f({0,0},{cf.ng,cf.ngi}),
            KOKKOS_LAMBDA (const int g, const int i){
                metrics(i,cf.ng-g-1,0,0) = metrics(i,cf.ng+g,0,0);
                metrics(i,cf.ng-g-1,0,1) = metrics(i,cf.ng+g,0,1);
                metrics(i,cf.ng-g-1,1,0) = metrics(i,cf.ng+g,1,0);
                metrics(i,cf.ng-g-1,1,1) = metrics(i,cf.ng+g,1,1);
        });
    }
    if (cf.xPlus < 0){
        Kokkos::parallel_for("metricbcr",policy_f({0,0},{cf.ng,cf.ngi}),
            KOKKOS_LAMBDA (const int g, const int i){
                metrics(i,cf.ngi-1-g,0,0) = metrics(i,cf.nci+g,0,0);
                metrics(i,cf.ngi-1-g,0,1) = metrics(i,cf.nci+g,0,1);
                metrics(i,cf.ngi-1-g,1,0) = metrics(i,cf.nci+g,1,0);
                metrics(i,cf.ngi-1-g,1,1) = metrics(i,cf.nci+g,1,1);
        });
    }
    if (cf.yMinus < 0){
        Kokkos::parallel_for("metricbcb",policy_f({0,0},{cf.ng,cf.ngj}),
            KOKKOS_LAMBDA (const int g, const int j){
                metrics(cf.ng-g-1,j,0,0) = metrics(cf.ng+g,j,0,0);
                metrics(cf.ng-g-1,j,0,1) = metrics(cf.ng+g,j,0,1);
                metrics(cf.ng-g-1,j,1,0) = metrics(cf.ng+g,j,1,0);
                metrics(cf.ng-g-1,j,1,1) = metrics(cf.ng+g,j,1,1);
        });
    }
    if (cf.yPlus < 0){
        Kokkos::parallel_for("metricbct",policy_f({0,0},{cf.ng,cf.ngj}),
            KOKKOS_LAMBDA (const int g, const int j){
                metrics(cf.ngj-1-g,j,0,0) = metrics(cf.ncj+g,j,0,0);
                metrics(cf.ngj-1-g,j,0,1) = metrics(cf.ncj+g,j,0,1);
                metrics(cf.ngj-1-g,j,1,0) = metrics(cf.ncj+g,j,1,0);
                metrics(cf.ngj-1-g,j,1,1) = metrics(cf.ncj+g,j,1,1);
        });
    }


    Kokkos::fence();
    timers["calcMetrics"].accumulate();

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
    Kokkos::parallel_for( ghost_pol, computeTransformedVelocity2D(var,metrics,rho,tvel) );
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
        Kokkos::parallel_for( cell_pol, advectGen2D(dvar,fluxx,fluxy,cd,v) );
        Kokkos::fence();
        timers["flux"].accumulate();
    }

    // Apply Pressure Gradient Term
    timers["pressgrad"].reset();
    Kokkos::parallel_for( cell_pol, applyGenPressure2D(dvar,metrics,p,cd) );
    Kokkos::fence();
    timers["pressgrad"].accumulate();
}
