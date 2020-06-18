#ifndef VISCOSOTY_HPP
#define VISCOSOTY_HPP

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
#endif
