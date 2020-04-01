struct maxWaveSpeed2D {
    
    FS4D var;
    FS2D p;
    FS2D rho;
    Kokkos::View<double*> cd;

    maxWaveSpeed2D (FS4D var_, FS2D p_, FS2D rho_, Kokkos::View<double*> cd_)
        : var(var_), p(p_), rho(rho_), cd(cd_) {}
    
    KOKKOS_INLINE_FUNCTION
    void operator()(const int i, const int j, double& lmax) const {
        int ns = (int)cd(0);
        double gamma, gammas, Rs;
        double Cp = 0;
        double Cv = 0;

        double a, s;

        for (int s=0; s<ns; ++s){
            gammas = cd(5+2*s);
            Rs = cd(5+2*s+1);

            Cp = Cp + (var(i,j,0,4+s)/rho(i,j))*( gammas*Rs/(gammas-1) );
            Cv = Cv + (var(i,j,0,4+s)/rho(i,j))*( Rs/(gammas-1) );
        }

        gamma = Cp/Cv;

        a = sqrt(gamma*p(i,j,0)/rho(i,j));
        s = a + sqrt(var(i,j,0,0)*var(i,j,0,0) + var(i,j,0,1)*var(i,j,0,1));

        if (s > lmax)
            lmax = s;
    }
};

struct calculateRhoGrad2D {

    FS4D var;
    FS2D rho;
    FS3D gradRho;
    Kokkos::View<double*> cd;

    calculateRhoGrad2D (FS4D var_, FS2D rho_, FS3D gradRho_, Kokkos::View<double*> cd_)
        : var(var_), rho(rho_), gradRho(gradRho_), cd(cd_) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const int i, const int j) const {

        int i1 = 0;
        int i2 = 0;

        double ux1 = var(i-2,j,0,0)/rho(i,j-2);
        double ux2 = var(i-1,j,0,0)/rho(i,j-1);
        double ux3 = var(i+1,j,0,0)/rho(i,j+1);
        double ux4 = var(i+2,j,0,0)/rho(i,j+2);

        double uy1 = var(i,j-2,0,0)/rho(i,j-2);
        double uy2 = var(i,j-1,0,0)/rho(i,j-1);
        double uy3 = var(i,j+1,0,0)/rho(i,j+1);
        double uy4 = var(i,j+2,0,0)/rho(i,j+2);

        double vx1 = var(i-2,j,0,1)/rho(i,j-2);
        double vx2 = var(i-1,j,0,1)/rho(i,j-1);
        double vx3 = var(i+1,j,0,1)/rho(i,j+1);
        double vx4 = var(i+2,j,0,1)/rho(i,j+2);

        double vy1 = var(i,j-2,0,1)/rho(i,j-2);
        double vy2 = var(i,j-1,0,1)/rho(i,j-1);
        double vy3 = var(i,j+1,0,1)/rho(i,j+1);
        double vy4 = var(i,j+2,0,1)/rho(i,j+2);

        double ex1 = var(i-2,j,0,2)/rho(i-2,j)-0.5*rho(i-2,j)*(ux1*ux1+vx1*vx1);
        double ex2 = var(i-1,j,0,2)/rho(i-1,j)-0.5*rho(i-1,j)*(ux2*ux2+vx2*vx2);
        double ex3 = var(i+1,j,0,2)/rho(i+1,j)-0.5*rho(i+1,j)*(ux3*ux3+vx3*vx3);
        double ex4 = var(i+2,j,0,2)/rho(i+2,j)-0.5*rho(i+2,j)*(ux4*ux4+vx4*vx4);

        double ey1 = var(i,j-2,0,2)/rho(i,j-2)-0.5*rho(i,j-2)*(uy1*uy1+vy1*vy1);
        double ey2 = var(i,j-1,0,2)/rho(i,j-1)-0.5*rho(i,j-1)*(uy2*uy2+vy2*vy2);
        double ey3 = var(i,j+1,0,2)/rho(i,j+1)-0.5*rho(i,j+1)*(uy3*uy3+vy3*vy3);
        double ey4 = var(i,j+2,0,2)/rho(i,j+2)-0.5*rho(i,j+2)*(uy4*uy4+vy4*vy4);

        double dxr = (rho(i-2,j) - 8.0*rho(i-1,j) + 8.0*rho(i+1,j) - rho(i+2,j))/(12.0*cd(1));
        double dyr = (rho(i,j-2) - 8.0*rho(i,j-1) + 8.0*rho(i,j+1) - rho(i,j+2))/(12.0*cd(2));

        double dxe = (ex1 - 8.0*ex2 + 8.0*ex3 - ex4)/(12.0*cd(1));
        double dye = (ey1 - 8.0*ey2 + 8.0*ey3 - ey4)/(12.0*cd(2));

        double dxu = (ux1 - 8.0*ux2 + 8.0*ux3 - ux4)/(12.0*cd(1));
        double dyv = (vy1 - 8.0*vy2 + 8.0*vy3 - vy4)/(12.0*cd(2));

        double n1 = dxr;
        double n2 = dyr;

        double rgrad = sqrt(dxr*dxr+dyr*dyr);
        double divu = dxu+dyv;

        double dnednr = (n1*dxe+n2*dye)*(n1*dxr+n2*dyr);

        //compression switch
        if (dnednr > 1.0e-6)
            i1 = 1;
        //else
        //    i1 = 0;
        
        // detect shock front (equation 5a)
        if (divu < -1.0e-6)
            i2 = 1;
            //gradRho(i,j,0) = (1.0-indicator)*rgrad;
            //gradRho(i,j,k,0) = (1-indicator)*divu*rgrad;
        //else
        //    i2 = 0;
            //gradRho(i,j,0) = 0.0;

        if (i1 == 0 && i2 == 1)
            gradRho(i,j,0) = rgrad;
        else
            gradRho(i,j,0) = 0.0;

        if (i1 == 1 && i2 == 0){
            gradRho(i,j,1) = rgrad;
            gradRho(i,j,2) =-dxr;
            gradRho(i,j,3) = dyr;
        }else{
            gradRho(i,j,1) = 0.0;
            gradRho(i,j,2) = 0.0;
            gradRho(i,j,3) = 0.0;
        }
            

        //gradRho(i,j,0) = (1.0-indicator)*indicator2*rgrad;
        //gradRho(i,j,0) = indicator2*rgrad;
        // detect contact surface
        //gradRho(i,j,1) = indicator*rgrad;
        //gradRho(i,j,1) = (1.0-indicator2)*indicator*rgrad;

        // gradient components
        //gradRho(i,j,2) = -i1*dxr;
        //gradRho(i,j,3) =  i1*dyr;

        //if (j == 10 && i == 27){
        //    std::cout << gradRho(i-1,j,0) << ", "
        //              << gradRho(i  ,j,0) << ", "
        //              << gradRho(i+1,j,0) << std::endl;
        //}
        //if (j == 10 && i < 50)
        //if (j == 10)
        //    std::cout << i << ":  " << dnednr << endl;
        //    std::cout << gradRho(i-1,j,0) << endl;
    }
};

struct updateCeq2D {

    FS4D dvar;
    FS4D var;
    FS3D gradRho;
    double maxS,kap,eps;
    Kokkos::View<double*> cd;

    updateCeq2D (FS4D dvar_, FS4D var_, FS3D gradRho_, double maxS_, Kokkos::View<double*> cd_,
                      double kap_, double eps_)
        : dvar(dvar_), var(var_), gradRho(gradRho_), maxS(maxS_), cd(cd_), kap(kap_), eps(eps_) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const int i, const int j) const {
        double dx = cd(1);
        double dy = cd(2);

        // calculate cequation variable indices based on number of species (c variables come after species densities)
        int nv = (int)cd(0) + 3;
        int nc  = nv;

        double dx_right,dx_left,dy_top,dy_bot;
        double lap;

        // average cell size
        //double dxmag = pow(dx*dy,1.0/2.0);
        double dxmag = sqrt(dx*dx + dy*dy);
          
        for (int n=0; n<4; ++n){
            dx_right = (var(i+1,j  ,0,nc+n) - var(i  ,j  ,0,nc+n))/dx;
            dx_left  = (var(i  ,j  ,0,nc+n) - var(i-1,j  ,0,nc+n))/dx;

            dy_top   = (var(i  ,j+1,0,nc+n) - var(i  ,j  ,0,nc+n))/dy;
            dy_bot   = (var(i  ,j  ,0,nc+n) - var(i  ,j-1,0,nc+n))/dy;

            // update ceq right hand side
            lap = 0.0;
            lap = (dx_right - dx_left)/dx + (dy_top - dy_bot)/dy;
            dvar(i,j,0,nc+n) = maxS/(eps*dxmag)*(gradRho(i,j,n) - var(i,j,0,nc+n)) + kap*maxS*dxmag*lap;
        }
    }
};
