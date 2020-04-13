struct maxCvar2D {
    
    FS4D dvar;
    int v;
    Kokkos::View<double*> cd;

    maxCvar2D (FS4D dvar_, int v_, Kokkos::View<double*> cd_) : dvar(dvar_), v(v_), cd(cd_){}
    
    KOKKOS_INLINE_FUNCTION
    void operator()(const int i, const int j, double& lmax) const {
        int ns = (int)cd(0);
        int nv = ns + 3;

        double s = dvar(i,j,0,nv+v);

        if (s > lmax)
            lmax = s;
    }
};

struct maxGradRho2D {
    
    FS3D gradRho;
    int v;

    maxGradRho2D (FS3D gradRho_, int v_) : gradRho(gradRho_), v(v_) {}
    
    KOKKOS_INLINE_FUNCTION
    void operator()(const int i, const int j, double& lmax) const {

        double s = gradRho(i,j,v);

        if (s > lmax)
            lmax = s;
    }
};

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

            Cp = Cp + (var(i,j,0,3+s)/rho(i,j))*( gammas*Rs/(gammas-1) );
            Cv = Cv + (var(i,j,0,3+s)/rho(i,j))*( Rs/(gammas-1) );
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

        double i1 = 0;
        double i2 = 0;

        double ux1 = var(i-2,j,0,0)/rho(i-2,j);
        double ux2 = var(i-1,j,0,0)/rho(i-1,j);
        double ux3 = var(i+1,j,0,0)/rho(i+1,j);
        double ux4 = var(i+2,j,0,0)/rho(i+2,j);

        double uy1 = var(i,j-2,0,0)/rho(i,j-2);
        double uy2 = var(i,j-1,0,0)/rho(i,j-1);
        double uy3 = var(i,j+1,0,0)/rho(i,j+1);
        double uy4 = var(i,j+2,0,0)/rho(i,j+2);

        double vx1 = var(i-2,j,0,1)/rho(i-2,j);
        double vx2 = var(i-1,j,0,1)/rho(i-1,j);
        double vx3 = var(i+1,j,0,1)/rho(i+1,j);
        double vx4 = var(i+2,j,0,1)/rho(i+2,j);

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

        //double dxr = (rho(i+1,j) - rho(i-1,j))/(2.0*cd(1));
        //double dyr = (rho(i,j+1) - rho(i,j-1))/(2.0*cd(2));

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
        if (dnednr < -1.0e-6)
            i1 = 1.0;
        //else
        //    i1 = 0;
        
        // detect shock front (equation 5a)
        if (divu < -1.0e-6){
            i2 = 1.0;
        }

        //i1 = 1.0;
        //i2 = 1.0;
            //gradRho(i,j,0) = (1.0-indicator)*rgrad;
            //gradRho(i,j,k,0) = (1-indicator)*divu*rgrad;
        //else
        //    i2 = 0;
            //gradRho(i,j,0) = 0.0;

        //if (i1 == 0 && i2 == 1){
        //    gradRho(i,j,0) = rgrad;
        ////    printf("### DEBUG ### i2 rgrad (%d,%d) %f\n",i,j,rgrad);
        //}else
        //    gradRho(i,j,0) = 0.0;
        //
        //if (i1 == 1 && i2 == 0){
        ////    printf("### DEBUG ### i1 rgrad (%d,%d) %f, %f, %f\n",i,j,rgrad,-dxr,dyr);
        //    gradRho(i,j,1) = rgrad;
        //    gradRho(i,j,2) =-dxr;
        //    gradRho(i,j,3) = dyr;
        //}else{
        //    gradRho(i,j,1) = 0.0;
        //    gradRho(i,j,2) = 0.0;
        //    gradRho(i,j,3) = 0.0;
        //}
        //

        gradRho(i,j,0) = (1-i1)*i2*rgrad;
        gradRho(i,j,1) = i1*rgrad;
        gradRho(i,j,2) =-i1*dxr;
        gradRho(i,j,3) = i1*dyr;


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
    double maxG,maxGh,maxTau1,maxTau2;

    updateCeq2D (FS4D dvar_, FS4D var_, FS3D gradRho_, double maxS_, Kokkos::View<double*> cd_,
                      double kap_, double eps_, double maxG_, double maxGh_, double maxTau1_, double maxTau2_)
        : dvar(dvar_), var(var_), gradRho(gradRho_), maxS(maxS_), cd(cd_), kap(kap_), eps(eps_), maxG(maxG_),
          maxGh(maxGh_),maxTau1(maxTau1_), maxTau2(maxTau2_) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const int i, const int j) const {
        double dx = cd(1);
        double dy = cd(2);
        double maxGr[5];

        // calculate cequation variable indices based on number of species (c variables come after species densities)
        int nv = (int)cd(0) + 3;
        int nc  = nv;

        double dx_right,dx_left,dy_top,dy_bot;
        double lap;

        // average cell size
        //double dxmag = pow(dx*dy,1.0/2.0);
        double dxmag = sqrt(dx*dx + dy*dy);

        maxGr[0] = maxG;
        maxGr[1] = maxGh;
        maxGr[2] = maxTau1;
        maxGr[3] = maxTau2;
        maxGr[4] = maxTau2;
          
        for (int n=0; n<4; ++n){
            dx_right = (var(i+1,j  ,0,nc+n) - var(i  ,j  ,0,nc+n))/dx;
            dx_left  = (var(i  ,j  ,0,nc+n) - var(i-1,j  ,0,nc+n))/dx;

            dy_top   = (var(i  ,j+1,0,nc+n) - var(i  ,j  ,0,nc+n))/dy;
            dy_bot   = (var(i  ,j  ,0,nc+n) - var(i  ,j-1,0,nc+n))/dy;


            // update ceq right hand side
            lap = 0.0;
            lap = dx_right - dx_left + dy_top - dy_bot;
            //lap = (dx_right - dx_left)/dx + (dy_top - dy_bot)/dy;
            //dvar(i,j,0,nc+n) = (maxS/(eps*dxmag))*(gradRho(i,j,n) - var(i,j,0,nc+n)) + kap*maxS*dxmag*lap;
            dvar(i,j,0,nc+n) = (maxS/(eps*dxmag))*(gradRho(i,j,n)/maxGr[n] - var(i,j,0,nc+n)) + kap*maxS*dxmag*lap;

            //if (n==0)
            //    printf("### DEBUG ### dvar,maxS,dxmag (%d,%d) %f, %f, %f, %f\n",i,j,dvar(i,j,0,nc+n),dxmag,maxS,lap);
        }
    }
};

struct applyCeq2D {
    FS4D dvar;
    FS4D var;
    FS2D rho;
    Kokkos::View<double*> cd;
    double betau, betae, maxC, maxCh, alpha, mu;

    applyCeq2D (FS4D dvar_, FS4D var_, FS2D rho_,
                double betau_, double betae_, double alpha_, double maxC_, double maxCh_, double mu_,
                Kokkos::View<double*> cd_)
               : dvar(dvar_), var(var_), rho(rho_),
                 betau(betau_), betae(betae_), alpha(alpha_), maxC(maxC_), maxCh(maxCh_), mu(mu_),
                 cd(cd_) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const int i, const int j) const {
        int nv = (int)cd(0) + 3;
        int nc  = nv;
        double dx = cd(1);
        double dy = cd(2);
        double dxmag = dx*dx+dy*dy;

        double bTildU = (dxmag/maxC)*betau;
        double bTildE = (dxmag/maxC)*betae;
        double aTilde = (dxmag/(mu*mu*maxCh))*alpha;

        double u  = var(i  ,j  ,0,0)/rho(i  ,j  ); // u-velocity in current cell
        double ur = var(i+1,j  ,0,0)/rho(i+1,j  ); // u-velocity in right cell
        double ul = var(i-1,j  ,0,0)/rho(i-1,j  ); // u-velocity in left cell
        double ut = var(i  ,j+1,0,0)/rho(i  ,j+1); // u-velocity in top cell
        double ub = var(i  ,j-1,0,0)/rho(i  ,j-1); // u-velocity in bottom cell

        double v  = var(i  ,j  ,0,1)/rho(i  ,j  ); // v-velocity
        double vr = var(i+1,j  ,0,1)/rho(i+1,j  );
        double vl = var(i-1,j  ,0,1)/rho(i-1,j  );
        double vt = var(i  ,j+1,0,1)/rho(i  ,j+1);
        double vb = var(i  ,j-1,0,1)/rho(i  ,j-1);

        double e  = var(i  ,j  ,0,2)/rho(i  ,j  ); // energy
        double er = var(i+1,j  ,0,2)/rho(i+1,j  );
        double el = var(i-1,j  ,0,2)/rho(i-1,j  );
        double et = var(i  ,j+1,0,2)/rho(i  ,j+1);
        double eb = var(i  ,j-1,0,2)/rho(i  ,j-1);

        double rr = (rho(i  ,j  ) + rho(i+1,j  ))/2.0; // density on right face
        double rl = (rho(i-1,j  ) + rho(i  ,j  ))/2.0; // density on left face
        double rt = (rho(i  ,j  ) + rho(i  ,j+1))/2.0; // density on top face
        double rb = (rho(i  ,j-1) + rho(i  ,j  ))/2.0; // density on bottom face

        double cr = (var(i  ,j  ,0,nc) + var(i+1,j  ,0,nc))/2.0; // iso C in right face
        double cl = (var(i-1,j  ,0,nc) + var(i  ,j  ,0,nc))/2.0; // iso C in left face
        double ct = (var(i  ,j  ,0,nc) + var(i  ,j+1,0,nc))/2.0; // iso C in top face
        double cb = (var(i  ,j-1,0,nc) + var(i  ,j  ,0,nc))/2.0; // iso C in bottom face

        double chr = (var(i  ,j  ,0,nc+1) + var(i+1,j  ,0,nc+1))/2.0; // iso C in right face
        double chl = (var(i-1,j  ,0,nc+1) + var(i  ,j  ,0,nc+1))/2.0; // iso C in left face
        double cht = (var(i  ,j  ,0,nc+1) + var(i  ,j+1,0,nc+1))/2.0; // iso C in top face
        double chb = (var(i  ,j-1,0,nc+1) + var(i  ,j  ,0,nc+1))/2.0; // iso C in bottom face

        double t1r = (var(i  ,j  ,0,nc+2) + var(i+1,j  ,0,nc+2))/2.0; // iso C in right face
        double t1l = (var(i-1,j  ,0,nc+2) + var(i  ,j  ,0,nc+2))/2.0; // iso C in left face
        double t1t = (var(i  ,j  ,0,nc+2) + var(i  ,j+1,0,nc+2))/2.0; // iso C in top face
        double t1b = (var(i  ,j-1,0,nc+2) + var(i  ,j  ,0,nc+2))/2.0; // iso C in bottom face

        double t2r = (var(i  ,j  ,0,nc+3) + var(i+1,j  ,0,nc+3))/2.0; // iso C in right face
        double t2l = (var(i-1,j  ,0,nc+3) + var(i  ,j  ,0,nc+3))/2.0; // iso C in left face
        double t2t = (var(i  ,j  ,0,nc+3) + var(i  ,j+1,0,nc+3))/2.0; // iso C in top face
        double t2b = (var(i  ,j-1,0,nc+3) + var(i  ,j  ,0,nc+3))/2.0; // iso C in bottom face

        double dur = (ur- u )/dx;  // u-velocity on right face
        double dul = (u - ul)/dx;  // u-velocity on left face
        double dut = (ut- u )/dy;  // u-velocity on top face
        double dub = (u - ub)/dy;  // u-velocity on bottom face

        double dvr = (vr- v )/dx;  // v-velocity on right face
        double dvl = (v - vl)/dx;  // v-velocity on left face
        double dvt = (vt- v )/dy;  // v-velocity on top face
        double dvb = (v - vb)/dy;  // v-velocity on bottom face

        double der = (er- e )/dx;  // energy on right face
        double del = (e - el)/dx;  // energy on left face
        double det = (et- e )/dy;  // energy on top face
        double deb = (e - eb)/dy;  // energy on bottom face

        // isotropic diffusion operator
        double diffu = (bTildU/dx)*(rr*cr*dur - rl*cl*dul) + 
                       (bTildU/dy)*(rt*ct*dut - rb*cb*dub);
        double diffv = (bTildU/dx)*(rr*cr*dvr - rl*cl*dvl) + 
                       (bTildU/dy)*(rt*ct*dvt - rb*cb*dvb);
        double diffe = (bTildU/dx)*(rr*cr*der - rl*cl*del) + 
                       (bTildU/dy)*(rt*ct*det - rb*cb*deb);

        dvar(i,j,0,0) += diffu;
        dvar(i,j,0,1) += diffv;
        dvar(i,j,0,2) += diffe;

        diffu = aTilde*( (rr*chr*t1r*t1r*dur - rl*chl*t1l*t1l*dul)/dx
                       +(rt*cht*t1t*t1t*dut - rb*chb*t1b*t1b*dub)/dy
                       +(rr*chr*t1r*t2r*dur - rl*chl*t1l*t2l*dul)/dx
                       +(rt*cht*t1t*t2t*dut - rb*chb*t1b*t2b*dub)/dy);

        diffv = aTilde*( (rr*chr*t1r*t1r*dvr - rl*chl*t1l*t1l*dvl)/dx
                       +(rt*cht*t1t*t1t*dvt - rb*chb*t1b*t1b*dvb)/dy
                       +(rr*chr*t1r*t2r*dvr - rl*chl*t1l*t2l*dvl)/dx
                       +(rt*cht*t1t*t2t*dvt - rb*chb*t1b*t2b*dvb)/dy);

        dvar(i,j,0,0) += diffu;
        dvar(i,j,0,1) += diffv;

        //if (diffu > 0.0 || diffv > 0.0)
        //    printf("### DEBUG ### %e: %e - %e\n",aTilde,diffu,diffv);
    }
};
