struct computeFluxWeno2D {
    FS4D var;
    FS2D p;
    FS2D rho;
    FS2D wenox;
    FS2D wenoy;
    Kokkos::View<double*> cd;
    int v;
    double eps = 0.000001;

    computeFluxWeno2D (FS4D var_, FS2D p_, FS2D rho_,
                         FS2D wenox_, FS2D wenoy_, Kokkos::View<double*> cd_, int v_)
                         : var(var_), p(p_), rho(rho_),
                           wenox(wenox_), wenoy(wenoy_), cd(cd_), v(v_) {}
    
    KOKKOS_INLINE_FUNCTION
    void operator()(const int i, const int j) const {

        int ns = (int)cd(0);
        double ur,vr,w,b1,b2,b3,w1,w2,w3,p1,p2,p3,f1,f2,f3,f4,f5;
        double dx = cd(1);
        double dy = cd(2);

        //calculate cell face velocities (in positive direction) with 4th order interpolation
        //velocity is momentum divided by total density
        ur = ( -     var(i+2,j,0,0)/rho(i+2,j) + 7.0*var(i+1,j,0,0)/rho(i+1,j)
               + 7.0*var(i  ,j,0,0)/rho(i  ,j) -     var(i-1,j,0,0)/rho(i-1,j) )/12.0;

        vr = ( -     var(i,j+2,0,1)/rho(i,j+2) + 7.0*var(i,j+1,0,1)/rho(i,j+1)
               + 7.0*var(i,j  ,0,1)/rho(i,j  ) -     var(i,j-1,0,1)/rho(i,j-1) )/12.0;

        //for each direction
        for (int idx=0; idx<2; ++idx){
            //get stencil data.  the flux for the energy equation includes pressure so only add 
            //pressure for the energy variable (index 2 for 2d problem)
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
                wenox(i,j) = ur*w/dx;
            }
            if (idx == 1) {
                wenoy(i,j) = vr*w/dy;
            }
        }
    }
};

struct computeFluxCentered2D {
    
    FS4D var;
    FS2D p;
    FS2D rho;
    FS2D fluxx;
    FS2D fluxy;
    Kokkos::View<double*> cd;
    int v;

    computeFluxCentered2D (FS4D var_, FS2D p_, FS2D rho_, FS2D fx_, FS2D fy_, Kokkos::View<double*> cd_, int v_)
                         : var(var_), p(p_), rho(rho_), fluxx(fx_), fluxy(fy_), cd(cd_), v(v_) {}
    
    KOKKOS_INLINE_FUNCTION
    void operator()(const int i, const int j) const {

        double dx = cd(1);
        double dy = cd(2);
        double ur,vr,x1,y1;
        double px,py;
        ur = 0.0;
        vr = 0.0;
        x1 = 0.0;
        y1 = 0.0;

        ur = ( 0.0 - var(i+2,j,0,0)/rho(i+2,j) + 7.0*var(i+1,j,0,0)/rho(i+1,j)
               + 7.0*var(i  ,j,0,0)/rho(i  ,j) - var(i-1,j,0,0)/rho(i-1,j) )/12.0;

        vr = ( 0.0 - var(i,j+2,0,1)/rho(i,j+2) + 7.0*var(i,j+1,0,1)/rho(i,j+1)
               + 7.0*var(i,j  ,0,1)/rho(i,j  ) - var(i,j-1,0,1)/rho(i,j-1) )/12.0;

        px =  (0.0 - p(i+2,j) + 7.0*p(i+1,j) + 7.0*p(i  ,j) - p(i-1,j) )/12.0;

        py =  (0.0 - p(i,j+2) + 7.0*p(i,j+1) + 7.0*p(i  ,j) - p(i,j-1) )/12.0;

        x1 = ( 0.0 - var(i+2,j,0,v) + 7.0*var(i+1,j,0,v)
               + 7.0*var(i  ,j,0,v) -     var(i-1,j,0,v) )/12.0;

        y1 = (0.0 - var(i,j+2,0,v) + 7.0*var(i,j+1,0,v)
               + 7.0*var(i  ,j,0,v) -     var(i,j-1,0,v) )/12.0;

        if (v == 2){
            x1 = x1+px;
            y1 = y1+py;
        }

        fluxx(i,j) = ur*x1/dx;
        fluxy(i,j) = vr*y1/dy;

    }
};
