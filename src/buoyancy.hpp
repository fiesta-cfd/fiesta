#ifndef bouyancy_hpp
#define bouyancy_hpp
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
#endif
