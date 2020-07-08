#ifndef bouyancy_hpp
#define bouyancy_hpp
struct computeBuoyancy {
    FS4D dvar;
    FS4D var;
    FS2D rho;

    computeBuoyancy (FS4D dvar_, FS4D var_, FS2D rho_)
         : dvar(dvar_), var(var_), rho(rho_) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const int i, const int j) const {

        double v = var(i,j,0,1)/rho(i,j);
        double rhop = rho(i,j) - 1.0;
        double eps = 1e-6;

        if (rhop >= eps || rhop <= -eps){
            double f = -9.81*rhop;;

            dvar(i,j,0,1) += f;
            dvar(i,j,0,2) += v*f;

//            printf("%f\n",f);
        }
    }
};
#endif
