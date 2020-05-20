struct detectNoise2D {
    FS4D var;
    FS2D noise;
    int v;
    double dh;
    Kokkos::View<double*> cd;

    detectNoise2D (FS4D var_, FS2D n_, double dh_, Kokkos::View<double*> cd_, int v_)
                         : var(var_), noise(n_), dh(dh_), cd(cd_), v(v_) {}
    
    KOKKOS_INLINE_FUNCTION
    void operator()(const int ii, const int jj) const {

        int ng = (int)cd(5);

        double dx = cd(1);
        double dy = cd(2);
        
        int i = 2*ii+ng;
        int j = 2*jj+ng;

        double a = sqrt(6.0*dx*dy);

        double c = (a/192.0)+(var(i-1,j+1,0,v)+var(i  ,j+1,0,v)+var(i+1,j+1,0,v)
                                          +var(i-1,j  ,0,v)+var(i  ,j  ,0,v)+var(i+1,j  ,0,v)
                                          +var(i-1,j-1,0,v)+var(i  ,j-1,0,v)+var(i+1,j-1,0,v));

        double cref = dh*a/16.0;

        if (abs(c) >= cref){
            noise(i-1,j+1) = 1;
            noise(i  ,j+1) = 1;
            noise(i+1,j+1) = 1;
            noise(i-1,j  ) = 1;
            noise(i  ,j  ) = 1;
            noise(i+1,j  ) = 1;
            noise(i-1,j-1) = 1;
            noise(i  ,j-1) = 1;
            noise(i+1,j-1) = 1;
        }

    }
};
