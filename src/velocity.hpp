#ifndef VELOCITY_HPP
#define VELOCITY_HPP

struct computeGenVelocity3D {
    FS4D var;
    FS5D metrics;
    FS3D rho;
    FS4D vel;

    computeGenVelocity3D (FS4D var_, FS5D m_, FS3D r_, FS4D v_)
         : var(var_), metrics(m_), rho(r_), vel(v_) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const int i, const int j, const int k) const {

        vel(i,j,k,0) = ( metrics(i,j,k,0,0)*var(i,j,k,0)/rho(i,j,k)
                       + metrics(i,j,k,0,1)*var(i,j,k,1)/rho(i,j,k)
                       + metrics(i,j,k,0,2)*var(i,j,k,2)/rho(i,j,k) );
        vel(i,j,k,1) = ( metrics(i,j,k,1,0)*var(i,j,k,0)/rho(i,j,k)
                       + metrics(i,j,k,1,1)*var(i,j,k,1)/rho(i,j,k)
                       + metrics(i,j,k,1,2)*var(i,j,k,2)/rho(i,j,k) );
        vel(i,j,k,2) = ( metrics(i,j,k,2,0)*var(i,j,k,0)/rho(i,j,k)
                       + metrics(i,j,k,2,1)*var(i,j,k,1)/rho(i,j,k)
                       + metrics(i,j,k,2,2)*var(i,j,k,2)/rho(i,j,k) );
//        if (k==25)
//            printf("%d, %d, %f, %f, %f\n",i,j,vel(i,j,k,0),vel(i,j,k,1),vel(i,j,k,2));
    }
};
struct computeGenVelocity2D {
    FS4D var;
    FS4D metrics;
    FS2D rho;
    FS3D vel;

    computeGenVelocity2D (FS4D var_, FS4D m_, FS2D r_, FS3D v_)
         : var(var_), metrics(m_), rho(r_), vel(v_) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const int i, const int j) const {

        vel(i,j,0) = ( metrics(i,j,0,0)*var(i,j,0,0)/rho(i,j) + metrics(i,j,0,1)*var(i,j,0,1)/rho(i,j) );
        vel(i,j,1) = ( metrics(i,j,1,0)*var(i,j,0,0)/rho(i,j) + metrics(i,j,1,1)*var(i,j,0,1)/rho(i,j) );
    }
};

struct computeVelocity2D {
    FS4D var;
    FS2D rho;
    FS3D vel;

    computeVelocity2D (FS4D var_, FS2D r_, FS3D v_)
         : var(var_), rho(r_), vel(v_) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const int i, const int j) const {

        vel(i,j,0) = var(i,j,0,0)/rho(i,j);
        vel(i,j,1) = var(i,j,0,1)/rho(i,j);
    }
};
#endif
