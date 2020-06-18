#ifndef VELOCITY_HPP
#define VELOCITY_HPP

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
