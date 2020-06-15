#include "fiesta.hpp"

struct particleStruct {
    int ci,cj;
    double x,y;
    int state;
};

typedef typename Kokkos::View<particleStruct*> FSP2D;

struct findInitialCell2D{
    FS4D grid;
    FSP2D particles;
    int p;

    findInitialCell2D(FS4D g_, FSP2D part_, int p_) : grid(g_), particles(part_), p(p_) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const int i, const int j) const {

        if (particles(p).x >= grid(i,j,0,0) && particles(p).x < grid(i+1,j,0,0)){
            if (particles(p).y >= grid(i,j,0,1) && particles(p).y < grid(i,j+1,0,1)){
                particles(p).ci = i;
                particles(p).cj = j;
            }
        }
    }
};

struct advectParticles2D{
    FS4D var;
    FS2D rho;
    FS4D grid;
    FSP2D particles;
    double dt;
    int imax,jmax;

    advectParticles2D(FS4D v_, FS2D r_, FS4D g_, FSP2D part_, double dt_, int im_, int jm_)
           : var(v_), rho(r_), grid(g_), particles(part_), dt(dt_), imax(im_), jmax(jm_) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const int p) const {
        if (particles(p).state){
            // get particle local cell indices
            int i = particles(p).ci;
            int j = particles(p).cj;

            // compute local cell average velocity
            double u = var(i,j,0,0)/rho(i,j);
            double v = var(i,j,0,1)/rho(i,j);
            

            // calculate position change
            double dx = u*dt;
            double dy = v*dt;

            // deactivate particle if it leaves domain
            if (particles(p).x+dx < grid(0,j,0,0))
                particles(p).state = 0;
            if (particles(p).y+dy < grid(i,0,0,1))
                particles(p).state = 0;
            if (particles(p).x+dx > grid(imax,j,0,0))
                particles(p).state = 0;
            if (particles(p).y+dy > grid(i,jmax,0,1))
                particles(p).state = 0;

            // move particle
            particles(p).x += dx;
            particles(p).y += dy;

            // check if particle moved to new cell
            if (particles(p).x < grid(i,j,0,0)){
                particles(p).ci -= 1;
                if (particles(p).y < grid(i,j,0,1)){
                    particles(p).cj -= 1;
                }else if (particles(p).y >= grid(i,j+1,0,1)){
                    particles(p).cj += 1;
                }
            }else if (particles(p).x >= grid(i+1,j,0,0)){
                particles(p).ci += 1;
                if (particles(p).y < grid(i,j,0,1)){
                    particles(p).cj -= 1;
                }else if (particles(p).y >= grid(i,j+1,0,1)){
                    particles(p).cj += 1;
                }
            }
        }
    }
};
