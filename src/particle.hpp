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
    int ng;

    findInitialCell2D(FS4D g_, FSP2D part_, int p_, int ng_) : grid(g_), particles(part_), p(p_), ng(ng_) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const int i, const int j) const {

        if (particles(p).x >= grid(i,j,0,0) && particles(p).x < grid(i+1,j,0,0)){
            if (particles(p).y >= grid(i,j,0,1) && particles(p).y < grid(i,j+1,0,1)){
                // convert grid index to cell index (take into account ghost cells)
                // cells own their left and bottom edges
                particles(p).ci = i+ng;
                particles(p).cj = j+ng;
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
    int ng;

    advectParticles2D(FS4D v_, FS2D r_, FS4D g_, FSP2D part_, double dt_, int im_, int jm_, int ng_)
           : var(v_), rho(r_), grid(g_), particles(part_), dt(dt_), imax(im_), jmax(jm_), ng(ng_) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const int p) const {
        if (particles(p).state){
            // get particle local cell indices
            int ci = particles(p).ci;
            int cj = particles(p).cj;
            // convert to grid indices
            int gi = ci-ng;
            int gj = cj-ng;

            int im = grid.extent(0);
            int jm = grid.extent(1);

            // compute local cell average velocity
            double u = var(ci,cj,0,0)/rho(ci,cj);
            double v = var(ci,cj,0,1)/rho(ci,cj);
            

            // calculate position change
            double dx = u*dt;
            double dy = v*dt;

            // deactivate particle if it leaves domain
            if (particles(p).x+dx < grid(0,gj,0,0))
                particles(p).state = 0;
            if (particles(p).y+dy < grid(gi,0,0,1))
                particles(p).state = 0;
            if (particles(p).x+dx > grid(im-1,gj,0,0))
                particles(p).state = 0;
            if (particles(p).y+dy > grid(gi,jm-1,0,1))
                particles(p).state = 0;

            // move particle
            particles(p).x += dx;
            particles(p).y += dy;

            // check if particle moved to new cell
            if (particles(p).x < grid(gi,gj,0,0)){
                particles(p).ci -= 1;
                if (particles(p).y < grid(gi,gj,0,1)){
                    particles(p).cj -= 1;
                }else if (particles(p).y >= grid(gi,gj+1,0,1)){
                    particles(p).cj += 1;
                }
            }else if (particles(p).x >= grid(gi+1,gj,0,0)){
                particles(p).ci += 1;
                if (particles(p).y < grid(gi,gj,0,1)){
                    particles(p).cj -= 1;
                }else if (particles(p).y >= grid(gi,gj+1,0,1)){
                    particles(p).cj += 1;
                }
            }
        }
    }
};
