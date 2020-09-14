#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include "fiesta.hpp"
#include "input.hpp"

struct particleStruct2D {
  int ci, cj;  // particle cell index
  double x, y; // particle position
  int state;   // particle state (1=enabled, 0=disabled)
  double u,v;
  double m;
  double r;
};

typedef typename Kokkos::View<particleStruct2D *>
    FSP2D; // kokkos typename for array or particle structs
typedef typename Kokkos::View<particleStruct2D *>::HostMirror
    FSP2DH; // kokkos typename for array or particle structs

struct findInitialCell2D {
  // functor for doing initial cell search for particles created by position
  FS4D grid;
  FSP2D particles;
  int p;
  int ng;

  findInitialCell2D(FS4D g_, FSP2D part_, int p_, int ng_)
      : grid(g_), particles(part_), p(p_), ng(ng_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, const int j) const {

    if (particles(p).x >= grid(i, j, 0, 0) &&
        particles(p).x < grid(i + 1, j, 0, 0)) {
      if (particles(p).y >= grid(i, j, 0, 1) &&
          particles(p).y < grid(i, j + 1, 0, 1)) {
        // convert grid index to cell index (take into account ghost cells)
        // cells own their left and bottom edges
        particles(p).ci = i + ng;
        particles(p).cj = j + ng;
//        printf("Found Particle in (%d, %d) with coordinates (%.2f, %.2f)\n",particles(p).ci,particles(p).cj,particles(p).x,particles(p).y);
      }
    }
  }
};

struct advectParticles2D {
  // functor for moving particles and finding their new cells and states
  FS4D var;
  FS2D rho;
  FS4D grid;
  FSP2D particles;
  double dt;
  int ng;

  advectParticles2D(FS4D v_, FS2D r_, FS4D g_, FSP2D part_, double dt_, int ng_)
      : var(v_), rho(r_), grid(g_), particles(part_), dt(dt_), ng(ng_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int p) const {
    if (particles(p).state) {
      // get particle local cell indices
      int ci = particles(p).ci;
      int cj = particles(p).cj;
      // convert to grid indices
      int gi = ci - ng;
      int gj = cj - ng;

      int im = grid.extent(0);
      int jm = grid.extent(1);

      // get velocities of current cell and surrounding cells
      double u[3][3];
      double v[3][3];
      int li, lj;
      for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
          li = ci - 1 + i;
          lj = cj - 1 + j;
          u[j][i] = var(li, lj, 0, 0) / rho(li, lj);
          v[j][i] = var(li, lj, 0, 1) / rho(li, lj);
        }
      }

      // average velocities to cell corners (interpolate to cell corners for
      // regular grid)
      double Qu11 = (u[0][0] + u[0][1] + u[1][0] + u[1][1]) / 4.0; // top left
      double Qv11 = (v[0][0] + v[0][1] + v[1][0] + v[1][1]) / 4.0;

      double Qu12 = (u[0][1] + u[0][2] + u[1][1] + u[1][2]) / 4.0; // top right
      double Qv12 = (v[0][1] + v[0][2] + v[1][1] + v[1][2]) / 4.0;

      double Qu21 =
          (u[1][0] + u[1][1] + u[2][0] + u[2][1]) / 4.0; // bottom left
      double Qv21 = (v[1][0] + v[1][1] + v[2][0] + v[2][1]) / 4.0;

      double Qu22 =
          (u[1][1] + u[1][2] + u[2][1] + u[2][2]) / 4.0; // bottom right
      double Qv22 = (v[1][1] + v[1][2] + v[2][1] + v[2][2]) / 4.0;

      // get cell grid points
      double x1 = grid(gi, gj, 0, 0);
      double x2 = grid(gi + 1, gj, 0, 0);
      double y1 = grid(gi, gj, 0, 1);
      double y2 = grid(gi, gj + 1, 0, 1);

      // get particle current location
      double x = particles(p).x;
      double y = particles(p).y;

      // interpolate to x on top and bottom edges
      double Ru1 = ((x2 - x) / (x2 - x1)) * Qu11 +
                   ((x - x1) / (x2 - x1)) * Qu12; // top edge u
      double Rv1 = ((x2 - x) / (x2 - x1)) * Qv11 +
                   ((x - x1) / (x2 - x1)) * Qv12; // top edge v
      double Rv2 = ((x2 - x) / (x2 - x1)) * Qv21 +
                   ((x - x1) / (x2 - x1)) * Qv22; // bottom edge u
      double Ru2 = ((x2 - x) / (x2 - x1)) * Qu21 +
                   ((x - x1) / (x2 - x1)) * Qu22; // bottom edge v
      // interpolate to y
      double lu = ((y2 - y) / (y2 - y1)) * Ru2 + ((y - y1) / (y2 - y1)) * Ru1;
      double lv = ((y2 - y) / (y2 - y1)) * Rv2 + ((y - y1) / (y2 - y1)) * Rv1;

      // calculate force F=6*pi*mu*r*v
      double Fx = 0.000552 * particles(p).r*(lu-particles(p).u);
      double Fy = 0.000552 * particles(p).r*(lv-particles(p).v);
      //double Fx = 0.000552 * 0.002 *(lu-particles(p).u);
      //double Fy = 0.000552 * 0.002 *(lv-particles(p).v);

      // calculate acceleration
      double ax = Fx/particles(p).m;
      double ay = Fy/particles(p).m;
      //double ax = Fx/0.00000001;
      //double ay = Fy/0.00000001;

      // calculate change in velocity
      particles(p).u += ax*dt;
      particles(p).v += ay*dt;

      // calculate change in position
      double dx = particles(p).u*dt + 0.5*ax*dt*dt;
      double dy = particles(p).v*dt + 0.5*ay*dt*dt;
 
      //double dx = lu*dt;
      //double dy = lv*dt;

      // deactivate particle if it will leave domain
      if (particles(p).x + dx < grid(0, gj, 0, 0))
        particles(p).state = 0;
      if (particles(p).y + dy < grid(gi, 0, 0, 1))
        particles(p).state = 0;
      if (particles(p).x + dx > grid(im - 1, gj, 0, 0))
        particles(p).state = 0;
      if (particles(p).y + dy > grid(gi, jm - 1, 0, 1))
        particles(p).state = 0;
      //            printf("%d, %d, %d, %d, %d\n",p,im,jm,gi,gj);

      // move particle if active
      if (particles(p).state) {
        particles(p).x += dx;
        particles(p).y += dy;
      }

      // check if particle moved to new cell
      if (particles(p).x < grid(gi, gj, 0, 0)) { // check if exited left
        particles(p).ci -= 1;
      } else if (particles(p).x >=
                 grid(gi + 1, gj, 0, 0)) { // check if exited right
        particles(p).ci += 1;
      }
      if (particles(p).y < grid(gi, gj, 0, 1)) { // check if exited bottom
        particles(p).cj -= 1;
      } else if (particles(p).y >=
                 grid(gi, gj + 1, 0, 1)) { // check if exited top
        particles(p).cj += 1;
      }
      // printf("%d, %d, %d, %f, %f, %d\n",p,ci,cj,x,y,particles(p).state);
    }
  }
};

// typedef typename Kokkos::View<particleStruct2D*>::HostMirror FSP2DH; //
// kokkos typename for array or particle structs inline void
// writeParticles(struct inputConfig cf,
// Kokkos::View<particleStruct2D*>::HostMirror particlesH){ inline void
// writeParticles(struct inputConfig cf,
// Kokkos::View<particleStruct2D*>::HostMirror pH){
//    using namespace std;
//    stringstream ss;
//    ss << "particle-" << setw(7) << setfill('0') << cf.t << ".vtk";
//    ofstream f;
//    //f.open("particle.vtk");
//    f.open(ss.str());
//    f << "# vtk DataFile Version 4.2" << endl;
//    f << "Test Particles" << endl;
//    f << "ASCII" << endl;
//    f << "DATASET POLYDATA" << endl;
//    f << "POINTS " << cf.p_np << " float" << endl;
//    for (int p=0; p<cf.p_np; ++p){
//        f << pH(p).x << " " << pH(p).y << " " << "0.0" << endl;
//    }
//    f << "VERTICES " << cf.p_np << " " << cf.p_np*2 <<endl;
//    for (int p=0; p<cf.p_np; ++p){
//        f << "1 " << p << endl;
//    }
//    f << "POINT_DATA " << cf.p_np << endl;
//    f << "SCALARS state float" << endl;
//    f << "LOOKUP_TABLE default" << endl;
//    for (int p=0; p<cf.p_np; ++p){
//        f << pH(p).state << endl;
//    }
//    f << "SCALARS ci float" << endl;
//    f << "LOOKUP_TABLE default" << endl;
//    for (int p=0; p<cf.p_np; ++p){
//        f << pH(p).ci << endl;
//    }
//    f << "SCALARS cj float" << endl;
//    f << "LOOKUP_TABLE default" << endl;
//    for (int p=0; p<cf.p_np; ++p){
//        f << pH(p).cj << endl;
//    }
//    f.flush();
//    f.close();
//}
#endif
