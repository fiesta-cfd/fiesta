/*
  Copyright 2019-2021 The University of New Mexico

  This file is part of FIESTA.
  
  FIESTA is free software: you can redistribute it and/or modify it under the
  terms of the GNU Lesser General Public License as published by the Free
  Software Foundation, either version 3 of the License, or (at your option) any
  later version.
  
  FIESTA is distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
  A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
  details.
  
  You should have received a copy of the GNU Lesser General Public License
  along with FIESTA.  If not, see <https://www.gnu.org/licenses/>.
*/

#ifndef CEQ_HPP
#define CEQ_HPP
#include <cstdio>
struct maxCvar3D {
  FS4D dvar;
  int v;

  maxCvar3D(FS4D dvar_, int v_)
      : dvar(dvar_), v(v_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, const int j, const int k, double &lmax) const {

    double s = abs(dvar(i, j, k, v));

    if (s > lmax)
      lmax = s;
  }
};

struct maxCvar2D {
  FS4D dvar;
  int v;
  Kokkos::View<double *> cd;

  maxCvar2D(FS4D dvar_, int v_, Kokkos::View<double *> cd_)
      : dvar(dvar_), v(v_), cd(cd_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, const int j, double &lmax) const {
    int ns = (int)cd(0);
    int nv = ns + 3;

    double s = abs(dvar(i, j, 0, nv + v));

    if (s > lmax)
      lmax = s;
  }
};

struct maxGradRho2D {
  FS3D gradRho;
  int v;

  maxGradRho2D(FS3D gradRho_, int v_) : gradRho(gradRho_), v(v_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, const int j, double &lmax) const {

    double s = abs(gradRho(i, j, v));

    if (s > lmax)
      lmax = s;
  }
};

struct maxWaveSpeed2D {
  FS4D var;
  FS2D p;
  FS2D rho;
  Kokkos::View<double *> cd;

  maxWaveSpeed2D(FS4D var_, FS2D p_, FS2D rho_, Kokkos::View<double *> cd_)
      : var(var_), p(p_), rho(rho_), cd(cd_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, const int j, double &lmax) const {
    int ns = (int)cd(0);
    double gamma, gammas, Rs;
    double Cp = 0;
    double Cv = 0;

    double a, s;

    for (int s = 0; s < ns; ++s) {
      gammas = cd(6 + 2 * s);
      Rs = cd(6 + 2 * s + 1);

      Cp = Cp + (var(i,j,0,3+s) / rho(i,j)) * (gammas * Rs / (gammas - 1));
      Cv = Cv + (var(i,j,0,3+s) / rho(i,j)) * (Rs / (gammas - 1));
    }

    gamma = Cp / Cv;

    a = sqrt(gamma * p(i, j) / rho(i, j));
    s = a + sqrt(var(i, j, 0, 0) * var(i, j, 0, 0) +
                 var(i, j, 0, 1) * var(i, j, 0, 1));

    if (s > lmax)
      lmax = s;
  }
};

struct calculateRhoGrad2D {
  FS4D var;
  FS3D vel;
  FS2D rho;
  FS3D gradRho;
  double dx,dy;

  calculateRhoGrad2D(FS4D var_, FS3D vel_, FS2D rho_, FS3D gradRho_, double dx_, double dy_)
      : var(var_), vel(vel_), rho(rho_), gradRho(gradRho_), dx(dx_), dy(dy_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, const int j) const {

    double i1 = 0;
    double i2 = 0;

    double ex1 = var(i-2,j,0,2)/rho(i-2,j) - 0.5*( vel(i-2,j,0)*vel(i-2,j,0) + vel(i-2,j,1)*vel(i-2,j,1) );
    double ex2 = var(i-1,j,0,2)/rho(i-1,j) - 0.5*( vel(i-1,j,0)*vel(i-1,j,0) + vel(i-1,j,1)*vel(i-1,j,1) );
    double ex3 = var(i+1,j,0,2)/rho(i+1,j) - 0.5*( vel(i+1,j,0)*vel(i+1,j,0) + vel(i+1,j,1)*vel(i+1,j,1) );
    double ex4 = var(i+2,j,0,2)/rho(i+2,j) - 0.5*( vel(i+2,j,0)*vel(i+2,j,0) + vel(i+2,j,1)*vel(i+2,j,1) );

    double ey1 = var(i,j-2,0,2)/rho(i,j-2) - 0.5*( vel(i,j-2,0)*vel(i,j-2,0) + vel(i,j-2,1)*vel(i,j-2,1) );
    double ey2 = var(i,j-1,0,2)/rho(i,j-1) - 0.5*( vel(i,j-1,0)*vel(i,j-1,0) + vel(i,j-1,1)*vel(i,j-1,1) );
    double ey3 = var(i,j+1,0,2)/rho(i,j+1) - 0.5*( vel(i,j+1,0)*vel(i,j+1,0) + vel(i,j+1,1)*vel(i,j+1,1) );
    double ey4 = var(i,j+2,0,2)/rho(i,j+2) - 0.5*( vel(i,j+2,0)*vel(i,j+2,0) + vel(i,j+2,1)*vel(i,j+2,1) );

    double dxr = (rho(i-2,j) - 8.0*rho(i-1,j) + 8.0*rho(i+1,j) - rho(i+2,j)) / (12.0*dx);
    double dyr = (rho(i,j-2) - 8.0*rho(i,j-1) + 8.0*rho(i,j+1) - rho(i,j+2)) / (12.0*dy);

    double dxe = (ex1 - 8.0 * ex2 + 8.0 * ex3 - ex4) / (12.0 * dx);
    double dye = (ey1 - 8.0 * ey2 + 8.0 * ey3 - ey4) / (12.0 * dy);

    double dxu = (vel(i-2,j,0) - 8.0*vel(i-1,j,0) + 8.0*vel(i+1,j,0) - vel(i+2,j,0)) / (12.0*dx);
    double dyv = (vel(i,j-2,1) - 8.0*vel(i,j-1,1) + 8.0*vel(i,j+1,1) - vel(i,j+2,1)) / (12.0*dy);

    double n1 = dxr;
    double n2 = dyr;

    double rgrad = sqrt(dxr * dxr + dyr * dyr);
    double divu = dxu + dyv;

    double dnednr = (n1 * dxe + n2 * dye) * (n1 * dxr + n2 * dyr);

    // compression switch
    if (dnednr < 0.0)
      i1 = 1.0;
    else
      i1 = 0.0;

    // detect shock front (equation 5a)
    if (divu < 0.0) {
      i2 = 1.0;
    }

    gradRho(i, j, 0) = (1 - i1) * i2 * rgrad;
    gradRho(i, j, 1) = i1 * rgrad;
    gradRho(i, j, 2) = -i1 * dyr;
    gradRho(i, j, 3) = i1 * dxr;
  }
};

struct updateCeq2D {
  FS4D dvar;
  FS4D var;
  FS3D gradRho;
  double maxS, kap, eps;
  Kokkos::View<double *> cd;

  updateCeq2D(FS4D dvar_, FS4D var_, FS3D gradRho_, double maxS_,
              Kokkos::View<double *> cd_, double kap_, double eps_)
      : dvar(dvar_), var(var_), gradRho(gradRho_), maxS(maxS_), cd(cd_),
        kap(kap_), eps(eps_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, const int j) const {
    double dx = cd(1);
    double dy = cd(2);

    // calculate cequation variable indices based on number of species (c
    // variables come after species densities)
    int nv = (int)cd(0) + 3;
    int nc = nv;

    double dx_right, dx_left, dy_top, dy_bot;
    double lap;
    //double normalizer=1.0;

    // average cell size
    double dxmag = sqrt(dx * dx + dy * dy);

    for (int n = 0; n < 4; ++n) {
      dx_right = (var(i + 1, j, 0, nc + n) - var(i, j, 0, nc + n)) / dx;
      dx_left = (var(i, j, 0, nc + n) - var(i - 1, j, 0, nc + n)) / dx;

      dy_top = (var(i, j + 1, 0, nc + n) - var(i, j, 0, nc + n)) / dy;
      dy_bot = (var(i, j, 0, nc + n) - var(i, j - 1, 0, nc + n)) / dy;

      // update ceq right hand side
      lap = dx_right - dx_left + dy_top - dy_bot;

      dvar(i,j,0,nc+n) =
          (maxS / (eps * dxmag)) * (gradRho(i,j,n) - var(i,j,0,nc+n)) + kap * maxS * dxmag * lap;
    }
  }
};

struct computeCeqFlux2D {
  FS4D var;
  FS5D m; // m(face,direction of derivative, i, j, velocity component)
  FS2D r;
  FS3D v;
  double alpha;
  int nv;

  computeCeqFlux2D(FS4D var_, FS5D m_, FS3D v_, FS2D r_, double a_, int nv_)
      : var(var_), m(m_), v(v_), r(r_), alpha(a_), nv(nv_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, const int j) const {
    int ih; // ihat - i-direction unit vector
    int jh; // jhat - j-direction unit vector

    // Compute alph*rho*C*Ctau1*Ctau2 on each face for each direction
    for (int w=0;w<2;++w){
      for (int f=0;f<2;++f){
        for (int d=0;d<2;++d){
          // set direction unit vector
          ih=0; jh=0;
          if (f==0) ih=1;
          if (f==1) jh=1;

          m(f,d,i,j,w) =(r(i,j)+r(i+ih,j+jh))/2.0;
          m(f,d,i,j,w)+=(var(i,j,0,nv+1)+var(i+ih,j+jh,0,nv+1))/2.0;
          m(f,d,i,j,w)+=(var(i,j,0,nv+2)+var(i+ih,j+jh,0,nv+2))/2.0;
          m(f,d,i,j,w)+=(var(i,j,0,nv+3)+var(i+ih,j+jh,0,nv+3))/2.0;
          m(f,d,i,j,w)*=alpha;
        }
      }
    }
  }
};

struct computeCeqFaces2D {
  FS5D m; // m(face,direction of derivative, i, j, velocity component)
  FS3D v;
  FS1D cd;

  computeCeqFaces2D(FS5D m_, FS3D v_, FS1D cd_) : m(m_), v(v_), cd(cd_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, const int j) const {
    double dx = cd(1);
    double dy = cd(2);

    // Compute dux and duy on each face for each velocity
    for (int w=0;w<2;++w){
      m(0,0,i,j,w) *= ( v(i+1,j,w)-v(i,j,w) )/dx;
      m(0,1,i,j,w) *= ( (v(i,j+1,w)+v(i+1,j+1,w))-(v(i,j-1,w)+v(i+1,j-1,w)) )/(4.0*dx);
      m(1,1,i,j,w) *= ( v(i,j+1,w)-v(i,j,w) )/dy;
      m(1,0,i,j,w) *= ( (v(i+1,j,w)+v(i+1,j+1,w))-(v(i-1,j,w)+v(i-1,j+1,w)) )/(4.0*dx);
    }
  }
};

struct applyCeq2D {
  FS5D m;
  FS4D dvar,varx;
  double dx,dy;

  applyCeq2D(FS4D dvar_, FS4D varx_, FS5D m_, double dx_, double dy_)
    : dvar(dvar_), varx(varx_), m(m_), dx(dx_), dy(dy_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, const int j) const {
    double diffu;

    // assembly artificial viscosity term
    for (int w=0;w<2;++w){
      diffu=0;
      diffu += (m(0,0,i+1,j,w)-m(0,0,i,j,w))/dx;
      diffu += (m(0,1,i+1,j,w)-m(0,1,i,j,w))/dx;
      diffu += (m(1,0,i,j+1,w)-m(1,0,i,j,w))/dy;
      diffu += (m(1,1,i,j+1,w)-m(1,1,i,j,w))/dy;

      dvar(i,j,0,w)+=diffu;
      varx(i,j,0,5+w)=diffu;
    }
  }
};

struct maxWaveSpeed {
  FS4D var;
  FS3D p;
  FS3D rho;
  Kokkos::View<double *> cd;

  maxWaveSpeed(FS4D var_, FS3D p_, FS3D rho_, Kokkos::View<double *> cd_)
      : var(var_), p(p_), rho(rho_), cd(cd_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, const int j, const int k, double &lmax) const {
    int ns = (int)cd(0);
    double gamma, gammas, Rs;
    double Cp = 0;
    double Cv = 0;

    double a, s;

    for (int s = 0; s < ns; ++s) {
      gammas = cd(6 + 3 * s);
      Rs = cd(6 + 3 * s + 1);

      Cp = Cp +
           (var(i, j, k, 4 + s) / rho(i, j, k)) * (gammas * Rs / (gammas - 1));
      Cv = Cv + (var(i, j, k, 4 + s) / rho(i, j, k)) * (Rs / (gammas - 1));
    }

    gamma = Cp / Cv;

    a = sqrt(gamma * p(i, j, k) / rho(i, j, k));
    s = a + sqrt(var(i, j, k, 0) * var(i, j, k, 0) +
                 var(i, j, k, 1) * var(i, j, k, 1) +
                 var(i, j, k, 2) * var(i, j, k, 2)) /
                rho(i, j, k);

    if (s > lmax)
      lmax = s;
  }
};

struct calculateRhoGrad {
  FS4D var;
  FS3D rho;
  FS4D gradRho;
  Kokkos::View<double *> cd;

  calculateRhoGrad(FS4D var_, FS3D rho_, FS4D gradRho_,
                   Kokkos::View<double *> cd_)
      : var(var_), rho(rho_), gradRho(gradRho_), cd(cd_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, const int j, const int k) const {

    int indicator = 0;

    double u1 = var(i - 2, j, k, 0) / rho(i - 2, j, k);
    double u2 = var(i - 1, j, k, 0) / rho(i - 1, j, k);
    double u3 = var(i + 1, j, k, 0) / rho(i + 1, j, k);
    double u4 = var(i + 2, j, k, 0) / rho(i + 2, j, k);

    double v1 = var(i, j - 2, k, 0) / rho(i, j - 2, k);
    double v2 = var(i, j - 1, k, 0) / rho(i, j - 1, k);
    double v3 = var(i, j + 1, k, 0) / rho(i, j + 1, k);
    double v4 = var(i, j + 2, k, 0) / rho(i, j + 2, k);

    double w1 = var(i, j, k - 2, 0) / rho(i, j, k - 2);
    double w2 = var(i, j, k - 1, 0) / rho(i, j, k - 1);
    double w3 = var(i, j, k + 1, 0) / rho(i, j, k + 1);
    double w4 = var(i, j, k + 2, 0) / rho(i, j, k + 2);

    double ex1 = var(i - 2, j, k, 3) / rho(i, j, k - 2) -
                 0.5 * (1 / (rho(i - 2, j, k) * rho(i - 2, j, k))) *
                     (var(i - 2, j, k, 0) * var(i - 2, j, k, 0) +
                      var(i - 2, j, k, 1) * var(i - 2, j, k, 1) +
                      var(i - 2, j, k, 2) * var(i - 2, j, k, 2));
    double ex2 = var(i - 1, j, k, 3) / rho(i, j, k - 1) -
                 0.5 * (1 / (rho(i - 1, j, k) * rho(i - 1, j, k))) *
                     (var(i - 1, j, k, 0) * var(i - 1, j, k, 0) +
                      var(i - 1, j, k, 1) * var(i - 1, j, k, 1) +
                      var(i - 1, j, k, 2) * var(i - 1, j, k, 2));
    double ex3 = var(i + 1, j, k, 3) / rho(i, j, k + 1) -
                 0.5 * (1 / (rho(i + 1, j, k) * rho(i + 1, j, k))) *
                     (var(i + 1, j, k, 0) * var(i + 1, j, k, 0) +
                      var(i + 1, j, k, 1) * var(i + 1, j, k, 1) +
                      var(i + 1, j, k, 2) * var(i + 1, j, k, 2));
    double ex4 = var(i + 2, j, k, 3) / rho(i, j, k + 2) -
                 0.5 * (1 / (rho(i + 2, j, k) * rho(i + 2, j, k))) *
                     (var(i + 2, j, k, 0) * var(i + 2, j, k, 0) +
                      var(i + 2, j, k, 1) * var(i + 2, j, k, 1) +
                      var(i + 2, j, k, 2) * var(i + 2, j, k, 2));

    double ey1 = var(i, j - 2, k, 3) / rho(i, j, k - 2) -
                 0.5 * (1 / (rho(i, j - 2, k) * rho(i, j - 2, k))) *
                     (var(i, j - 2, k, 0) * var(i, j - 2, k, 0) +
                      var(i, j - 2, k, 1) * var(i, j - 2, k, 1) +
                      var(i, j - 2, k, 2) * var(i, j - 2, k, 2));
    double ey2 = var(i, j - 1, k, 3) / rho(i, j, k - 1) -
                 0.5 * (1 / (rho(i, j - 1, k) * rho(i, j - 1, k))) *
                     (var(i, j - 1, k, 0) * var(i, j - 1, k, 0) +
                      var(i, j - 1, k, 1) * var(i, j - 1, k, 1) +
                      var(i, j - 1, k, 2) * var(i, j - 1, k, 2));
    double ey3 = var(i, j + 1, k, 3) / rho(i, j, k + 1) -
                 0.5 * (1 / (rho(i, j + 1, k) * rho(i, j + 1, k))) *
                     (var(i, j + 1, k, 0) * var(i, j + 1, k, 0) +
                      var(i, j + 1, k, 1) * var(i, j + 1, k, 1) +
                      var(i, j + 1, k, 2) * var(i, j + 1, k, 2));
    double ey4 = var(i, j + 2, k, 3) / rho(i, j, k + 2) -
                 0.5 * (1 / (rho(i, j + 2, k) * rho(i, j + 2, k))) *
                     (var(i, j + 2, k, 0) * var(i, j + 2, k, 0) +
                      var(i, j + 2, k, 1) * var(i, j + 2, k, 1) +
                      var(i, j + 2, k, 2) * var(i, j + 2, k, 2));

    double ez1 = var(i, j, k - 2, 3) / rho(i, j, k - 2) -
                 0.5 * (1 / (rho(i, j, k - 2) * rho(i, j, k - 2))) *
                     (var(i, j, k - 2, 0) * var(i, j, k - 2, 0) +
                      var(i, j, k - 2, 1) * var(i, j, k - 2, 1) +
                      var(i, j, k - 2, 2) * var(i, j, k - 2, 2));
    double ez2 = var(i, j, k - 1, 3) / rho(i, j, k - 1) -
                 0.5 * (1 / (rho(i, j, k - 1) * rho(i, j, k - 1))) *
                     (var(i, j, k - 1, 0) * var(i, j, k - 1, 0) +
                      var(i, j, k - 1, 1) * var(i, j, k - 1, 1) +
                      var(i, j, k - 1, 2) * var(i, j, k - 1, 2));
    double ez3 = var(i, j, k + 1, 3) / rho(i, j, k + 1) -
                 0.5 * (1 / (rho(i, j, k + 1) * rho(i, j, k + 1))) *
                     (var(i, j, k + 1, 0) * var(i, j, k + 1, 0) +
                      var(i, j, k + 1, 1) * var(i, j, k + 1, 1) +
                      var(i, j, k + 1, 2) * var(i, j, k + 1, 2));
    double ez4 = var(i, j, k + 2, 3) / rho(i, j, k + 2) -
                 0.5 * (1 / (rho(i, j, k + 2) * rho(i, j, k + 2))) *
                     (var(i, j, k + 2, 0) * var(i, j, k + 2, 0) +
                      var(i, j, k + 2, 1) * var(i, j, k + 2, 1) +
                      var(i, j, k + 2, 2) * var(i, j, k + 2, 2));

    double dxr = (rho(i-2,j,k) - 8.0*rho(i-1,j,k) + 8.0*rho(i+1,j,k) - rho(i+2,j,k)) / (12.0*cd(1));
    double dyr = (rho(i,j-2,k) - 8.0*rho(i,j-1,k) + 8.0*rho(i,j+1,k) - rho(i,j+2,k)) / (12.0*cd(2));
    double dzr = (rho(i,j,k-2) - 8.0*rho(i,j,k-1) + 8.0*rho(i,j,k+1) - rho(i,j,k+2)) / (12.0*cd(3));

    double dxe = (ex1 - 8.0 * ex2 + 8.0 * ex3 - ex4) / (12.0 * cd(1));
    double dye = (ey1 - 8.0 * ey2 + 8.0 * ey3 - ey4) / (12.0 * cd(2));
    double dze = (ez1 - 8.0 * ez2 + 8.0 * ez3 - ez4) / (12.0 * cd(3));

    double dxu = (u1 - 8.0 * u2 + 8.0 * u3 - u4) / (12.0 * cd(1));
    double dyv = (v1 - 8.0 * v2 + 8.0 * v3 - v4) / (12.0 * cd(2));
    double dzw = (w1 - 8.0 * w2 + 8.0 * w3 - w4) / (12.0 * cd(3));

    double n1 = dxr;
    double n2 = dyr;
    double n3 = dzr;

    double rgrad = sqrt(dxr * dxr + dyr * dyr + dzr * dzr);
    double divu = dxu + dyv + dzw;

    double dnednr =
        (n1 * dxe + n2 * dye + n3 * dze) * (n1 * dxr + n2 * dyr + n3 * dzr);

    // compression switch
    if (dnednr <= 0)
      indicator = 1;
    else
      indicator = 0;

    // detect shock front (equation 5a)
    if (divu <= 0)
      gradRho(i, j, k, 0) = (1 - indicator) * rgrad;
    // gradRho(i,j,k,0) = (1-indicator)*divu*rgrad;
    else
      gradRho(i, j, k, 0) = 0;

    // detect contact surface
    gradRho(i, j, k, 1) = indicator * rgrad;

    // gradient components
    gradRho(i, j, k, 2) = dxr;
    gradRho(i, j, k, 3) = dyr;
    gradRho(i, j, k, 4) = dzr;
  }
};

struct maxGradFunctor {
  FS4D gradRho;
  int n;

  maxGradFunctor(FS4D gradRho_, int n_) : gradRho(gradRho_), n(n_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, const int j, const int k, double &lmax) const {
    if (abs(gradRho(i, j, k, n)) > lmax)
      lmax = abs(gradRho(i, j, k, n));
  }
};

struct updateCeq {
  FS4D dvar;
  FS4D var;
  FS4D varx;
  FS4D gradRho;
  double maxS, kap, eps;
  Kokkos::View<double *> cd;

  updateCeq(FS4D dvar_, FS4D var_, FS4D varx_, FS4D gradRho_, double maxS_,
            Kokkos::View<double *> cd_, double kap_, double eps_)
      : dvar(dvar_), var(var_), varx(varx_), gradRho(gradRho_), maxS(maxS_), cd(cd_),
        kap(kap_), eps(eps_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, const int j, const int k) const {
    double dx = cd(1);
    double dy = cd(2);
    double dz = cd(3);

    int nv = (int)cd(0) + 4;
    int nc = nv;

    double lap;

    // average cell size
    double dxmag = sqrt(dx*dx + dy*dy + dz*dz);

    for (int n = 0; n < 5; ++n) {
      lap = (var(i-1,j,k,nc+n)-2*var(i,j,k,nc+n)+var(i+1,j,k,nc+n))/dx
          + (var(i,j-1,k,nc+n)-2*var(i,j,k,nc+n)+var(i,j+1,k,nc+n))/dy
          + (var(i,j,k-1,nc+n)-2*var(i,j,k,nc+n)+var(i,j,k+1,nc+n))/dz;

      dvar(i,j,k,nc+n) = (maxS/(eps*dxmag))*(gradRho(i,j,k,n)-var(i,j,k,nc+n)) + kap*maxS*dxmag*lap;
    }
  }
};

struct calculateCeqFlux {
  FS4D var;
  FS3D rho;
  FS6D mFlux; //(m,n,i,j,k,dir)
  FS4D cFlux; //(i,j,k,dir)
  Kokkos::View<double *> cd;

  calculateCeqFlux(FS4D var_, FS3D rho_, FS6D mFlux_, FS4D cFlux_,
                   Kokkos::View<double *> cd_)
      : var(var_), rho(rho_), mFlux(mFlux_), cFlux(cFlux_), cd(cd_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, const int j, const int k) const {

    int nv = (int)cd(0) + 4;
    int nc = nv;
    int nch = nv + 1;
    int nc1 = nv + 2;

    int ip = 0;
    int jp = 0;
    int kp = 0;

    double c_left;
    double c_right;
    double ch_left;
    double ch_right;

    double cn_left[3];
    double cn_right[3];
    double cmag_left = 0;
    double cmag_right = 0;

    double m_left;
    double m_right;

    double rho_left;
    double rho_right;

    // for each direction (i=0, j=1, k=2)
    for (int dir = 0; dir < 3; ++dir) {
      ip = 0;
      jp = 0;
      kp = 0;
      if (dir == 0)
        ip = 1;
      if (dir == 1)
        jp = 1;
      if (dir == 2)
        kp = 1;

      // left is current cell, right is cell in positive direction

      // get right and left components of isotropic C
      c_left = var(i, j, k, nc);
      c_right = var(i + ip, j + jp, k + kp, nc);

      // get right and left components of anisotropic C
      ch_left = var(i, j, k, nch);
      ch_right = var(i + ip, j + jp, k + kp, nch);

      rho_left = rho(i, j, k);
      rho_right = rho(i + ip, j + jp, k + kp);

      // get right and left values of each directional C
      for (int idx = 0; idx < 3; ++idx) {
        cn_left[idx] = var(i, j, k, nc1 + idx);
        cn_right[idx] = var(i + ip, j + jp, k + kp, nc1 + idx);
      }

      // calculate magnitude of directional C
      for (int idx = 0; idx < 3; ++idx)
        cmag_left += cn_left[idx] * cn_left[idx];
      // cmag_left  = sqrt(cmag_left);
      for (int idx = 0; idx < 3; ++idx)
        cmag_right += cn_right[idx] * cn_right[idx];
      // cmag_right = sqrt(cmag_right);

      if (cmag_left <= 0.000001 || cmag_right <= 0.000001) {
        for (int m = 0; m < 3; ++m) {
          for (int n = 0; n < 3; ++n) {
            mFlux(m, n, i, j, k, dir) = 0.0;
            cFlux(i, j, k, dir) = 0.0;
          }
        }
      } else {
        // tensor components
        for (int m = 0; m < 3; ++m) {
          for (int n = 0; n < 3; ++n) {

            // dirac delta
            double d = 0;
            if (m == n)
              d = 1;

            // calculate right and left tensor components
            m_left = d - (cn_left[m] * cn_left[n] / cmag_left);
            m_right = d - (cn_right[m] * cn_right[n] / cmag_right);
            //m_left = d - (cn_left[m] * cn_left[n]);
            //m_right = d - (cn_right[m] * cn_right[n]);

            // include isotropic c
            m_left = m_left * ch_left;
            m_right = m_right * ch_right;

            // include density
            m_left = m_left * rho_left;
            m_right = m_right * rho_right;

            // find flux
            mFlux(m, n, i, j, k, dir) = (m_right - m_left) / 2.0;
          }
        }

        // calcualte isotropic C flux
        cFlux(i, j, k, dir) = (c_right * rho_right - c_left * rho_left) / 2.0;

      }
    }
  }
};

struct applyCeq {
  FS4D dvar;
  FS4D var;
  FS4D varx;
  FS4D vel;
  FS3D rho;
  FS6D mFlux; //(m,n,i,j,k,dir)
  FS4D cFlux; //(i,j,k,dir)
  double alpha, beta, betae;
  Kokkos::View<double *> cd;

  applyCeq(FS4D dvar_, FS4D var_, FS4D varx_, FS4D vel_, FS3D rho_, FS6D mFlux_, FS4D cFlux_,
           Kokkos::View<double *> cd_, double alpha_, double beta_, double betae_)
      : dvar(dvar_), var(var_), varx(varx_), vel(vel_), rho(rho_), mFlux(mFlux_), cFlux(cFlux_),
        cd(cd_), alpha(alpha_), beta(beta_), betae(betae_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, const int j, const int k) const {

    int nv = (int)cd(0) + 4;
    double diffu;

    double du_right[3][3]; // face, direction
    double du_left[3][3];

    double dx[3];

    dx[0] = cd(1);
    dx[1] = cd(2);
    dx[2] = cd(3);

    double an; // anisotropic part
    double is; // isotropic part

    int ip, jp, kp;

    // for each velocity component and energy
    for (int n = 0; n < 3; ++n) {
      // left face
      du_left[0][0] = (vel(i,j,k,n) - vel(i-1,j,k,n)) / cd(1);
      du_left[0][1] = ( (vel(i-1,j+1,k,n) + vel(i,j+1,k,n))
                       -(vel(i-1,j-1,k,n) + vel(i,j-1,k,n)) ) / (4*cd(2));
      du_left[0][2] = ( (vel(i-1,j,k+1,n) + vel(i,j,k+1,n))
                       -(vel(i-1,j,k-1,n) + vel(i,j,k-1,n)) ) / (4*cd(3));

      // bottom face
      du_left[1][0] = ( (vel(i+1,j,k,n) + vel(i+1,j-1,k,n))
                       -(vel(i-1,j,k,n) + vel(i-1,j-1,k,n)) ) / (4*cd(1));
      du_left[1][1] = (vel(i,j,k,n) - vel(i,j-1,k,n)) / cd(2);
      du_left[1][2] = ( (vel(i,j,k+1,n) + vel(i,j-1,k+1,n))
                       -(vel(i,j,k-1,n) + vel(i,j-1,k-1,n)) ) / (4*cd(3));

      // back (hind) face
      du_left[2][0] = ( (vel(i+1,j,k,n) + vel(i+1,j,k-1,n))
                       -(vel(i-1,j,k,n) + vel(i-1,j,k-1,n)) ) / (4*cd(1));
      du_left[2][1] = ( (vel(i,j+1,k,n) + vel(i,j+1,k-1,n))
                       -(vel(i,j-1,k,n) + vel(i,j-1,k-1,n)) ) / (4*cd(2));
      du_left[2][2] = (vel(i,j,k,n) - vel(i,j,k-1,n)) / cd(3);

      // right face
      du_right[0][0] = (vel(i+1,j,k,n) - vel(i,j,k,n)) / cd(1);
      du_right[0][1] = ( (vel(i,j+1,k,n) + vel(i+1,j+1,k,n))
                       -(vel(i,j-1,k,n) + vel(i+1,j-1,k,n)) ) / (4*cd(2));
      du_right[0][2] = ( (vel(i,j,k+1,n) + vel(i+1,j,k+1,n))
                       -(vel(i,j,k-1,n) + vel(i+1,j,k-1,n)) ) / (4*cd(3));

      // top face
      du_right[1][0] = ( (vel(i+1,j+1,k,n) + vel(i+1,j,k,n))
                       -(vel(i-1,j+1,k,n) + vel(i-1,j,k,n)) ) / (4*cd(1));
      du_right[1][1] = (vel(i,j,k,n) - vel(i,j-1,k,n)) / cd(2);
      du_right[1][2] = ( (vel(i,j+1,k+1,n) + vel(i,j,k+1,n))
                       -(vel(i,j+1,k-1,n) + vel(i,j,k-1,n)) ) / (4*cd(3));

      // front face
      du_right[2][0] = ( (vel(i+1,j,k+1,n) + vel(i+1,j,k,n))
                       -(vel(i-1,j,k+1,n) + vel(i-1,j,k,n)) ) / (4*cd(1));
      du_right[2][1] = ( (vel(i,j+1,k+1,n) + vel(i,j+1,k,n))
                       -(vel(i,j-1,k+1,n) + vel(i,j-1,k,n)) ) / (4*cd(2));
      du_right[2][2] = (vel(i,j,k,n) - vel(i,j,k-1,n)) / cd(3);

      an = 0;
      is = 0;
      ip = 0;
      jp = 0;
      kp = 0;

      //if (n < 3) {
        diffu = 0;
        for (int d = 0; d < 3; ++d) { // each direction for divergence
          ip = 0; jp = 0; kp = 0;
          if (d == 0) ip = 1;
          if (d == 1) jp = 1;
          if (d == 2) kp = 1;

          for (int f = 0; f < 3; ++f) { // each direction for gradient
            an = an + (mFlux(d,f,i,j,k,d) * du_right[d][f] -
                       mFlux(d,f,i-ip,j-jp,k-kp,d) * du_left[d][f]) / dx[d];
            is = is + (cFlux(i,j,k,d) * du_right[d][f] -
                       cFlux(i-ip,j-jp,k-kp,d) * du_left[d][f]) / dx[d];
          }

        }

        diffu = alpha*an + beta*is;

        dvar(i,j,k,n) += diffu;
        //if (n==2) varx(i,j,k,9) = dvar(i,j,k,n);
        //if (n==1) varx(i,j,k,9) = 1.0;
        if (n==1) varx(i,j,k,9) = dvar(i,j,k,n);
        //if (an!=0) printf("(%d,%d,%d) %d: %f  %f  %f %f\n",i,j,k,n,diffu,alpha,an,is);
        //if (an!=0) printf("(%d,%d,%d) %d: %f  %f  %f\n",i,j,k,n,diffu,alpha,an);
      //} else {
      //  diffu = 0;
      //  for (int d = 0; d < 3; ++d) { // each direction for divergence
      //    ip = 0; jp = 0; kp = 0;
      //    if (d == 0) ip = 1;
      //    if (d == 1) jp = 1;
      //    if (d == 2) kp = 1;

      //    for (int f = 0; f < 3; ++f) { // each direction for gradient
      //      diffu = diffu + (cFlux(i,j,k,d) * du_right[d][f] -
      //                       cFlux(i-ip,j-jp,k-kp,d) * du_left[d][f]) / dx[d];
      //    }

      //  }
      //  dvar(i,j,k,n) += betae*diffu;
      //}
    }
  }
};
#endif
