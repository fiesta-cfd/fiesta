
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

#ifndef CEQ3D_HPP
#define CEQ3D_HPP
#include <cstdio>

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

struct calculateRhoGrad {
  FS4D var,vel;
  FS3D rho;
  FS4D gradRho;
  double dx,dy,dz;

  calculateRhoGrad(FS4D var_, FS4D vel_, FS3D rho_, FS4D gradRho_,
                   double dx_, double dy_, double dz_)
      : var(var_), vel(vel_), rho(rho_), gradRho(gradRho_), dx(dx_), dy(dy_), dz(dz_) {}

  // central difference scheme for 1st derivative in 2d on three index variable 
  KOKKOS_INLINE_FUNCTION
  //double deriv24d(const FS4D &var,
  double derivVel(
      const int i, const int j, const int k,
      const int ih, const int jh, const int kh,
      const int v, const double d) const {
    return (vel(i-2*ih,j-2*jh,k-2*kh,v) - 8.0*vel(i-ih,j-jh,k-kh,v)
        + 8.0*vel(i+ih,j+jh,k+kh,v) - vel(i+2*ih,j+2*jh,k+2*kh,v)) / (12.0*d);
  }

  // central difference scheme for 1st derivative in 2d on two index variable 
  KOKKOS_INLINE_FUNCTION
  double derivRho(
      const int i, const int j, const int k,
      const int ih, const int jh, const int kh,
      const double d) const {
    return (rho(i-2*ih,j-2*jh,k-2*kh) - 8.0*rho(i-ih,j-jh,k-kh)
        + 8.0*rho(i+ih,j+jh,k+kh) - rho(i+2*ih,j+2*jh,k+2*kh)) / (12.0*d);
  }

  // central difference scheme for 1st derivative of specific internal energy
  KOKKOS_INLINE_FUNCTION
  double derivEnergy(
      const int i, const int j, const int k,
      const int ih, const int jh, const int kh,
      const double d) const {
    return (nrg(i-2*ih,j-2*jh,k-2*kh) - 8.0*nrg(i-ih,j-jh,k-kh)
        + 8.0*nrg(i+ih,j+jh,k+kh) - nrg(i+2*ih,j+2*jh,k+2*kh)) / (12.0*d);
  }

  // convert total energy to specific internal energy (TE-KE)/rho
  KOKKOS_INLINE_FUNCTION
  double nrg( const int i, const int j, const int k) const {
    return var(i,j,k,3)/rho(i,j,k)
      -0.5*rho(i,j,k)*(vel(i,j,k,0)*vel(i,j,k,0) + vel(i,j,k,1)*vel(i,j,k,1) + vel(i,j,k,2)*vel(i,j,k,2));
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, const int j, const int k) const {
    double dxr = derivRho(i,j,k,1,0,0,dx);
    double dyr = derivRho(i,j,k,0,1,0,dy);
    double dzr = derivRho(i,j,k,0,0,1,dz);

    double dxe = derivEnergy(i,j,k,1,0,0,dx);
    double dye = derivEnergy(i,j,k,0,1,0,dy);
    double dze = derivEnergy(i,j,k,0,0,1,dz);

    double dxu = derivVel(i,j,k,1,0,0,0,dx);
    double dyv = derivVel(i,j,k,0,1,0,1,dy);
    double dzw = derivVel(i,j,k,0,0,1,2,dz);

    double n1 = dxr;
    double n2 = dyr;
    double n3 = dzr;

    double rgrad = sqrt(dxr*dxr + dyr*dyr + dzr*dzr);
    double divu = dxu + dyv + dzw;

    double dnednr = (n1*dxe + n2*dye + n3*dze)*(n1*dxr + n2*dyr + n3*dzr);

    // compression switch
    int indicator = 0;
    if (dnednr < 0.0)
      indicator = 1;
    else
      indicator = 0;

    // detect shock front (equation 5a)
    if (divu < 0)
      gradRho(i, j, k, 0) = (1 - indicator) * rgrad;
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

struct calculateCeqFaces {
  FS4D var,varx;
  FS3D rho;
  FS6D mFlux;
  double alpha;
  int nv;

  calculateCeqFaces(FS4D var_, FS4D varx_, FS3D rho_, FS6D mFlux_, double a_, int nv_)
      : var(var_), varx(varx_), rho(rho_), mFlux(mFlux_), alpha(a_), nv(nv_){}

  KOKKOS_INLINE_FUNCTION
  double interpolateRho(const int i, const int j, const int k, const int ih, const int jh, const int kh) const {
    return ( rho(i,j,k) + rho(i+ih,j+jh,k+kh) )/2.0;
  }

  KOKKOS_INLINE_FUNCTION
  double interpolateC(const int i, const int j, const int k, const int ih, const int jh, const int kh, const int v) const {
    return ( var(i,j,k,nv+v) + var(i+ih,j+jh,k+kh,nv+v) )/2.0;
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, const int j, const int k) const {
    int ip = 0;
    int jp = 0;
    int kp = 0;

    double r,chat,cn[3];

    double M,cmag;

    // for each direction (i=0, j=1, k=2)
    for (int face = 0; face < 3; ++face) {
      ip = 0;
      jp = 0;
      kp = 0;
      if (face == 0)
        ip = 1;
      if (face == 1)
        jp = 1;
      if (face == 2)
        kp = 1;

      r     = interpolateRho(i,j,k,ip,jp,kp);
      chat  = interpolateC(i,j,k,ip,jp,kp,1);
      cn[0] = interpolateC(i,j,k,ip,jp,kp,2);
      cn[1] = interpolateC(i,j,k,ip,jp,kp,3);
      cn[2] = interpolateC(i,j,k,ip,jp,kp,4);

      cmag  = 1.0e-6;
      cmag += cn[0]*cn[0];
      cmag += cn[1]*cn[1];
      cmag += cn[2]*cn[2];

      for (int dir = 0; dir < 3; ++dir) {
        double dirac = 0.0;
        if (face == dir)
          dirac = 1.0;
        M=dirac - (cn[face]*cn[dir])/cmag;
        for(int w=0; w<3; ++w){
          mFlux(face,dir,i,j,k,w) = alpha*r*chat*M;
        }
      }
    }
    //varx(i,j,k,12)=mFlux(0,0,i,j,k,2);
    //varx(i,j,k,13)=mFlux(1,1,i,j,k,2);
    //varx(i,j,k,14)=mFlux(2,2,i,j,k,2);
    //varx(i,j,k,15)=mFlux(1,2,i,j,k,2);
    //varx(i,j,k,16)=mFlux(0,2,i,j,k,2);
    //varx(i,j,k,17)=mFlux(0,1,i,j,k,2);
  }
};

struct calculateCeqGrads {
  FS4D vel;
  FS6D mFlux;
  double dx,dy,dz;

  calculateCeqGrads(FS4D vel_, FS6D mFlux_, double dx_, double dy_, double dz_)
      : vel(vel_), mFlux(mFlux_), dx(dx_), dy(dy_), dz(dz_){}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, const int j, const int k) const {

    for (int w = 0; w < 3; ++w) {
      // right face
      mFlux(0,0,i,j,k,w) *= (vel(i+1,j,k,w) - vel(i,j,k,w)) / dx;
      mFlux(0,1,i,j,k,w) *= ( (vel(i,j+1,k,w) + vel(i+1,j+1,k,w))
                             -(vel(i,j-1,k,w) + vel(i+1,j-1,k,w)) ) / (4*dy);
      mFlux(0,2,i,j,k,w) *= ( (vel(i,j,k+1,w) + vel(i+1,j,k+1,w))
                             -(vel(i,j,k-1,w) + vel(i+1,j,k-1,w)) ) / (4*dz);

      // top face
      mFlux(1,0,i,j,k,w) *= ( (vel(i+1,j+1,k,w) + vel(i+1,j,k,w))
                             -(vel(i-1,j+1,k,w) + vel(i-1,j,k,w)) ) / (4*dx);
      mFlux(1,1,i,j,k,w) *= (vel(i,j+1,k,w) - vel(i,j,k,w)) / dy;
      mFlux(1,2,i,j,k,w) *= ( (vel(i,j+1,k+1,w) + vel(i,j,k+1,w))
                             -(vel(i,j+1,k-1,w) + vel(i,j,k-1,w)) ) / (4*dz);

      // front face
      mFlux(2,0,i,j,k,w) *= ( (vel(i+1,j,k+1,w) + vel(i+1,j,k,w))
                             -(vel(i-1,j,k+1,w) + vel(i-1,j,k,w)) ) / (4*dx);
      mFlux(2,1,i,j,k,w) *= ( (vel(i,j+1,k+1,w) + vel(i,j+1,k,w))
                             -(vel(i,j-1,k+1,w) + vel(i,j-1,k,w)) ) / (4*dy);
      mFlux(2,2,i,j,k,w) *= (vel(i,j,k+1,w) - vel(i,j,k,w)) / dz;
    }
  }
};

struct applyCeq {
  FS4D dvar,varx;
  FS6D mFlux;
  double dx,dy,dz;

  applyCeq(FS4D dvar_, FS4D varx_, FS6D mFlux_, double dx_, double dy_, double dz_)
      : dvar(dvar_), varx(varx_), mFlux(mFlux_), dx(dx_), dy(dy_), dz(dz_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, const int j, const int k) const {
    double diffu;
    int ih, jh, kh;
    double d[3];
    d[0]=dx;
    d[1]=dy;
    d[2]=dz;

    for (int w = 0; w < 3; ++w) {
      diffu = 0.0;
      for (int face=0; face<3; ++face){
        ih = 0; jh = 0; kh = 0;
        if (face == 0) ih = 1;
        if (face == 1) jh = 1;
        if (face == 2) kh = 1;
        for (int dir=0; dir<3; ++dir){
          diffu += (mFlux(face,dir,i,j,k,w)-mFlux(face,dir,i-ih,j-jh,k-kh,w))/d[face];
        }
      }
      dvar(i,j,k,w) += diffu;
      //varx(i,j,k,9+w) = diffu;
    }
  }
};
#endif
