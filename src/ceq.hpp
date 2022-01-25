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

struct maxCvar2D {
  FS4D var;
  int v;
  Kokkos::View<FSCAL *> cd;

  maxCvar2D(FS4D var_, int v_, Kokkos::View<FSCAL *> cd_)
      : var(var_), v(v_), cd(cd_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, const int j, FSCAL &lmax) const {
    int ns = (int)cd(0);
    int nv = ns + 3;

    FSCAL s = abs(var(i, j, 0, nv + v));

    if (s > lmax)
      lmax = s;
  }
};

struct maxWaveSpeed2D {
  FS4D var;
  FS2D p;
  FS2D rho;
  Kokkos::View<FSCAL *> cd;

  maxWaveSpeed2D(FS4D var_, FS2D p_, FS2D rho_, Kokkos::View<FSCAL *> cd_)
      : var(var_), p(p_), rho(rho_), cd(cd_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, const int j, FSCAL &lmax) const {
    int ns = (int)cd(0);
    FSCAL gamma, gammas, Rs;
    FSCAL Cp = 0;
    FSCAL Cv = 0;

    FSCAL a, s;

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
  FSCAL dx,dy;

  calculateRhoGrad2D(FS4D var_, FS3D vel_, FS2D rho_, FS3D gradRho_, FSCAL dx_, FSCAL dy_)
      : var(var_), vel(vel_), rho(rho_), gradRho(gradRho_), dx(dx_), dy(dy_) {}

  // central difference scheme for 1st derivative in 2d on two index variable 
  KOKKOS_INLINE_FUNCTION
  FSCAL rhoDerivative(const int i, const int j, const int ih, const int jh, const FSCAL dx) const {
    return (rho(i-2*ih,j-2*jh) - 8.0*rho(i-ih,j-jh) + 8.0*rho(i+ih,j+jh) - rho(i+2*ih,j+2*jh)) / (12.0*dx);
  }

  // central difference scheme for 1st derivative of specific internal energy
  KOKKOS_INLINE_FUNCTION
  FSCAL energyDerivative(const int i, const int j, const int ih, const int jh, const FSCAL dx) const {
    return (nrg(i-2*ih,j-2*jh) - 8.0*nrg(i-ih,j-jh) + 8.0*nrg(i+ih,j+jh) - nrg(i+2*ih,j+2*jh)) / (12.0*dx);
  }

  // convert total energy to specific internal energy (TE-KE)/rho
  KOKKOS_INLINE_FUNCTION
  FSCAL nrg(const int i, const int j) const {
    return var(i,j,0,2)/rho(i,j)-0.5*rho(i,j)*(vel(i,j,0)*vel(i,j,0)+vel(i,j,1)*vel(i,j,1));
  }

  // central difference scheme for 1st derivative in 2d on three index variable 
  KOKKOS_INLINE_FUNCTION
  FSCAL velocityDerivative(const int i, const int j, const int ih, const int jh, const int v, const FSCAL dx) const {
    return (vel(i-2*ih,j-2*jh,v) - 8.0*vel(i-ih,j-jh,v) + 8.0*vel(i+ih,j+jh,v) - vel(i+2*ih,j+2*jh,v)) / (12.0*dx);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, const int j) const {
    FSCAL dxr = rhoDerivative(i,j,1,0,dx);
    FSCAL dyr = rhoDerivative(i,j,0,1,dy);

    FSCAL dxe = energyDerivative(i,j,1,0,dx);
    FSCAL dye = energyDerivative(i,j,0,1,dy);

    FSCAL dxu = velocityDerivative(i,j,1,0,0,dx);
    FSCAL dyv = velocityDerivative(i,j,0,1,1,dy);

    FSCAL n1 = dxr;
    FSCAL n2 = dyr;

    FSCAL rgrad = sqrt(dxr*dxr + dyr*dyr);
    FSCAL divu = dxu + dyv;

    FSCAL dnednr = (n1*dxe + n2*dye)*(n1*dxr + n2*dyr);

    // compression switch
    FSCAL i1 = 0;
    if (dnednr < 0.0)
      i1 = 1.0;
    else
      i1 = 0.0;

    // detect shock front (equation 5a)
    FSCAL i2 = 0;
    if (divu < 0.0) {
      i2 = 1.0;
    }

    gradRho(i, j, 0) = (1 - i1)*i2*rgrad;
    gradRho(i, j, 1) = i1*rgrad;
    gradRho(i, j, 2) = -dyr;
    gradRho(i, j, 3) = dxr;
  }
};

struct updateCeq2D {
  FS4D dvar;
  FS4D var;
  FS3D gradRho;
  FSCAL maxS, kap, eps;
  Kokkos::View<FSCAL *> cd;

  updateCeq2D(FS4D dvar_, FS4D var_, FS3D gradRho_, FSCAL maxS_,
              Kokkos::View<FSCAL *> cd_, FSCAL kap_, FSCAL eps_)
      : dvar(dvar_), var(var_), gradRho(gradRho_), maxS(maxS_), cd(cd_),
        kap(kap_), eps(eps_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, const int j) const {
    FSCAL dx = cd(1);
    FSCAL dy = cd(2);

    // calculate cequation variable indices based on number of species (c
    // variables come after species densities)
    int nv = (int)cd(0) + 3;
    int nc = nv;

    FSCAL dx_right, dx_left, dy_top, dy_bot;
    FSCAL lap;

    // average cell size
    FSCAL dxmag = sqrt(dx * dx + dy * dy);

    for (int n = 0; n < 4; ++n) {
      dx_right = ( var(i+1,j,0,nc+n) - var(i,j,0,nc+n) )/dx;
      dx_left  = ( var(i,j,0,nc+n) - var(i-1,j,0,nc+n) )/dx;

      dy_top = ( var(i,j+1,0,nc+n) - var(i,j,0,nc+n) )/dy;
      dy_bot = ( var(i,j,0,nc+n) - var(i,j-1,0,nc+n) )/dy;

      // update ceq right hand side
      lap = dx_right - dx_left + dy_top - dy_bot;

      dvar(i,j,0,nc+n) =
          (maxS/(eps*dxmag))*( gradRho(i,j,n) - var(i,j,0,nc+n) ) + kap*maxS*dxmag*lap;
    }
  }
};

struct computeCeqFlux2D {
  FS4D var;
  FS5D m; // m(face,direction of derivative, i, j, velocity component)
  FS2D rho;
  FSCAL alpha;
  int nv;
  FSCAL maxCh;

  computeCeqFlux2D(FS4D var_, FS5D m_, FS2D rho_, FSCAL a_, int nv_, FSCAL maxCh_)
      : var(var_), m(m_), rho(rho_), alpha(a_), nv(nv_), maxCh(maxCh_) {}

  KOKKOS_INLINE_FUNCTION
  FSCAL interpolateRho(const int i, const int j, const int ih, const int jh) const {
    return ( rho(i,j) + rho(i+ih,j+jh) )/2.0;
  }

  KOKKOS_INLINE_FUNCTION
  FSCAL interpolateC(const int i, const int j, const int ih, const int jh, const int v) const {
    return ( var(i,j,0,nv+v) + var(i+ih,j+jh,0,nv+v) )/2.0;
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, const int j) const {
    int ih; // ihat - i-direction unit vector
    int jh; // jhat - j-direction unit vector
    FSCAL chat,r,ctmag;
    FSCAL ct[2];

    // Compute alph*rho*C*Ctau1*Ctau2 on each face for each direction
    for (int f=0;f<2;++f){
      // set direction unit vector
      ih=0; jh=0;
      if (f==0) ih=1;
      if (f==1) jh=1;

      r     = interpolateRho(i,j,ih,jh);
      chat  = interpolateC(i,j,ih,jh,1);
      ct[0] = interpolateC(i,j,ih,jh,2);
      ct[1] = interpolateC(i,j,ih,jh,3);
      ctmag = (ct[0]*ct[0] + ct[1]*ct[1])+1.0e-6;

      for (int d=0;d<2;++d){
        for (int w=0;w<2;++w){
          m(f,d,i,j,w) = alpha*r*chat*(ct[f]*ct[d])/ctmag;
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
    FSCAL dx = cd(1);
    FSCAL dy = cd(2);

    // Compute dux and duy on each face for each velocity
    for (int w=0;w<2;++w){
      m(0,0,i,j,w) *= (  v(i+1,j,w)              - v(i,j,w) )/dx;
      m(0,1,i,j,w) *= ( (v(i,j+1,w)+v(i+1,j+1,w))-(v(i,j-1,w)+v(i+1,j-1,w)) )/(4.0*dy);
      m(1,0,i,j,w) *= ( (v(i+1,j,w)+v(i+1,j+1,w))-(v(i-1,j,w)+v(i-1,j+1,w)) )/(4.0*dx);
      m(1,1,i,j,w) *= (  v(i,j+1,w)              - v(i,j,w) )/dy;
    }
  }
};

struct applyCeq2D {
  FS5D m;
  FS4D dvar,varx;
  FSCAL dx,dy;

  applyCeq2D(FS4D dvar_, FS4D varx_, FS5D m_, FSCAL dx_, FSCAL dy_)
    : dvar(dvar_), varx(varx_), m(m_), dx(dx_), dy(dy_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, const int j) const {
    FSCAL diffu;

    // assembly artificial viscosity term
    for (int w=0;w<2;++w){
      diffu=0;
      diffu += (m(0,0,i,j,w)-m(0,0,i-1,j,w))/dx;
      diffu += (m(0,1,i,j,w)-m(0,1,i-1,j,w))/dx;
      diffu += (m(1,0,i,j,w)-m(1,0,i,j-1,w))/dy;
      diffu += (m(1,1,i,j,w)-m(1,1,i,j-1,w))/dy;

      dvar(i,j,0,w)+=diffu;
      //varx(i,j,0,5+w)=diffu;
    }
  }
};
#endif
