
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
      const int v, const double dx) const {
    return (vel(i-2*ih,j-2*jh,k-2*kh,v) - 8.0*vel(i-ih,j-jh,j-kh,v)
        + 8.0*vel(i+ih,j+jh,k+kh,v) - vel(i+2*ih,j+2*jh,k+2*kh,v)) / (12.0*dx);
  }

  // central difference scheme for 1st derivative in 2d on two index variable 
  KOKKOS_INLINE_FUNCTION
  double derivRho(
      const int i, const int j, const int k,
      const int ih, const int jh, const int kh,
      const double dx) const {
    return (rho(i-2*ih,j-2*jh,k-2*kh) - 8.0*rho(i-ih,j-jh,k-kh)
        + 8.0*rho(i+ih,j+jh,k+kh) - rho(i+2*ih,j+2*jh,k+2*kh)) / (12.0*dx);
  }

  // central difference scheme for 1st derivative of specific internal energy
  KOKKOS_INLINE_FUNCTION
  double derivEnergy(
      const int i, const int j, const int k,
      const int ih, const int jh, const int kh,
      const double dx) const {
    return (nrg(i-2*ih,j-2*jh,k-2*kh) - 8.0*nrg(i-ih,j-jh,k-kh)
        + 8.0*nrg(i+ih,j+jh,k+kh) - nrg(i+2*ih,j+2*jh,k+2*kh)) / (12.0*dx);
  }

  // convert total energy to specific internal energy (TE-KE)/rho
  KOKKOS_INLINE_FUNCTION
  double nrg(
      const int i, const int j, const int k) const {
    return var(i,j,k,2)/rho(i,j,k)
      -0.5*(vel(i,j,k,0)*vel(i,j,k,0)+vel(i,j,k,1)*vel(i,j,k,1)+vel(i,j,k,2)*vel(i,j,k,2));
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

    double rgrad = sqrt(dxr * dxr + dyr * dyr + dzr * dzr);
    double divu = dxu + dyv + dzw;

    double dnednr =
        (n1 * dxe + n2 * dye + n3 * dze) * (n1 * dxr + n2 * dyr + n3 * dzr);

    // compression switch
    int indicator = 0;
    if (dnednr < 0)
      indicator = 1;
    else
      indicator = 0;

    // detect shock front (equation 5a)
    if (divu < 0)
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

      //if (cmag_left <= 0.00000001 || cmag_right <= 0.00000001) {
      if (cmag_left <= 1.0e-10 || cmag_right <= 1.0e-10) {
        double dirac;
        for (int m = 0; m < 3; ++m) {
          for (int n = 0; n < 3; ++n) {
            dirac=0.0;
            if (m==n)
              dirac=1.0;
            mFlux(m, n, i, j, k, dir) = dirac;
            cFlux(i, j, k, dir) = 0.0;
          }
        }
      } else {
        // tensor components
        //for (int m = 0; m < 3; ++m) {
          for (int n = 0; n < 3; ++n) {

            // dirac delta
            double d = 0;
            if (dir == n)
              d = 1;

            // calculate right and left tensor components
            m_left = d - (cn_left[dir] * cn_left[n] / cmag_left);
            m_right = d - (cn_right[dir] * cn_right[n] / cmag_right);
            //m_left = d - (cn_left[m] * cn_left[n]);
            //m_right = d - (cn_right[m] * cn_right[n]);

            // include isotropic c
            m_left = m_left * ch_left;
            m_right = m_right * ch_right;

            // include density
            m_left = m_left * rho_left;
            m_right = m_right * rho_right;

            // find flux
            mFlux(dir, n, i, j, k, 0) = (m_right + m_left) / 2.0;
          }
        //}

        // calcualte isotropic C flux
        cFlux(i, j, k, dir) = (c_right * rho_right + c_left * rho_left) / 2.0;
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

      an = 0.0;
      is = 0.0;
      ip = 0.0;
      jp = 0.0;
      kp = 0.0;

      //if (n < 3) {
        diffu = 0.0;
        for (int d = 0; d < 3; ++d) { // each direction for divergence
          ip = 0; jp = 0; kp = 0;
          if (d == 0) ip = 1;
          if (d == 1) jp = 1;
          if (d == 2) kp = 1;

          for (int f = 0; f < 3; ++f) { // each direction for gradient
            an = an + (mFlux(d,f,i,j,k,0) * du_right[d][f] -
                       mFlux(d,f,i-ip,j-jp,k-kp,0) * du_left[d][f]) / dx[d];
            is = is + (cFlux(i,j,k,d) * du_right[d][f] -
                       cFlux(i-ip,j-jp,k-kp,d) * du_left[d][f]) / dx[d];
          }

        }

        diffu = alpha*an + beta*is;

        dvar(i,j,k,n) += diffu;
        //if (n==2) varx(i,j,k,9) = dvar(i,j,k,n);
        //if (n==1) varx(i,j,k,9) = 1.0;
        if (n==2) varx(i,j,k,9) = diffu;
        //if (n==2) varx(i,j,k,9) = dvar(i,j,k,n);
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
