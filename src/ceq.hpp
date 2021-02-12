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

  FS4D dvar;
  int v;
  Kokkos::View<double *> cd;

  maxCvar2D(FS4D dvar_, int v_, Kokkos::View<double *> cd_)
      : dvar(dvar_), v(v_), cd(cd_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, const int j, double &lmax) const {
    int ns = (int)cd(0);
    int nv = ns + 3;

    double s = dvar(i, j, 0, nv + v);

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

    double s = gradRho(i, j, v);

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

      Cp =
          Cp + (var(i, j, 0, 3 + s) / rho(i, j)) * (gammas * Rs / (gammas - 1));
      Cv = Cv + (var(i, j, 0, 3 + s) / rho(i, j)) * (Rs / (gammas - 1));
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
  FS2D rho;
  FS3D gradRho;
  Kokkos::View<double *> cd;

  calculateRhoGrad2D(FS4D var_, FS2D rho_, FS3D gradRho_,
                     Kokkos::View<double *> cd_)
      : var(var_), rho(rho_), gradRho(gradRho_), cd(cd_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, const int j) const {

    double i1 = 0;
    double i2 = 0;

    double ux1 = var(i - 2, j, 0, 0) / rho(i - 2, j);
    double ux2 = var(i - 1, j, 0, 0) / rho(i - 1, j);
    double ux3 = var(i + 1, j, 0, 0) / rho(i + 1, j);
    double ux4 = var(i + 2, j, 0, 0) / rho(i + 2, j);

    double uy1 = var(i, j - 2, 0, 0) / rho(i, j - 2);
    double uy2 = var(i, j - 1, 0, 0) / rho(i, j - 1);
    double uy3 = var(i, j + 1, 0, 0) / rho(i, j + 1);
    double uy4 = var(i, j + 2, 0, 0) / rho(i, j + 2);

    double vx1 = var(i - 2, j, 0, 1) / rho(i - 2, j);
    double vx2 = var(i - 1, j, 0, 1) / rho(i - 1, j);
    double vx3 = var(i + 1, j, 0, 1) / rho(i + 1, j);
    double vx4 = var(i + 2, j, 0, 1) / rho(i + 2, j);

    double vy1 = var(i, j - 2, 0, 1) / rho(i, j - 2);
    double vy2 = var(i, j - 1, 0, 1) / rho(i, j - 1);
    double vy3 = var(i, j + 1, 0, 1) / rho(i, j + 1);
    double vy4 = var(i, j + 2, 0, 1) / rho(i, j + 2);

    double ex1 = var(i - 2, j, 0, 2) / rho(i - 2, j) -
                 0.5 * rho(i - 2, j) * (ux1 * ux1 + vx1 * vx1);
    double ex2 = var(i - 1, j, 0, 2) / rho(i - 1, j) -
                 0.5 * rho(i - 1, j) * (ux2 * ux2 + vx2 * vx2);
    double ex3 = var(i + 1, j, 0, 2) / rho(i + 1, j) -
                 0.5 * rho(i + 1, j) * (ux3 * ux3 + vx3 * vx3);
    double ex4 = var(i + 2, j, 0, 2) / rho(i + 2, j) -
                 0.5 * rho(i + 2, j) * (ux4 * ux4 + vx4 * vx4);

    double ey1 = var(i, j - 2, 0, 2) / rho(i, j - 2) -
                 0.5 * rho(i, j - 2) * (uy1 * uy1 + vy1 * vy1);
    double ey2 = var(i, j - 1, 0, 2) / rho(i, j - 1) -
                 0.5 * rho(i, j - 1) * (uy2 * uy2 + vy2 * vy2);
    double ey3 = var(i, j + 1, 0, 2) / rho(i, j + 1) -
                 0.5 * rho(i, j + 1) * (uy3 * uy3 + vy3 * vy3);
    double ey4 = var(i, j + 2, 0, 2) / rho(i, j + 2) -
                 0.5 * rho(i, j + 2) * (uy4 * uy4 + vy4 * vy4);

    double dxr = (rho(i - 2, j) - 8.0 * rho(i - 1, j) + 8.0 * rho(i + 1, j) -
                  rho(i + 2, j)) /
                 (12.0 * cd(1));
    double dyr = (rho(i, j - 2) - 8.0 * rho(i, j - 1) + 8.0 * rho(i, j + 1) -
                  rho(i, j + 2)) /
                 (12.0 * cd(2));

    // double dxr = (rho(i+1,j) - rho(i-1,j))/(2.0*cd(1));
    // double dyr = (rho(i,j+1) - rho(i,j-1))/(2.0*cd(2));

    double dxe = (ex1 - 8.0 * ex2 + 8.0 * ex3 - ex4) / (12.0 * cd(1));
    double dye = (ey1 - 8.0 * ey2 + 8.0 * ey3 - ey4) / (12.0 * cd(2));

    double dxu = (ux1 - 8.0 * ux2 + 8.0 * ux3 - ux4) / (12.0 * cd(1));
    double dyv = (vy1 - 8.0 * vy2 + 8.0 * vy3 - vy4) / (12.0 * cd(2));

    double n1 = dxr;
    double n2 = dyr;

    double rgrad = sqrt(dxr * dxr + dyr * dyr);
    double divu = dxu + dyv;

    double dnednr = (n1 * dxe + n2 * dye) * (n1 * dxr + n2 * dyr);

    // compression switch
    if (dnednr < -1.0e-6)
      i1 = 1.0;
    // else
    //    i1 = 0;

    // detect shock front (equation 5a)
    if (divu < -1.0e-6) {
      i2 = 1.0;
    }

    // i1 = 1.0;
    // i2 = 1.0;
    // gradRho(i,j,0) = (1.0-indicator)*rgrad;
    // gradRho(i,j,k,0) = (1-indicator)*divu*rgrad;
    // else
    //    i2 = 0;
    // gradRho(i,j,0) = 0.0;

    // if (i1 == 0 && i2 == 1){
    //    gradRho(i,j,0) = rgrad;
    ////    printf("### DEBUG ### i2 rgrad (%d,%d) %f\n",i,j,rgrad);
    //}else
    //    gradRho(i,j,0) = 0.0;
    //
    // if (i1 == 1 && i2 == 0){
    ////    printf("### DEBUG ### i1 rgrad (%d,%d) %f, %f,
    ///%f\n",i,j,rgrad,-dxr,dyr);
    //    gradRho(i,j,1) = rgrad;
    //    gradRho(i,j,2) =-dxr;
    //    gradRho(i,j,3) = dyr;
    //}else{
    //    gradRho(i,j,1) = 0.0;
    //    gradRho(i,j,2) = 0.0;
    //    gradRho(i,j,3) = 0.0;
    //}
    //

    gradRho(i, j, 0) = (1 - i1) * i2 * rgrad;
    gradRho(i, j, 1) = i1 * rgrad;
    gradRho(i, j, 2) = -i1 * dyr;
    gradRho(i, j, 3) = i1 * dxr;

    // gradRho(i,j,0) = (1.0-indicator)*indicator2*rgrad;
    // gradRho(i,j,0) = indicator2*rgrad;
    // detect contact surface
    // gradRho(i,j,1) = indicator*rgrad;
    // gradRho(i,j,1) = (1.0-indicator2)*indicator*rgrad;

    // gradient components
    // gradRho(i,j,2) = -i1*dxr;
    // gradRho(i,j,3) =  i1*dyr;

    // if (j == 10 && i == 27){
    //    std::cout << gradRho(i-1,j,0) << ", "
    //              << gradRho(i  ,j,0) << ", "
    //              << gradRho(i+1,j,0) << std::endl;
    //}
    // if (j == 10 && i < 50)
    // if (j == 10)
    //    std::cout << i << ":  " << dnednr << endl;
    //    std::cout << gradRho(i-1,j,0) << endl;
  }
};

struct updateCeq2D {

  FS4D dvar;
  FS4D var;
  FS3D gradRho;
  double maxS, kap, eps;
  Kokkos::View<double *> cd;
  double maxG, maxGh, maxTau1, maxTau2;

  updateCeq2D(FS4D dvar_, FS4D var_, FS3D gradRho_, double maxS_,
              Kokkos::View<double *> cd_, double kap_, double eps_,
              double maxG_, double maxGh_, double maxTau1_, double maxTau2_)
      : dvar(dvar_), var(var_), gradRho(gradRho_), maxS(maxS_), cd(cd_),
        kap(kap_), eps(eps_), maxG(maxG_), maxGh(maxGh_), maxTau1(maxTau1_),
        maxTau2(maxTau2_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, const int j) const {
    double dx = cd(1);
    double dy = cd(2);
    double maxGr[5];

    // calculate cequation variable indices based on number of species (c
    // variables come after species densities)
    int nv = (int)cd(0) + 3;
    int nc = nv;

    double dx_right, dx_left, dy_top, dy_bot;
    double lap;

    // average cell size
    // double dxmag = pow(dx*dy,1.0/2.0);
    double dxmag = sqrt(dx * dx + dy * dy);

    maxGr[0] = maxG;
    maxGr[1] = maxGh;
    maxGr[2] = maxTau1;
    maxGr[3] = maxTau2;
    maxGr[4] = maxTau2;

    for (int n = 0; n < 4; ++n) {
      dx_right = (var(i + 1, j, 0, nc + n) - var(i, j, 0, nc + n)) / dx;
      dx_left = (var(i, j, 0, nc + n) - var(i - 1, j, 0, nc + n)) / dx;

      dy_top = (var(i, j + 1, 0, nc + n) - var(i, j, 0, nc + n)) / dy;
      dy_bot = (var(i, j, 0, nc + n) - var(i, j - 1, 0, nc + n)) / dy;

      // update ceq right hand side
      lap = 0.0;
      lap = dx_right - dx_left + dy_top - dy_bot;
      // lap = (dx_right - dx_left)/dx + (dy_top - dy_bot)/dy;
      // dvar(i,j,0,nc+n) = (maxS/(eps*dxmag))*(gradRho(i,j,n) -
      // var(i,j,0,nc+n)) + kap*maxS*dxmag*lap;
      dvar(i, j, 0, nc + n) =
          (maxS / (eps * dxmag)) *
              (gradRho(i, j, n) / maxGr[n] - var(i, j, 0, nc + n)) +
          kap * maxS * dxmag * lap;

      // if (n==0)
      //    printf("### DEBUG ### dvar,maxS,dxmag (%d,%d) %f, %f, %f,
      //    %f\n",i,j,dvar(i,j,0,nc+n),dxmag,maxS,lap);
    }
  }
};

struct applyCeq2D {
  FS4D dvar;
  FS4D var;
  FS2D rho;
  Kokkos::View<double *> cd;
  double betau, betae, maxC, maxCh, alpha, mu;

  applyCeq2D(FS4D dvar_, FS4D var_, FS2D rho_, double betau_, double betae_,
             double alpha_, double maxC_, double maxCh_, double mu_,
             Kokkos::View<double *> cd_)
      : dvar(dvar_), var(var_), rho(rho_), betau(betau_), betae(betae_),
        alpha(alpha_), maxC(maxC_), maxCh(maxCh_), mu(mu_), cd(cd_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, const int j) const {
    int nv = (int)cd(0) + 3;
    int nc = nv;
    double dx = cd(1);
    double dy = cd(2);
    double dxmag = dx * dx + dy * dy;

    double bTildU = (dxmag / maxC) * betau;
    double bTildE = (dxmag / maxC) * betae;
    double aTilde = (dxmag / (mu * mu * maxCh)) * alpha;

    double u = var(i, j, 0, 0) / rho(i, j); // u-velocity in current cell
    double ur = var(i + 1, j, 0, 0) / rho(i + 1, j); // u-velocity in right cell
    double ul = var(i - 1, j, 0, 0) / rho(i - 1, j); // u-velocity in left cell
    double ut = var(i, j + 1, 0, 0) / rho(i, j + 1); // u-velocity in top cell
    double ub =
        var(i, j - 1, 0, 0) / rho(i, j - 1); // u-velocity in bottom cell

    double v = var(i, j, 0, 1) / rho(i, j); // v-velocity
    double vr = var(i + 1, j, 0, 1) / rho(i + 1, j);
    double vl = var(i - 1, j, 0, 1) / rho(i - 1, j);
    double vt = var(i, j + 1, 0, 1) / rho(i, j + 1);
    double vb = var(i, j - 1, 0, 1) / rho(i, j - 1);

    double e = var(i, j, 0, 2) / rho(i, j); // energy
    double er = var(i + 1, j, 0, 2) / rho(i + 1, j);
    double el = var(i - 1, j, 0, 2) / rho(i - 1, j);
    double et = var(i, j + 1, 0, 2) / rho(i, j + 1);
    double eb = var(i, j - 1, 0, 2) / rho(i, j - 1);

    double rr = (rho(i, j) + rho(i + 1, j)) / 2.0; // density on right face
    double rl = (rho(i - 1, j) + rho(i, j)) / 2.0; // density on left face
    double rt = (rho(i, j) + rho(i, j + 1)) / 2.0; // density on top face
    double rb = (rho(i, j - 1) + rho(i, j)) / 2.0; // density on bottom face

    double cr =
        (var(i, j, 0, nc) + var(i + 1, j, 0, nc)) / 2.0; // iso C in right face
    double cl =
        (var(i - 1, j, 0, nc) + var(i, j, 0, nc)) / 2.0; // iso C in left face
    double ct =
        (var(i, j, 0, nc) + var(i, j + 1, 0, nc)) / 2.0; // iso C in top face
    double cb =
        (var(i, j - 1, 0, nc) + var(i, j, 0, nc)) / 2.0; // iso C in bottom face

    double chr = (var(i, j, 0, nc + 1) + var(i + 1, j, 0, nc + 1)) /
                 2.0; // iso C in right face
    double chl = (var(i - 1, j, 0, nc + 1) + var(i, j, 0, nc + 1)) /
                 2.0; // iso C in left face
    double cht = (var(i, j, 0, nc + 1) + var(i, j + 1, 0, nc + 1)) /
                 2.0; // iso C in top face
    double chb = (var(i, j - 1, 0, nc + 1) + var(i, j, 0, nc + 1)) /
                 2.0; // iso C in bottom face

    double t1r = (var(i, j, 0, nc + 2) + var(i + 1, j, 0, nc + 2)) /
                 2.0; // iso C in right face
    double t1l = (var(i - 1, j, 0, nc + 2) + var(i, j, 0, nc + 2)) /
                 2.0; // iso C in left face
    double t1t = (var(i, j, 0, nc + 2) + var(i, j + 1, 0, nc + 2)) /
                 2.0; // iso C in top face
    double t1b = (var(i, j - 1, 0, nc + 2) + var(i, j, 0, nc + 2)) /
                 2.0; // iso C in bottom face

    double t2r = (var(i, j, 0, nc + 3) + var(i + 1, j, 0, nc + 3)) /
                 2.0; // iso C in right face
    double t2l = (var(i - 1, j, 0, nc + 3) + var(i, j, 0, nc + 3)) /
                 2.0; // iso C in left face
    double t2t = (var(i, j, 0, nc + 3) + var(i, j + 1, 0, nc + 3)) /
                 2.0; // iso C in top face
    double t2b = (var(i, j - 1, 0, nc + 3) + var(i, j, 0, nc + 3)) /
                 2.0; // iso C in bottom face

    double dur = (ur - u) / dx; // u-velocity on right face
    double dul = (u - ul) / dx; // u-velocity on left face
    double dut = (ut - u) / dy; // u-velocity on top face
    double dub = (u - ub) / dy; // u-velocity on bottom face

    double dvr = (vr - v) / dx; // v-velocity on right face
    double dvl = (v - vl) / dx; // v-velocity on left face
    double dvt = (vt - v) / dy; // v-velocity on top face
    double dvb = (v - vb) / dy; // v-velocity on bottom face

    double der = (er - e) / dx; // energy on right face
    double del = (e - el) / dx; // energy on left face
    double det = (et - e) / dy; // energy on top face
    double deb = (e - eb) / dy; // energy on bottom face

    // isotropic diffusion operator
    double diffu = (bTildU / dx) * (rr * cr * dur - rl * cl * dul) +
                   (bTildU / dy) * (rt * ct * dut - rb * cb * dub);
    double diffv = (bTildU / dx) * (rr * cr * dvr - rl * cl * dvl) +
                   (bTildU / dy) * (rt * ct * dvt - rb * cb * dvb);
    double diffe = (bTildE / dx) * (rr * cr * der - rl * cl * del) +
                   (bTildE / dy) * (rt * ct * det - rb * cb * deb);

    dvar(i, j, 0, 0) += diffu;
    dvar(i, j, 0, 1) += diffv;
    dvar(i, j, 0, 2) += diffe;

    diffu = aTilde *
            ((rr * chr * t1r * t1r * dur - rl * chl * t1l * t1l * dul) / dx +
             (rt * cht * t1t * t1t * dut - rb * chb * t1b * t1b * dub) / dy +
             (rr * chr * t1r * t2r * dur - rl * chl * t1l * t2l * dul) / dx +
             (rt * cht * t1t * t2t * dut - rb * chb * t1b * t2b * dub) / dy);

    diffv = aTilde *
            ((rr * chr * t1r * t1r * dvr - rl * chl * t1l * t1l * dvl) / dx +
             (rt * cht * t1t * t1t * dvt - rb * chb * t1b * t1b * dvb) / dy +
             (rr * chr * t1r * t2r * dvr - rl * chl * t1l * t2l * dvl) / dx +
             (rt * cht * t1t * t2t * dvt - rb * chb * t1b * t2b * dvb) / dy);

    dvar(i, j, 0, 0) += diffu;
    dvar(i, j, 0, 1) += diffv;

    // if (diffu > 0.0 || diffv > 0.0)
    //    printf("### DEBUG ### %e: %e - %e\n",aTilde,diffu,diffv);
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
      gammas = cd(6 + 2 * s);
      Rs = cd(6 + 2 * s + 1);

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

    double dxr = (rho(i - 2, j, k) - 8.0 * rho(i - 1, j, k) +
                  8.0 * rho(i + 1, j, k) - rho(i + 2, j, k)) /
                 (12.0 * cd(1));
    double dyr = (rho(i, j - 2, k) - 8.0 * rho(i, j - 1, k) +
                  8.0 * rho(i, j + 1, k) - rho(i, j + 2, k)) /
                 (12.0 * cd(2));
    double dzr = (rho(i, j, k - 2) - 8.0 * rho(i, j, k - 1) +
                  8.0 * rho(i, j, k + 1) - rho(i, j, k + 2)) /
                 (12.0 * cd(3));

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
  FS4D gradRho;
  double maxS, kap, eps;
  Kokkos::View<double *> cd;

  updateCeq(FS4D dvar_, FS4D var_, FS4D gradRho_, double maxS_,
            Kokkos::View<double *> cd_, double kap_, double eps_)
      : dvar(dvar_), var(var_), gradRho(gradRho_), maxS(maxS_), cd(cd_),
        kap(kap_), eps(eps_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, const int j, const int k) const {
    double dx[3];
    dx[0] = cd(1);
    dx[1] = cd(2);
    dx[2] = cd(3);

    // calculate cequation variable indices based on number of species (c
    // variables come after species densities)
    int nv = (int)cd(0) + 4;
    int nc = nv;

    double dc_left[3][3];
    double dc_right[3][3];
    double lap;

    // average cell size
    double dxmag = pow(dx[0] * dx[1] * dx[2], 1.0 / 3.0);
    // double dxmag = sqrt(dx*dx + dy*dy + dz*dz);

    for (int n = 0; n < 5; ++n) {
      // left face
      dc_left[0][0] = (var(i, j, k, nc + n) - var(i - 1, j, k, nc + n)) / dx[0];
      dc_left[0][1] =
          (var(i - 1, j + 1, k, nc + n) + var(i, j + 1, k, nc + n) +
           var(i - 1, j - 1, k, nc + n) + var(i, j - 1, k, nc + n)) /
          (4 * dx[1]);
      dc_left[0][2] =
          (var(i - 1, j, k + 1, nc + n) + var(i, j, k + 1, nc + n) +
           var(i - 1, j, k - 1, nc + n) + var(i, j, k - 1, nc + n)) /
          (4 * dx[2]);

      // bottom face
      dc_left[1][0] =
          (var(i - 1, j, k, nc + n) + var(i + 1, j, k, nc + n) +
           var(i - 1, j - 1, k, nc + n) + var(i + 1, j - 1, k, nc + n)) /
          (4 * dx[0]);
      dc_left[1][1] = (var(i, j, k, nc + n) - var(i, j - 1, k, nc + n)) / dx[1];
      dc_left[1][2] =
          (var(i, j, k - 1, nc + n) + var(i, j, k + 1, nc + n) +
           var(i, j - 1, k - 1, nc + n) + var(i, j - 1, k + 1, nc + n)) /
          (4 * dx[2]);

      // back (hind) face
      dc_left[2][0] =
          (var(i - 1, j, k, nc + n) + var(i + 1, j, k, nc + n) +
           var(i - 1, j, k - 1, nc + n) + var(i + 1, j, k - 1, nc + n)) /
          (4 * dx[0]);
      dc_left[2][1] =
          (var(i, j - 1, k, nc + n) + var(i, j + 1, k, nc + n) +
           var(i, j - 1, k - 1, nc + n) + var(i, j + 1, k - 1, nc + n)) /
          (4 * dx[1]);
      dc_left[2][2] = (var(i, j, k, nc + n) - var(i, j, k - 1, nc + n)) / dx[2];

      // right face
      dc_right[0][0] =
          (var(i + 1, j, k, nc + n) - var(i, j, k, nc + n)) / dx[0];
      dc_right[0][1] =
          (var(i, j + 1, k, nc + n) + var(i + 1, j + 1, k, nc + n) +
           var(i, j - 1, k, nc + n) + var(i + 1, j - 1, k, nc + n)) /
          (4 * dx[1]);
      dc_right[0][2] =
          (var(i, j, k + 1, nc + n) + var(i + 1, j, k + 1, nc + n) +
           var(i, j, k - 1, nc + n) + var(i + 1, j, k - 1, nc + n)) /
          (4 * dx[2]);

      // top face
      dc_right[1][0] =
          (var(i - 1, j + 1, k, nc + n) + var(i + 1, j + 1, k, nc + n) +
           var(i - 1, j, k, nc + n) + var(i + 1, j, k, nc + n)) /
          (4 * dx[0]);
      dc_right[1][1] =
          (var(i, j + 1, k, nc + n) - var(i, j, k, nc + n)) / dx[1];
      dc_right[1][2] =
          (var(i, j + 1, k - 1, nc + n) + var(i, j + 1, k + 1, nc + n) +
           var(i, j, k - 1, nc + n) + var(i, j, k + 1, nc + n)) /
          (4 * dx[2]);

      // front face
      dc_right[2][0] =
          (var(i - 1, j, k + 1, nc + n) + var(i + 1, j, k + 1, nc + n) +
           var(i - 1, j, k, nc + n) + var(i + 1, j, k, nc + n)) /
          (4 * dx[0]);
      dc_right[2][1] =
          (var(i, j - 1, k + 1, nc + n) + var(i, j + 1, k + 1, nc + n) +
           var(i, j - 1, k, nc + n) + var(i, j + 1, k, nc + n)) /
          (4 * dx[1]);
      dc_right[2][2] =
          (var(i, j, k + 1, nc + n) - var(i, j, k, nc + n)) / dx[2];

      // update ceq right hand side
      lap = 0;
      for (int d = 0; d < 3; ++d) {
        for (int f = 0; f < 3; ++f) {
          lap = lap + (dc_right[d][f] - dc_left[d][f]) / dx[d];
        }
      }
      dvar(i, j, k, nc + n) =
          maxS / (eps * dxmag) * (gradRho(i, j, k, n) - var(i, j, k, nc + n)) +
          kap * maxS * dxmag * lap;
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

        ip = 0;
        jp = 0;
        kp = 0;
      }
    }
  }
};

struct applyCeq {

  FS4D dvar;
  FS4D var;
  FS3D rho;
  FS6D mFlux; //(m,n,i,j,k,dir)
  FS4D cFlux; //(i,j,k,dir)
  double alpha, beta, betae;
  Kokkos::View<double *> cd;

  applyCeq(FS4D dvar_, FS4D var_, FS3D rho_, FS6D mFlux_, FS4D cFlux_,
           Kokkos::View<double *> cd_, double alpha_, double beta_,
           double betae_)
      : dvar(dvar_), var(var_), rho(rho_), mFlux(mFlux_), cFlux(cFlux_),
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
    for (int n = 0; n < 4; ++n) {
      // left face
      du_left[0][0] = (var(i, j, k, n) / rho(i, j, k) -
                       var(i - 1, j, k, n) / rho(i - 1, j, k)) /
                      cd(1);
      du_left[0][1] = (var(i - 1, j + 1, k, n) / rho(i - 1, j + 1, k) +
                       var(i, j + 1, k, n) / rho(i, j + 1, k) +
                       var(i - 1, j - 1, k, n) / rho(i - 1, j - 1, k) +
                       var(i, j - 1, k, n) / rho(i, j - 1, k)) /
                      (4 * cd(2));
      du_left[0][2] = (var(i - 1, j, k + 1, n) / rho(i - 1, j, k + 1) +
                       var(i, j, k + 1, n) / rho(i, j, k + 1) +
                       var(i - 1, j, k - 1, n) / rho(i - 1, j, k - 1) +
                       var(i, j, k - 1, n) / rho(i, j, k - 1)) /
                      (4 * cd(3));

      // bottom face
      du_left[1][0] = (var(i - 1, j, k, n) / rho(i - 1, j, k) +
                       var(i + 1, j, k, n) / rho(i + 1, j, k) +
                       var(i - 1, j - 1, k, n) / rho(i - 1, j - 1, k) +
                       var(i + 1, j - 1, k, n) / rho(i + 1, j - 1, k)) /
                      (4 * cd(1));
      du_left[1][1] = (var(i, j, k, n) / rho(i, j, k) -
                       var(i, j - 1, k, n) / rho(i, j - 1, k)) /
                      cd(2);
      du_left[1][2] = (var(i, j, k - 1, n) / rho(i, j, k - 1) +
                       var(i, j, k + 1, n) / rho(i, j, k + 1) +
                       var(i, j - 1, k - 1, n) / rho(i, j - 1, k - 1) +
                       var(i, j - 1, k + 1, n) / rho(i, j - 1, k + 1)) /
                      (4 * cd(3));

      // back (hind) face
      du_left[2][0] = (var(i - 1, j, k, n) / rho(i - 1, j, k) +
                       var(i + 1, j, k, n) / rho(i + 1, j, k) +
                       var(i - 1, j, k - 1, n) / rho(i - 1, j, k - 1) +
                       var(i + 1, j, k - 1, n) / rho(i + 1, j, k - 1)) /
                      (4 * cd(1));
      du_left[2][1] = (var(i, j - 1, k, n) / rho(i, j - 1, k) +
                       var(i, j + 1, k, n) / rho(i, j + 1, k) +
                       var(i, j - 1, k - 1, n) / rho(i, j - 1, k - 1) +
                       var(i, j + 1, k - 1, n) / rho(i, j + 1, k - 1)) /
                      (4 * cd(2));
      du_left[2][2] = (var(i, j, k, n) / rho(i, j, k) -
                       var(i, j, k - 1, n) / rho(i, j, k - 1)) /
                      cd(3);

      // right face
      du_right[0][0] = (var(i + 1, j, k, n) / rho(i + 1, j, k) -
                        var(i, j, k, n) / rho(i, j, k)) /
                       cd(1);
      du_right[0][1] = (var(i, j + 1, k, n) / rho(i, j + 1, k) +
                        var(i + 1, j + 1, k, n) / rho(i + 1, j + 1, k) +
                        var(i, j - 1, k, n) / rho(i, j - 1, k) +
                        var(i + 1, j - 1, k, n) / rho(i + 1, j - 1, k)) /
                       (4 * cd(2));
      du_right[0][2] = (var(i, j, k + 1, n) / rho(i, j, k + 1) +
                        var(i + 1, j, k + 1, n) / rho(i + 1, j, k + 1) +
                        var(i, j, k - 1, n) / rho(i, j, k - 1) +
                        var(i + 1, j, k - 1, n) / rho(i + 1, j, k - 1)) /
                       (4 * cd(3));

      // top face
      du_right[1][0] = (var(i - 1, j + 1, k, n) / rho(i - 1, j + 1, k) +
                        var(i + 1, j + 1, k, n) / rho(i + 1, j + 1, k) +
                        var(i - 1, j, k, n) / rho(i - 1, j, k) +
                        var(i + 1, j, k, n) / rho(i + 1, j, k)) /
                       (4 * cd(1));
      du_right[1][1] = (var(i, j + 1, k, n) / rho(i, j + 1, k) -
                        var(i, j, k, n) / rho(i, j, k)) /
                       cd(2);
      du_right[1][2] = (var(i, j + 1, k - 1, n) / rho(i, j + 1, k - 1) +
                        var(i, j + 1, k + 1, n) / rho(i, j + 1, k + 1) +
                        var(i, j, k - 1, n) / rho(i, j, k - 1) +
                        var(i, j, k + 1, n) / rho(i, j, k + 1)) /
                       (4 * cd(3));

      // front face
      du_right[2][0] = (var(i - 1, j, k + 1, n) / rho(i - 1, j, k + 1) +
                        var(i + 1, j, k + 1, n) / rho(i + 1, j, k + 1) +
                        var(i - 1, j, k, n) / rho(i - 1, j, k) +
                        var(i + 1, j, k, n) / rho(i + 1, j, k)) /
                       (4 * cd(1));
      du_right[2][1] = (var(i, j - 1, k + 1, n) / rho(i, j - 1, k + 1) +
                        var(i, j + 1, k + 1, n) / rho(i, j + 1, k + 1) +
                        var(i, j - 1, k, n) / rho(i, j - 1, k) +
                        var(i, j + 1, k, n) / rho(i, j + 1, k)) /
                       (4 * cd(2));
      du_right[2][2] = (var(i, j, k + 1, n) / rho(i, j, k + 1) -
                        var(i, j, k, n) / rho(i, j, k)) /
                       cd(3);

      an = 0;
      is = 0;
      ip = 0;
      jp = 0;
      kp = 0;

      if (n < 3) {
        diffu = 0;
        for (int d = 0; d < 3; ++d) { // each direction for divergence
          if (d == 0)
            ip = 1;
          if (d == 1)
            jp = 1;
          if (d == 2)
            kp = 1;

          for (int f = 0; f < 3; ++f) { // each direction for gradient
            an = an + (mFlux(d, f, i, j, k, d) * du_right[d][f] -
                       mFlux(d, f, i - ip, j - jp, k - kp, d) * du_left[d][f]) /
                          dx[d];
            is = is + (cFlux(i, j, k, d) * du_right[d][f] -
                       cFlux(i - ip, j - jp, k - kp, d) * du_left[d][f]) /
                          dx[d];
          }

          ip = 0;
          jp = 0;
          kp = 0;
        }

        diffu = alpha * an + beta * is;

        dvar(i, j, k, n) = dvar(i, j, k, n) + diffu;
      } else {
        diffu = 0;
        for (int d = 0; d < 3; ++d) { // each direction for divergence
          if (d == 0)
            ip = 1;
          if (d == 1)
            jp = 1;
          if (d == 2)
            kp = 1;

          for (int f = 0; f < 3; ++f) { // each direction for gradient
            diffu = diffu + (cFlux(i, j, k, d) * du_right[d][f] -
                             cFlux(i - ip, j - jp, k - kp, d) * du_left[d][f]) /
                                dx[d];
          }

          ip = 0;
          jp = 0;
          kp = 0;
        }
        dvar(i, j, k, n) = dvar(i, j, k, n) + betae * diffu;
      }
    }
  }
};
#endif
