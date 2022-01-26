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

#include "kokkosTypes.hpp"

#ifndef VISCOSOTY_HPP
#define VISCOSOTY_HPP

struct calculateStressTensor2dv {
  FS4D var;
  FS2D rho;
  FS4D stressx;
  FS4D stressy;
  FS3D vel;
  Kokkos::View<FSCAL *> cd;

  calculateStressTensor2dv(FS4D var_, FS2D rho_, FS3D v_, FS4D strx_,
                           FS4D stry_, Kokkos::View<FSCAL *> cd_)
      : var(var_), rho(rho_), vel(v_), stressx(strx_), stressy(stry_), cd(cd_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, const int j) const {
    int ns = (int)cd(0);
    FSCAL dx = cd(1);
    FSCAL dy = cd(2);
    // FSCAL mu1 = 2.928e-5;
    // FSCAL mu2 = 1.610e-5;
    FSCAL dudx, dvdy, dudy, dvdx;

    //  FSCAL muij = 0.0;
    //  FSCAL muip = 0.0;
    //  FSCAL mujp = 0.0;

    // for (int s = 0; s < ns; ++s) {
    //   muij += var(i, j, 0, 3 + s) * cd(6 + 3 * s + 2) / rho(i, j);
    //   muip += var(i + 1, j, 0, 3 + s) * cd(6 + 3 * s + 2) / rho(i + 1, j);
    //   mujp += var(i, j + 1, 0, 3 + s) * cd(6 + 3 * s + 2) / rho(i, j + 1);
    // }

    // // FSCAL muij = (var(i,j,0,3)*mu1 + var(i,j,0,4)*mu2)/rho(i,j);
    // // FSCAL muip = (var(i+1,j,0,3)*mu1 + var(i+1,j,0,4)*mu2)/rho(i+1,j);
    // // FSCAL mujp = (var(i,j+1,0,3)*mu1 + var(i,j+1,0,4)*mu2)/rho(i,j+1);

    // FSCAL mur = (muip + muij) / 2;
    // FSCAL mut = (mujp + muij) / 2;
    FSCAL mur = 2.928e-5;
    FSCAL mut = 2.928e-5;

    // xface
    dudx = (vel(i+1,j,0) -vel(i,j,0))/dx;
    dvdx = (vel(i+1,j,1) -vel(i,j,1))/dx;
    dudy =( (vel(i,j+1,0)+vel(i+1,j+1,0)) - (vel(i,j-1,1)+vel(i+1,j-1,0)) )/(4.0*dy);
    dvdy =( (vel(i,j+1,1)+vel(i+1,j+1,1)) - (vel(i,j-1,1)+vel(i+1,j-1,1)) )/(4.0*dy);

    // if (dudx <= 1e-8)
    //   dudx = 0.00;
    // if (dvdx <= 1e-8)
    //   dvdx = 0.00;
    // if (dudy <= 1e-8)
    //   dudy = 0.00;
    // if (dvdy <= 1e-8)
    //   dvdy = 0.00;

    stressx(i, j, 0, 0) = (2.0 / 3.0) * mur * (2.0 * dudx - dvdy);
    stressx(i, j, 1, 1) = (2.0 / 3.0) * mur * (2.0 * dvdy - dudx);
    stressx(i, j, 0, 1) = mur * (dudy + dvdx);
    stressx(i, j, 1, 0) = stressx(i, j, 0, 1);

    // yface
    dudx = ( (vel(i+1,j,0)+vel(i+1,j+1,0)) - (vel(i-1,j,0)+vel(i-1,j+1,0)) )/(4.0*dx);
    dvdx = ( (vel(i+1,j,1)+vel(i+1,j+1,1)) - (vel(i-1,j,1)+vel(i-1,j+1,1)) )/(4.0*dx);
    dudy = ( vel(i,j+1,0)-vel(i,j,0) )/dy;
    dvdy = ( vel(i,j+1,1)-vel(i,j,1) )/dy;

    // if (dudx <= 1e-8)
    //   dudx = 0.00;
    // if (dvdx <= 1e-8)
    //   dvdx = 0.00;
    // if (dudy <= 1e-8)
    //   dudy = 0.00;
    // if (dvdy <= 1e-8)
    //   dvdy = 0.00;

    stressy(i, j, 0, 0) = (2.0 / 3.0) * mut * (2.0 * dudx - dvdy);
    stressy(i, j, 1, 1) = (2.0 / 3.0) * mut * (2.0 * dvdy - dudx);
    stressy(i, j, 0, 1) = mut * (dudy + dvdx);
    stressy(i, j, 1, 0) = stressy(i, j, 0, 1);

    // printf("(%2d,%2d) %.1e, %.1e, %.1e == %.1e, %.1e, %.1e, \n",i,j,
    //     stressx(i,j,0,0),stressx(i,j,1,1),stressx(i,j,0,1),
    //     stressy(i,j,0,0),stressy(i,j,1,1),stressy(i,j,0,1));
  }
};

struct calculateHeatFlux2dv {
  FS4D var;
  FS2D rho;
  FS2D qx;
  FS2D qy;
  FS2D T;
  Kokkos::View<FSCAL *> cd;

  FSCAL k = 0.026;

  calculateHeatFlux2dv(FS4D var_, FS2D rho_, FS2D T_, FS2D qx_, FS2D qy_,
                       Kokkos::View<FSCAL *> cd_)
      : var(var_), rho(rho_), T(T_), qx(qx_), qy(qy_), cd(cd_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, const int j) const {
    FSCAL dx = cd(1);
    FSCAL dy = cd(2);

    qx(i, j) = -k * (T(i + 1, j) - T(i, j)) / dx;
    qy(i, j) = -k * (T(i, j + 1) - T(i, j)) / dy;

    // if (i==2002 && j==3)
    //    printf("%f\n",qx(2002,3));
  }
};

struct applyViscousTerm2dv {
  FS4D dvar;
  FS4D var;
  FS2D rho;
  FS3D vel;
  FS2D qx;
  FS2D qy;
  FS4D stressx;
  FS4D stressy;
  Kokkos::View<FSCAL *> cd;

  applyViscousTerm2dv(FS4D dvar_, FS4D var_, FS2D rho_, FS3D vel_, FS4D strx_, FS4D stry_,
                      FS2D qx_, FS2D qy_, Kokkos::View<FSCAL *> cd_)
      : dvar(dvar_), var(var_), rho(rho_), vel(vel_), stressx(strx_), stressy(stry_),
        qx(qx_), qy(qy_), cd(cd_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, const int j) const {
    FSCAL dx = cd(1);
    FSCAL dy = cd(2);
    FSCAL a, b, c1, c2;

    FSCAL ur = (vel(i+1,j,0)+vel(i,j,0))/2.0;
    FSCAL ul = (vel(i-1,j,0)+vel(i,j,0))/2.0;
    FSCAL vr = (vel(i+1,j,1)+vel(i,j,1))/2.0;
    FSCAL vl = (vel(i-1,j,1)+vel(i,j,1))/2.0;

    FSCAL ut = (vel(i,j+1,0)+vel(i,j,0))/2.0;
    FSCAL ub = (vel(i,j-1,0)+vel(i,j,0))/2.0;
    FSCAL vt = (vel(i,j+1,1)+vel(i,j,1))/2.0;
    FSCAL vb = (vel(i,j-1,1)+vel(i,j,1))/2.0;

    a = (stressx(i,j,0,0) - stressx(i-1,j,0,0)) / dx +
        (stressy(i,j,0,1) - stressy(i,j-1,0,1)) / dy;

    b = (stressx(i,j,1,0) - stressx(i-1,j,1,0)) / dx +
        (stressy(i,j,1,1) - stressy(i,j-1,1,1)) / dy;

    c1 = (ur*stressx(i,j,0,0) - ul*stressx(i-1,j,0,0)) / dx +
         (vr*stressx(i,j,0,1) - vl*stressx(i-1,j,0,1)) / dx +

         (ut*stressy(i,j,1,0) - ub*stressy(i,j-1,1,0)) / dy +
         (vt*stressy(i,j,1,1) - vb*stressy(i,j-1,1,1)) / dy;

    c2 = ( qx(i,j)-qx(i-1,j) )/dx + ( qy(i,j)-qy(i,j-1) )/dy;

    // if (i==2002 && j==3)
    //    printf("div %f: %f, %f\n",a1,stressx(i-1,j,0,0),stressx(i,j,0,0));
    dvar(i, j, 0, 0) = dvar(i, j, 0, 0) + a;
    dvar(i, j, 0, 1) = dvar(i, j, 0, 1) + b;
    dvar(i, j, 0, 2) = dvar(i, j, 0, 2) + c1 - c2;

    //printf("(%2d,%2d) %.1e, %.1e, %.1e == %.1e, %.1e, %.1e, \n",i,j,
    //    stressx(i-1,j,0,0),stressx(i-1,j,1,1),stressx(i-1,j,0,1),
    //    stressy(i,j-1,0,0),stressy(i,j-1,1,1),stressy(i,j-1,0,1));

    // printf("(%2d,%2d) %.1e, %.1e, %.1e, %.1e\n",i,j,a,b,c1,c2);

    // if (i == 2000)
    //    printf("%f, %f, %f\n",a1,a2,a3);
  }
};

struct calculateStressTensor3dv {
  FS4D var;
  FS3D rho;
  FS5D stressx;
  FS5D stressy;
  FS5D stressz;
  FS4D vel;
  Kokkos::View<FSCAL *> cd;

  calculateStressTensor3dv(FS4D var_, FS3D rho_, FS4D v_, FS5D strx_,
                           FS5D stry_, FS5D strz_, Kokkos::View<FSCAL *> cd_)
      : var(var_), rho(rho_), vel(v_), stressx(strx_), stressy(stry_), stressz(strz_), cd(cd_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, const int j, const int k) const {
    int ns = (int)cd(0);
    FSCAL dx = cd(1);
    FSCAL dy = cd(2);
    FSCAL dz = cd(3);
    FSCAL dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz;

    FSCAL mu = 2.928e-5;

    // xface
    dudx = (vel(i+1,j,k,0) -vel(i,j,k,0))/dx;
    dvdx = (vel(i+1,j,k,1) -vel(i,j,k,1))/dx;
    dwdx = (vel(i+1,j,k,2) -vel(i,j,k,2))/dx;

    dudy = ( (vel(i+1,j+1,k,0)+vel(i,j+1,k,0)) - (vel(i+1,j-1,k,0)+vel(i,j-1,k,0)) ) / (4*dy);
    dvdy = ( (vel(i+1,j+1,k,1)+vel(i,j+1,k,1)) - (vel(i+1,j-1,k,1)+vel(i,j-1,k,1)) ) / (4*dy);
    dwdy = ( (vel(i+1,j+1,k,2)+vel(i,j+1,k,2)) - (vel(i+1,j-1,k,2)+vel(i,j-1,k,2)) ) / (4*dy);

    dudz = ( (vel(i+1,j,k+1,0)+vel(i,j,k+1,0)) - (vel(i+1,j,k-1,0)+vel(i,j,k-1,0)) ) / (4*dz);
    dvdz = ( (vel(i+1,j,k+1,1)+vel(i,j,k+1,1)) - (vel(i+1,j,k-1,1)+vel(i,j,k-1,1)) ) / (4*dz);
    dwdz = ( (vel(i+1,j,k+1,2)+vel(i,j,k+1,2)) - (vel(i+1,j,k-1,2)+vel(i,j,k-1,2)) ) / (4*dz);

    stressx(i,j,k,0,0) = (2.0 / 3.0) * mu * (2.0 * dudx - dvdy - dwdz);
    stressx(i,j,k,1,1) = (2.0 / 3.0) * mu * (2.0 * dvdy - dudx - dwdz);
    stressx(i,j,k,2,2) = (2.0 / 3.0) * mu * (2.0 * dwdz - dudx - dvdy);
    stressx(i,j,k,0,1) = mu*(dudy + dvdx);
    stressx(i,j,k,0,2) = mu*(dudz + dwdx);
    stressx(i,j,k,1,2) = mu*(dvdz + dwdy);
    stressx(i,j,k,1,0) = stressx(i,j,k,0,1);
    stressx(i,j,k,2,0) = stressx(i,j,k,0,1);
    stressx(i,j,k,2,1) = stressx(i,j,k,1,2);

    // yface
    dudx = ( (vel(i+1,j+1,k,0)+vel(i+1,j,k,0)) - (vel(i-1,j+1,k,0)+vel(i-1,j,k,0)) ) / (4*dx);
    dvdx = ( (vel(i+1,j+1,k,1)+vel(i+1,j,k,1)) - (vel(i-1,j+1,k,1)+vel(i-1,j,k,1)) ) / (4*dx);
    dwdx = ( (vel(i+1,j+1,k,2)+vel(i+1,j,k,2)) - (vel(i-1,j+1,k,2)+vel(i-1,j,k,2)) ) / (4*dx);

    dudy = (vel(i,j+1,k,0) -vel(i,j,k,0))/dy;
    dvdy = (vel(i,j+1,k,1) -vel(i,j,k,1))/dy;
    dwdy = (vel(i,j+1,k,2) -vel(i,j,k,2))/dy;

    dudz = ( (vel(i,j+1,k+1,0)+vel(i,j,k+1,0)) - (vel(i,j+1,k-1,0)+vel(i,j,k-1,0)) ) / (4*dz);
    dvdz = ( (vel(i,j+1,k+1,1)+vel(i,j,k+1,1)) - (vel(i,j+1,k-1,1)+vel(i,j,k-1,1)) ) / (4*dz);
    dwdz = ( (vel(i,j+1,k+1,2)+vel(i,j,k+1,2)) - (vel(i,j+1,k-1,2)+vel(i,j,k-1,2)) ) / (4*dz);

    stressy(i,j,k,0,0) = (2.0 / 3.0) * mu * (2.0 * dudx - dvdy - dwdz);
    stressy(i,j,k,1,1) = (2.0 / 3.0) * mu * (2.0 * dvdy - dudx - dwdz);
    stressy(i,j,k,2,2) = (2.0 / 3.0) * mu * (2.0 * dwdz - dudx - dvdy);
    stressy(i,j,k,0,1) = mu*(dudy + dvdx);
    stressy(i,j,k,0,2) = mu*(dudz + dwdx);
    stressy(i,j,k,1,2) = mu*(dvdz + dwdy);
    stressy(i,j,k,1,0) = stressy(i,j,k,0,1);
    stressy(i,j,k,2,0) = stressy(i,j,k,0,1);
    stressy(i,j,k,2,1) = stressy(i,j,k,1,2);

    // zface
    dudx = ( (vel(i+1,j,k+1,0)+vel(i+1,j,k,0)) - (vel(i-1,j,k+1,0)+vel(i-1,j,k,0)) ) / (4*dx);
    dvdx = ( (vel(i+1,j,k+1,1)+vel(i+1,j,k,1)) - (vel(i-1,j,k+1,1)+vel(i-1,j,k,1)) ) / (4*dx);
    dwdx = ( (vel(i+1,j,k+1,2)+vel(i+1,j,k,2)) - (vel(i-1,j,k+1,2)+vel(i-1,j,k,2)) ) / (4*dx);

    dudy = ( (vel(i,j+1,k+1,0)+vel(i,j+1,k,0)) - (vel(i,j-1,k+1,0)+vel(i,j-1,k,0)) ) / (4*dy);
    dvdy = ( (vel(i,j+1,k+1,1)+vel(i,j+1,k,1)) - (vel(i,j-1,k+1,1)+vel(i,j-1,k,1)) ) / (4*dy);
    dwdy = ( (vel(i,j+1,k+1,2)+vel(i,j+1,k,2)) - (vel(i,j-1,k+1,2)+vel(i,j-1,k,2)) ) / (4*dy);

    dudz = (vel(i,j,k+1,0) -vel(i,j,k,0))/dz;
    dvdz = (vel(i,j,k+1,1) -vel(i,j,k,1))/dz;
    dwdz = (vel(i,j,k+1,2) -vel(i,j,k,2))/dz;
    
    stressz(i,j,k,0,0) = (2.0 / 3.0) * mu * (2.0 * dudx - dvdy - dwdz);
    stressz(i,j,k,1,1) = (2.0 / 3.0) * mu * (2.0 * dvdy - dudx - dwdz);
    stressz(i,j,k,2,2) = (2.0 / 3.0) * mu * (2.0 * dwdz - dudx - dvdy);
    stressz(i,j,k,0,1) = mu*(dudy + dvdx);
    stressz(i,j,k,0,2) = mu*(dudz + dwdx);
    stressz(i,j,k,1,2) = mu*(dvdz + dwdy);
    stressz(i,j,k,1,0) = stressz(i,j,k,0,1);
    stressz(i,j,k,2,0) = stressz(i,j,k,0,1);
    stressz(i,j,k,2,1) = stressz(i,j,k,1,2);
  }
};
struct calculateHeatFlux3dv {
  FS4D var;
  FS3D rho;
  FS3D qx;
  FS3D qy;
  FS3D qz;
  FS3D T;
  Kokkos::View<FSCAL *> cd;

  FSCAL k = 0.026;

  calculateHeatFlux3dv(FS4D var_, FS3D rho_, FS3D T_, FS3D qx_, FS3D qy_, FS3D qz_,
                       Kokkos::View<FSCAL *> cd_)
      : var(var_), rho(rho_), T(T_), qx(qx_), qy(qy_), qz(qy_), cd(cd_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, const int j, const int k) const {
    FSCAL dx = cd(1);
    FSCAL dy = cd(2);
    FSCAL dz = cd(3);

    qx(i,j,k) = -k*(T(i+1,j,k) - T(i,j,k)) /dx;
    qy(i,j,k) = -k*(T(i,j+1,k) - T(i,j,k)) /dy;
    qz(i,j,k) = -k*(T(i,j,k+1) - T(i,j,k)) /dz;
  }
};

struct applyViscousTerm3dv {
  FS4D dvar;
  FS4D var;
  FS3D rho;
  FS4D vel;
  FS3D qx;
  FS3D qy;
  FS3D qz;
  FS5D stressx;
  FS5D stressy;
  FS5D stressz;
  Kokkos::View<FSCAL *> cd;

  applyViscousTerm3dv(FS4D dvar_, FS4D var_, FS3D rho_, FS4D vel_, FS5D strx_, FS5D stry_, FS5D strz_,
                      FS3D qx_, FS3D qy_, FS3D qz_, Kokkos::View<FSCAL *> cd_)
      : dvar(dvar_), var(var_), rho(rho_), vel(vel_), stressx(strx_), stressy(stry_), stressz(strz_),
        qx(qx_), qy(qy_), qz(qz_), cd(cd_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, const int j, const int k) const {
    FSCAL dx = cd(1);
    FSCAL dy = cd(2);
    FSCAL dz = cd(3);
    FSCAL a, b, c, d1, d2;

    FSCAL ur = (vel(i+1,j,k,0)+vel(i,j,k,0))/2.0;
    FSCAL ul = (vel(i-1,j,k,0)+vel(i,j,k,0))/2.0;
    FSCAL vr = (vel(i+1,j,k,1)+vel(i,j,k,1))/2.0;
    FSCAL vl = (vel(i-1,j,k,1)+vel(i,j,k,1))/2.0;
    FSCAL wr = (vel(i+1,j,k,2)+vel(i,j,k,2))/2.0;
    FSCAL wl = (vel(i-1,j,k,2)+vel(i,j,k,2))/2.0;

    FSCAL ut = (vel(i,j+1,k,0)+vel(i,j,k,0))/2.0;
    FSCAL ub = (vel(i,j-1,k,0)+vel(i,j,k,0))/2.0;
    FSCAL vt = (vel(i,j+1,k,1)+vel(i,j,k,1))/2.0;
    FSCAL vb = (vel(i,j-1,k,1)+vel(i,j,k,1))/2.0;
    FSCAL wt = (vel(i,j+1,k,2)+vel(i,j,k,2))/2.0;
    FSCAL wb = (vel(i,j-1,k,2)+vel(i,j,k,2))/2.0;

    FSCAL uf = (vel(i,j,k+1,0)+vel(i,j,k,0))/2.0;
    FSCAL uh = (vel(i,j,k-1,0)+vel(i,j,k,0))/2.0;
    FSCAL vf = (vel(i,j,k+1,1)+vel(i,j,k,1))/2.0;
    FSCAL vh = (vel(i,j,k-1,1)+vel(i,j,k,1))/2.0;
    FSCAL wf = (vel(i,j,k+1,2)+vel(i,j,k,2))/2.0;
    FSCAL wh = (vel(i,j,k-1,2)+vel(i,j,k,2))/2.0;

    a = (stressx(i,j,k,0,0) - stressx(i-1,j,k,0,0)) / dx +
        (stressy(i,j,k,1,0) - stressy(i,j-1,k,1,0)) / dy +
        (stressz(i,j,k,2,0) - stressz(i,j,k-1,2,0)) / dz;

    b = (stressx(i,j,k,0,1) - stressx(i-1,j,k,0,1)) / dx +
        (stressy(i,j,k,1,1) - stressy(i,j-1,k,1,1)) / dy +
        (stressz(i,j,k,2,1) - stressz(i,j,k-1,2,1)) / dz;

    c = (stressx(i,j,k,0,2) - stressx(i-1,j,k,0,2)) / dx +
        (stressy(i,j,k,1,2) - stressy(i,j-1,k,1,2)) / dy +
        (stressz(i,j,k,2,2) - stressz(i,j,k-1,2,2)) / dz;

    d1 = (ur*stressx(i,j,k,0,0) - ul*stressx(i-1,j,k,0,0)) / dx +
         (vr*stressx(i,j,k,0,1) - vl*stressx(i-1,j,k,0,1)) / dx +
         (wr*stressx(i,j,k,0,2) - wl*stressx(i-1,j,k,0,2)) / dx +

         (ut*stressy(i,j,k,1,0) - ub*stressy(i,j-1,k,1,0)) / dy +
         (vt*stressy(i,j,k,1,1) - vb*stressy(i,j-1,k,1,1)) / dy +
         (wt*stressy(i,j,k,1,2) - wb*stressy(i,j-1,k,1,2)) / dy +

         (uf*stressz(i,j,k,2,0) - uh*stressz(i,j,k-1,2,0)) / dz +
         (vf*stressz(i,j,k,2,1) - vh*stressz(i,j,k-1,2,1)) / dz +
         (wf*stressz(i,j,k,2,2) - wh*stressz(i,j,k-1,2,2)) / dz;

    d2 = ( qx(i,j,k)-qx(i-1,j,k) )/dx +
         ( qy(i,j,k)-qy(i,j-1,k) )/dy +
         ( qz(i,j,k)-qz(i,j,k-1) )/dz;

    dvar(i,j,k,0) = dvar(i,j,k,0) + a;
    dvar(i,j,k,1) = dvar(i,j,k,1) + b;
    dvar(i,j,k,2) = dvar(i,j,k,2) + c;
    dvar(i,j,k,3) = dvar(i,j,k,3) + d1 - d2;
  }
};
#endif
