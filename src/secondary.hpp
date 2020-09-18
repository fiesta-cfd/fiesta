#ifndef SECONDARY_H
#define SECONDARY_H

struct calculateRhoPT2D {
  //Kokkos::View<double ****> var;
  //Kokkos::View<double **> rho,p,T;
  //Kokkos::View<double *> cd;
  FS4D var;
  FS2D rho, p, T;
  FS1D cd;

  calculateRhoPT2D(FS4D var_, FS2D p_, FS2D rho_, FS2D T_, FS1D cd_)
      : var(var_), p(p_), rho(rho_), T(T_), cd(cd_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, const int j) const {

    int ns = (int)cd(0);
    int nv = (int)cd(4);
    double gamma, gammas, Rs;
    double Cp = 0;
    double Cv = 0;

    rho(i, j) = 0.0;

    // Total Density for this cell
    for (int s = 0; s < ns; ++s) {
      rho(i, j) = rho(i, j) + var(i, j, 0, 3 + s);
    }

    // Calculate mixture ratio of specific heats
    for (int s = 0; s < ns; ++s) {
      gammas = cd(6 + 3 * s);
      Rs = cd(6 + 3 * s + 1);

      // accumulate mixture heat capacity by mass fraction weights
      Cp =
          Cp + (var(i, j, 0, 3 + s) / rho(i, j)) * (gammas * Rs / (gammas - 1));
      Cv = Cv + (var(i, j, 0, 3 + s) / rho(i, j)) * (Rs / (gammas - 1));
    }
    gamma = Cp / Cv;

    // calculate pressure assuming perfect gas
    p(i, j) =
        (gamma - 1) * (var(i, j, 0, 2) -
                       (0.5 / rho(i, j)) * (var(i, j, 0, 0) * var(i, j, 0, 0) +
                                            var(i, j, 0, 1) * var(i, j, 0, 1)));

    // calcualte cell temperature
    T(i, j) = p(i, j) / ((Cp - Cv) * rho(i, j));
  }
};

struct calculateRhoPT3D {
  FS4D var;
  FS3D p;
  FS3D rho;
  FS1D cd;

  calculateRhoPT3D(FS4D var_, FS3D p_, FS3D rho_, FS1D cd_)
      : var(var_), p(p_), rho(rho_), cd(cd_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, const int j, const int k) const {

    int ns = (int)cd(0);
    int nv = (int)cd(4);
    double gamma, gammas, Rs;
    double Cp = 0;
    double Cv = 0;

    rho(i, j, k) = 0.0;

    for (int s = 0; s < ns; ++s) {
      rho(i, j, k) = rho(i, j, k) + var(i, j, k, 4 + s);
    }

    for (int s = 0; s < ns; ++s) {
      gammas = cd(6 + 3 * s);
      Rs = cd(6 + 3 * s + 1);

      Cp = Cp +
           (var(i, j, k, 4 + s) / rho(i, j, k)) * (gammas * Rs / (gammas - 1));
      Cv = Cv + (var(i, j, k, 4 + s) / rho(i, j, k)) * (Rs / (gammas - 1));
    }

    gamma = Cp / Cv;

    p(i, j, k) = (gamma - 1) *
                 (var(i, j, k, 3) -
                  (0.5 / rho(i, j, k)) * (var(i, j, k, 0) * var(i, j, k, 0) +
                                          var(i, j, k, 1) * var(i, j, k, 1) +
                                          var(i, j, k, 2) * var(i, j, k, 2)));
  }
};

#endif
