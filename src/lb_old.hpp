#ifndef SRC_LB_OLD
#define SRC_LB_OLD

#include "utils/scalar.hpp"
#include "utils/sqr.hpp"


namespace LB_OLD {

/** Data structure holding the parameters for the Lattice Boltzmann system. */
struct LB_Parameters {
  /** number density (LJ units) */
  double rho;

  /** kinematic viscosity (LJ units) */
  double viscosity;

  /** bulk viscosity (LJ units) */
  double bulk_viscosity;

  /** lattice spacing (LJ units) */
  double agrid;

  /** time step for fluid propagation (LJ units)
   *  Note: Has to be larger than MD time step! */
  double tau;

  /** friction coefficient for viscous coupling (LJ units)
   * Note that the friction coefficient is quite high and may
   * lead to numerical artifacts with low order integrators */
  double friction;

  /** external force density applied to the fluid at each lattice site (MD
   * units) */
  double ext_force_density[3]; /* Open question: Do we want a local force or
                          global force? */
  double rho_lb_units;
  /** relaxation of the odd kinetic modes */
  double gamma_odd;
  /** relaxation of the even kinetic modes */
  double gamma_even;
  /** relaxation rate of shear modes */
  double gamma_shear;
  /** relaxation rate of bulk modes */
  double gamma_bulk;

  /** Flag determining whether lbpar.gamma_shear, gamma_odd, and gamma_even are
   * calculated
   *  from lbpar.gamma_shear in such a way to yield a TRT LB with minimized slip
   * at
   *  bounce-back boundaries */
  bool is_TRT;

  int resend_halo;

  /** \name Derived parameters */
  /** Flag indicating whether fluctuations are present. */
  int fluct;
  /** amplitudes of the fluctuations of the modes */
  double phi[19];
};

LB_Parameters lbpar = {
    // rho
    1.0,
    // viscosity
    1.0,
    // bulk_viscosity
    1.0,
    // agrid
    1.0,
    // tau
    .01,
    // friction
    0.0,
    // ext_force_density
    {0.0, 0.0, 0.0},
    // rho_lb_units
    0.,
    // gamma_odd
    0.,
    // gamma_even
    0.,
    // gamma_shear
    0.,
    // gamma_bulk
    0.,
    // is_TRT
    false,
    // resend_halo
    0};

std::array<double, 19>
lb_calc_m_from_n(const std::array<double, 19> &populations) {
  std::array<double, 19> mode{};
  double n0, n1p, n1m, n2p, n2m, n3p, n3m, n4p, n4m, n5p, n5m, n6p, n6m, n7p,
      n7m, n8p, n8m, n9p, n9m;

  n0 = populations[0];
  n1p = populations[1] + populations[2];
  n1m = populations[1] - populations[2];
  n2p = populations[3] + populations[4];
  n2m = populations[3] - populations[4];
  n3p = populations[5] + populations[6];
  n3m = populations[5] - populations[6];
  n4p = populations[7] + populations[8];
  n4m = populations[7] - populations[8];
  n5p = populations[9] + populations[10];
  n5m = populations[9] - populations[10];
  n6p = populations[11] + populations[12];
  n6m = populations[11] - populations[12];
  n7p = populations[13] + populations[14];
  n7m = populations[13] - populations[14];
  n8p = populations[15] + populations[16];
  n8m = populations[15] - populations[16];
  n9p = populations[17] + populations[18];
  n9m = populations[17] - populations[18];

  /* mass mode */
  mode[0] = n0 + n1p + n2p + n3p + n4p + n5p + n6p + n7p + n8p + n9p;

  /* momentum modes */
  mode[1] = n1m + n4m + n5m + n6m + n7m;
  mode[2] = n2m + n4m - n5m + n8m + n9m;
  mode[3] = n3m + n6m - n7m + n8m - n9m;

  /* stress modes */
  mode[4] = -n0 + n4p + n5p + n6p + n7p + n8p + n9p;
  mode[5] = n1p - n2p + n6p + n7p - n8p - n9p;
  mode[6] = n1p + n2p - n6p - n7p - n8p - n9p - 2. * (n3p - n4p - n5p);
  mode[7] = n4p - n5p;
  mode[8] = n6p - n7p;
  mode[9] = n8p - n9p;

  /* kinetic modes */
  mode[10] = -2. * n1m + n4m + n5m + n6m + n7m;
  mode[11] = -2. * n2m + n4m - n5m + n8m + n9m;
  mode[12] = -2. * n3m + n6m - n7m + n8m - n9m;
  mode[13] = n4m + n5m - n6m - n7m;
  mode[14] = n4m - n5m - n8m - n9m;
  mode[15] = n6m - n7m - n8m + n9m;
  mode[16] = n0 + n4p + n5p + n6p + n7p + n8p + n9p - 2. * (n1p + n2p + n3p);
  mode[17] = -n1p + n2p + n6p + n7p - n8p - n9p;
  mode[18] = -n1p - n2p - n6p - n7p - n8p - n9p + 2. * (n3p + n4p + n5p);
  return mode;
}

std::vector<double> force_density{{0.1, 0.2, 0.3}};

inline void lb_relax_modes(double* mode) {
  double rho, j[3], pi_eq[6];

  /* re-construct the real density
   * remember that the populations are stored as differences to their
   * equilibrium value */
  rho = mode[0] + lbpar.rho * lbpar.agrid * lbpar.agrid * lbpar.agrid;

  j[0] = mode[1] + 0.5 * force_density[0];
  j[1] = mode[2] + 0.5 * force_density[1];
  j[2] = mode[3] + 0.5 * force_density[2];

  /* equilibrium part of the stress modes */
  pi_eq[0] = Utils::scalar(j, j) / rho;
  pi_eq[1] = (Utils::sqr(j[0]) - Utils::sqr(j[1])) / rho;
  pi_eq[2] = (Utils::scalar(j, j) - 3.0 * Utils::sqr(j[2])) / rho;
  pi_eq[3] = j[0] * j[1] / rho;
  pi_eq[4] = j[0] * j[2] / rho;
  pi_eq[5] = j[1] * j[2] / rho;

  /* relax the stress modes */
  mode[4] = pi_eq[0] + lbpar.gamma_bulk * (mode[4] - pi_eq[0]);
  mode[5] = pi_eq[1] + lbpar.gamma_shear * (mode[5] - pi_eq[1]);
  mode[6] = pi_eq[2] + lbpar.gamma_shear * (mode[6] - pi_eq[2]);
  mode[7] = pi_eq[3] + lbpar.gamma_shear * (mode[7] - pi_eq[3]);
  mode[8] = pi_eq[4] + lbpar.gamma_shear * (mode[8] - pi_eq[4]);
  mode[9] = pi_eq[5] + lbpar.gamma_shear * (mode[9] - pi_eq[5]);

  /* relax the ghost modes (project them out) */
  /* ghost modes have no equilibrium part due to orthogonality */
  mode[10] = lbpar.gamma_odd * mode[10];
  mode[11] = lbpar.gamma_odd * mode[11];
  mode[12] = lbpar.gamma_odd * mode[12];
  mode[13] = lbpar.gamma_odd * mode[13];
  mode[14] = lbpar.gamma_odd * mode[14];
  mode[15] = lbpar.gamma_odd * mode[15];
  mode[16] = lbpar.gamma_even * mode[16];
  mode[17] = lbpar.gamma_even * mode[17];
  mode[18] = lbpar.gamma_even * mode[18];
}

inline void lb_apply_forces(double *mode) {

  double rho, *f, u[3], C[6];

  f = force_density.data();

  rho = mode[0] + lbpar.rho * lbpar.agrid * lbpar.agrid * lbpar.agrid;

  /* hydrodynamic momentum density is redefined when external forces present */
  u[0] = (mode[1] + 0.5 * f[0]) / rho;
  u[1] = (mode[2] + 0.5 * f[1]) / rho;
  u[2] = (mode[3] + 0.5 * f[2]) / rho;

  C[0] = (1. + lbpar.gamma_bulk) * u[0] * f[0] +
         1. / 3. * (lbpar.gamma_bulk - lbpar.gamma_shear) * Utils::scalar(u, f);
  C[2] = (1. + lbpar.gamma_bulk) * u[1] * f[1] +
         1. / 3. * (lbpar.gamma_bulk - lbpar.gamma_shear) * Utils::scalar(u, f);
  C[5] = (1. + lbpar.gamma_bulk) * u[2] * f[2] +
         1. / 3. * (lbpar.gamma_bulk - lbpar.gamma_shear) * Utils::scalar(u, f);
  C[1] = 1. / 2. * (1. + lbpar.gamma_shear) * (u[0] * f[1] + u[1] * f[0]);
  C[3] = 1. / 2. * (1. + lbpar.gamma_shear) * (u[0] * f[2] + u[2] * f[0]);
  C[4] = 1. / 2. * (1. + lbpar.gamma_shear) * (u[1] * f[2] + u[2] * f[1]);

  /* update momentum modes */
  mode[1] += f[0];
  mode[2] += f[1];
  mode[3] += f[2];

  /* update stress modes */
  mode[4] += C[0] + C[2] + C[5];
  mode[5] += C[0] - C[2];
  mode[6] += C[0] + C[2] - 2. * C[5];
  mode[7] += C[1];
  mode[8] += C[3];
  mode[9] += C[4];
}



}

#endif
