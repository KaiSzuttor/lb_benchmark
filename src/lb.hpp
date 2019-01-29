#ifndef LB_HPP
#define LB_HPP

#include "lb-d3q19.hpp"
#include "utils/scalar.hpp"
#include "utils/sqr.hpp"
#include "utils/matrix_vector_product.hpp"
#include "utils/Span.hpp"

namespace LB {

/** Data structure holding the parameters for the Lattice Boltzmann system. */
struct LB_Parameters {
  /** number density (LB units) */
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

  /** friction coefficient for viscous coupling (LJ units) */
  double friction;

  /** external force density applied to the fluid at each lattice site (MD
   * units) */
  std::array<double, 3> ext_force_density;

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

  /** \name Derived parameters */
  /** Flag indicating whether fluctuations are present. */
  int fluct;
  /** amplitudes of the fluctuations of the modes */
  std::array<double, 19> phi;
};

/** Struct holding the Lattice Boltzmann parameters */
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
    0.01,
    // friction
    1.0,
    // ext_force_density
    {0.0, 0.0, 0.0},
    // gamma_odd
    .1,
    // gamma_even
    .1,
    // gamma_shear
    .1,
    // gamma_bulk
    .1,
    // is_TRT
    false,
    // resend_halo
    0};


template <typename T>
inline std::array<T, 19> normalize_modes(const std::array<T, 19> &modes) {
  auto normalized_modes = modes;
  for (int i = 0; i < modes.size(); i++) {
    normalized_modes[i] /= ::D3Q19::w_k[i];
  }
  return normalized_modes;
}

template <typename T, std::size_t N>
inline std::array<T, N> lb_calc_n_from_m(const std::array<T, N> &modes) {
  auto const normalized_modes = normalize_modes(modes);
  auto ret = Utils::matrix_vector_product<T, N, ::D3Q19::e_ki_transposed>(
      normalized_modes);
  std::transform(ret.begin(), ret.end(), ::D3Q19::w.begin(), ret.begin(),
                 std::multiplies<T>());
  return ret;
}

template <typename T, std::size_t N>
inline std::array<T, N> lb_calc_m_from_n(const std::array<T, N>& n) {
  return Utils::matrix_vector_product<T, N, ::D3Q19::e_ki>(n);
}

/** Calculation of hydrodynamic modes */
inline std::array<double, 19> lb_calc_modes(std::array<Utils::Span<double>, 19>& lbfluid, std::size_t index) {
  std::array<double, 19> n;
  for (int i = 0; i < 19; i++) {
    n[i] = lbfluid[i][index];
  }
  return lb_calc_m_from_n(n);
}


std::vector<double> force_density{{0.1, 0.2, 0.3}};

template <typename T>
inline std::array<T, 19>
lb_relax_modes(const std::array<T, 19> &modes) {
  auto relaxed_modes = modes;
  double rho, j[3], pi_eq[6];

  /* re-construct the real density
   * remember that the populations are stored as differences to their
   * equilibrium value */
  rho = modes[0] + lbpar.rho;

  j[0] = modes[1] + 0.5 * force_density[0];
  j[1] = modes[2] + 0.5 * force_density[1];
  j[2] = modes[3] + 0.5 * force_density[2];

  /* equilibrium part of the stress modes */
  pi_eq[0] = Utils::scalar(j, j) / rho;
  pi_eq[1] = (Utils::sqr(j[0]) - Utils::sqr(j[1])) / rho;
  pi_eq[2] = (Utils::scalar(j, j) - 3.0 * Utils::sqr(j[2])) / rho;
  pi_eq[3] = j[0] * j[1] / rho;
  pi_eq[4] = j[0] * j[2] / rho;
  pi_eq[5] = j[1] * j[2] / rho;

  /* relax the stress modes */
  relaxed_modes[4] = pi_eq[0] + lbpar.gamma_bulk * (modes[4] - pi_eq[0]);
  relaxed_modes[5] = pi_eq[1] + lbpar.gamma_shear * (modes[5] - pi_eq[1]);
  relaxed_modes[6] = pi_eq[2] + lbpar.gamma_shear * (modes[6] - pi_eq[2]);
  relaxed_modes[7] = pi_eq[3] + lbpar.gamma_shear * (modes[7] - pi_eq[3]);
  relaxed_modes[8] = pi_eq[4] + lbpar.gamma_shear * (modes[8] - pi_eq[4]);
  relaxed_modes[9] = pi_eq[5] + lbpar.gamma_shear * (modes[9] - pi_eq[5]);

  /* relax the ghost modes (project them out) */
  /* ghost modes have no equilibrium part due to orthogonality */
  relaxed_modes[10] = lbpar.gamma_odd * modes[10];
  relaxed_modes[11] = lbpar.gamma_odd * modes[11];
  relaxed_modes[12] = lbpar.gamma_odd * modes[12];
  relaxed_modes[13] = lbpar.gamma_odd * modes[13];
  relaxed_modes[14] = lbpar.gamma_odd * modes[14];
  relaxed_modes[15] = lbpar.gamma_odd * modes[15];
  relaxed_modes[16] = lbpar.gamma_even * modes[16];
  relaxed_modes[17] = lbpar.gamma_even * modes[17];
  relaxed_modes[18] = lbpar.gamma_even * modes[18];
  return relaxed_modes;
}

template <typename T>
std::array<T, 19> lb_apply_forces(const std::array<T, 19> &modes) {
  std::array<T, 19> modes_with_forces = modes;
  T rho, u[3], C[6];

  const auto &f = force_density;

  rho = modes[0] + lbpar.rho;

  /* hydrodynamic momentum density is redefined when external forces present */
  u[0] = (modes[1] + 0.5 * f[0]) / rho;
  u[1] = (modes[2] + 0.5 * f[1]) / rho;
  u[2] = (modes[3] + 0.5 * f[2]) / rho;

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
  modes_with_forces[1] += f[0];
  modes_with_forces[2] += f[1];
  modes_with_forces[3] += f[2];

  /* update stress modes */
  modes_with_forces[4] += C[0] + C[2] + C[5];
  modes_with_forces[5] += C[0] - C[2];
  modes_with_forces[6] += C[0] + C[2] - 2. * C[5];
  modes_with_forces[7] += C[1];
  modes_with_forces[8] += C[3];
  modes_with_forces[9] += C[4];
  return modes_with_forces;
}

} // namespace LB

#endif
