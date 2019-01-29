#include <benchmark/benchmark.h>
#include <random>
#include <boost/multi_array.hpp>
#include <boost/range/numeric.hpp>

#include "lb.hpp"
#include "lb_old.hpp"
#include "lb-d3q19.hpp"
#include "utils/Span.hpp"

namespace LB {

struct Lattice {
  std::array<int, 3> grid;
  std::array<int, 3> halo_grid;
  int halo_grid_volume;
  int halo_offset;

  Lattice(double agrid, double box_l) : grid({{static_cast<int>(box_l/agrid),
                                              static_cast<int>(box_l/agrid),
                                              static_cast<int>(box_l/agrid)}}),
                                        halo_grid({{grid[0] + 2, grid[1] + 2, grid[2] + 2}}),
                                        halo_grid_volume(halo_grid[0] * halo_grid[1] * halo_grid[2]),
                                        halo_offset(1 + halo_grid[0] * (1 + halo_grid[1] * 1)) {};
};


Lattice lblattice(0.5, 10);


using LB_FluidData = boost::multi_array<double, 2>;
static LB_FluidData lbfluid_a;
static LB_FluidData lbfluid_b;
using LB_Fluid = std::array<Utils::Span<double>, 19>;
LB_Fluid lbfluid;
LB_Fluid lbfluid_post;

template <typename T>
inline std::array<T, 19> lb_relax_modes(std::size_t index,
                                        const std::array<T, 19> &modes) {
  T rho, j[3], pi_eq[6];

  /* re-construct the real density
   * remember that the populations are stored as differences to their
   * equilibrium value */
  rho = modes[0] + lbpar.rho;

  j[0] = modes[1];
  j[1] = modes[2];
  j[2] = modes[3];

  /* equilibrium part of the stress modes */
  pi_eq[0] = Utils::scalar(j, j) / rho;
  pi_eq[1] = (Utils::sqr(j[0]) - Utils::sqr(j[1])) / rho;
  pi_eq[2] = (Utils::scalar(j, j) - 3.0 * Utils::sqr(j[2])) / rho;
  pi_eq[3] = j[0] * j[1] / rho;
  pi_eq[4] = j[0] * j[2] / rho;
  pi_eq[5] = j[1] * j[2] / rho;

  return {{modes[0], modes[1], modes[2], modes[3],
           /* relax the stress modes */
           pi_eq[0] + lbpar.gamma_bulk * (modes[4] - pi_eq[0]),
           pi_eq[1] + lbpar.gamma_shear * (modes[5] - pi_eq[1]),
           pi_eq[2] + lbpar.gamma_shear * (modes[6] - pi_eq[2]),
           pi_eq[3] + lbpar.gamma_shear * (modes[7] - pi_eq[3]),
           pi_eq[4] + lbpar.gamma_shear * (modes[8] - pi_eq[4]),
           pi_eq[5] + lbpar.gamma_shear * (modes[9] - pi_eq[5]),
           /* relax the ghost modes (project them out) */
           /* ghost modes have no equilibrium part due to orthogonality */
           lbpar.gamma_odd * modes[10], lbpar.gamma_odd * modes[11],
           lbpar.gamma_odd * modes[12], lbpar.gamma_odd * modes[13],
           lbpar.gamma_odd * modes[14], lbpar.gamma_odd * modes[15],
           lbpar.gamma_even * modes[16], lbpar.gamma_even * modes[17],
           lbpar.gamma_even * modes[18]}};
}


template <typename T>
std::array<T, 19> lb_apply_forces(std::size_t index,
                                  const std::array<T, 19> &modes) {
  T rho, u[3], C[6];

  const auto &f = std::array<T, 3>{{0.1, 0.0, 0.0}};

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

  return {{modes[0],
           /* update momentum modes */
           modes[1] + f[0], modes[2] + f[1], modes[3] + f[2],
           /* update stress modes */
           modes[4] + C[0] + C[2] + C[5], modes[5] + C[0] - C[2],
           modes[6] + C[0] + C[2] - 2. * C[5], modes[7] + C[1], modes[8] + C[3],
           modes[9] + C[4], modes[10], modes[11], modes[12], modes[13],
           modes[14], modes[15], modes[16], modes[17], modes[18]}};
}


template <typename T>
inline void lb_calc_n_from_modes_push(LB_Fluid &lbfluid, std::size_t index,
                                      const std::array<T, 19> m) {
  const std::array<int, 3> period = {
      {1, lblattice.halo_grid[0],
       lblattice.halo_grid[0] * lblattice.halo_grid[1]}};
  auto const f = lb_calc_n_from_m(m);
  for (int i = 0; i < 19; i++) {
    auto const next = index + boost::inner_product(period, D3Q19::c[i], 0);
    lbfluid[i][next] = f[i];
  }
}

}


static void BM_LB_hot_loop_new(benchmark::State &state) {
  const std::array<int, 2> size = {19, LB::lblattice.halo_grid_volume};
  LB::lbpar.gamma_shear = 1. - 2. / (6. * LB::lbpar.viscosity + 1.);
  LB::lbpar.gamma_bulk = 1. - 2. / (9. * LB::lbpar.bulk_viscosity + 1.);
  LB::lbfluid_a.resize(size);
  LB::lbfluid_b.resize(size);
  for (int i=0; i < size[0]; i++) {
      LB::lbfluid[i] = Utils::Span<double>(LB::lbfluid_a[i].origin(), size[1]);
      LB::lbfluid_post[i] = Utils::Span<double>(LB::lbfluid_b[i].origin(), size[1]);
  }
  std::array<double, 3> j = {{0., 0., 0.}};
  std::array<double, 6> pi = {{0., 0., 0., 0., 0., 0.}};
  for (std::size_t index = 0; index < LB::lblattice.halo_grid_volume; ++index) {
      LB::lbfluid[0][index] = 0.0;
      LB::lbfluid[1][index] = 0.0;
      LB::lbfluid[2][index] = 0.0;
      LB::lbfluid[3][index] = 0.0;
      LB::lbfluid[4][index] = 0.0;
      LB::lbfluid[5][index] = 0.0;
      LB::lbfluid[6][index] = 0.0;
      LB::lbfluid[7][index] = 0.0;
      LB::lbfluid[8][index] = 0.0;
      LB::lbfluid[9][index] = 0.0;
      LB::lbfluid[10][index] = 0.0;
      LB::lbfluid[11][index] = 0.0;
      LB::lbfluid[12][index] = 0.0;
      LB::lbfluid[13][index] = 0.0;
      LB::lbfluid[14][index] = 0.0;
      LB::lbfluid[15][index] = 0.0;
      LB::lbfluid[16][index] = 0.0;
      LB::lbfluid[17][index] = 0.0;
      LB::lbfluid[18][index] = 0.0;
  }
  for (auto _ : state) {
    std::size_t index = LB::lblattice.halo_offset;
    for (int z = 1; z <= LB::lblattice.grid[2]; z++) {
      for (int y = 1; y <= LB::lblattice.grid[1]; y++) {
        for (int x = 1; x <= LB::lblattice.grid[0]; x++) {
          auto const modes = LB::lb_calc_modes(LB::lbfluid, index);
          auto const relaxed_modes = LB::lb_relax_modes(index, modes);
          auto const modes_with_forces = LB::lb_apply_forces(index, relaxed_modes);
          LB::lb_calc_n_from_modes_push(LB::lbfluid_post, index, modes);
          ++index;
        }
        index += 2;
      }
      index += 2 * LB::lblattice.halo_grid[0];
    }
  }
}

namespace LB_Old {
/** Calculation of hydrodynamic modes */
void lb_calc_modes(std::size_t index, double *mode) {
  double n0, n1p, n1m, n2p, n2m, n3p, n3m, n4p, n4m, n5p, n5m, n6p, n6m, n7p,
      n7m, n8p, n8m, n9p, n9m;
  n0 = LB::lbfluid[0][index];
  n1p = LB::lbfluid[1][index] + LB::lbfluid[2][index];
  n1m = LB::lbfluid[1][index] - LB::lbfluid[2][index];
  n2p = LB::lbfluid[3][index] + LB::lbfluid[4][index];
  n2m = LB::lbfluid[3][index] - LB::lbfluid[4][index];
  n3p = LB::lbfluid[5][index] + LB::lbfluid[6][index];
  n3m = LB::lbfluid[5][index] - LB::lbfluid[6][index];
  n4p = LB::lbfluid[7][index] + LB::lbfluid[8][index];
  n4m = LB::lbfluid[7][index] - LB::lbfluid[8][index];
  n5p = LB::lbfluid[9][index] + LB::lbfluid[10][index];
  n5m = LB::lbfluid[9][index] - LB::lbfluid[10][index];
  n6p = LB::lbfluid[11][index] + LB::lbfluid[12][index];
  n6m = LB::lbfluid[11][index] - LB::lbfluid[12][index];
  n7p = LB::lbfluid[13][index] + LB::lbfluid[14][index];
  n7m = LB::lbfluid[13][index] - LB::lbfluid[14][index];
  n8p = LB::lbfluid[15][index] + LB::lbfluid[16][index];
  n8m = LB::lbfluid[15][index] - LB::lbfluid[16][index];
  n9p = LB::lbfluid[17][index] + LB::lbfluid[18][index];
  n9m = LB::lbfluid[17][index] - LB::lbfluid[18][index];

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
}

inline void lb_relax_modes(std::size_t index, double *mode) {
  double rho, j[3], pi_eq[6];
  using namespace LB;
  /* re-construct the real density
   * remember that the populations are stored as differences to their
   * equilibrium value */
  rho = mode[0] + lbpar.rho * lbpar.agrid * lbpar.agrid * lbpar.agrid;

  j[0] = mode[1];
  j[1] = mode[2];
  j[2] = mode[3];

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

inline void lb_apply_forces(std::size_t index, double *mode) {
  using namespace LB;
  double rho, u[3], C[6];

  const auto &f = std::array<double, 3>{{0.1, 0.0, 0.0}};

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

inline void lb_calc_n_from_modes_push(LB::LB_Fluid &lbfluid, std::size_t index, double *m) {
  using namespace LB;
  int yperiod = lblattice.halo_grid[0];
  int zperiod = lblattice.halo_grid[0] * lblattice.halo_grid[1];
  std::size_t next[19];
  next[0] = index;
  next[1] = index + 1;
  next[2] = index - 1;
  next[3] = index + yperiod;
  next[4] = index - yperiod;
  next[5] = index + zperiod;
  next[6] = index - zperiod;
  next[7] = index + (1 + yperiod);
  next[8] = index - (1 + yperiod);
  next[9] = index + (1 - yperiod);
  next[10] = index - (1 - yperiod);
  next[11] = index + (1 + zperiod);
  next[12] = index - (1 + zperiod);
  next[13] = index + (1 - zperiod);
  next[14] = index - (1 - zperiod);
  next[15] = index + (yperiod + zperiod);
  next[16] = index - (yperiod + zperiod);
  next[17] = index + (yperiod - zperiod);
  next[18] = index - (yperiod - zperiod);

  /* normalization factors enter in the back transformation */
  for (int i = 0; i < 19; i++)
    m[i] = (1. / D3Q19::w[i]) * m[i];

  lbfluid[0][next[0]] = m[0] - m[4] + m[16];
  lbfluid[1][next[1]] =
      m[0] + m[1] + m[5] + m[6] - m[17] - m[18] - 2. * (m[10] + m[16]);
  lbfluid[2][next[2]] =
      m[0] - m[1] + m[5] + m[6] - m[17] - m[18] + 2. * (m[10] - m[16]);
  lbfluid[3][next[3]] =
      m[0] + m[2] - m[5] + m[6] + m[17] - m[18] - 2. * (m[11] + m[16]);
  lbfluid[4][next[4]] =
      m[0] - m[2] - m[5] + m[6] + m[17] - m[18] + 2. * (m[11] - m[16]);
  lbfluid[5][next[5]] = m[0] + m[3] - 2. * (m[6] + m[12] + m[16] - m[18]);
  lbfluid[6][next[6]] = m[0] - m[3] - 2. * (m[6] - m[12] + m[16] - m[18]);
  lbfluid[7][next[7]] = m[0] + m[1] + m[2] + m[4] + 2. * m[6] + m[7] + m[10] +
                        m[11] + m[13] + m[14] + m[16] + 2. * m[18];
  lbfluid[8][next[8]] = m[0] - m[1] - m[2] + m[4] + 2. * m[6] + m[7] - m[10] -
                        m[11] - m[13] - m[14] + m[16] + 2. * m[18];
  lbfluid[9][next[9]] = m[0] + m[1] - m[2] + m[4] + 2. * m[6] - m[7] + m[10] -
                        m[11] + m[13] - m[14] + m[16] + 2. * m[18];
  lbfluid[10][next[10]] = m[0] - m[1] + m[2] + m[4] + 2. * m[6] - m[7] - m[10] +
                          m[11] - m[13] + m[14] + m[16] + 2. * m[18];
  lbfluid[11][next[11]] = m[0] + m[1] + m[3] + m[4] + m[5] - m[6] + m[8] +
                          m[10] + m[12] - m[13] + m[15] + m[16] + m[17] - m[18];
  lbfluid[12][next[12]] = m[0] - m[1] - m[3] + m[4] + m[5] - m[6] + m[8] -
                          m[10] - m[12] + m[13] - m[15] + m[16] + m[17] - m[18];
  lbfluid[13][next[13]] = m[0] + m[1] - m[3] + m[4] + m[5] - m[6] - m[8] +
                          m[10] - m[12] - m[13] - m[15] + m[16] + m[17] - m[18];
  lbfluid[14][next[14]] = m[0] - m[1] + m[3] + m[4] + m[5] - m[6] - m[8] -
                          m[10] + m[12] + m[13] + m[15] + m[16] + m[17] - m[18];
  lbfluid[15][next[15]] = m[0] + m[2] + m[3] + m[4] - m[5] - m[6] + m[9] +
                          m[11] + m[12] - m[14] - m[15] + m[16] - m[17] - m[18];
  lbfluid[16][next[16]] = m[0] - m[2] - m[3] + m[4] - m[5] - m[6] + m[9] -
                          m[11] - m[12] + m[14] + m[15] + m[16] - m[17] - m[18];
  lbfluid[17][next[17]] = m[0] + m[2] - m[3] + m[4] - m[5] - m[6] - m[9] +
                          m[11] - m[12] - m[14] + m[15] + m[16] - m[17] - m[18];
  lbfluid[18][next[18]] = m[0] - m[2] + m[3] + m[4] - m[5] - m[6] - m[9] -
                          m[11] + m[12] + m[14] - m[15] + m[16] - m[17] - m[18];

  /* weights enter in the back transformation */
  for (int i = 0; i < 19; i++)
    lbfluid[i][next[i]] *= D3Q19::w_k[i];
}

}

static void BM_LB_hot_loop_old(benchmark::State &state) {
  const std::array<int, 2> size = {19, LB::lblattice.halo_grid_volume};
  LB::lbpar.gamma_shear = 1. - 2. / (6. * LB::lbpar.viscosity + 1.);
  LB::lbpar.gamma_bulk = 1. - 2. / (9. * LB::lbpar.bulk_viscosity + 1.);
  LB::lbfluid_a.resize(size);
  LB::lbfluid_b.resize(size);
  for (int i=0; i < size[0]; i++) {
      LB::lbfluid[i] = Utils::Span<double>(LB::lbfluid_a[i].origin(), size[1]);
      LB::lbfluid_post[i] = Utils::Span<double>(LB::lbfluid_b[i].origin(), size[1]);
  }
  std::array<double, 3> j = {{0., 0., 0.}};
  std::array<double, 6> pi = {{0., 0., 0., 0., 0., 0.}};
  for (std::size_t index = 0; index < LB::lblattice.halo_grid_volume; ++index) {
      LB::lbfluid[0][index] = 0.0;
      LB::lbfluid[1][index] = 0.0;
      LB::lbfluid[2][index] = 0.0;
      LB::lbfluid[3][index] = 0.0;
      LB::lbfluid[4][index] = 0.0;
      LB::lbfluid[5][index] = 0.0;
      LB::lbfluid[6][index] = 0.0;
      LB::lbfluid[7][index] = 0.0;
      LB::lbfluid[8][index] = 0.0;
      LB::lbfluid[9][index] = 0.0;
      LB::lbfluid[10][index] = 0.0;
      LB::lbfluid[11][index] = 0.0;
      LB::lbfluid[12][index] = 0.0;
      LB::lbfluid[13][index] = 0.0;
      LB::lbfluid[14][index] = 0.0;
      LB::lbfluid[15][index] = 0.0;
      LB::lbfluid[16][index] = 0.0;
      LB::lbfluid[17][index] = 0.0;
      LB::lbfluid[18][index] = 0.0;
  }
  for (auto _ : state) {
    std::size_t index = LB::lblattice.halo_offset;
    for (int z = 1; z <= LB::lblattice.grid[2]; z++) {
      for (int y = 1; y <= LB::lblattice.grid[1]; y++) {
        for (int x = 1; x <= LB::lblattice.grid[0]; x++) {
          double modes[19];
          LB_Old::lb_calc_modes(index, modes);
          LB_Old::lb_relax_modes(index, modes);
          LB_Old::lb_apply_forces(index, modes);
          LB_Old::lb_calc_n_from_modes_push(LB::lbfluid_post, index, modes);
          ++index;
        }
        index += 2;
      }
      index += 2 * LB::lblattice.halo_grid[0];
    }
  }
}


BENCHMARK(BM_LB_hot_loop_new);
BENCHMARK(BM_LB_hot_loop_old);
BENCHMARK_MAIN();
