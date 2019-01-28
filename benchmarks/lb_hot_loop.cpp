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

  const auto &f = std::array<T, 3>{{0.0, 0.0, 0.0}};

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
  LB::lbpar.fluct = 0;
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
          LB::lb_calc_n_from_modes_push(LB::lbfluid_post, index, modes_with_forces);
          ++index;
        }
        index += 2;
      }
      index += 2 * LB::lblattice.halo_grid[0];
    }
  }
}


BENCHMARK(BM_LB_hot_loop_new);
BENCHMARK_MAIN();
