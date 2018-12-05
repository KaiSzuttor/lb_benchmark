#include <random>
#include <benchmark/benchmark.h>

#include "lb.hpp"

std::mt19937 mt(10);
std::uniform_real_distribution<double> dist(1.0, 10.0);

static void BM_LB_Modes_Calculation_New(benchmark::State& state) {
  for (auto _ : state) {
    state.PauseTiming();
    std::array<double, 19> n;
    std::generate(n.begin(), n.end(), [](){return dist(mt);});
    state.ResumeTiming();
    auto const modes = LB::lb_calc_m_from_n(n);
  }
}

void lb_calc_m_from_n_old(const std::array<double, 19>& populations, double *mode) {
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
}

static void BM_LB_Modes_Calculation_Old(benchmark::State& state) {
  for (auto _ : state) {
    state.PauseTiming();
    std::array<double, 19> n;
    std::array<double, 19> m;
    std::generate(n.begin(), n.end(), [](){return dist(mt);});
    state.ResumeTiming();
    lb_calc_m_from_n_old(n, m.data());
  }
}

BENCHMARK(BM_LB_Modes_Calculation_Old);
BENCHMARK(BM_LB_Modes_Calculation_New);
BENCHMARK_MAIN();
