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

BENCHMARK(BM_LB_Modes_Calculation_New);
BENCHMARK_MAIN();
