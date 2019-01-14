#include <benchmark/benchmark.h>
#include <random>

#include "lb.hpp"
#include "lb_old.hpp"

std::mt19937 mt(10);
std::uniform_real_distribution<double> dist(1.0, 10.0);


static void BM_LB_mode_transformation_Old(benchmark::State &state) {
  for (auto _ : state) {
    state.PauseTiming();
    std::array<double, 19> n;
    std::generate(n.begin(), n.end(), []() { return dist(mt); });
    state.ResumeTiming();
    benchmark::DoNotOptimize(LB_OLD::lb_calc_m_from_n(n));
  }
}

static void BM_LB_mode_transformation_New(benchmark::State &state) {
  for (auto _ : state) {
    state.PauseTiming();
    std::array<double, 19> n;
    std::generate(n.begin(), n.end(), []() { return dist(mt); });
    state.ResumeTiming();
    benchmark::DoNotOptimize(LB::lb_calc_m_from_n(n));
  }
}

BENCHMARK(BM_LB_mode_transformation_Old);
BENCHMARK(BM_LB_mode_transformation_New);
BENCHMARK_MAIN();
