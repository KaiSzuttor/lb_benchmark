#include <benchmark/benchmark.h>
#include <random>

#include "lb.hpp"
#include "lb_old.hpp"

std::mt19937 mt(10);
std::uniform_real_distribution<double> dist(1.0, 10.0);

static void BM_LB_mode_relaxation_old(benchmark::State &state) {
  for (auto _ : state) {
    state.PauseTiming();
    std::array<double, 19> n;
    std::generate(n.begin(), n.end(), []() { return dist(mt); });
    state.ResumeTiming();
    LB_OLD::lb_relax_modes(n.data());
    benchmark::DoNotOptimize(n);
  }
}

static void BM_LB_mode_relaxation_new(benchmark::State &state) {
  for (auto _ : state) {
    state.PauseTiming();
    std::array<double, 19> n;
    std::generate(n.begin(), n.end(), []() { return dist(mt); });
    state.ResumeTiming();
    benchmark::DoNotOptimize(LB::lb_relax_modes(n));
  }
}

BENCHMARK(BM_LB_mode_relaxation_old);
BENCHMARK(BM_LB_mode_relaxation_new);
BENCHMARK_MAIN();
