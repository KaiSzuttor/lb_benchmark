#ifndef D3Q19_H
#define D3Q19_H

#include <array>

namespace D3Q19 {

/** Velocity sub-lattice of the D3Q19 model */
static constexpr const std::array<std::array<double, 3>, 19> c = {
    {{{0., 0., 0.}},
     {{1., 0., 0.}},
     {{-1., 0., 0.}},
     {{0., 1., 0.}},
     {{0., -1., 0.}},
     {{0., 0., 1.}},
     {{0., 0., -1.}},
     {{1., 1., 0.}},
     {{-1., -1., 0.}},
     {{1., -1., 0.}},
     {{-1., 1., 0.}},
     {{1., 0., 1.}},
     {{-1., 0., -1.}},
     {{1., 0., -1.}},
     {{-1., 0., 1.}},
     {{0., 1., 1.}},
     {{0., -1., -1.}},
     {{0., 1., -1.}},
     {{0., -1., 1.}}}};

/** Coefficients for pseudo-equilibrium distribution of the D3Q19 model */
static constexpr const std::array<std::array<double, 4>, 19> coefficients = {
    {{{1. / 3., 1., 3. / 2., -1. / 2.}},
     {{1. / 18., 1. / 6., 1. / 4., -1. / 12.}},
     {{1. / 18., 1. / 6., 1. / 4., -1. / 12.}},
     {{1. / 18., 1. / 6., 1. / 4., -1. / 12.}},
     {{1. / 18., 1. / 6., 1. / 4., -1. / 12.}},
     {{1. / 18., 1. / 6., 1. / 4., -1. / 12.}},
     {{1. / 18., 1. / 6., 1. / 4., -1. / 12.}},
     {{1. / 36., 1. / 12., 1. / 8., -1. / 24.}},
     {{1. / 36., 1. / 12., 1. / 8., -1. / 24.}},
     {{1. / 36., 1. / 12., 1. / 8., -1. / 24.}},
     {{1. / 36., 1. / 12., 1. / 8., -1. / 24.}},
     {{1. / 36., 1. / 12., 1. / 8., -1. / 24.}},
     {{1. / 36., 1. / 12., 1. / 8., -1. / 24.}},
     {{1. / 36., 1. / 12., 1. / 8., -1. / 24.}},
     {{1. / 36., 1. / 12., 1. / 8., -1. / 24.}},
     {{1. / 36., 1. / 12., 1. / 8., -1. / 24.}},
     {{1. / 36., 1. / 12., 1. / 8., -1. / 24.}},
     {{1. / 36., 1. / 12., 1. / 8., -1. / 24.}},
     {{1. / 36., 1. / 12., 1. / 8., -1. / 24.}}}};

/** Coefficients in the functional for the equilibrium distribution */
static constexpr const std::array<double, 19> w = {
    {1. / 3., 1. / 18., 1. / 18., 1. / 18., 1. / 18., 1. / 18., 1. / 18.,
     1. / 36., 1. / 36., 1. / 36., 1. / 36., 1. / 36., 1. / 36., 1. / 36.,
     1. / 36., 1. / 36., 1. / 36., 1. / 36., 1. / 36.}};

/** Basis of the mode space as described in [Duenweg, Schiller, Ladd] */
static constexpr const std::array<std::array<int, 19>, 19> e_ki = {
    {{{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}},
     {{0, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1, 0, 0, 0, 0}},
     {{0, 0, 0, 1, -1, 0, 0, 1, -1, -1, 1, 0, 0, 0, 0, 1, -1, 1, -1}},
     {{0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 1, -1, -1, 1, 1, -1, -1, 1}},
     {{-1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}},
     {{0, 1, 1, -1, -1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, -1, -1, -1, -1}},
     {{-0, 1, 1, 1, 1, -2, -2, 2, 2, 2, 2, -1, -1, -1, -1, -1, -1, -1, -1}},
     {{0, 0, 0, 0, 0, 0, 0, 1, 1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0}},
     {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, -1, -1, 0, 0, 0, 0}},
     {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, -1, -1}},
     {{0, -2, 2, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1, 0, 0, 0, 0}},
     {{0, 0, 0, -2, 2, 0, 0, 1, -1, -1, 1, 0, 0, 0, 0, 1, -1, 1, -1}},
     {{0, 0, 0, 0, 0, -2, 2, 0, 0, 0, 0, 1, -1, -1, 1, 1, -1, -1, 1}},
     {{0, -0, 0, 0, 0, 0, 0, 1, -1, 1, -1, -1, 1, -1, 1, 0, 0, 0, 0}},
     {{0, 0, 0, -0, 0, 0, 0, 1, -1, -1, 1, 0, 0, 0, 0, -1, 1, -1, 1}},
     {{0, 0, 0, 0, 0, -0, 0, 0, 0, 0, 0, 1, -1, -1, 1, -1, 1, 1, -1}},
     {{1, -2, -2, -2, -2, -2, -2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}},
     {{0, -1, -1, 1, 1, -0, -0, 0, 0, 0, 0, 1, 1, 1, 1, -1, -1, -1, -1}},
     {{0, -1, -1, -1, -1, 2, 2, 2, 2, 2, 2, -1, -1, -1, -1, -1, -1, -1, -1}}}};

static constexpr std::array<int, 19> e_ki_0 = {
    {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}};
static constexpr std::array<int, 19> e_ki_1 = {
    {0, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1, 0, 0, 0, 0}};
static constexpr std::array<int, 19> e_ki_2 = {
    {0, 0, 0, 1, -1, 0, 0, 1, -1, -1, 1, 0, 0, 0, 0, 1, -1, 1, -1}};
static constexpr std::array<int, 19> e_ki_3 = {
    {0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 1, -1, -1, 1, 1, -1, -1, 1}};
static constexpr std::array<int, 19> e_ki_4 = {
    {-1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}};
static constexpr std::array<int, 19> e_ki_5 = {
    {0, 1, 1, -1, -1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, -1, -1, -1, -1}};
static constexpr std::array<int, 19> e_ki_6 = {
    {-0, 1, 1, 1, 1, -2, -2, 2, 2, 2, 2, -1, -1, -1, -1, -1, -1, -1, -1}};
static constexpr std::array<int, 19> e_ki_7 = {
    {0, 0, 0, 0, 0, 0, 0, 1, 1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0}};
static constexpr std::array<int, 19> e_ki_8 = {
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, -1, -1, 0, 0, 0, 0}};
static constexpr std::array<int, 19> e_ki_9 = {
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, -1, -1}};
static constexpr std::array<int, 19> e_ki_10 = {
    {0, -2, 2, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1, 0, 0, 0, 0}};
static constexpr std::array<int, 19> e_ki_11 = {
    {0, 0, 0, -2, 2, 0, 0, 1, -1, -1, 1, 0, 0, 0, 0, 1, -1, 1, -1}};
static constexpr std::array<int, 19> e_ki_12 = {
    {0, 0, 0, 0, 0, -2, 2, 0, 0, 0, 0, 1, -1, -1, 1, 1, -1, -1, 1}};
static constexpr std::array<int, 19> e_ki_13 = {
    {0, -0, 0, 0, 0, 0, 0, 1, -1, 1, -1, -1, 1, -1, 1, 0, 0, 0, 0}};
static constexpr std::array<int, 19> e_ki_14 = {
    {0, 0, 0, -0, 0, 0, 0, 1, -1, -1, 1, 0, 0, 0, 0, -1, 1, -1, 1}};
static constexpr std::array<int, 19> e_ki_15 = {
    {0, 0, 0, 0, 0, -0, 0, 0, 0, 0, 0, 1, -1, -1, 1, -1, 1, 1, -1}};
static constexpr std::array<int, 19> e_ki_16 = {
    {1, -2, -2, -2, -2, -2, -2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}};
static constexpr std::array<int, 19> e_ki_17 = {
    {0, -1, -1, 1, 1, -0, -0, 0, 0, 0, 0, 1, 1, 1, 1, -1, -1, -1, -1}};
static constexpr std::array<int, 19> e_ki_18 = {
    {0, -1, -1, -1, -1, 2, 2, 2, 2, 2, 2, -1, -1, -1, -1, -1, -1, -1, -1}};

/* the following values are the (weighted) lengths of the vectors */
static constexpr const std::array<double, 19> w_k = {
    {1.0, 1. / 3., 1. / 3., 1. / 3., 2. / 3., 4. / 9., 4. / 3., 1. / 9.,
     1. / 9., 1. / 9., 2. / 3., 2. / 3., 2. / 3., 2. / 9., 2. / 9., 2. / 9.,
     2.0, 4. / 9., 4. / 3.}};

static constexpr const std::array<std::array<int, 19>, 19> e_ki_transposed = {
    {{{1, 0, 0, 0, -1, 0, -0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0}},
     {{1, 1, 0, 0, 0, 1, 1, 0, 0, 0, -2, 0, 0, -0, 0, 0, -2, -1, -1}},
     {{1, -1, 0, 0, 0, 1, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, -2, -1, -1}},
     {{1, 0, 1, 0, 0, -1, 1, 0, 0, 0, 0, -2, 0, 0, -0, 0, -2, 1, -1}},
     {{1, 0, -1, 0, 0, -1, 1, 0, 0, 0, 0, 2, 0, 0, 0, 0, -2, 1, -1}},
     {{1, 0, 0, 1, 0, 0, -2, 0, 0, 0, 0, 0, -2, 0, 0, -0, -2, -0, 2}},
     {{1, 0, 0, -1, 0, 0, -2, 0, 0, 0, 0, 0, 2, 0, 0, 0, -2, -0, 2}},
     {{1, 1, 1, 0, 1, 0, 2, 1, 0, 0, 1, 1, 0, 1, 1, 0, 1, 0, 2}},
     {{1, -1, -1, 0, 1, 0, 2, 1, 0, 0, -1, -1, 0, -1, -1, 0, 1, 0, 2}},
     {{1, 1, -1, 0, 1, 0, 2, -1, 0, 0, 1, -1, 0, 1, -1, 0, 1, 0, 2}},
     {{1, -1, 1, 0, 1, 0, 2, -1, 0, 0, -1, 1, 0, -1, 1, 0, 1, 0, 2}},
     {{1, 1, 0, 1, 1, 1, -1, 0, 1, 0, 1, 0, 1, -1, 0, 1, 1, 1, -1}},
     {{1, -1, 0, -1, 1, 1, -1, 0, 1, 0, -1, 0, -1, 1, 0, -1, 1, 1, -1}},
     {{1, 1, 0, -1, 1, 1, -1, 0, -1, 0, 1, 0, -1, -1, 0, -1, 1, 1, -1}},
     {{1, -1, 0, 1, 1, 1, -1, 0, -1, 0, -1, 0, 1, 1, 0, 1, 1, 1, -1}},
     {{1, 0, 1, 1, 1, -1, -1, 0, 0, 1, 0, 1, 1, 0, -1, -1, 1, -1, -1}},
     {{1, 0, -1, -1, 1, -1, -1, 0, 0, 1, 0, -1, -1, 0, 1, 1, 1, -1, -1}},
     {{1, 0, 1, -1, 1, -1, -1, 0, 0, -1, 0, 1, -1, 0, -1, 1, 1, -1, -1}},
     {{1, 0, -1, 1, 1, -1, -1, 0, 0, -1, 0, -1, 1, 0, 1, -1, 1, -1, -1}}}};

static constexpr const double c_sound_sq = 1. / 3.;

} // namespace D3Q19

#endif /* D3Q19_H */
