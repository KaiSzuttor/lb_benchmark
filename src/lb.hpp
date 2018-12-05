#ifndef LB_HPP
#define LB_HPP

#include "lb-d3q19.hpp"
#include "utils/inner_product.hpp"

namespace LB {

template <typename T>
std::array<T, 19> normalize_modes(const std::array<T, 19> &modes) {
  auto normalized_modes = modes;
  for (int i = 0; i < modes.size(); i++) {
    normalized_modes[i] /= ::D3Q19::w_k[i];
  }
  return normalized_modes;
}

template <typename T>
std::array<T, 19> lb_calc_n_from_m(const std::array<T, 19> &modes) {
  std::array<T, 19> ret;
  const auto normalized_modes = normalize_modes(modes);

  for (int i = 0; i < 19; i++) {
    ret[i] =
        Utils::inner_product(::D3Q19::e_ki_transposed[i], normalized_modes) *
        ::D3Q19::w[i];
  }
  return ret;
}

template <typename T>
std::array<T, 19> lb_calc_m_from_n(const std::array<T, 19> &n) {
  std::array<T, 19> m;
  for (int i = 0; i < 19; i++) {
    m[i] = Utils::inner_product(::D3Q19::e_ki[i], n);
  }
  return m;
}

} // namespace LB

#endif
