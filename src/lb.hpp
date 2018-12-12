#ifndef LB_HPP
#define LB_HPP

#include "lb-d3q19.hpp"
#include "utils/inner_product.hpp"
#include "utils/inner_product_template.hpp"

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

template <typename T>
std::array<T, 19> lb_calc_m_from_n_template(const std::array<T, 19> &n) {
  std::array<T, 19> m;
  m[0] = Utils::inner_product_template<T, 19, ::D3Q19::e_ki, 0>(n);
  m[1] = Utils::inner_product_template<T, 19, ::D3Q19::e_ki, 1>(n);
  m[2] = Utils::inner_product_template<T, 19, ::D3Q19::e_ki, 2>(n);
  m[3] = Utils::inner_product_template<T, 19, ::D3Q19::e_ki, 3>(n);
  m[4] = Utils::inner_product_template<T, 19, ::D3Q19::e_ki, 4>(n);
  m[5] = Utils::inner_product_template<T, 19, ::D3Q19::e_ki, 5>(n);
  m[6] = Utils::inner_product_template<T, 19, ::D3Q19::e_ki, 6>(n);
  m[7] = Utils::inner_product_template<T, 19, ::D3Q19::e_ki, 7>(n);
  m[8] = Utils::inner_product_template<T, 19, ::D3Q19::e_ki, 8>(n);
  m[9] = Utils::inner_product_template<T, 19, ::D3Q19::e_ki, 9>(n);
  m[10] = Utils::inner_product_template<T, 19, ::D3Q19::e_ki, 10>(n);
  m[11] = Utils::inner_product_template<T, 19, ::D3Q19::e_ki, 11>(n);
  m[12] = Utils::inner_product_template<T, 19, ::D3Q19::e_ki, 12>(n);
  m[13] = Utils::inner_product_template<T, 19, ::D3Q19::e_ki, 13>(n);
  m[14] = Utils::inner_product_template<T, 19, ::D3Q19::e_ki, 14>(n);
  m[15] = Utils::inner_product_template<T, 19, ::D3Q19::e_ki, 15>(n);
  m[16] = Utils::inner_product_template<T, 19, ::D3Q19::e_ki, 16>(n);
  m[17] = Utils::inner_product_template<T, 19, ::D3Q19::e_ki, 17>(n);
  m[18] = Utils::inner_product_template<T, 19, ::D3Q19::e_ki, 18>(n);
  return m;
}

} // namespace LB

#endif
