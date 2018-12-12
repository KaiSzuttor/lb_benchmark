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
  m[0] = Utils::inner_product_template<T, 19, ::D3Q19::e_ki_0>(n);
  m[1] = Utils::inner_product_template<T, 19, ::D3Q19::e_ki_1>(n);
  m[2] = Utils::inner_product_template<T, 19, ::D3Q19::e_ki_2>(n);
  m[3] = Utils::inner_product_template<T, 19, ::D3Q19::e_ki_3>(n);
  m[4] = Utils::inner_product_template<T, 19, ::D3Q19::e_ki_4>(n);
  m[5] = Utils::inner_product_template<T, 19, ::D3Q19::e_ki_5>(n);
  m[6] = Utils::inner_product_template<T, 19, ::D3Q19::e_ki_6>(n);
  m[7] = Utils::inner_product_template<T, 19, ::D3Q19::e_ki_7>(n);
  m[8] = Utils::inner_product_template<T, 19, ::D3Q19::e_ki_8>(n);
  m[9] = Utils::inner_product_template<T, 19, ::D3Q19::e_ki_9>(n);
  m[10] = Utils::inner_product_template<T, 19, ::D3Q19::e_ki_10>(n);
  m[11] = Utils::inner_product_template<T, 19, ::D3Q19::e_ki_11>(n);
  m[12] = Utils::inner_product_template<T, 19, ::D3Q19::e_ki_12>(n);
  m[13] = Utils::inner_product_template<T, 19, ::D3Q19::e_ki_13>(n);
  m[14] = Utils::inner_product_template<T, 19, ::D3Q19::e_ki_14>(n);
  m[15] = Utils::inner_product_template<T, 19, ::D3Q19::e_ki_15>(n);
  m[16] = Utils::inner_product_template<T, 19, ::D3Q19::e_ki_16>(n);
  m[17] = Utils::inner_product_template<T, 19, ::D3Q19::e_ki_17>(n);
  m[18] = Utils::inner_product_template<T, 19, ::D3Q19::e_ki_18>(n);
  return m;
}

} // namespace LB

#endif
