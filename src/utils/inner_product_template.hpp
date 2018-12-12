#ifndef UTILS_INNER_PRODUCT_TEMPLATE_HPP
#define UTILS_INNER_PRODUCT_TEMPLATE_HPP

#include <array>
#include <utility>

namespace Utils {

using Array = std::array<double, 19>;

template <int c> constexpr double mul(double ai) {
  return (c == 1) ? ai : c * ai;
}

template <int I, int c, int... cs> struct inner_product_template_impl {
  constexpr double operator()(Array const &a) const {
    return (c) ? mul<c>(std::get<I>(a)) +
                     inner_product_template_impl<I + 1, cs...>{}(a)
               : inner_product_template_impl<I + 1, cs...>{}(a);
  }
};

template <int I, int c> struct inner_product_template_impl<I, c> {
  constexpr double operator()(Array const &a) const {
    return (c) ? mul<c>(std::get<I>(a)) : 0.;
  }
};

template <typename T, std::size_t N, const std::array<int, N> &vec,
          std::size_t... I>
constexpr double inner_product_helper(const std::array<T, N> &n,
                                      std::index_sequence<I...>) {
  return inner_product_template_impl<0, std::get<I>(vec)...>{}(n);
}

template <typename T, std::size_t N, const std::array<int, N> &vec,
          typename Indices = std::make_index_sequence<N>>
constexpr double inner_product_template(const std::array<T, N> &n) {
  return inner_product_helper<T, N, vec>(n, Indices{});
}

} // namespace Utils

#endif
