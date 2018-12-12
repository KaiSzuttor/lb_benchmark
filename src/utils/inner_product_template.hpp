#ifndef UTILS_INNER_PRODUCT_TEMPLATE_HPP
#define UTILS_INNER_PRODUCT_TEMPLATE_HPP

#include <array>
#include <utility>

namespace Utils {

template <int c> constexpr double mul(double ai) {
  return (c == 1) ? ai : c * ai;
}

template <int I, typename T, std::size_t N, int c, int... cs>
struct inner_product_template_impl {
  constexpr T operator()(std::array<T, N> const &a) const {
    return (c) ? mul<c>(std::get<I>(a)) +
                     inner_product_template_impl<I + 1, T, N, cs...>{}(a)
               : inner_product_template_impl<I + 1, T, N, cs...>{}(a);
  }
};

template <int I, typename T, std::size_t N, int c>
struct inner_product_template_impl<I, T, N, c> {
  constexpr T operator()(std::array<T, N> const &a) const {
    return (c) ? mul<c>(std::get<I>(a)) : 0.;
  }
};

template <typename T, std::size_t N,
          const std::array<std::array<int, N>, N> &vec, std::size_t index,
          std::size_t... I>
constexpr T inner_product_helper(const std::array<T, N> &n,
                                 std::index_sequence<I...>) {
  return inner_product_template_impl<0, T, N, vec[index][I]...>{}(n);
}

template <typename T, std::size_t N,
          const std::array<std::array<int, N>, N> &vec, std::size_t index,
          typename Indices = std::make_index_sequence<N>>
constexpr T inner_product_template(const std::array<T, N> &n) {
  return inner_product_helper<T, N, vec, index>(n, Indices{});
}

} // namespace Utils

#endif
