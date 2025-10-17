//           Copyright Maksym Zhelyenzyakov 2025-2026.
// Distributed under the Boost Software License, Version 1.0.
//      (See accompanying file LICENSE_1_0.txt or copy at
//           https://www.boost.org/LICENSE_1_0.txt)
#ifndef DIFFERENTIABLE_OPT_UTILITIES_HPP
#define DIFFERENTIABLE_OPT_UTILITIES_HPP
#include <boost/math/differentiation/autodiff_reverse.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <cmath>
#include <random>
#include <type_traits>
#include <vector>

namespace boost {
namespace math {
namespace optimization {

namespace rdiff = boost::math::differentiation::reverse_mode;

/** @brief> helper to get the underlying realtype from
 *  update policy
 *  */
template<typename UpdPol>
struct update_policy_real_type;

template<template<typename> class UpdPol, typename RealType>
struct update_policy_real_type<UpdPol<RealType>>
{
  using type = RealType;
};

template<typename UpdPol>
using update_policy_real_type_t =
  typename update_policy_real_type<typename std::decay<UpdPol>::type>::type;

/** @brief> get realtype from argument container
 *  */
template<class Container>
struct argument_container_t;

template<template<typename, typename...> class Container,
         typename ValueType,
         typename... Args>
struct argument_container_t<Container<ValueType, Args...>>
{
  using type = ValueType;
};
template<template<typename, typename...> class Container,
         typename RealType,
         int N,
         typename... Args>
struct argument_container_t<Container<rdiff::rvar<RealType, N>, Args...>>
{
  using type = RealType;
};
template<typename ValueType, std::size_t N>
struct argument_container_t<std::array<ValueType, N>>
{
  using type = ValueType;
};

template<typename RealType, int M, std::size_t N>
struct argument_container_t<std::array<rdiff::rvar<RealType, M>, N>>
{
  using type = RealType;
};
/******************************************************************************/
/** @brief simple blas helpers
 * may optimize later if benchmarks show its needed, or just switch to Eigen
 */
template<typename Container>
auto
dot(const Container& x, const Container& y) -> typename Container::value_type
{
  using T = typename Container::value_type;
  BOOST_MATH_ASSERT(x.size() == y.size());
  return std::inner_product(x.begin(), x.end(), y.begin(), T(0));
}

template<typename Container>
auto
norm_2(const Container& x) -> typename Container::value_type
{
  return sqrt(dot(x, x));
}

template<typename Container>
auto
norm_1(const Container& x) -> typename Container::value_type
{
  using T = typename Container::value_type;
  T ret{ 0 };
  for (auto& xi : x) {
    ret += abs(xi);
  }
  return ret;
}

template<typename T>
T
norm_inf(const std::vector<T>& x)
{
  assert(!x.empty());

  T max_val = std::abs(x[0]);
  const std::size_t n = x.size();

  for (std::size_t i = 1; i < n; ++i) {
    const T abs_val = std::abs(x[i]);
    if (abs_val > max_val)
      max_val = abs_val;
  }
  return max_val;
}
/** @brief alpha*x (alpha is scalar, x is vector */
template<typename Container, typename RealType>
void
scale(Container& x, const RealType& alpha)
{
  for (auto& xi : x) {
    xi *= alpha;
  }
}

/** @brief y += alpha * x
 */
template<typename ContainerX, typename ContainerY, typename RealType>
void
axpy(RealType alpha, const ContainerX& x, ContainerY& y)
{
  BOOST_MATH_ASSERT(x.size() == y.size());
  const size_t n = x.size();
  for (size_t i = 0; i < n; ++i) {
    y[i] += alpha * x[i];
  }
}
/******************************************************************************/
template<typename RealType>
std::vector<RealType>
random_vector(size_t n)
{
  /** @brief> generates a random std::vector<RealType> of size n
   * using mt19937 algorithm
   */

  /** TODO: these may need to be marked thread local
   * in the future
   *
   * TODO: benchmark.
   */
  static boost::random::mt19937 rng{ std::random_device{}() };
  static boost::random::uniform_real_distribution<RealType> dist(0.0, 1.0);

  std::vector<RealType> result(n);
  std::generate(result.begin(), result.end(), [&] { return dist(rng); });
  return result;
}

} // namespace optimization
} // namespace math
} // namespace boost
#endif
