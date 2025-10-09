#include "test_autodiff_reverse.hpp" // reuse for same test infra
#include "test_functions_for_optimization.hpp"
#include <boost/math/differentiation/autodiff_reverse.hpp>
#include <boost/math/optimization/lbfgs.hpp>
#include <boost/math/optimization/minimizer.hpp>

namespace rdiff = boost::math::differentiation::reverse_mode;
namespace bopt = boost::math::optimization;

BOOST_AUTO_TEST_SUITE(basic_lbfgs)

BOOST_AUTO_TEST_CASE(default_lbfgs_test) {
  using T = double;
  constexpr size_t NITER = 10;
  constexpr size_t M = 10;
  const T eps = T{1e-8};

  RandomSample<T> rng{T(-10), T(10)};
  std::array<rdiff::rvar<T, 1>, 2> x;
  x[0] = rng.next();
  x[1] = rng.next();

  auto opt = bopt::make_lbfgs<decltype(&rosenbrock_saddle<rdiff::rvar<T, 1>>),
                              std::array<rdiff::rvar<T, 1>, 2>, T>(
      &rosenbrock_saddle<rdiff::rvar<T, 1>>, x, M);

  auto result = minimize(opt);
  std::cout << result << std::endl;
  for (auto &xi : x) {
    BOOST_REQUIRE_CLOSE(xi, T{1.0}, eps);
  }
}

BOOST_AUTO_TEST_SUITE_END()
