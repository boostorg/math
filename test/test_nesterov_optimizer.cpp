//           Copyright Maksym Zhelyenzyakov 2025-2026.
// Distributed under the Boost Software License, Version 1.0.
//      (See accompanying file LICENSE_1_0.txt or copy at
//           https://www.boost.org/LICENSE_1_0.txt)
#include "test_autodiff_reverse.hpp" // reuse for some basic options
#include "test_functions_for_optimization.hpp"
#include <boost/math/differentiation/autodiff_reverse.hpp>
#include <boost/math/optimization/minimizer.hpp>
#include <boost/math/optimization/nesterov.hpp>
namespace rdiff = boost::math::differentiation::reverse_mode;
namespace bopt = boost::math::optimization;
BOOST_AUTO_TEST_SUITE(nesterov_descent)

BOOST_AUTO_TEST_CASE_TEMPLATE(default_nesterov_test, T, all_float_types)
{
  T lr = T{ 1e-5 };
  T mu = T{ 0.95 };
  RandomSample<T> rng{ T(-10), (10) };
  std::vector<rdiff::rvar<T, 1>> x;
  x.push_back(rng.next());
  x.push_back(rng.next());
  T eps = T{ 1e-8 };
  auto nag =
    bopt::make_nag(&quadratic_high_cond_2D<rdiff::rvar<T, 1>>, x, lr, mu);
  auto constraint = bopt::unconstrained_policy<std::vector<rdiff::rvar<T, 1>>>{};
  auto convergence_policy = bopt::gradient_norm_convergence_policy<T>(T{ 1e-8 });

  auto z = minimize(nag, constraint, convergence_policy);
  for (auto& xi : x) {
    BOOST_REQUIRE_SMALL(xi.item(), eps);
  }
}
BOOST_AUTO_TEST_SUITE_END()
