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

template <class T>
void default_nesterov_test()
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
    CHECK_ABSOLUTE_ERROR(0,xi.item(), eps);
  }
}

int main()
{
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { default_nesterov_test<typename decltype(t)::type>(); });
   return boost::math::test::report_errors();
}
