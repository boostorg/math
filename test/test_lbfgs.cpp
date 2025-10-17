//           Copyright Maksym Zhelyenzyakov 2025-2026.
// Distributed under the Boost Software License, Version 1.0.
//      (See accompanying file LICENSE_1_0.txt or copy at
//           https://www.boost.org/LICENSE_1_0.txt)
#include "test_autodiff_reverse.hpp"
#include "test_functions_for_optimization.hpp"
#include <boost/math/differentiation/autodiff_reverse.hpp>
#include <boost/math/optimization/lbfgs.hpp>
#include <boost/math/optimization/minimizer.hpp>

namespace rdiff = boost::math::differentiation::reverse_mode;
namespace bopt = boost::math::optimization;

BOOST_AUTO_TEST_SUITE(basic_lbfgs)

BOOST_AUTO_TEST_CASE_TEMPLATE(default_lbfgs_test, T, all_float_types)
{
  constexpr size_t NITER = 10;
  constexpr size_t M = 10;
  const T eps = T{ 1e-8 };

  RandomSample<T> rng{ T(-10), T(10) };
  std::array<rdiff::rvar<T, 1>, 2> x;
  x[0] = rng.next();
  x[1] = rng.next();

  auto opt = bopt::make_lbfgs(&rosenbrock_saddle<rdiff::rvar<T, 1>>, x, M);

  auto result = minimize(opt);
  for (auto& xi : x) {
    BOOST_REQUIRE_CLOSE(xi, T{ 1.0 }, eps);
  }
}

// Custom initialization policy that zeros out the parameters
template<typename RealType>
struct zero_init_policy
{
  void operator()(std::vector<RealType>& x) const noexcept
  {
    std::fill(x.begin(), x.end(), RealType{ 0 });
  }
};

template<typename RealType>
struct analytic_objective_eval_pol
{
  template<typename Objective, typename ArgumentContainer>
  RealType operator()(Objective&& objective, ArgumentContainer& x)
  {
    return objective(x);
  }
};
template<typename RealType>
struct analytic_gradient_eval_pol
{
  template<class Objective,
           class ArgumentContainer,
           class FunctionEvaluationPolicy>
  void operator()(Objective&& obj_f,
                  ArgumentContainer& x,
                  FunctionEvaluationPolicy&& f_eval_pol,
                  RealType& obj_v,
                  std::vector<RealType>& grad_container)
  {
    RealType v = f_eval_pol(obj_f, x);
    obj_v = v;
    grad_container.resize(x.size());
    for (size_t i = 0; i < x.size(); ++i) {
      grad_container[i] = 2 * x[i];
    }
  }
};

// -- Test L-BFGS with custom initialization policy (zero_init_policy)
BOOST_AUTO_TEST_CASE_TEMPLATE(custom_init_lbfgs_test, T, all_float_types)
{
  constexpr size_t M = 8;
  const T eps = T{ 1e-6 };

  RandomSample<T> rng{ T(-5), T(5) };
  std::array<rdiff::rvar<T, 1>, 2> x;
  x[0] = rng.next();
  x[1] = rng.next();

  auto opt = bopt::make_lbfgs(&rosenbrock_saddle<rdiff::rvar<T, 1>>,
                              x,
                              M,
                              bopt::costant_initializer_rvar<T>(0.0));
  auto result = minimize(opt);

  for (auto& xi : x) {
    BOOST_REQUIRE_CLOSE(xi, T{ 1.0 }, eps);
  }
}

// // -- Test L-BFGS with analytic derivative policies
// BOOST_AUTO_TEST_CASE_TEMPLATE(analytic_lbfgs_test, T, all_float_types)
// {
//   constexpr size_t M = 10;
//   const T eps = T{ 1e-3 };

//   RandomSample<T> rng{ T(-5), T(5) };
//   std::vector<T> x(3);
//   for (auto& xi : x)
//     xi = rng.next();

//   auto opt = bopt::make_lbfgs(
//     &quadratic<rdiff::rvar<T, 1>>,              // Objective
//     x,                                          // Arguments
//     M,                                          // History size
//     bopt::random_uniform_initializer_rvar<T>{}, // Initialization
//     analytic_objective_eval_pol<T>{},           // Function eval
//     analytic_gradient_eval_pol<T>{}             // Gradient eval
//   );

//   // Run optimization (manual loop or minimize wrapper)
//   auto result = minimize(opt);

//   for (auto& xi : x) {
//     BOOST_REQUIRE_SMALL(xi, eps);
//   }
// }

// // -- Test L-BFGS with analytic policies and custom line search
// // (strong_wolfe_line_search_policy)
// BOOST_AUTO_TEST_CASE_TEMPLATE(analytic_lbfgs_strong_wolfe_test,
//                               T,
//                               all_float_types)
// {
//   constexpr size_t M = 12;
//   const T eps = T{ 1e-4 };

//   RandomSample<T> rng{ T(-10), T(10) };
//   std::vector<T> x(5);
//   for (auto& xi : x)
//     xi = rng.next();

//   auto opt = bopt::make_lbfgs(
//     &quadratic<T>,
//     x,
//     M,
//     boost::math::optimization::random_uniform_initializer_rvar<T>{},
//     analytic_objective_eval_pol<T>{},
//     analytic_gradient_eval_pol<T>{},
//     boost::math::optimization::armijo_line_search_policy<T>{});

//   auto result = minimize(opt);

//   for (auto& xi : x) {
//     BOOST_REQUIRE_SMALL(xi, eps);
//   }
// }

// // -- Test L-BFGS with random init policy (demonstrates flexible
// initialization) template<typename RealType> struct random_init_policy
// {
//   RandomSample<RealType> rng{ RealType(-1), RealType(1) };
//   void operator()(std::vector<RealType>& x) const noexcept
//   {
//     for (auto& xi : x)
//       xi = rng.next();
//   }
// };

// BOOST_AUTO_TEST_CASE_TEMPLATE(random_init_lbfgs_test, T, all_float_types)
// {
//   constexpr size_t M = 6;
//   const T eps = T{ 1e-4 };

//   std::vector<T> x(4);
//   random_init_policy<T>{}(x); // Apply initialization manually

//   auto opt = bopt::make_lbfgs(&quadratic<T>,
//                               x,
//                               M,
//                               random_init_policy<T>{},
//                               analytic_objective_eval_pol<T>{},
//                               analytic_gradient_eval_pol<T>{});

//   auto result = minimize(opt);

//   for (auto& xi : x) {
//     BOOST_REQUIRE_SMALL(xi, eps);
//   }
// }
BOOST_AUTO_TEST_SUITE_END()
