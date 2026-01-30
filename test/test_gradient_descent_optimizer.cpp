//           Copyright Maksym Zhelyenzyakov 2025-2026.
// Distributed under the Boost Software License, Version 1.0.
//      (See accompanying file LICENSE_1_0.txt or copy at
//           https://www.boost.org/LICENSE_1_0.txt)
#include "test_autodiff_reverse.hpp" // reuse for some basic options
#include "test_functions_for_optimization.hpp"
#include <boost/math/differentiation/autodiff_reverse.hpp>
#include <boost/math/optimization/gradient_descent.hpp>
#include <boost/math/optimization/minimizer.hpp>
namespace rdiff = boost::math::differentiation::reverse_mode;
namespace bopt = boost::math::optimization;
BOOST_AUTO_TEST_SUITE(basic_gradient_descent)

BOOST_AUTO_TEST_CASE_TEMPLATE(default_gd_test, T, all_float_types)
{
  size_t NITER = 2000;
  size_t N = 15;
  T lr = T{ 1e-2 };
  RandomSample<T> rng{ T(-100), (100) };
  std::vector<rdiff::rvar<T, 1>> x_ad;
  T eps = T{ 1e-3 };
  for (size_t i = 0; i < N; ++i) {
    x_ad.push_back(rng.next());
  }
  auto gdopt =
    bopt::make_gradient_descent(&quadratic<rdiff::rvar<T, 1>>, x_ad, lr);
  for (size_t i = 0; i < NITER; ++i) {
    gdopt.step();
  }
  for (auto& x : x_ad) {
    BOOST_REQUIRE_SMALL(x.item(), eps);
  }
}
BOOST_AUTO_TEST_CASE_TEMPLATE(test_minimize, T, all_float_types)
{
  size_t NITER = 2000;
  size_t N = 15;
  T lr = T{ 1e-2 };
  RandomSample<T> rng{ T(-100), (100) };
  std::vector<rdiff::rvar<T, 1>> x_ad;
  T eps = T{ 1e-3 };
  for (size_t i = 0; i < N; ++i) {
    x_ad.push_back(rng.next());
  }
  auto gdopt =
    bopt::make_gradient_descent(&quadratic<rdiff::rvar<T, 1>>, x_ad, lr);
  auto z = minimize(gdopt);
  for (auto& x : x_ad) {
    BOOST_REQUIRE_SMALL(x.item(), eps);
  }
}
BOOST_AUTO_TEST_CASE_TEMPLATE(random_initializer_test, T, all_float_types)
{
  size_t N = 10;
  T lr = T{ 1e-2 };
  std::vector<rdiff::rvar<T, 1>> x(N);

  auto gdopt =
    bopt::make_gradient_descent(&quadratic<rdiff::rvar<T, 1>>,
                                x,
                                lr,
                                bopt::random_uniform_initializer_rvar<T>(
                                  T(-2.0), T(2.0), 1234));
  for (auto& xi : x) {
    T v = xi.item();
    BOOST_TEST(v >= -2);
    BOOST_TEST(v <= 2);
  }
  gdopt.step();
}
BOOST_AUTO_TEST_CASE_TEMPLATE(const_initializer_test, T, all_float_types)
{
  size_t N = 10;
  T lr = T{ 1e-2 };
  std::vector<rdiff::rvar<T, 1>> x(N);

  auto gdopt = bopt::make_gradient_descent(
    &quadratic<rdiff::rvar<T, 1>>,
    x,
    lr,
    bopt::costant_initializer_rvar<T>(T{ 5.0 })); // all initialized to 5

  for (auto& xi : x) {
    T v = xi.item();
    BOOST_REQUIRE_CLOSE(v, T{ 5.0 }, T{ 1e-3 });
  }
  gdopt.step();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(box_constraint_test, T, all_float_types)
{
  size_t N = 5;
  T lr = T{ 1e-2 };
  std::vector<rdiff::rvar<T, 1>> x(N, T{ 10 });

  auto gdopt =
    bopt::make_gradient_descent(&quadratic<rdiff::rvar<T, 1>>, x, lr);

  auto res = bopt::minimize(
    gdopt, bopt::box_constraints<std::vector<rdiff::rvar<T, 1>>, T>(-1.0, 1.0));

  for (auto& xi : x) {
    BOOST_TEST(xi.item() >= -1.0);
    BOOST_TEST(xi.item() <= 1.0);
  }
}
BOOST_AUTO_TEST_CASE_TEMPLATE(max_iter_test, T, all_float_types)
{
  size_t N = 2;
  T lr = T{ 1e-6 }; // very slow learning
  std::vector<rdiff::rvar<T, 1>> x = { T{ 5 }, T{ 5 } };

  auto gdopt =
    bopt::make_gradient_descent(&quadratic<rdiff::rvar<T, 1>>, x, lr);

  size_t max_iter = 50;
  auto res =
    bopt::minimize(gdopt,
                   bopt::unconstrained_policy<std::vector<rdiff::rvar<T, 1>>>{},
                   bopt::gradient_norm_convergence_policy<T>(T{ 1e-20 }),
                   bopt::max_iter_termination_policy(max_iter));

  BOOST_TEST(!res.converged); // should not converge with tiny lr
  BOOST_REQUIRE_EQUAL(res.num_iter, max_iter);
}
BOOST_AUTO_TEST_CASE_TEMPLATE(history_tracking_test, T, all_float_types)
{
  size_t N = 3;
  T lr = T{ 1e-2 };
  std::vector<rdiff::rvar<T, 1>> x = { T{ 3 }, T{ -4 }, T{ 5 } };

  auto gdopt =
    bopt::make_gradient_descent(&quadratic<rdiff::rvar<T, 1>>, x, lr);

  auto res =
    bopt::minimize(gdopt,
                   bopt::unconstrained_policy<std::vector<rdiff::rvar<T, 1>>>{},
                   bopt::gradient_norm_convergence_policy<T>(T{ 1e-6 }),
                   bopt::max_iter_termination_policy(1000),
                   true); // enable history

  BOOST_TEST(!res.objective_history.empty());
  BOOST_TEST(res.objective_history.front() > res.objective_history.back());
}
BOOST_AUTO_TEST_CASE_TEMPLATE(rosenbrock_test, T, all_float_types)
{
  std::array<rdiff::rvar<T, 1>, 2> x = { T{ -1.2 }, T{ 1.0 } }; // bad start
  T lr = T{ 1e-3 };

  auto gdopt =
    bopt::make_gradient_descent(&rosenbrock_saddle<rdiff::rvar<T, 1>>, x, lr);

  auto res = bopt::minimize(
    gdopt,
    bopt::unconstrained_policy<std::array<rdiff::rvar<T, 1>, 2>>{},
    bopt::gradient_norm_convergence_policy<T>(T{ 1e-4 }),
    bopt::max_iter_termination_policy(50000));
  BOOST_TEST(res.converged);
  BOOST_REQUIRE_CLOSE(x[0].item(), T{ 1.0 }, T{ 1e-1 });
  BOOST_REQUIRE_CLOSE(x[1].item(), T{ 1.0 }, T{ 1e-1 });
}
BOOST_AUTO_TEST_CASE_TEMPLATE(objective_tol_convergence_test,
                              T,
                              all_float_types)
{
  using policy_t = bopt::objective_tol_convergence_policy<T>;
  policy_t pol(1e-3);
  std::vector<T> dummy_grad;

  BOOST_TEST(!pol(dummy_grad, T{100.0}));
  BOOST_TEST(!pol(dummy_grad, T{99.0}));
  BOOST_TEST(pol(dummy_grad, T{99.0005}));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(relative_objective_tol_test, T, all_float_types)
{
  using policy_t = bopt::relative_objective_tol_policy<T>;
  policy_t pol(1e-3);

  std::vector<T> dummy_grad;
  BOOST_TEST(!pol(dummy_grad, T{1000.0}));
  BOOST_TEST(!pol(dummy_grad, T{1010.0}));
  BOOST_TEST(pol(dummy_grad, T{1010.5}));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(combined_policy_test, T, all_float_types)
{
  using pol_abs = bopt::objective_tol_convergence_policy<T>;
  using pol_rel = bopt::relative_objective_tol_policy<T>;
  using pol_comb = bopt::combined_convergence_policy<pol_abs, pol_rel>;

  pol_abs abs_pol(1e-6);
  pol_rel rel_pol(1e-3);
  pol_comb comb(abs_pol, rel_pol);

  std::vector<T> dummy_grad;

  BOOST_TEST(!comb(dummy_grad, T{100.0}));
  BOOST_TEST(!comb(dummy_grad, T{110.0}));
  BOOST_TEST(comb(dummy_grad, T{110.1}));
  BOOST_TEST(comb(dummy_grad, T{110.1000001}));
}
BOOST_AUTO_TEST_CASE_TEMPLATE(nonnegativity_constraint_test, T, all_float_types)
{
    std::vector<T> x = {T{1.0}, T{-2.0}, T{3.0}, T{-4.0}};

    bopt::nonnegativity_constraint<std::vector<T>, T> proj;
    proj(x);

    for (auto& xi : x)
        BOOST_TEST(xi >= 0.0);
    BOOST_TEST(x == std::vector<T>({T{1.0}, T{0.0}, T{3.0}, T{0.0}}));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(l2_ball_constraint_test, T, all_float_types)
{
    std::vector<T>                              x = {T{3.0}, T{4.0}}; // norm = 5
    bopt::l2_ball_constraint<std::vector<T>, T> proj(1.0);
    proj(x);

    T norm = sqrt(x[0] * x[0] + x[1] * x[1]);
    BOOST_TEST(abs(norm - T{1.0}) < T{1e-12}); // projected to unit circle
}

BOOST_AUTO_TEST_CASE_TEMPLATE(l1_ball_constraint_test, T, all_float_types)
{
    std::vector<T> x = {T{3.0}, T{4.0}}; // L1 norm = 7

    bopt::l1_ball_constraint<std::vector<T>, T> proj(2.0);
    proj(x);

    T norm1 = abs(x[0]) + abs(x[1]);
    BOOST_TEST(abs(norm1 - T{2.0}) < T{1e-12});
}

BOOST_AUTO_TEST_CASE_TEMPLATE(simplex_constraint_test, T, all_float_types)
{
    std::vector<T> x = {T{-1.0}, T{2.0}, T{3.0}}; // has negative and sum != 1

    bopt::simplex_constraint<std::vector<T>, T> proj;
    proj(x);

    T sum = 0.0;
    for (auto& xi : x) {
        BOOST_TEST(xi >= 0.0); // all nonnegative
        sum += xi;
    }
    BOOST_TEST(abs(sum - T{1.0}) < T{1e-12}); // normalized to sum=1
}

BOOST_AUTO_TEST_CASE_TEMPLATE(unit_sphere_constraint_test, T, all_float_types)
{
    std::vector<T>                                  x = {T{0.3}, T{0.4}}; // norm = 0.5
    bopt::unit_sphere_constraint<std::vector<T>, T> proj;
    proj(x);

    T norm = sqrt(x[0] * x[0] + x[1] * x[1]);
    BOOST_TEST(abs(norm - T{1.0}) < T{1e-12}); // always projected to sphere
}

BOOST_AUTO_TEST_CASE_TEMPLATE(function_constraint_test, T, all_float_types)
{
  auto clip_to_half = [](std::vector<T>& v) {
    for (auto& xi : v)
      if (xi > 0.5)
        xi = 0.5;
  };

  bopt::function_constraint<std::vector<T>> proj(clip_to_half);
  std::vector<T>                            x = {T{0.2}, T{0.7}, T{1.5}};
  proj(x);

  BOOST_TEST(x == std::vector<T>({T{0.2}, T{0.5}, T{0.5}}));
}

template<typename RealType>
struct no_init_policy
{
  void operator()(std::vector<RealType>& x) const noexcept {}
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

BOOST_AUTO_TEST_CASE_TEMPLATE(analytic_derivative_policies, T, all_float_types)
{
  std::vector<T> x;
  size_t NITER = 2000;
  size_t N = 15;
  T lr = T{ 1e-2 };
  RandomSample<T> rng{ T(-100), (100) };
  T eps = T{ 1e-3 };
  for (size_t i = 0; i < N; ++i) {
    x.push_back(rng.next());
  }

  auto gdopt = bopt::make_gradient_descent(&quadratic<T>,
                                           x,
                                           lr,
                                           no_init_policy<T>{},
                                           analytic_objective_eval_pol<T>{},
                                           analytic_gradient_eval_pol<T>{});

  for (size_t i = 0; i < NITER; ++i) {
    gdopt.step();
  }
  for (auto& xi : x) {
    BOOST_REQUIRE_SMALL(xi, eps);
  }
}
BOOST_AUTO_TEST_SUITE_END()
