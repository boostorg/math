/*
 * Copyright Nick Thompson, 2026
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */

#ifndef BOOST_MATH_BUILD_MODULE
#include <boost/math/optimization/turbo.hpp>
#else
import boost.math;
#endif

#include "math_unit_test.hpp"
#include "test_functions_for_optimization.hpp"
#ifndef BOOST_MATH_BUILD_MODULE
#include <array>
#include <atomic>
#include <cmath>
#include <limits>
#include <random>
#include <stdexcept>
#include <vector>
#endif
using boost::math::optimization::turbo;
using boost::math::optimization::turbo_parameters;
using boost::math::optimization::validate_turbo_parameters;

// The random number generators of the standard library are fixed, but its distributions are not, so the same seed
// gives different runs on different standard libraries. The tolerances below hold for any seed we tried.

template <class ArgType, class Real>
void check_in_bounds(std::vector<std::pair<ArgType, Real>> const &queries, ArgType const &lo, ArgType const &hi) {
  bool inside = true;
  for (auto const &q : queries) {
    for (size_t i = 0; i < lo.size(); ++i) {
      inside = inside && lo[i] <= q.first[i] && q.first[i] <= hi[i];
    }
  }
  CHECK_TRUE(inside);
}

template <class Real> Real branin(std::array<Real, 2> const &v) {
  using boost::math::constants::pi;
  using std::cos;
  Real x = v[0];
  Real y = v[1];
  Real t = y - Real(5.1) / (4 * pi<Real>() * pi<Real>()) * x * x + 5 * x / pi<Real>() - 6;
  return t * t + 10 * (1 - 1 / (8 * pi<Real>())) * cos(x) + 10;
}

template <class Real> void test_branin() {
  std::cout << "Testing TuRBO on the Branin function . . .\n";
  using ArgType = std::array<Real, 2>;
  turbo_parameters<ArgType> params;
  params.lower_bounds = {-5, 0};
  params.upper_bounds = {10, 15};
  params.max_function_calls = 200;
  std::mt19937_64 gen(12345);
  std::vector<std::pair<ArgType, Real>> queries;
  std::atomic<Real> current_minimum_cost = std::numeric_limits<Real>::quiet_NaN();
  auto x = turbo(branin<Real>, params, gen, std::numeric_limits<Real>::quiet_NaN(), nullptr, &current_minimum_cost,
                 &queries);
  // The global minimum value, attained at (-pi, 12.275), (pi, 2.275) and (9.42478, 2.475):
  CHECK_ABSOLUTE_ERROR(Real(0.39788735772973816), branin(x), Real(1e-2));
  CHECK_EQUAL(queries.size(), size_t(200));
  check_in_bounds(queries, params.lower_bounds, params.upper_bounds);
  // The returned point and the current minimum are those of the best query:
  Real best = std::numeric_limits<Real>::infinity();
  for (auto const &q : queries) {
    CHECK_EQUAL(branin(q.first), q.second);
    best = (std::min)(best, q.second);
  }
  CHECK_EQUAL(branin(x), best);
  CHECK_EQUAL(current_minimum_cost.load(), best);
}

template <class Real> void test_six_hump_camel() {
  std::cout << "Testing TuRBO on the six-hump camel function . . .\n";
  using ArgType = std::vector<Real>;
  auto camel = [](ArgType const &v) {
    Real x = v[0];
    Real y = v[1];
    return (4 - Real(21) / 10 * x * x + x * x * x * x / 3) * x * x + x * y + (-4 + 4 * y * y) * y * y;
  };
  turbo_parameters<ArgType> params;
  params.lower_bounds = {-3, -2};
  params.upper_bounds = {3, 2};
  params.max_function_calls = 200;
  std::mt19937_64 gen(2);
  auto x = turbo(camel, params, gen);
  CHECK_EQUAL(x.size(), size_t(2));
  CHECK_ABSOLUTE_ERROR(Real(-1.0316284534898774), camel(x), Real(1e-2));
}

void test_hartmann6() {
  std::cout << "Testing TuRBO on the Hartmann 6 function . . .\n";
  using ArgType = std::vector<double>;
  auto hartmann = [](ArgType const &x) {
    static const double A[4][6] = {{10, 3, 17, 3.5, 1.7, 8}, {0.05, 10, 17, 0.1, 8, 14}, {3, 3.5, 1.7, 10, 17, 8},
                                   {17, 8, 0.05, 10, 0.1, 14}};
    static const double P[4][6] = {{0.1312, 0.1696, 0.5569, 0.0124, 0.8283, 0.5886},
                                   {0.2329, 0.4135, 0.8307, 0.3736, 0.1004, 0.9991},
                                   {0.2348, 0.1451, 0.3522, 0.2883, 0.3047, 0.6650},
                                   {0.4047, 0.8828, 0.8732, 0.5743, 0.1091, 0.0381}};
    static const double c[4] = {1, 1.2, 3, 3.2};
    double s = 0;
    for (int i = 0; i < 4; ++i) {
      double e = 0;
      for (int j = 0; j < 6; ++j) {
        e += A[i][j] * (x[j] - P[i][j]) * (x[j] - P[i][j]);
      }
      s -= c[i] * std::exp(-e);
    }
    return s;
  };
  turbo_parameters<ArgType> params;
  params.lower_bounds.resize(6, 0);
  params.upper_bounds.resize(6, 1);
  params.max_function_calls = 300;
  // Fewer candidates than the default 600, to keep the test fast:
  params.candidates = 300;
  std::mt19937_64 gen(3);
  auto x = turbo(hartmann, params, gen);
  // The global minimum is -3.32237, and the next best local minimum -3.2032. TuRBO-1 is a local method with restarts,
  // and about one run in five ends in the second basin after 300 calls (one in twenty after 500).
  CHECK_LE(hartmann(x), -3.19);
}

template <class Real> void test_rosenbrock() {
  std::cout << "Testing TuRBO on the Rosenbrock function . . .\n";
  using ArgType = std::array<Real, 2>;
  turbo_parameters<ArgType> params;
  params.lower_bounds = {-1.5, -0.5};
  params.upper_bounds = {2.5, 3};
  params.max_function_calls = 400;
  std::mt19937_64 gen(4);
  auto x = turbo(rosenbrock_saddle<Real>, params, gen);
  // The function varies over five orders of magnitude on this box, which is hard for the Gaussian process. Over 20
  // seeds, the median is 0.006 and the worst 0.12.
  CHECK_LE(rosenbrock_saddle(x), Real(0.25));
}

// The sphere of test_functions_for_optimization.hpp is NaN outside the unit ball, which covers most of these boxes.
// It also tests a return type different from the argument type.
void test_sphere() {
  std::cout << "Testing TuRBO on the sphere, on asymmetric boxes . . .\n";
  using ArgType = std::vector<float>;
  turbo_parameters<ArgType> params;
  params.lower_bounds = {-1, -0.4f, -0.7f};
  params.upper_bounds = {0.6f, 1, 0.3f};
  params.max_function_calls = 150;
  std::mt19937_64 gen(5);
  std::vector<std::pair<ArgType, double>> queries;
  auto x = turbo(sphere, params, gen, std::numeric_limits<double>::quiet_NaN(), nullptr, nullptr, &queries);
  CHECK_LE(sphere(x), 1e-3);
  check_in_bounds(queries, params.lower_bounds, params.upper_bounds);

  // A box that excludes the minimizer, which is then at the corner (0.2, 0.1):
  params.lower_bounds = {0.2f, 0.1f};
  params.upper_bounds = {0.9f, 0.5f};
  queries.clear();
  x = turbo(sphere, params, gen, std::numeric_limits<double>::quiet_NaN(), nullptr, nullptr, &queries);
  CHECK_LE(sphere(x) - 0.05, 2e-3);
  check_in_bounds(queries, params.lower_bounds, params.upper_bounds);
}

void test_termination() {
  std::cout << "Testing TuRBO termination: budget, cancellation, target value, batches . . .\n";
  using ArgType = std::array<double, 2>;
  turbo_parameters<ArgType> params;
  params.lower_bounds = {0.5, 0.5};
  params.upper_bounds = {2.048, 2.048};
  params.max_function_calls = 37;
  std::mt19937_64 gen(6);
  std::vector<std::pair<ArgType, double>> queries;
  auto const nan = std::numeric_limits<double>::quiet_NaN();
  turbo(rosenbrock_saddle<double>, params, gen, nan, nullptr, nullptr, &queries);
  CHECK_EQUAL(queries.size(), size_t(37));

  // The budget is exact even when it ends a batch early:
  params.batch_size = 3;
  queries.clear();
  turbo(rosenbrock_saddle<double>, params, gen, nan, nullptr, nullptr, &queries);
  CHECK_EQUAL(queries.size(), size_t(37));
  check_in_bounds(queries, params.lower_bounds, params.upper_bounds);
  params.batch_size = 1;

  // The default budget is 100 n:
  params.max_function_calls = 0;
  queries.clear();
  turbo(rosenbrock_saddle<double>, params, gen, nan, nullptr, nullptr, &queries);
  CHECK_EQUAL(queries.size(), size_t(200));

  // Cancellation from inside the cost function stops it after exactly that call:
  std::atomic<bool> cancel = false;
  size_t calls = 0;
  auto f = [&](ArgType const &x) {
    if (++calls == 10) {
      cancel = true;
    }
    return rosenbrock_saddle<double>(x);
  };
  queries.clear();
  turbo(f, params, gen, nan, &cancel, nullptr, &queries);
  CHECK_EQUAL(queries.size(), size_t(10));
  // Cancelled before the first call: no queries, and the midpoint of the box comes back.
  queries.clear();
  auto x = turbo(f, params, gen, nan, &cancel, nullptr, &queries);
  CHECK_EQUAL(queries.size(), size_t(0));
  CHECK_EQUAL(x[0], (params.lower_bounds[0] + params.upper_bounds[0]) / 2);
  CHECK_EQUAL(x[1], (params.lower_bounds[1] + params.upper_bounds[1]) / 2);

  // The run stops at the first call that reaches the target value:
  double const target = 0.05;
  queries.clear();
  x = turbo(rosenbrock_saddle<double>, params, gen, target, nullptr, nullptr, &queries);
  CHECK_LE(rosenbrock_saddle(x), target);
  CHECK_LE(queries.back().second, target);
  for (size_t k = 0; k + 1 < queries.size(); ++k) {
    CHECK_GE(queries[k].second, target);
  }
  CHECK_LE(queries.size(), size_t(199));

  // The initial guess is the first call; if it reaches the target, it is the only one:
  ArgType guess{1, 1};
  params.initial_guess = &guess;
  queries.clear();
  x = turbo(rosenbrock_saddle<double>, params, gen, 0.0, nullptr, nullptr, &queries);
  CHECK_EQUAL(queries.size(), size_t(1));
  CHECK_EQUAL(x[0], 1.0);
  CHECK_EQUAL(x[1], 1.0);
}

void test_validation() {
  std::cout << "Testing TuRBO parameter validation . . .\n";
  using ArgType = std::array<double, 2>;
  turbo_parameters<ArgType> params;
  params.lower_bounds = {-1, -1};
  params.upper_bounds = {1, 1};
  validate_turbo_parameters(params);
  std::mt19937_64 gen(7);

  auto bad = params;
  bad.upper_bounds[1] = -1;
  CHECK_THROW(validate_turbo_parameters(bad), std::domain_error);
  auto run = [&]() { return turbo(rosenbrock_saddle<double>, bad, gen); };
  CHECK_THROW(run(), std::domain_error);

  // TuRBO scales the box to the unit cube, so it must have positive width:
  bad = params;
  bad.upper_bounds[0] = bad.lower_bounds[0];
  CHECK_THROW(validate_turbo_parameters(bad), std::domain_error);

  bad = params;
  bad.upper_bounds[0] = std::numeric_limits<double>::infinity();
  CHECK_THROW(validate_turbo_parameters(bad), std::domain_error);

  ArgType outside{2, 0};
  bad = params;
  bad.initial_guess = &outside;
  CHECK_THROW(validate_turbo_parameters(bad), std::domain_error);

  bad = params;
  bad.batch_size = 0;
  CHECK_THROW(validate_turbo_parameters(bad), std::invalid_argument);

  bad = params;
  bad.batch_size = 4;
  bad.candidates = 3;
  CHECK_THROW(validate_turbo_parameters(bad), std::invalid_argument);

  bad = params;
  bad.max_training_points = 0;
  CHECK_THROW(validate_turbo_parameters(bad), std::invalid_argument);
}

int main() {
  test_branin<float>();
  test_branin<double>();
  test_six_hump_camel<float>();
  test_six_hump_camel<double>();
  test_hartmann6();
  test_rosenbrock<float>();
  test_rosenbrock<double>();
  test_sphere();
  test_termination();
  test_validation();
  return boost::math::test::report_errors();
}
