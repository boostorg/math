/*
 * Copyright Nick Thompson, 2026
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */

#ifndef BOOST_MATH_BUILD_MODULE
#include <boost/math/optimization/dycors.hpp>
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
using boost::math::optimization::dycors;
using boost::math::optimization::dycors_parameters;
using boost::math::optimization::validate_dycors_parameters;

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

template <class Real> void test_ackley() {
  std::cout << "Testing DYCORS on the Ackley function . . .\n";
  using ArgType = std::array<Real, 2>;
  dycors_parameters<ArgType> params;
  params.lower_bounds = {-5, -3};
  params.upper_bounds = {4, 6};
  params.max_function_calls = 200;
  std::mt19937_64 gen(12345);
  std::vector<std::pair<ArgType, Real>> queries;
  std::atomic<Real> current_minimum_cost = std::numeric_limits<Real>::quiet_NaN();
  auto x = dycors(ackley<Real>, params, gen, std::numeric_limits<Real>::quiet_NaN(), nullptr, &current_minimum_cost,
                  &queries);
  CHECK_LE(ackley(x), Real(1e-2));
  CHECK_EQUAL(queries.size(), size_t(200));
  check_in_bounds(queries, params.lower_bounds, params.upper_bounds);
  // The returned point and the current minimum are those of the best query:
  Real best = std::numeric_limits<Real>::infinity();
  for (auto const &q : queries) {
    CHECK_EQUAL(ackley(q.first), q.second);
    best = (std::min)(best, q.second);
  }
  CHECK_EQUAL(ackley(x), best);
  CHECK_EQUAL(current_minimum_cost.load(), best);
}

template <class Real> void test_rosenbrock() {
  std::cout << "Testing DYCORS on the Rosenbrock function . . .\n";
  using ArgType = std::array<Real, 2>;
  dycors_parameters<ArgType> params;
  params.lower_bounds = {-1.5, -0.5};
  params.upper_bounds = {2.5, 3};
  params.max_function_calls = 500;
  std::mt19937_64 gen(4);
  auto x = dycors(rosenbrock_saddle<Real>, params, gen);
  // DYCORS is meant for many dimensions, and is slow in the curved valley of this one. Over 60 seeds, the median is
  // 0.008 and the worst 0.18.
  CHECK_LE(rosenbrock_saddle(x), Real(0.25));
}

// A ten-dimensional shifted sphere, with the minimizer off the centre of the box.
template <class Real> void test_sphere() {
  std::cout << "Testing DYCORS on a shifted sphere in 10 dimensions . . .\n";
  using ArgType = std::vector<Real>;
  auto f = [](ArgType const &x) {
    Real s = 0;
    for (size_t i = 0; i < x.size(); ++i) {
      Real const t = x[i] - Real(i) / 10;
      s += t * t;
    }
    return s;
  };
  dycors_parameters<ArgType> params;
  params.lower_bounds.assign(10, -3);
  params.upper_bounds.assign(10, 4);
  params.max_function_calls = 300;
  std::mt19937_64 gen(5);
  std::vector<std::pair<ArgType, Real>> queries;
  auto x = dycors(f, params, gen, std::numeric_limits<Real>::quiet_NaN(), nullptr, nullptr, &queries);
  CHECK_EQUAL(x.size(), size_t(10));
  CHECK_LE(f(x), Real(1e-3));
  check_in_bounds(queries, params.lower_bounds, params.upper_bounds);
}

// The sphere of test_functions_for_optimization.hpp is NaN outside the unit ball, which covers most of these boxes.
// It also tests a return type different from the argument type.
void test_nan_and_corner() {
  std::cout << "Testing DYCORS on the sphere, with NaNs and a minimizer at a corner . . .\n";
  using ArgType = std::vector<float>;
  dycors_parameters<ArgType> params;
  params.lower_bounds = {-1, -0.4f, -0.7f};
  params.upper_bounds = {0.6f, 1, 0.3f};
  params.max_function_calls = 150;
  std::mt19937_64 gen(5);
  std::vector<std::pair<ArgType, double>> queries;
  auto x = dycors(sphere, params, gen, std::numeric_limits<double>::quiet_NaN(), nullptr, nullptr, &queries);
  CHECK_LE(sphere(x), 1e-3);
  check_in_bounds(queries, params.lower_bounds, params.upper_bounds);

  // A box that excludes the minimizer, which is then at the corner (0.2, 0.1):
  params.lower_bounds = {0.2f, 0.1f};
  params.upper_bounds = {0.9f, 0.5f};
  queries.clear();
  x = dycors(sphere, params, gen, std::numeric_limits<double>::quiet_NaN(), nullptr, nullptr, &queries);
  CHECK_LE(sphere(x) - 0.05, 1e-4);
  check_in_bounds(queries, params.lower_bounds, params.upper_bounds);
}

void test_determinism() {
  std::cout << "Testing that DYCORS is deterministic for a fixed seed . . .\n";
  using ArgType = std::vector<double>;
  dycors_parameters<ArgType> params;
  params.lower_bounds.assign(5, -5.12);
  params.upper_bounds.assign(5, 5.12);
  params.max_function_calls = 100;
  std::vector<std::pair<ArgType, double>> q1, q2;
  std::mt19937_64 gen1(17);
  auto x1 = dycors(rastrigin<double>, params, gen1, std::numeric_limits<double>::quiet_NaN(), nullptr, nullptr, &q1);
  std::mt19937_64 gen2(17);
  auto x2 = dycors(rastrigin<double>, params, gen2, std::numeric_limits<double>::quiet_NaN(), nullptr, nullptr, &q2);
  CHECK_TRUE(x1 == x2);
  CHECK_TRUE(q1 == q2);
  // The initial design is a symmetric Latin hypercube of 2(d + 1) points: x_{11 - k} = -x_k on this box.
  bool symmetric = true;
  for (size_t k = 0; k < 12; ++k) {
    for (size_t i = 0; i < 5; ++i) {
      symmetric = symmetric && std::abs(q1[k].first[i] + q1[11 - k].first[i]) < 1e-12;
    }
  }
  CHECK_TRUE(symmetric);
}

void test_termination() {
  std::cout << "Testing DYCORS termination: budget, cancellation, target value . . .\n";
  using ArgType = std::array<double, 2>;
  dycors_parameters<ArgType> params;
  params.lower_bounds = {0.5, 0.5};
  params.upper_bounds = {2.048, 2.048};
  params.max_function_calls = 37;
  std::mt19937_64 gen(6);
  std::vector<std::pair<ArgType, double>> queries;
  auto const nan = std::numeric_limits<double>::quiet_NaN();
  dycors(rosenbrock_saddle<double>, params, gen, nan, nullptr, nullptr, &queries);
  CHECK_EQUAL(queries.size(), size_t(37));

  // The default budget is 100 d, and the budget can be just the initial design of 2(d + 1) points:
  params.max_function_calls = 0;
  queries.clear();
  dycors(rosenbrock_saddle<double>, params, gen, nan, nullptr, nullptr, &queries);
  CHECK_EQUAL(queries.size(), size_t(200));
  params.max_function_calls = 6;
  queries.clear();
  dycors(rosenbrock_saddle<double>, params, gen, nan, nullptr, nullptr, &queries);
  CHECK_EQUAL(queries.size(), size_t(6));
  params.max_function_calls = 0;

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
  dycors(f, params, gen, nan, &cancel, nullptr, &queries);
  CHECK_EQUAL(queries.size(), size_t(10));
  // Cancelled before the first call: no queries, and the midpoint of the box comes back.
  queries.clear();
  auto x = dycors(f, params, gen, nan, &cancel, nullptr, &queries);
  CHECK_EQUAL(queries.size(), size_t(0));
  CHECK_EQUAL(x[0], (params.lower_bounds[0] + params.upper_bounds[0]) / 2);
  CHECK_EQUAL(x[1], (params.lower_bounds[1] + params.upper_bounds[1]) / 2);

  // The run stops at the first call that reaches the target value:
  double const target = 0.5;
  queries.clear();
  x = dycors(rosenbrock_saddle<double>, params, gen, target, nullptr, nullptr, &queries);
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
  x = dycors(rosenbrock_saddle<double>, params, gen, 0.0, nullptr, nullptr, &queries);
  CHECK_EQUAL(queries.size(), size_t(1));
  CHECK_EQUAL(x[0], 1.0);
  CHECK_EQUAL(x[1], 1.0);

  // Otherwise the guess, exactly as given, is followed by the whole symmetric design, here of 6 points:
  guess = {0.7, 1.9};
  params.max_function_calls = 7;
  queries.clear();
  dycors(rosenbrock_saddle<double>, params, gen, nan, nullptr, nullptr, &queries);
  CHECK_EQUAL(queries.size(), size_t(7));
  CHECK_TRUE(queries[0].first == guess);
  bool symmetric = true;
  for (size_t k = 1; k < 7; ++k) {
    for (size_t i = 0; i < 2; ++i) {
      double const mid = (params.lower_bounds[i] + params.upper_bounds[i]) / 2;
      symmetric = symmetric && std::abs(queries[k].first[i] + queries[7 - k].first[i] - 2 * mid) < 1e-12;
    }
  }
  CHECK_TRUE(symmetric);
}

void test_validation() {
  std::cout << "Testing DYCORS parameter validation . . .\n";
  using ArgType = std::array<double, 2>;
  dycors_parameters<ArgType> params;
  params.lower_bounds = {-1, -1};
  params.upper_bounds = {1, 1};
  validate_dycors_parameters(params);
  std::mt19937_64 gen(7);

  auto bad = params;
  bad.upper_bounds[1] = -1;
  CHECK_THROW(validate_dycors_parameters(bad), std::domain_error);
  auto run = [&]() { return dycors(rosenbrock_saddle<double>, bad, gen); };
  CHECK_THROW(run(), std::domain_error);

  // DYCORS scales the box to the unit cube, so it must have positive width:
  bad = params;
  bad.upper_bounds[0] = bad.lower_bounds[0];
  CHECK_THROW(validate_dycors_parameters(bad), std::domain_error);

  bad = params;
  bad.upper_bounds[0] = std::numeric_limits<double>::infinity();
  CHECK_THROW(validate_dycors_parameters(bad), std::domain_error);

  ArgType outside{2, 0};
  bad = params;
  bad.initial_guess = &outside;
  CHECK_THROW(validate_dycors_parameters(bad), std::domain_error);

  // The linear tail needs at least 2d = 4 points of a symmetric design:
  bad = params;
  bad.initial_points = 3;
  CHECK_THROW(validate_dycors_parameters(bad), std::invalid_argument);
  bad.initial_points = 4;
  validate_dycors_parameters(bad);
  // The budget must cover the initial design, of 4 points here and 2(d + 1) = 6 by default:
  bad.max_function_calls = 3;
  CHECK_THROW(validate_dycors_parameters(bad), std::invalid_argument);
  bad.initial_points = 0;
  bad.max_function_calls = 5;
  CHECK_THROW(validate_dycors_parameters(bad), std::invalid_argument);
  bad.initial_points = 4;
  bad.max_function_calls = 0;
  CHECK_EQUAL(dycors(rosenbrock_saddle<double>, bad, gen).size(), size_t(2));
}

int main() {
  test_ackley<float>();
  test_ackley<double>();
  test_rosenbrock<float>();
  test_rosenbrock<double>();
  test_sphere<float>();
  test_sphere<double>();
  test_nan_and_corner();
  test_determinism();
  test_termination();
  test_validation();
  return boost::math::test::report_errors();
}
