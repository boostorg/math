/*
 * Copyright Nick Thompson, 2024
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */

#ifndef BOOST_MATH_BUILD_MODULE
#include <cmath>
#include <boost/math/optimization/jso.hpp>
#else
import boost.math;
#endif

#include "math_unit_test.hpp"
#include "test_functions_for_optimization.hpp"
#ifndef BOOST_MATH_BUILD_MODULE
#include <random>
#include <limits>
#endif

using boost::math::optimization::jso;
using boost::math::optimization::jso_parameters;
using boost::math::optimization::detail::weighted_lehmer_mean;

void test_default_population_with_many_threads() {
  using ArgType = std::array<double, 2>;
  auto jso_params = jso_parameters<ArgType>();
  jso_params.lower_bounds = {-5, -5};
  jso_params.upper_bounds = {5, 5};
  jso_params.threads = 128;

  boost::math::optimization::validate_jso_parameters(jso_params);
  CHECK_GE(jso_params.initial_population_size, jso_params.threads);
}

// A budget smaller than the initial population used to end in a remainder by zero.
void test_budget_smaller_than_population() {
  using ArgType = std::array<double, 2>;
  auto jso_params = jso_parameters<ArgType>();
  jso_params.lower_bounds = {-5, -5};
  jso_params.upper_bounds = {5, 5};
  jso_params.threads = 1;
  jso_params.initial_population_size = 30;
  jso_params.max_function_evaluations = 10;
  std::mt19937_64 gen(12345);
  auto f = [](ArgType const &x) { return x[0] * x[0] + x[1] * x[1]; };
  CHECK_THROW(jso(f, jso_params, gen), std::invalid_argument);
  // The default population of two variables is 39, and the default budget is larger:
  jso_params.initial_population_size = 0;
  jso_params.max_function_evaluations = 38;
  CHECK_THROW(jso(f, jso_params, gen), std::invalid_argument);
  jso_params.initial_population_size = 50;
  jso_params.max_function_evaluations = 0;
  CHECK_NO_THROW(boost::math::optimization::validate_jso_parameters(jso_params));
  // One generation after the initial population:
  jso_params.max_function_evaluations = 50;
  std::vector<std::pair<ArgType, double>> queries;
  auto x = jso(f, jso_params, gen, std::numeric_limits<double>::quiet_NaN(), nullptr, nullptr, &queries);
  CHECK_LE(f(x), 50.0);
  CHECK_LE(queries.size(), size_t(100));
}

void test_weighted_lehmer_mean() {
  size_t n = 50;
  std::vector<double> weights(n, 1.0);
  std::vector<double> values(n, 2.5);
  // Technically, this is not a fully general weighted Lehmer mean,
  // but just a weighted contraharmonic mean.
  // So we have a few more invariants available to us:
  CHECK_ULP_CLOSE(2.5, weighted_lehmer_mean(values, weights), n);
  std::mt19937_64 gen(12345);
  std::uniform_real_distribution<double> unif(std::numeric_limits<double>::epsilon(),1);
  for (size_t i = 0; i < n; ++i) {
    weights[i] = unif(gen);
    values[i] = unif(gen);
  }
  auto mean = weighted_lehmer_mean(values, weights);
  CHECK_LE(mean, 1.0);
  CHECK_LE(std::numeric_limits<double>::epsilon(), mean);
}

template <class Real> void test_ackley() {
  std::cout << "Testing jSO on Ackley function . . .\n";
  using ArgType = std::array<Real, 2>;
  auto jso_params = jso_parameters<ArgType>();
  jso_params.lower_bounds = {-5, -5};
  jso_params.upper_bounds = {5, 5};

  std::mt19937_64 gen(12345);
  auto local_minima = jso(ackley<Real>, jso_params, gen);
  CHECK_LE(std::abs(local_minima[0]), std::sqrt(std::numeric_limits<Real>::epsilon()));
  CHECK_LE(std::abs(local_minima[1]), std::sqrt(std::numeric_limits<Real>::epsilon()));

  // Does it work with a lambda?
  auto ack = [](std::array<Real, 2> const &x) { return ackley<Real>(x); };
  local_minima = jso(ack, jso_params, gen);
  CHECK_LE(std::abs(local_minima[0]), std::sqrt(std::numeric_limits<Real>::epsilon()));
  CHECK_LE(std::abs(local_minima[1]), std::sqrt(std::numeric_limits<Real>::epsilon()));

  // Test that if an intial guess is the exact solution, the returned solution is the exact solution:
  std::array<Real, 2> initial_guess{0, 0};
  jso_params.initial_guess = &initial_guess;
  local_minima = jso(ack, jso_params, gen);

  using std::fabs;
  CHECK_LE(fabs(local_minima[0]), 128 * boost::math::tools::epsilon<Real>());
  CHECK_LE(fabs(local_minima[1]), 128 * boost::math::tools::epsilon<Real>());
}

template <class Real> void test_rosenbrock_saddle() {
  std::cout << "Testing jSO on Rosenbrock saddle . . .\n";
  using ArgType = std::array<Real, 2>;
  jso_parameters<ArgType> jso_params { };
  jso_params.lower_bounds = {0.5, 0.5};
  jso_params.upper_bounds = {2.048, 2.048};
  std::mt19937_64 gen(234568);
  auto local_minima = jso(rosenbrock_saddle<Real>, jso_params, gen);

  CHECK_ABSOLUTE_ERROR(Real(1), local_minima[0], std::sqrt(std::numeric_limits<Real>::epsilon()));
  CHECK_ABSOLUTE_ERROR(Real(1), local_minima[1], std::sqrt(std::numeric_limits<Real>::epsilon()));

  // Does cancellation work?
  std::atomic<bool> cancel = true;
  gen.seed(12345);
  local_minima =
      jso(rosenbrock_saddle<Real>, jso_params, gen, std::numeric_limits<Real>::quiet_NaN(), &cancel);
  CHECK_GE(std::abs(local_minima[0] - Real(1)), std::sqrt(std::numeric_limits<Real>::epsilon()));
}



template <class Real> void test_rastrigin() {
  std::cout << "Testing jSO on Rastrigin function (global minimum = (0,0,...,0))\n";
  using ArgType = std::vector<Real>;
  auto jso_params = jso_parameters<ArgType>();
  jso_params.lower_bounds.resize(3, static_cast<Real>(-5.12));
  jso_params.upper_bounds.resize(3, static_cast<Real>(5.12));
  jso_params.initial_population_size = 5000;
  jso_params.max_function_evaluations = 1000000;
  std::mt19937_64 gen(34567);

  // By definition, the value of the function which a target value is provided must be <= target_value.
  Real target_value = 1e-3;
  auto local_minima = jso(rastrigin<Real>, jso_params, gen, target_value);
  CHECK_LE(rastrigin(local_minima), target_value);
}

// Tests NaN return types and return type != input type:
void test_sphere() {
  std::cout << "Testing jSO on sphere . . .\n";
  using ArgType = std::vector<float>;
  auto jso_params = jso_parameters<ArgType>();
  jso_params.lower_bounds.resize(8, -1);
  jso_params.upper_bounds.resize(8, 1);
  std::mt19937_64 gen(56789);
  auto local_minima = jso(sphere, jso_params, gen);
  for (auto x : local_minima) {
    CHECK_ABSOLUTE_ERROR(0.0f, x, 2e-4f);
  }
}

template<typename Real>
void test_three_hump_camel() {
  std::cout << "Testing jSO on three hump camel . . .\n";
  using ArgType = std::array<Real, 2>;
  auto jso_params = jso_parameters<ArgType>();
  jso_params.lower_bounds[0] = -5.0;
  jso_params.lower_bounds[1] = -5.0;
  jso_params.upper_bounds[0] = 5.0;
  jso_params.upper_bounds[1] = 5.0;
  std::mt19937_64 gen(56789);
  auto local_minima = jso(three_hump_camel<Real>, jso_params, gen);
  for (auto x : local_minima) {
    CHECK_ABSOLUTE_ERROR(0.0f, x, 2e-4f);
  }
}

template<typename Real>
void test_beale() {
  std::cout << "Testing jSO on the Beale function . . .\n";
  using ArgType = std::array<Real, 2>;
  auto jso_params = jso_parameters<ArgType>();
  jso_params.lower_bounds[0] = -5.0;
  jso_params.lower_bounds[1] = -5.0;
  jso_params.upper_bounds[0]= 5.0;
  jso_params.upper_bounds[1]= 5.0;
  std::mt19937_64 gen(56789);
  auto local_minima = jso(beale<Real>, jso_params, gen);
  CHECK_ABSOLUTE_ERROR(Real(3), local_minima[0], Real(2e-4));
  CHECK_ABSOLUTE_ERROR(Real(1)/Real(2), local_minima[1], Real(2e-4));
}

int main() {
  test_default_population_with_many_threads();
  test_ackley<float>();
  test_ackley<double>();
  test_rosenbrock_saddle<double>();
  test_rastrigin<double>();
  test_three_hump_camel<float>();
  test_beale<double>();
  test_sphere();
  test_weighted_lehmer_mean();
  test_budget_smaller_than_population();
  return boost::math::test::report_errors();
}
