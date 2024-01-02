/*
 * Copyright Nick Thompson, 2023
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */

#include "math_unit_test.hpp"
#include <boost/math/constants/constants.hpp>
#include <boost/math/tools/differential_evolution.hpp>
#include <random>

using boost::math::constants::e;
using boost::math::constants::two_pi;
using boost::math::tools::differential_evolution;
using boost::math::tools::differential_evolution_parameters;
using std::cbrt;
using std::cos;
using std::exp;
using std::sqrt;

// Taken from: https://en.wikipedia.org/wiki/Test_functions_for_optimization
template <typename Real> Real ackley(std::array<Real, 2> const &v) {
  Real x = v[0];
  Real y = v[1];
  Real arg1 = -sqrt((x * x + y * y) / 2) / 5;
  Real arg2 = cos(two_pi<Real>() * x) + cos(two_pi<Real>() * y);
  return -20 * exp(arg1) - exp(arg2 / 2) + 20 + e<Real>();
}

template <class Real> void test_ackley() {
  using ArgType = std::array<Real, 2>;
  auto de_params = differential_evolution_parameters<ArgType>();
  de_params.lower_bounds = {-5, -5};
  de_params.upper_bounds = {5, 5};

  std::mt19937_64 gen(12345);
  auto local_minima = differential_evolution(ackley<Real>, de_params, gen);
  CHECK_LE(std::abs(local_minima[0]), 10 * std::numeric_limits<Real>::epsilon());
  CHECK_LE(std::abs(local_minima[1]), 10 * std::numeric_limits<Real>::epsilon());

  // Does it work with a lambda?
  auto ack = [](std::array<Real, 2> const &x) { return ackley<Real>(x); };
  local_minima = differential_evolution(ack, de_params, gen);
  CHECK_LE(std::abs(local_minima[0]), 10 * std::numeric_limits<Real>::epsilon());
  CHECK_LE(std::abs(local_minima[1]), 10 * std::numeric_limits<Real>::epsilon());

  // Test that if an intial guess is the exact solution, the returned solution is the exact solution:
  std::array<Real, 2> initial_guess{0, 0};
  de_params.initial_guess = &initial_guess;
  local_minima = differential_evolution(ack, de_params, gen);
  CHECK_EQUAL(local_minima[0], Real(0));
  CHECK_EQUAL(local_minima[1], Real(0));
}

template <typename Real> auto rosenbrock_saddle(std::array<Real, 2> const &v) {
  auto x = v[0];
  auto y = v[1];
  return 100 * (x * x - y) * (x * x - y) + (1 - x) * (1 - x);
}

template <class Real> void test_rosenbrock_saddle() {
  using ArgType = std::array<Real, 2>;
  auto de_params = differential_evolution_parameters<ArgType>();
  de_params.lower_bounds = {0.5, 0.5};
  de_params.upper_bounds = {2.048, 2.048};
  std::mt19937_64 gen(234568);
  auto local_minima = differential_evolution(rosenbrock_saddle<Real>, de_params, gen);

  CHECK_ABSOLUTE_ERROR(Real(1), local_minima[0], 10 * std::numeric_limits<Real>::epsilon());
  CHECK_ABSOLUTE_ERROR(Real(1), local_minima[1], 10 * std::numeric_limits<Real>::epsilon());

  // Does cancellation work?
  std::atomic<bool> cancel = true;
  gen.seed(12345);
  local_minima =
      differential_evolution(rosenbrock_saddle<Real>, de_params, gen, std::numeric_limits<Real>::quiet_NaN(), &cancel);
  CHECK_GE(std::abs(local_minima[0] - Real(1)), std::sqrt(std::numeric_limits<Real>::epsilon()));
}

template <class Real> Real rastrigin(std::vector<Real> const &v) {
  Real A = 10;
  Real y = 10 * v.size();
  for (auto x : v) {
    y += x * x - A * cos(two_pi<Real>() * x);
  }
  return y;
}

template <class Real> void test_rastrigin() {
  using ArgType = std::vector<Real>;
  auto de_params = differential_evolution_parameters<ArgType>();
  de_params.lower_bounds.resize(8, static_cast<Real>(-5.12));
  de_params.upper_bounds.resize(8, static_cast<Real>(5.12));
  std::mt19937_64 gen(34567);
  auto local_minima = differential_evolution(rastrigin<Real>, de_params, gen);
  for (auto x : local_minima) {
    CHECK_ABSOLUTE_ERROR(x, Real(0), Real(2e-4));
  }

  // By definition, the value of the function which a target value is provided must be <= target_value.
  Real target_value = 1e-3;
  local_minima = differential_evolution(rastrigin<Real>, de_params, gen, target_value);
  CHECK_LE(rastrigin(local_minima), target_value);
}

double sphere(std::vector<float> const &v) {
  double r = 0.0;
  for (auto x : v) {
    double x_ = static_cast<double>(x);
    r += x_ * x_;
  }
  if (r >= 1) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  return r;
}

// Tests NaN return types and return type != input type:
void test_sphere() {
  using ArgType = std::vector<float>;
  auto de_params = differential_evolution_parameters<ArgType>();
  de_params.lower_bounds.resize(8, -1);
  de_params.upper_bounds.resize(8, 1);
  de_params.NP *= 10;
  de_params.max_generations *= 10;
  std::mt19937_64 gen(56789);
  auto local_minima = differential_evolution(sphere, de_params, gen);
  for (auto x : local_minima) {
    CHECK_ABSOLUTE_ERROR(0.0f, x, 2e-4f);
  }
}

#define GCC_COMPILER (defined(__GNUC__) && !defined(__clang__))

int main() {
  // GCC<=8 rejects the function call syntax we use here.
  // Just do a workaround:
#if !defined(GCC_COMPILER) || (defined(GCC_COMPILER) && __GNUC__ >= 9)
  test_ackley<float>();
  test_ackley<double>();
  test_rosenbrock_saddle<double>();
  test_rastrigin<float>();
#endif
  test_sphere();
  return boost::math::test::report_errors();
}
