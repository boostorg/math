/*
 * Copyright Nick Thompson, 2026
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */

#ifndef BOOST_MATH_BUILD_MODULE
#include <boost/math/optimization/multilevel_coordinate_search.hpp>
#else
import boost.math;
#endif

#include "math_unit_test.hpp"
#include "test_functions_for_optimization.hpp"
#ifndef BOOST_MATH_BUILD_MODULE
#include <limits>
#include <stdexcept>
#endif
using boost::math::optimization::multilevel_coordinate_search;
using boost::math::optimization::multilevel_coordinate_search_parameters;
using boost::math::optimization::validate_multilevel_coordinate_search_parameters;


// The bounds are asymmetric throughout, so that the midpoint of the box (the default base point) is not the minimizer.
template <class Real> void test_ackley() {
  std::cout << "Testing multilevel coordinate search on Ackley function . . .\n";
  using ArgType = std::array<Real, 2>;
  auto mcs_params = multilevel_coordinate_search_parameters<ArgType>();
  mcs_params.lower_bounds = {-5, -4};
  mcs_params.upper_bounds = {3, 5};
  auto local_minima = multilevel_coordinate_search(ackley<Real>, mcs_params);
  CHECK_LE(std::abs(local_minima[0]), std::sqrt(std::numeric_limits<Real>::epsilon()));
  CHECK_LE(std::abs(local_minima[1]), std::sqrt(std::numeric_limits<Real>::epsilon()));

  // Does it work with a lambda?
  auto ack = [](std::array<Real, 2> const &x) { return ackley<Real>(x); };
  local_minima = multilevel_coordinate_search(ack, mcs_params);
  CHECK_LE(std::abs(local_minima[0]), std::sqrt(std::numeric_limits<Real>::epsilon()));
  CHECK_LE(std::abs(local_minima[1]), std::sqrt(std::numeric_limits<Real>::epsilon()));

  // Test that if an intial guess is the exact solution, the returned solution is the exact solution:
  std::array<Real, 2> initial_guess{0, 0};
  mcs_params.initial_guess = &initial_guess;
  local_minima = multilevel_coordinate_search(ack, mcs_params);
  CHECK_EQUAL(local_minima[0], Real(0));
  CHECK_EQUAL(local_minima[1], Real(0));

  std::atomic<bool> cancel = false;
  Real target_value = ackley<Real>(initial_guess);
  std::atomic<Real> current_minimum_cost = std::numeric_limits<Real>::quiet_NaN();
  // Test query storage:
  std::vector<std::pair<ArgType, Real>> queries;
  local_minima = multilevel_coordinate_search(ack, mcs_params, target_value, &cancel, &current_minimum_cost, &queries);
  CHECK_EQUAL(local_minima[0], Real(0));
  CHECK_EQUAL(local_minima[1], Real(0));
  CHECK_EQUAL(current_minimum_cost.load(), target_value);
  // The initial guess reaches the target, so it is the only call:
  CHECK_EQUAL(queries.size(), size_t(1));
  for (auto const & q : queries) {
    CHECK_EQUAL(ackley<Real>(q.first), q.second);
  }
}


template <class Real> void test_rosenbrock_saddle() {
  std::cout << "Testing multilevel coordinate search on Rosenbrock saddle . . .\n";
  using ArgType = std::array<Real, 2>;
  auto mcs_params = multilevel_coordinate_search_parameters<ArgType>();
  mcs_params.lower_bounds = {0.5, 0.5};
  mcs_params.upper_bounds = {2.048, 2.048};
  auto local_minima = multilevel_coordinate_search(rosenbrock_saddle<Real>, mcs_params);
  CHECK_ABSOLUTE_ERROR(Real(1), local_minima[0], Real(1e-7));
  CHECK_ABSOLUTE_ERROR(Real(1), local_minima[1], Real(1e-7));

  // Does cancellation work? Cancel from inside the cost function after a few calls:
  std::atomic<bool> cancel = false;
  size_t calls = 0;
  auto f = [&](ArgType const &x) {
    if (++calls == 10) {
      cancel = true;
    }
    return rosenbrock_saddle<Real>(x);
  };
  std::vector<std::pair<ArgType, Real>> queries;
  local_minima = multilevel_coordinate_search(f, mcs_params, std::numeric_limits<Real>::quiet_NaN(), &cancel,
                                              nullptr, &queries);
  CHECK_EQUAL(queries.size(), size_t(10));
  CHECK_GE(std::abs(local_minima[0] - Real(1)), std::sqrt(std::numeric_limits<Real>::epsilon()));

  // Cancelled before the first call: no queries, and the midpoint of the box comes back.
  queries.clear();
  local_minima = multilevel_coordinate_search(f, mcs_params, std::numeric_limits<Real>::quiet_NaN(), &cancel,
                                              nullptr, &queries);
  CHECK_EQUAL(queries.size(), size_t(0));
  CHECK_EQUAL(local_minima[0], (mcs_params.lower_bounds[0] + mcs_params.upper_bounds[0]) / 2);

  // max_function_calls is respected:
  mcs_params.max_function_calls = 50;
  queries.clear();
  multilevel_coordinate_search(rosenbrock_saddle<Real>, mcs_params, std::numeric_limits<Real>::quiet_NaN(), nullptr,
                               nullptr, &queries);
  CHECK_LE(queries.size(), size_t(50));
}


template <class Real> void test_rastrigin() {
  std::cout << "Testing multilevel coordinate search on Rastrigin function (global minimum = (0,0,...,0))\n";
  using ArgType = std::vector<Real>;
  auto mcs_params = multilevel_coordinate_search_parameters<ArgType>();
  mcs_params.lower_bounds = {-5.12, -4, -3};
  mcs_params.upper_bounds = {4, 5.12, 5};
  std::vector<std::pair<ArgType, Real>> queries;
  auto local_minima = multilevel_coordinate_search(rastrigin<Real>, mcs_params, std::numeric_limits<Real>::quiet_NaN(),
                                                   nullptr, nullptr, &queries);
  CHECK_EQUAL(local_minima.size(), size_t(3));
  CHECK_LE(rastrigin(local_minima), Real(1e-8));
  size_t const calls = queries.size();

  // By definition, the value of the function which a target value is provided must be <= target_value.
  Real target_value = 2.0;
  queries.clear();
  local_minima = multilevel_coordinate_search(rastrigin<Real>, mcs_params, target_value, nullptr, nullptr, &queries);
  CHECK_LE(rastrigin(local_minima), target_value);
  CHECK_LE(queries.size() + 1, calls);
}


// Tests NaN return types and return type != input type:
void test_sphere() {
  std::cout << "Testing multilevel coordinate search on sphere . . .\n";
  using ArgType = std::vector<float>;
  auto mcs_params = multilevel_coordinate_search_parameters<ArgType>();
  // The sphere function is NaN outside the unit ball, which includes the corners and most of the faces of this box:
  mcs_params.lower_bounds.resize(4, -1);
  mcs_params.upper_bounds.resize(4, 0.6f);
  auto local_minima = multilevel_coordinate_search(sphere, mcs_params);
  for (auto x : local_minima) {
    CHECK_ABSOLUTE_ERROR(0.0f, x, 0.001f);
  }
  // An easy function converges without the local search, too:
  mcs_params.local_search = false;
  local_minima = multilevel_coordinate_search(sphere, mcs_params);
  for (auto x : local_minima) {
    CHECK_ABSOLUTE_ERROR(0.0f, x, 0.02f);
  }
}


template<typename Real>
void test_beale() {
  std::cout << "Testing multilevel coordinate search on the Beale function . . .\n";
  using ArgType = std::array<Real, 2>;
  auto mcs_params = multilevel_coordinate_search_parameters<ArgType>();
  mcs_params.lower_bounds = {-5, -5};
  mcs_params.upper_bounds = {5, 5};
  auto local_minima = multilevel_coordinate_search(beale<Real>, mcs_params);
  CHECK_ABSOLUTE_ERROR(Real(3), local_minima[0], Real(1e-8));
  CHECK_ABSOLUTE_ERROR(Real(1)/Real(2), local_minima[1], Real(1e-8));
}


// Two of the Dixon-Szego test functions of Section 7 of the paper, on the boxes used there.
template<typename Real>
void test_six_hump_camel() {
  std::cout << "Testing multilevel coordinate search on the six-hump camel function . . .\n";
  using ArgType = std::array<Real, 2>;
  auto camel = [](ArgType const &v) {
    Real x = v[0];
    Real y = v[1];
    return (4 - Real(21) / 10 * x * x + x * x * x * x / 3) * x * x + x * y + (-4 + 4 * y * y) * y * y;
  };
  auto mcs_params = multilevel_coordinate_search_parameters<ArgType>();
  mcs_params.lower_bounds = {-3, -2};
  mcs_params.upper_bounds = {3, 2};
  auto local_minima = multilevel_coordinate_search(camel, mcs_params);
  CHECK_ABSOLUTE_ERROR(Real(-1.0316284534898774), camel(local_minima), Real(1e-12));
}


template<typename Real>
void test_branin() {
  std::cout << "Testing multilevel coordinate search on the Branin function . . .\n";
  using ArgType = std::array<Real, 2>;
  auto branin = [](ArgType const &v) {
    using boost::math::constants::pi;
    using std::cos;
    Real x = v[0];
    Real y = v[1];
    Real t = y - Real(5.1) / (4 * pi<Real>() * pi<Real>()) * x * x + 5 * x / pi<Real>() - 6;
    return t * t + 10 * (1 - 1 / (8 * pi<Real>())) * cos(x) + 10;
  };
  auto mcs_params = multilevel_coordinate_search_parameters<ArgType>();
  mcs_params.lower_bounds = {-5, 0};
  mcs_params.upper_bounds = {10, 15};
  auto local_minima = multilevel_coordinate_search(branin, mcs_params);
  CHECK_ABSOLUTE_ERROR(Real(0.39788735772973816), branin(local_minima), Real(1e-12));
}


void test_validation() {
  std::cout << "Testing multilevel coordinate search parameter validation . . .\n";
  using ArgType = std::array<double, 2>;
  auto mcs_params = multilevel_coordinate_search_parameters<ArgType>();
  mcs_params.lower_bounds = {-1, -1};
  mcs_params.upper_bounds = {1, 1};
  validate_multilevel_coordinate_search_parameters(mcs_params);

  auto bad_params = mcs_params;
  bad_params.upper_bounds[1] = -1;
  CHECK_THROW(validate_multilevel_coordinate_search_parameters(bad_params), std::domain_error);
  auto run = [&]() { return multilevel_coordinate_search(rosenbrock_saddle<double>, bad_params); };
  CHECK_THROW(run(), std::domain_error);

  bad_params = mcs_params;
  bad_params.line_search_points = 2;
  CHECK_THROW(validate_multilevel_coordinate_search_parameters(bad_params), std::invalid_argument);

  bad_params = mcs_params;
  bad_params.max_level = 1;
  CHECK_THROW(validate_multilevel_coordinate_search_parameters(bad_params), std::invalid_argument);
}

int main() {
  test_ackley<float>();
  test_ackley<double>();
  test_rosenbrock_saddle<double>();
  test_rastrigin<double>();
  test_sphere();
  test_beale<double>();
  test_six_hump_camel<double>();
  test_branin<double>();
  test_validation();
  return boost::math::test::report_errors();
}
