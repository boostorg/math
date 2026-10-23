/*
 * Copyright Nick Thompson, 2024
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */

#ifndef BOOST_MATH_BUILD_MODULE
#include <boost/math/optimization/cma_es.hpp>
#else
import boost.math;
#endif

#include "math_unit_test.hpp"
#include "test_functions_for_optimization.hpp"
#ifndef BOOST_MATH_BUILD_MODULE
#include <array>
#include <algorithm>
#include <cmath>
#include <random>
#include <limits>
#include <vector>
#endif

using std::abs;
using boost::math::optimization::cma_es;
using boost::math::optimization::cma_es_parameters;

template <class Real> void test_ackley() {
  std::cout << "Testing CMA-ES on Ackley function . . .\n";
  using ArgType = std::array<Real, 2>;
  auto params = cma_es_parameters<ArgType>();
  params.lower_bounds = {-5, -5};
  params.upper_bounds = {5, 5};

  std::mt19937_64 gen(12345);
  auto local_minima = cma_es(ackley<Real>, params, gen);
  CHECK_LE(std::abs(local_minima[0]), Real(0.1));
  CHECK_LE(std::abs(local_minima[1]), Real(0.1));

  // Does it work with a lambda?
  auto ack = [](std::array<Real, 2> const &x) { return ackley<Real>(x); };
  local_minima = cma_es(ack, params, gen);
  CHECK_LE(std::abs(local_minima[0]), Real(0.1));
  CHECK_LE(std::abs(local_minima[1]), Real(0.1));

  // Test that if an intial guess is the exact solution, the returned solution is the exact solution:
  std::array<Real, 2> initial_guess{0, 0};
  params.initial_guess = &initial_guess;
  local_minima = cma_es(ack, params, gen);
  CHECK_EQUAL(local_minima[0], Real(0));
  CHECK_EQUAL(local_minima[1], Real(0));

  std::atomic<bool> cancel = false;
  Real target_value = 0.0;
  std::atomic<Real> current_minimum_cost = std::numeric_limits<Real>::quiet_NaN();
  // Test query storage:
  std::vector<std::pair<ArgType, Real>> queries;
  local_minima = cma_es(ack, params, gen, target_value, &cancel, &current_minimum_cost, &queries);
  CHECK_EQUAL(local_minima[0], Real(0));
  CHECK_EQUAL(local_minima[1], Real(0));
  CHECK_LE(size_t(1), queries.size());
  for (auto const & q : queries) {
    // Not exact: with -ffp-contract=fast, ackley can be fused differently where it is inlined.
    auto expected = ackley<Real>(q.first);
    CHECK_ULP_CLOSE(expected, q.second, 4);
  }
}


template <class Real> void test_rosenbrock_saddle() {
  std::cout << "Testing CMA-ES on Rosenbrock saddle . . .\n";
  using ArgType = std::array<Real, 2>;
  auto params = cma_es_parameters<ArgType>();
  params.lower_bounds = {0.5, 0.5};
  params.upper_bounds = {2.048, 2.048};
  params.max_generations = 2000;
  std::mt19937_64 gen(234568);
  auto local_minima = cma_es(rosenbrock_saddle<Real>, params, gen);

  CHECK_ABSOLUTE_ERROR(Real(1), local_minima[0], Real(0.05));
  CHECK_ABSOLUTE_ERROR(Real(1), local_minima[1], Real(0.05));

  // Does cancellation work?
  std::atomic<bool> cancel = true;
  gen.seed(12345);
  local_minima =
      cma_es(rosenbrock_saddle<Real>, params, gen, std::numeric_limits<Real>::quiet_NaN(), &cancel);
  CHECK_GE(std::abs(local_minima[0] - Real(1)), std::sqrt(std::numeric_limits<Real>::epsilon()));
}


template <class Real> void test_rastrigin() {
  std::cout << "Testing CMA-ES on Rastrigin function (global minimum = (0,0,...,0))\n";
  using ArgType = std::vector<Real>;
  auto params = cma_es_parameters<ArgType>();
  params.lower_bounds.resize(3, static_cast<Real>(-5.12));
  params.upper_bounds.resize(3, static_cast<Real>(5.12));
  params.max_generations = 1000000;
  params.population_size = 100;
  std::mt19937_64 gen(34567);

  // By definition, the value of the function which a target value is provided must be <= target_value.
  Real target_value = 2.0;
  auto local_minima = cma_es(rastrigin<Real>, params, gen, target_value);
  CHECK_LE(rastrigin(local_minima), target_value);
}


// Tests NaN return types and return type != input type:
void test_sphere() {
  std::cout << "Testing CMA-ES on sphere . . .\n";
  using ArgType = std::vector<float>;
  auto params = cma_es_parameters<ArgType>();
  params.lower_bounds.resize(4, -1);
  params.upper_bounds.resize(4, 1);
  params.max_generations = 100000;
  std::mt19937_64 gen(56789);
  auto local_minima = cma_es(sphere, params, gen, 1e-6f);
  for (auto x : local_minima) {
    CHECK_ABSOLUTE_ERROR(0.0f, x, 0.5f);
  }
}


template<typename Real>
void test_three_hump_camel() {
  std::cout << "Testing CMA-ES on three hump camel . . .\n";
  using ArgType = std::array<Real, 2>;
  auto params = cma_es_parameters<ArgType>();
  params.lower_bounds[0] = -5.0;
  params.lower_bounds[1] = -5.0;
  params.upper_bounds[0] = 5.0;
  params.upper_bounds[1] = 5.0;
  std::mt19937_64 gen(56789);
  auto local_minima = cma_es(three_hump_camel<Real>, params, gen);
  for (auto x : local_minima) {
    CHECK_ABSOLUTE_ERROR(0.0f, x, 0.2f);
  }
}


template<typename Real>
void test_beale() {
  std::cout << "Testing CMA-ES on the Beale function . . .\n";
  using ArgType = std::array<Real, 2>;
  auto params = cma_es_parameters<ArgType>();
  params.lower_bounds[0] = -5.0;
  params.lower_bounds[1] = -5.0;
  params.upper_bounds[0]= 5.0;
  params.upper_bounds[1]= 5.0;
  std::mt19937_64 gen(56789);
  auto local_minima = cma_es(beale<Real>, params, gen);
  CHECK_ABSOLUTE_ERROR(Real(3), local_minima[0], Real(0.1));
  CHECK_ABSOLUTE_ERROR(Real(1)/Real(2), local_minima[1], Real(0.1));
}

// The unconstrained minimum of the quadratic is outside the bounds,
// and the constrained minimum on [1, 2]^4 is at (1, 1, 1, 1):
template<typename Real>
void test_bounds() {
  std::cout << "Testing that CMA-ES respects the bounds . . .\n";
  using ArgType = std::vector<Real>;
  auto params = cma_es_parameters<ArgType>();
  params.lower_bounds.resize(4, Real(1));
  params.upper_bounds.resize(4, Real(2));
  params.max_generations = 200;
  std::mt19937_64 gen(12345);
  size_t calls = 0;
  size_t out_of_bounds_calls = 0;
  auto f = [&](ArgType const & x) {
    ++calls;
    for (size_t i = 0; i < x.size(); ++i) {
      if (!(params.lower_bounds[i] <= x[i] && x[i] <= params.upper_bounds[i])) {
        ++out_of_bounds_calls;
        break;
      }
    }
    auto y = x;
    return quadratic(y);
  };
  auto local_minima = cma_es(f, params, gen);
  CHECK_LE(size_t(1), calls);
  CHECK_EQUAL(out_of_bounds_calls, size_t(0));
  for (auto x : local_minima) {
    CHECK_ABSOLUTE_ERROR(Real(1), x, Real(0.01));
  }
}

#ifndef BOOST_MATH_BUILD_MODULE
// detail::symmetric_eigendecomposition is not exported from the module.
template<typename Real>
void check_symmetric_eigendecomposition(std::vector<Real> const & C, size_t n) {
  using boost::math::optimization::detail::symmetric_eigendecomposition;
  std::vector<Real> B;
  std::vector<Real> D;
  CHECK_EQUAL(symmetric_eigendecomposition(C, n, B, D), true);
  // The residuals are computed in long double, so that they measure the error of the decomposition:
  long double C_norm_sq = 0;
  long double residual_sq = 0;
  long double orthogonality_sq = 0;
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < n; ++j) {
      C_norm_sq += static_cast<long double>(C[i*n + j])*C[i*n + j];
      long double CB = -static_cast<long double>(B[i*n + j])*D[j];
      long double BTB = (i == j) ? -1 : 0;
      for (size_t k = 0; k < n; ++k) {
        CB += static_cast<long double>(C[i*n + k])*B[k*n + j];
        BTB += static_cast<long double>(B[k*n + i])*B[k*n + j];
      }
      residual_sq += CB*CB;
      orthogonality_sq += BTB*BTB;
    }
  }
  const long double tol = (2*n + 8)*static_cast<long double>(std::numeric_limits<Real>::epsilon());
  CHECK_LE(std::sqrt(residual_sq/C_norm_sq), tol);
  CHECK_LE(std::sqrt(orthogonality_sq), tol);
}

template<typename Real>
void test_symmetric_eigendecomposition() {
  std::cout << "Testing the symmetric eigendecomposition used by CMA-ES . . .\n";
  std::mt19937_64 gen(12345);
  std::normal_distribution<Real> dis(0, 1);
  for (size_t n : {1, 2, 3, 10, 50}) {
    // G G^T + nI is symmetric positive definite:
    std::vector<Real> G(n*n);
    for (auto & g : G) {
      g = dis(gen);
    }
    std::vector<Real> C(n*n, Real(0));
    for (size_t i = 0; i < n; ++i) {
      for (size_t j = 0; j < n; ++j) {
        for (size_t k = 0; k < n; ++k) {
          C[i*n + j] += G[i*n + k]*G[j*n + k];
        }
      }
      C[i*n + i] += static_cast<Real>(n);
    }
    check_symmetric_eigendecomposition(C, n);
    // A diagonal matrix:
    std::vector<Real> diagonal(n*n, Real(0));
    for (size_t i = 0; i < n; ++i) {
      diagonal[i*n + i] = static_cast<Real>(i + 1);
    }
    check_symmetric_eigendecomposition(diagonal, n);
    // I + uu^T has the eigenvalue 1 repeated n - 1 times:
    std::vector<Real> u(n);
    for (auto & ui : u) {
      ui = dis(gen);
    }
    std::vector<Real> repeated(n*n);
    for (size_t i = 0; i < n; ++i) {
      for (size_t j = 0; j < n; ++j) {
        repeated[i*n + j] = (i == j ? Real(1) : Real(0)) + u[i]*u[j];
      }
    }
    check_symmetric_eigendecomposition(repeated, n);
  }
}
#endif

int main() {
#if (defined(__clang__) || defined(_MSC_VER))
  test_ackley<float>();
  test_ackley<double>();
  test_rosenbrock_saddle<double>();
  test_rastrigin<double>();
  test_three_hump_camel<float>();
  test_beale<double>();
#endif
  test_sphere();
  test_bounds<float>();
  test_bounds<double>();
#ifndef BOOST_MATH_BUILD_MODULE
  test_symmetric_eigendecomposition<float>();
  test_symmetric_eigendecomposition<double>();
#endif
  return boost::math::test::report_errors();
}
