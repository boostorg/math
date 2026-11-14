/*
 * Copyright Nick Thompson, 2026
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */

#ifndef BOOST_MATH_BUILD_MODULE
#include <boost/math/optimization/detail/bounded_lbfgs.hpp>
#else
import boost.math;
#endif

#include "math_unit_test.hpp"
#ifndef BOOST_MATH_BUILD_MODULE
#include <cmath>
#include <limits>
#include <random>
#include <vector>
#endif
using boost::math::optimization::detail::bounded_lbfgs;
using boost::math::optimization::detail::bounded_lbfgs_parameters;

double const inf = std::numeric_limits<double>::infinity();

double rosenbrock(std::vector<double> const &x, std::vector<double> &g) {
  double const a = x[1] - x[0] * x[0];
  g[0] = -400 * a * x[0] - 2 * (1 - x[0]);
  g[1] = 200 * a;
  return 100 * a * a + (1 - x[0]) * (1 - x[0]);
}

// With either line search.
void test_rosenbrock(bool strong_wolfe) {
  bounded_lbfgs_parameters params;
  params.strong_wolfe = strong_wolfe;
  params.function_tolerance = 0;
  std::vector<double> x{-1.2, 1};
  auto result = bounded_lbfgs(rosenbrock, x, {-inf, -inf}, {inf, inf}, params);
  CHECK_ABSOLUTE_ERROR(1.0, x[0], 1e-4);
  CHECK_ABSOLUTE_ERROR(1.0, x[1], 1e-4);
  CHECK_LE(result.value, 1e-9);
  CHECK_LE(result.projected_gradient, params.gradient_tolerance);
  CHECK_LE(result.iterations, size_t(100));

  // With x_0 <= 1/2 the minimizer is (1/2, 1/4), where the value is 1/4 and the gradient (-1, 0):
  x = {-1.2, 1};
  result = bounded_lbfgs(rosenbrock, x, {-1.5, -1.5}, {0.5, 2}, params);
  CHECK_EQUAL(0.5, x[0]);
  CHECK_ABSOLUTE_ERROR(0.25, x[1], 1e-6);
  CHECK_ABSOLUTE_ERROR(0.25, result.value, 1e-10);
  CHECK_LE(result.projected_gradient, params.gradient_tolerance);

  // Both bounds active: the minimizer is the corner (3/2, 2), and the start is outside the box.
  x = {4, -3};
  result = bounded_lbfgs(rosenbrock, x, {1.5, -1}, {3, 2}, params);
  CHECK_EQUAL(1.5, x[0]);
  CHECK_EQUAL(2.0, x[1]);
  CHECK_EQUAL(0.0, result.projected_gradient);

  // Maximizing the negative is the same problem:
  params.maximize = true;
  x = {-1.2, 1};
  auto negative = [](std::vector<double> const &p, std::vector<double> &g) {
    double const v = -rosenbrock(p, g);
    g[0] = -g[0];
    g[1] = -g[1];
    return v;
  };
  result = bounded_lbfgs(negative, x, {-1.5, -1.5}, {0.5, 2}, params);
  CHECK_EQUAL(0.5, x[0]);
  CHECK_ABSOLUTE_ERROR(0.25, x[1], 1e-6);
  CHECK_ABSOLUTE_ERROR(-0.25, result.value, 1e-10);
}

// x^T A x/2 - b^T x on the unit cube, with A = G^T G/n + I/10 dense: half the bounds are active at the minimizer, the
// situation of an acquisition function, and those variables have to be found and left on their bounds while the
// others converge at the step length of the face.
void test_quadratic(size_t n, double spread, bool start_on_bound, bool strong_wolfe = false) {
  std::mt19937_64 gen(n);
  std::normal_distribution<double> normal;
  std::vector<double> G(n * n), A(n * n), b(n);
  for (auto &v : G) {
    v = normal(gen);
  }
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < n; ++j) {
      double s = i == j ? 0.1 : 0.0;
      for (size_t k = 0; k < n; ++k) {
        s += G[k * n + i] * G[k * n + j] / static_cast<double>(n);
      }
      A[i * n + j] = s;
    }
    b[i] = spread * normal(gen);
  }
  auto f = [&](std::vector<double> const &x, std::vector<double> &g) {
    double v = 0;
    for (size_t i = 0; i < n; ++i) {
      double s = 0;
      for (size_t j = 0; j < n; ++j) {
        s += A[i * n + j] * x[j];
      }
      g[i] = s - b[i];
      v += x[i] * (s / 2 - b[i]);
    }
    return v;
  };
  std::vector<double> x(n, start_on_bound ? 0.0 : 0.5), lower(n, 0.0), upper(n, 1.0), g(n);
  bounded_lbfgs_parameters params;
  params.strong_wolfe = strong_wolfe;
  params.gradient_tolerance = 1e-7;
  params.function_tolerance = 0;
  auto result = bounded_lbfgs(f, x, lower, upper, params);
  CHECK_LE(result.projected_gradient, params.gradient_tolerance);
  CHECK_LE(result.evaluations, size_t(150));
  // The Karush-Kuhn-Tucker conditions: the gradient vanishes in the free variables and points out of the box in those
  // on a bound.
  CHECK_ABSOLUTE_ERROR(f(x, g), result.value, 1e-14);
  size_t active = 0;
  for (size_t i = 0; i < n; ++i) {
    if (x[i] == 0) {
      CHECK_LE(0.0, g[i]);
      ++active;
    }
    else if (x[i] == 1) {
      CHECK_LE(g[i], 0.0);
      ++active;
    }
    else {
      CHECK_TRUE(0 < x[i] && x[i] < 1);
      CHECK_ABSOLUTE_ERROR(0.0, g[i], 1e-7);
    }
  }
  std::cout << "n = " << n << ", spread " << spread << ": " << active << " active bounds, " << result.iterations
            << " iterations, " << result.evaluations << " evaluations\n";
  CHECK_LE(n / 3, active);
  CHECK_LE(active, 2 * n / 3);
}

// -sum log(x_i) + x_i, which is NaN for negative x_i and infinite at 0: the minimizer is 1.
void test_nan() {
  size_t const n = 5;
  auto f = [](std::vector<double> const &x, std::vector<double> &g) {
    double v = 0;
    for (size_t i = 0; i < x.size(); ++i) {
      v += x[i] - std::log(x[i]);
      g[i] = 1 - 1 / x[i];
    }
    return v;
  };
  std::vector<double> x{0.001, 5, 0.01, 9, 2}, lower(n, -1.0), upper(n, 10.0);
  upper[4] = 0.5;
  auto result = bounded_lbfgs(f, x, lower, upper);
  for (size_t i = 0; i + 1 < n; ++i) {
    CHECK_ABSOLUTE_ERROR(1.0, x[i], 1e-4);
  }
  CHECK_EQUAL(0.5, x[4]);
  CHECK_ABSOLUTE_ERROR(4 + 0.5 - std::log(0.5), result.value, 1e-8);
  // Not finite at the start: nothing can be done, and x is only projected.
  x.assign(n, -2.0);
  result = bounded_lbfgs(f, x, lower, upper);
  CHECK_TRUE(!std::isfinite(result.value));
  CHECK_EQUAL(-1.0, x[0]);
  CHECK_EQUAL(size_t(1), result.evaluations);
}

// From a corner of the box, where the gradient points inwards in some variables and outwards in the others.
void test_start_on_bound() {
  auto f = [](std::vector<double> const &x, std::vector<double> &g) {
    double v = 0;
    for (size_t i = 0; i < x.size(); ++i) {
      double const c = i % 2 == 0 ? 0.3 + 0.01 * static_cast<double>(i) : -1.0;
      g[i] = 2 * (x[i] - c) * static_cast<double>(i + 1);
      v += (x[i] - c) * (x[i] - c) * static_cast<double>(i + 1);
    }
    return v;
  };
  size_t const n = 20;
  std::vector<double> x(n, 0.0), lower(n, 0.0), upper(n, 1.0);
  bounded_lbfgs_parameters params;
  params.function_tolerance = 0;
  auto result = bounded_lbfgs(f, x, lower, upper, params);
  for (size_t i = 0; i < n; ++i) {
    if (i % 2 == 0) {
      CHECK_ABSOLUTE_ERROR(0.3 + 0.01 * static_cast<double>(i), x[i], 1e-5);
    }
    else {
      CHECK_EQUAL(0.0, x[i]);
    }
  }
  CHECK_LE(result.projected_gradient, 1e-5);
  // Already optimal, at a vertex: one evaluation.
  std::vector<double> c{0.0, 1.0};
  auto linear = [](std::vector<double> const &p, std::vector<double> &g) {
    g[0] = 1;
    g[1] = -2;
    return p[0] - 2 * p[1];
  };
  result = bounded_lbfgs(linear, c, {0.0, 0.0}, {1.0, 1.0});
  CHECK_EQUAL(size_t(1), result.evaluations);
  CHECK_EQUAL(size_t(0), result.iterations);
  CHECK_EQUAL(-2.0, result.value);
  // And from the opposite vertex it is reached exactly:
  c = {1.0, 0.0};
  result = bounded_lbfgs(linear, c, {0.0, 0.0}, {1.0, 1.0});
  CHECK_EQUAL(0.0, c[0]);
  CHECK_EQUAL(1.0, c[1]);
  CHECK_LE(result.evaluations, size_t(4));
}

int main() {
  test_rosenbrock(false);
  test_rosenbrock(true);
  test_quadratic(50, 0.2, false, true);
  test_quadratic(50, 0.2, true, true);
  test_quadratic(50, 0.2, false);
  test_quadratic(50, 0.2, true);
  test_quadratic(100, 0.2, false);
  test_nan();
  test_start_on_bound();
  return boost::math::test::report_errors();
}
