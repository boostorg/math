/*
 * Copyright Nick Thompson, 2026
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */

#include "math_unit_test.hpp"
#include "test_functions_for_optimization.hpp"
#include <boost/math/optimization/detail/rdiff_optimization_policies.hpp>
#ifndef BOOST_MATH_STANDALONE
#include <boost/multiprecision/cpp_bin_float.hpp>
#endif
#include <array>
#include <cmath>
#include <limits>
#include <vector>

namespace rdiff = boost::math::differentiation::reverse_mode;
using boost::math::optimization::reverse_mode_jacobian_evaluation_policy;

// Three outputs of two variables, in a std::array.
struct tall {
  template <class C> std::array<typename C::value_type, 3> operator()(C const &x) const {
    using std::exp;
    using std::sin;
    return {x[0] * x[1], sin(x[0]) + x[1] * x[1], exp(x[0] - x[1])};
  }
};

// Two outputs of four variables, in a std::vector; the first is an argument, as it is in the ZDT problems.
struct wide {
  template <class C> std::vector<typename C::value_type> operator()(C const &x) const {
    using std::log;
    return {x[0], x[0] * x[1] * x[2] + log(x[3]) * x[0]};
  }
};

// The policy the caller writes when the Jacobian is known.
template <class Real> struct analytic_jacobian_eval_pol {
  template <class Objective, class ArgumentContainer>
  void operator()(Objective &&obj_f, ArgumentContainer &x, std::vector<Real> &obj_v, std::vector<Real> &J) {
    using std::cos;
    using std::exp;
    auto y = obj_f(x);
    obj_v.assign(y.begin(), y.end());
    const Real e = exp(x[0] - x[1]);
    J = {x[1], x[0], cos(x[0]), 2 * x[1], e, -e};
  }
};

template <class Real> void test_analytic() {
  using std::log;
  const Real tol = 8 * std::numeric_limits<Real>::epsilon();
  auto &tape = rdiff::get_active_tape<Real, 1>();
  std::vector<Real> y, J, y_expected, J_expected;
  {
    std::array<Real, 2> x0{Real(3) / 4, Real(-5) / 4};
    std::array<rdiff::rvar<Real, 1>, 2> x{x0[0], x0[1]};
    tape.add_checkpoint();
    reverse_mode_jacobian_evaluation_policy<Real>{}(tall{}, x, y, J);
    analytic_jacobian_eval_pol<Real>{}(tall{}, x0, y_expected, J_expected);
    CHECK_EQUAL(y.size(), size_t(3));
    CHECK_EQUAL(J.size(), size_t(6));
    for (size_t i = 0; i < J.size(); ++i) {
      CHECK_RELATIVE_ERROR(J_expected[i], J[i], tol);
    }
    for (size_t i = 0; i < y.size(); ++i) {
      CHECK_RELATIVE_ERROR(y_expected[i], y[i], tol);
    }
  }
  {
    std::vector<Real> x0{Real(3) / 2, Real(-1) / 2, Real(2), Real(5) / 4};
    std::vector<rdiff::rvar<Real, 1>> x(x0.begin(), x0.end());
    tape.add_checkpoint();
    reverse_mode_jacobian_evaluation_policy<Real>{}(wide{}, x, y, J);
    J_expected = {1, 0, 0, 0, x0[1] * x0[2] + log(x0[3]), x0[0] * x0[2], x0[0] * x0[1], x0[0] / x0[3]};
    CHECK_EQUAL(y.size(), size_t(2));
    CHECK_EQUAL(J.size(), size_t(8));
    for (size_t i = 0; i < J.size(); ++i) {
      CHECK_RELATIVE_ERROR(J_expected[i], J[i], tol);
    }
    CHECK_EQUAL(y[0], x0[0]);
  }
  tape.clear();
}

template <class F> void test_central_differences(F f, size_t n, size_t m) {
  auto &tape = rdiff::get_active_tape<double, 1>();
  std::vector<double> x0(n);
  for (size_t j = 0; j < n; ++j) {
    x0[j] = 0.15 + 0.7 * static_cast<double>(j) / static_cast<double>(n);
  }
  std::vector<rdiff::rvar<double, 1>> x(x0.begin(), x0.end());
  tape.add_checkpoint();
  std::vector<double> y, J;
  reverse_mode_jacobian_evaluation_policy<double>{}(f, x, y, J);
  CHECK_EQUAL(y.size(), m);
  CHECK_EQUAL(J.size(), m * n);
  const double h = 1e-6;
  for (size_t j = 0; j < n; ++j) {
    auto a = x0, b = x0;
    a[j] -= h;
    b[j] += h;
    auto fa = f(a), fb = f(b);
    for (size_t i = 0; i < m; ++i) {
      CHECK_ABSOLUTE_ERROR((fb[i] - fa[i]) / (2 * h), J[i * n + j], 1e-7);
    }
  }
  tape.clear();
}

// A recording made before the checkpoint is not disturbed, and the calls reuse the tape after it.
void test_tape() {
  auto &tape = rdiff::get_active_tape<double, 1>();
  rdiff::rvar<double, 1> a = 2.0;
  rdiff::rvar<double, 1> z = a * a * a;
  std::array<rdiff::rvar<double, 1>, 2> x{0.75, -1.25};
  tape.add_checkpoint();
  std::vector<double> y, J, first;
  reverse_mode_jacobian_evaluation_policy<double> jacobian;
  jacobian(tall{}, x, y, first);
  const auto size = tape.end() - tape.begin();
  for (int i = 0; i < 3; ++i) {
    jacobian(tall{}, x, y, J);
    CHECK_EQUAL(tape.end() - tape.begin(), size);
    CHECK_TRUE(J == first);
  }
  tape.zero_grad();
  z.backward();
  CHECK_EQUAL(a.adjoint(), 12.0);
  // And the other way around: the Jacobian after a sweep of the other recording.
  jacobian(tall{}, x, y, J);
  CHECK_TRUE(J == first);
  CHECK_EQUAL(a.adjoint(), 12.0);
  // A checkpoint which is popped takes with it what was recorded after it, and the one before it is the last again.
  tape.rewind_to_last_checkpoint();
  const auto before = tape.end() - tape.begin();
  tape.add_checkpoint();
  rdiff::rvar<double, 1> more = z * a;
  CHECK_EQUAL(more.item(), 16.0);
  tape.add_checkpoint();
  tape.pop_checkpoint();
  tape.pop_checkpoint();
  CHECK_EQUAL(tape.end() - tape.begin(), before);
  jacobian(tall{}, x, y, J);
  CHECK_TRUE(J == first);
  CHECK_EQUAL(tape.end() - tape.begin(), size);
  tape.clear();
}

// Derivatives which are not finite arrive as they are, and stay in their own row.
void test_non_finite() {
  auto &tape = rdiff::get_active_tape<double, 1>();
  std::array<rdiff::rvar<double, 1>, 2> x{0.0, 3.0};
  tape.add_checkpoint();
  std::vector<double> y, J;
  auto f = [](std::array<rdiff::rvar<double, 1>, 2> const &v) {
    return std::array<rdiff::rvar<double, 1>, 2>{sqrt(v[0]), v[0] * v[1]};
  };
  reverse_mode_jacobian_evaluation_policy<double>{}(f, x, y, J);
  CHECK_TRUE(J == std::vector<double>({std::numeric_limits<double>::infinity(), 0.0, 3.0, 0.0}));
  tape.clear();
  // ZDT1 at x[0] = 0, where f2 has the factor 1 - sqrt(x[0]/g).
  std::vector<rdiff::rvar<double, 1>> u{0.0, 0.25, 0.25};
  tape.add_checkpoint();
  reverse_mode_jacobian_evaluation_policy<double>{}([](auto const &v) { return zdt(v, 1); }, u, y, J);
  CHECK_TRUE(std::vector<double>(J.begin(), J.begin() + 4) ==
             std::vector<double>({1.0, 0.0, 0.0, -std::numeric_limits<double>::infinity()}));
  CHECK_NAN(J[4]);
  tape.clear();
}

int main() {
  test_analytic<float>();
  test_analytic<double>();
#ifndef BOOST_MATH_STANDALONE
  test_analytic<boost::multiprecision::cpp_bin_float_50>();
#endif
  test_central_differences([](auto const &x) { return zdt(x, 1); }, 6, 2);
  test_central_differences([](auto const &x) { return zdt(x, 3); }, 6, 2);
  test_central_differences([](auto const &x) { return dtlz2(x); }, 7, 3);
  test_central_differences(tall{}, 2, 3);
  test_central_differences(wide{}, 4, 2);
  test_tape();
  test_non_finite();
  return boost::math::test::report_errors();
}
