//  (C) Copyright Nick Thompson 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#include "math_unit_test.hpp"
#include <boost/math/differentiation/autodiff_reverse.hpp>
#include <boost/math/optimization/minimizer.hpp>
#include <boost/math/optimization/trust_region_newton.hpp>
#include <stdexcept>
#include <vector>

namespace rdiff = boost::math::differentiation::reverse_mode;
namespace bopt = boost::math::optimization;
using boost::math::test::check_ulp_close;

template<class Real>
Real rosenbrock(std::vector<Real> const& x)
{
  Real f = Real(0);
  for (size_t i = 0; i + 1 < x.size(); ++i) {
    f += 100 * (x[i + 1] - x[i] * x[i]) * (x[i + 1] - x[i] * x[i]) + (1 - x[i]) * (1 - x[i]);
  }
  return f;
}

// Minimum 0 at (i + 1)/4.
template<class Real>
Real quadratic(std::vector<Real> const& x)
{
  Real f = Real(0);
  for (size_t i = 0; i < x.size(); ++i) {
    Real t = (i + 1) * (x[i] - Real(i + 1) / 4);
    f += t * t;
  }
  return f;
}

template<class Real>
Real himmelblau(std::vector<Real> const& v)
{
  Real a = v[0] * v[0] + v[1] - 11;
  Real b = v[0] + v[1] * v[1] - 7;
  return a * a + b * b;
}

template<class T>
using rvars = std::vector<rdiff::rvar<T, 2>>;

template<class T>
rvars<T> make_x(std::vector<T> const& x0)
{
  return rvars<T>(x0.begin(), x0.end());
}

template<class T>
void test_rosenbrock(size_t n)
{
  auto x = make_x<T>(std::vector<T>(n, T(-1.2)));
  auto opt = bopt::make_trust_region_newton(&rosenbrock<rdiff::rvar<T, 2>>, x);
  auto r = bopt::minimize(opt, bopt::unconstrained_policy<rvars<T>>{},
                          bopt::gradient_norm_convergence_policy<T>(T(16) * std::numeric_limits<T>::epsilon()),
                          bopt::max_iter_termination_policy(500));
  CHECK_EQUAL(r.converged, true);
  for (auto& xi : x) {
    CHECK_ABSOLUTE_ERROR(T(1), xi.item(), T(64) * std::sqrt(std::numeric_limits<T>::epsilon()));
  }
}

template<class T>
void test_quadratic()
{
  // The truncated CG step is an inexact Newton step, so a handful of steps
  // solve the problem; each is accepted and the radius never shrinks.
  const size_t n = 20;
  auto x = make_x<T>(std::vector<T>(n, T(0)));
  auto opt = bopt::make_trust_region_newton(&quadratic<rdiff::rvar<T, 2>>, x, T(100));
  auto r = bopt::minimize(opt, bopt::unconstrained_policy<rvars<T>>{},
                          bopt::gradient_norm_convergence_policy<T>(T(32) * std::numeric_limits<T>::epsilon()),
                          bopt::max_iter_termination_policy(100));
  CHECK_EQUAL(r.converged, true);
  CHECK_LE(r.num_iter, size_t(12));
  CHECK_EQUAL(opt.radius(), T(100));
  for (size_t i = 0; i < n; ++i) {
    // |x_i - c_i| <= |g_i| / 2
    CHECK_ABSOLUTE_ERROR(T(i + 1) / 4, x[i].item(), T(16) * std::numeric_limits<T>::epsilon());
  }
}

template<class T>
void test_himmelblau_escapes_maximum()
{
  // Near the local maximum the Hessian is negative definite; Steihaug-CG
  // follows the negative curvature to the boundary.
  auto x = make_x<T>({ T(-0.270845), T(-0.923039) });
  auto opt = bopt::make_trust_region_newton(&himmelblau<rdiff::rvar<T, 2>>, x, T(0.5));
  auto r = bopt::minimize(opt, bopt::unconstrained_policy<rvars<T>>{},
                          bopt::gradient_norm_convergence_policy<T>(T(100) * std::numeric_limits<T>::epsilon()),
                          bopt::max_iter_termination_policy(200));
  CHECK_EQUAL(r.converged, true);
  CHECK_LE(r.objective_value, T(64) * std::numeric_limits<T>::epsilon());
  // The objective decreases, except for rejected steps where it is unchanged.
  for (size_t i = 1; i < r.objective_history.size(); ++i) {
    CHECK_LE(r.objective_history[i], r.objective_history[i - 1]);
  }
}

template<class T>
void test_saddle()
{
  // f = x^2 - y^2 + y^4/4 has a saddle at the origin and minima at (0, +-sqrt(2)).
  auto f = [](std::vector<rdiff::rvar<T, 2>> const& v) {
    return v[0] * v[0] - v[1] * v[1] + v[1] * v[1] * v[1] * v[1] / 4;
  };
  auto x = make_x<T>({ T(1), T(1) / 1024 });
  auto opt = bopt::make_trust_region_newton(std::move(f), x);
  auto r = bopt::minimize(opt, bopt::unconstrained_policy<rvars<T>>{},
                          bopt::gradient_norm_convergence_policy<T>(T(16) * std::numeric_limits<T>::epsilon()),
                          bopt::max_iter_termination_policy(200));
  CHECK_EQUAL(r.converged, true);
  CHECK_ABSOLUTE_ERROR(T(0), x[0].item(), T(16) * std::numeric_limits<T>::epsilon());
  CHECK_ULP_CLOSE(std::sqrt(T(2)), x[1].item(), 16);
}

// The Hessian-vector product agrees with the analytic Hessian, at successive
// points and with seeds below epsilon.
template<class T>
void test_hessian_vector_product()
{
  auto x = make_x<T>({ T(0), T(0) });
  bopt::tape_initializer_rvar<T, 2>{}(x);
  bopt::reverse_mode_function_eval_policy<T, 2> fe;
  bopt::reverse_mode_gradient_evaluation_policy<T, 2> ge;
  bopt::reverse_mode_hessian_vector_product_policy<T> hvp;
  auto f = &rosenbrock<rdiff::rvar<T, 2>>;
  const T pts[3][2] = { { T(-1.2), T(1) }, { T(1), T(1) }, { T(2), T(-1) } };
  for (auto& pt : pts) {
    x[0].get_value().get_value() = pt[0];
    x[1].get_value().get_value() = pt[1];
    T v;
    std::vector<T> g, Hv;
    ge(f, x, fe, v, g);
    CHECK_ULP_CLOSE(rosenbrock(std::vector<T>{ pt[0], pt[1] }), v, 4);
    CHECK_ULP_CLOSE(-400 * pt[0] * (pt[1] - pt[0] * pt[0]) - 2 * (1 - pt[0]), g[0], 8);
    for (T scale : { T(1), std::numeric_limits<T>::epsilon() / 8 }) {
      hvp(f, x, std::vector<T>{ scale, 2 * scale }, Hv);
      const T H00 = 1200 * pt[0] * pt[0] - 400 * pt[1] + 2, H01 = -400 * pt[0];
      CHECK_ULP_CLOSE(scale * (H00 + 2 * H01), Hv[0], 8);
      CHECK_ULP_CLOSE(scale * (H01 + 400), Hv[1], 8);
    }
  }
}

// At x = 1 the adjoint of t is 2t = 0, but it still carries f'' = 12x^2 - 4.
template<class T>
void test_zero_adjoint_hessian()
{
  auto x = make_x<T>({ T(1) });
  bopt::tape_initializer_rvar<T, 2>{}(x);
  bopt::reverse_mode_function_eval_policy<T, 2> fe;
  bopt::reverse_mode_gradient_evaluation_policy<T, 2> ge;
  auto f = [](rvars<T> const& v) {
    rdiff::rvar<T, 2> t = v[0] * v[0] - 1;
    return rdiff::rvar<T, 2>(t * t + v[0]);
  };
  T value;
  std::vector<T> g, Hv;
  ge(f, x, fe, value, g);
  CHECK_EQUAL(g[0], T(1));
  bopt::reverse_mode_hessian_vector_product_policy<T>{}(f, x, std::vector<T>{ T(3) }, Hv);
  CHECK_EQUAL(Hv[0], T(24));
}

template<class T>
void test_box_constraints()
{
  // The minimum, at 2, is outside the box: minimize clips every step to 1.
  auto x = make_x<T>(std::vector<T>(3, T(0)));
  auto opt = bopt::make_trust_region_newton(
    [](rvars<T> const& v) {
      rdiff::rvar<T, 2> f = T(0);
      for (auto& vi : v) {
        f += (vi - 2) * (vi - 2);
      }
      return f;
    },
    x);
  bopt::minimize(opt, bopt::box_constraints<rvars<T>, T>(T(-1), T(1)),
                 bopt::gradient_norm_convergence_policy<T>(T(1e-6)), bopt::max_iter_termination_policy(10));
  for (auto& xi : x) {
    CHECK_EQUAL(xi.item(), T(1));
  }
  // Each step evaluates at the clipped point, with the variables' own nodes:
  // once the constraint is lifted, one Newton step solves the problem.
  CHECK_LE(T(2), opt.radius());
  opt.step();
  for (auto& xi : x) {
    CHECK_ULP_CLOSE(T(2), xi.item(), 4);
  }
}

void test_no_uphill_steps()
{
  // The model is linear at 0; a boundary step to -1 raises f by 16, which is measurable.
  auto x = make_x<double>({ 0.0 });
  auto opt = bopt::make_trust_region_newton(
    [](rvars<double> const& v) { return rdiff::rvar<double, 2>(1e16 + v[0] + 16 * v[0] * v[0] * v[0] * v[0]); }, x);
  for (int k = 0; k < 50; ++k) {
    opt.step();
    CHECK_LE(opt.objective_value(), 1e16);
  }
}

template<class T>
void test_invalid_radius()
{
  for (T r : { T(0), T(-1), std::numeric_limits<T>::infinity(), std::numeric_limits<T>::quiet_NaN() }) {
    auto x = make_x<T>({ T(0) });
    bool thrown = false;
    try {
      bopt::make_trust_region_newton(&quadratic<rdiff::rvar<T, 2>>, x, r);
    } catch (std::domain_error const&) {
      thrown = true;
    }
    CHECK_EQUAL(thrown, true);
  }
}

int main()
{
  test_zero_adjoint_hessian<float>();
  test_zero_adjoint_hessian<double>();
  test_box_constraints<float>();
  test_box_constraints<double>();
  test_no_uphill_steps();
  test_invalid_radius<double>();
  test_hessian_vector_product<float>();
  test_hessian_vector_product<double>();
  test_rosenbrock<float>(2);
  test_rosenbrock<double>(2);
  test_rosenbrock<double>(10);
  test_quadratic<float>();
  test_quadratic<double>();
  test_himmelblau_escapes_maximum<float>();
  test_himmelblau_escapes_maximum<double>();
  test_saddle<float>();
  test_saddle<double>();
  return boost::math::test::report_errors();
}
