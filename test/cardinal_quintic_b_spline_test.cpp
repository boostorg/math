/*
 * Copyright Nick Thompson, 2019
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */

#include "math_unit_test.hpp"
#include <numeric>
#include <utility>
#include <boost/math/interpolators/cardinal_quintic_b_spline.hpp>
using boost::math::interpolators::cardinal_quintic_b_spline;

template<class Real>
void test_constant()
{
    Real c = 7.2;
    Real t0 = 0;
    Real h = Real(1)/Real(16);
    size_t n = 513;
    std::vector<Real> v(n, c);
    std::pair<Real, Real> left_endpoint_derivatives{0, 0};
    std::pair<Real, Real> right_endpoint_derivatives{0, 0};
    auto qbs = cardinal_quintic_b_spline<Real>(v.data(), v.size(), t0, h, left_endpoint_derivatives, right_endpoint_derivatives);

    size_t i = 0;
    while (i < n) {
      Real t = t0 + i*h;
      CHECK_ULP_CLOSE(c, qbs(t), 3);
      //CHECK_MOLLIFIED_CLOSE(0, qbs.prime(t), 100*std::numeric_limits<Real>::epsilon());
      ++i;
    }

    i = 0;
    while (i < n - 1) {
      Real t = t0 + i*h + h/2;
      CHECK_ULP_CLOSE(c, qbs(t), 4);
      //CHECK_MOLLIFIED_CLOSE(0, qbs.prime(t), 300*std::numeric_limits<Real>::epsilon());
      t = t0 + i*h + h/4;
      CHECK_ULP_CLOSE(c, qbs(t), 4);
      //CHECK_MOLLIFIED_CLOSE(0, qbs.prime(t), 150*std::numeric_limits<Real>::epsilon());
      ++i;
    }
}


template<class Real>
void test_linear()
{
    Real m = 8.3;
    Real b = 7.2;
    Real t0 = 0;
    Real h = Real(1)/Real(16);
    size_t n = 512;
    std::vector<Real> y(n);
    for (size_t i = 0; i < n; ++i) {
      Real t = i*h;
      y[i] = m*t + b;
    }
    std::pair<Real, Real> left_endpoint_derivatives{m, 0};
    std::pair<Real, Real> right_endpoint_derivatives{m, 0};
    auto qbs = cardinal_quintic_b_spline<Real>(y.data(), y.size(), t0, h, left_endpoint_derivatives, right_endpoint_derivatives);

    size_t i = 0;
    while (i < n) {
      Real t = t0 + i*h;
      if (!CHECK_ULP_CLOSE(m*t+b, qbs(t), 3)) {
          std::cerr << "  Problem at t = " << t << "\n";
      }
      //CHECK_ULP_CLOSE(m, qbs.prime(t), 820);
      ++i;
    }

    i = 0;
    while (i < n - 1) {
      Real t = t0 + i*h + h/2;
      if(!CHECK_ULP_CLOSE(m*t+b, qbs(t), 4)) {
          std::cerr << "  Problem at t = " << t << "\n";
      }
      //CHECK_MOLLIFIED_CLOSE(m, qbs.prime(t), 1500*std::numeric_limits<Real>::epsilon());
      t = t0 + i*h + h/4;
      if(!CHECK_ULP_CLOSE(m*t+b, qbs(t), 4)) {
          std::cerr << "  Problem at t = " << t << "\n";
      }
      //CHECK_MOLLIFIED_CLOSE(m, qbs.prime(t), 1500*std::numeric_limits<Real>::epsilon());
      ++i;
    }
}


template<class Real>
void test_quadratic()
{
    Real a = Real(1)/Real(16);
    Real b = -3.5;
    Real c = -9;
    Real t0 = 0;
    Real h = Real(1)/Real(16);
    size_t n = 513;
    std::vector<Real> y(n);
    for (size_t i = 0; i < n; ++i) {
      Real t = i*h;
      y[i] = a*t*t + b*t + c;
    }
    Real t_max = t0 + (n-1)*h;
    std::pair<Real, Real> left_endpoint_derivatives{b, 2*a};
    std::pair<Real, Real> right_endpoint_derivatives{2*a*t_max + b, 2*a};

    auto qbs = cardinal_quintic_b_spline<Real>(y, t0, h, left_endpoint_derivatives, right_endpoint_derivatives);

    size_t i = 0;
    while (i < n) {
      Real t = t0 + i*h;
      CHECK_ULP_CLOSE(a*t*t + b*t + c, qbs(t), 3);
      ++i;
    }

    i = 0;
    while (i < n -1) {
      Real t = t0 + i*h + h/2;
      if(!CHECK_ULP_CLOSE(a*t*t + b*t + c, qbs(t), 5)) {
          std::cerr << "  Problem at abscissa t = " << t << "\n";
      }

      t = t0 + i*h + h/4;
      if (!CHECK_ULP_CLOSE(a*t*t + b*t + c, qbs(t), 5)) {
          std::cerr << "  Problem abscissa t = " << t << "\n";
      }
      ++i;
    }
}

int main()
{
    test_constant<float>();
    test_constant<double>();
    test_constant<long double>();

    test_linear<float>();
    test_linear<double>();
    test_linear<long double>();

    test_quadratic<long double>();
    //test_quadratic<long double>();

    return boost::math::test::report_errors();
}
