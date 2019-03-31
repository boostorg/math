/*
 * Copyright Nick Thompson, 2019
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */

#include "math_unit_test.hpp"
#include <numeric>
#include <utility>
#include <boost/core/demangle.hpp>
#include <boost/math/interpolators/cardinal_quadratic_b_spline.hpp>
#ifdef BOOST_HAS_FLOAT128
#include <boost/multiprecision/float128.hpp>
using boost::multiprecision::float128;
#endif

using boost::math::interpolators::cardinal_quadratic_b_spline;

template<class Real>
void test_constant()
{
    Real c = 7.2;
    Real t0 = 0;
    Real h = Real(1)/Real(16);
    size_t n = 512;
    std::vector<Real> v(n, c);
    auto qbs = cardinal_quadratic_b_spline(v.data(), v.size(), t0, h);
  
    size_t i = 0;
    while (i < n) {
      Real t = t0 + i*h;
      CHECK_ULP_CLOSE(c, qbs(t), 2);
      ++i;
    }

    i = 0;
    while (i < n) {
      Real t = t0 + i*h + h/2;
      CHECK_ULP_CLOSE(c, qbs(t), 2);

      t = t0 + i*h + h/4;
      CHECK_ULP_CLOSE(c, qbs(t), 2);
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
    auto qbs = cardinal_quadratic_b_spline(y.data(), y.size(), t0, h);
  
    size_t i = 0;
    while (i < n) {
      Real t = t0 + i*h;
      CHECK_ULP_CLOSE(m*t+b, qbs(t), 2);
      ++i;
    }

    i = 0;
    while (i < n) {
      Real t = t0 + i*h + h/2;
      CHECK_ULP_CLOSE(m*t+b, qbs(t), 2);

      t = t0 + i*h + h/4;
      CHECK_ULP_CLOSE(m*t+b, qbs(t), 2);
      ++i;
    }
}

template<class Real>
void test_quadratic()
{
    Real a = 8.2;
    Real b = 7.2;
    Real c = -9.2;
    Real t0 = 0;
    Real h = Real(1)/Real(16);
    size_t n = 513;
    std::vector<Real> y(n);
    for (size_t i = 0; i < n; ++i) {
      Real t = i*h;
      y[i] = a*t*t + b*t + c;
    }
    Real t_max = t0 + (n-1)*h;
    auto qbs = cardinal_quadratic_b_spline(y.data(), y.size(), t0, h, b, 2*a*t_max + b);
  
    size_t i = 0;
    while (i < n) {
      Real t = t0 + i*h;
      CHECK_ULP_CLOSE(a*t*t + b*t + c, qbs(t), 2);
      ++i;
    }

    i = 0;
    while (i < n) {
      Real t = t0 + i*h + h/2;
      CHECK_ULP_CLOSE(a*t*t + b*t + c, qbs(t), 8);

      t = t0 + i*h + h/4;
      if (!CHECK_ULP_CLOSE(a*t*t + b*t + c, qbs(t), 12)) {
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
#ifdef BOOST_HAS_FLOAT128
    test_constant<float128>(); 
#endif

    test_linear<float>();
    test_linear<double>();
    test_linear<long double>();
#ifdef BOOST_HAS_FLOAT128
    test_linear<float128>(); 
#endif

    test_quadratic<double>();
    test_quadratic<long double>();
#ifdef BOOST_HAS_FLOAT128
    test_quadratic<float128>(); 
#endif


    return boost::math::test::report_errors();
}
