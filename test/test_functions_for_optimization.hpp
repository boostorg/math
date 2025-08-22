/*
 * Copyright Nick Thompson, 2024
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef TEST_FUNCTIONS_FOR_OPTIMIZATION_HPP
#define TEST_FUNCTIONS_FOR_OPTIMIZATION_HPP

#include <boost/math/constants/constants.hpp>

#include <array>
#include <vector>

// Taken from: https://en.wikipedia.org/wiki/Test_functions_for_optimization
template <typename Real> Real ackley(std::array<Real, 2> const &v) {
  using std::sqrt;
  using std::cos;
  using std::exp;
  using boost::math::constants::two_pi;
  using boost::math::constants::e;
  Real x = v[0];
  Real y = v[1];
  Real arg1 = -sqrt((x * x + y * y) / 2) / 5;
  Real arg2 = cos(two_pi<Real>() * x) + cos(two_pi<Real>() * y);
  return -20 * exp(arg1) - exp(arg2 / 2) + 20 + e<Real>();
}

template <typename Real> auto rosenbrock_saddle(std::array<Real, 2> const &v) -> Real {
  Real x { v[0] };
  Real y { v[1] };
  return static_cast<Real>(100 * (x * x - y) * (x * x - y) + (1 - x) * (1 - x));
}


template <class Real> Real rastrigin(std::vector<Real> const &v) {
  using std::cos;
  using boost::math::constants::two_pi;
  auto A = static_cast<Real>(10);
  auto y = static_cast<Real>(10 * v.size());
  for (auto x : v) {
    y += x * x - A * cos(two_pi<Real>() * x);
  }
  return y;
}

// Useful for testing return-type != scalar argument type,
// and robustness to NaNs:
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

template<typename Real>
Real three_hump_camel(std::array<Real, 2> const & v) {
  Real x = v[0];
  Real y = v[1];
  auto xsq = x*x;
  return 2*xsq - (1 + Real(1)/Real(20))*xsq*xsq  + xsq*xsq*xsq/6 + x*y + y*y;
}

// Minima occurs at (3, 1/2) with value 0:
template<typename Real>
Real beale(std::array<Real, 2> const & v) {
  Real x = v[0];
  Real y = v[1];
  Real t1 = Real(3)/Real(2) -x + x*y;
  Real t2 = Real(9)/Real(4) -x  + x*y*y;
  Real t3 = Real(21)/Real(8) -x  + x*y*y*y;
  return t1*t1 + t2*t2 + t3*t3;
}


#endif
