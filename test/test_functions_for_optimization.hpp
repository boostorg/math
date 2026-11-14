/*
 * Copyright Nick Thompson, 2024
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef TEST_FUNCTIONS_FOR_OPTIMIZATION_HPP
#define TEST_FUNCTIONS_FOR_OPTIMIZATION_HPP

#ifndef BOOST_MATH_BUILD_MODULE
#include <boost/math/constants/constants.hpp>
#include <array>
#include <vector>
#endif

/* simple n-d quadratic function */
template<typename RealType>
RealType
quadratic(std::vector<RealType>& x)
{
  RealType res{ 0.0 };
  for (auto& item : x) {
    res += item * item;
  }
  return res;
}

template<typename RealType>
RealType
quadratic_high_cond_2D(std::vector<RealType>& x)
{
  return 1000 * x[0] * x[0] + x[1] * x[1];
}

// Taken from: https://en.wikipedia.org/wiki/Test_functions_for_optimization
template<typename Real>
Real
ackley(std::array<Real, 2> const& v)
{
  using boost::math::constants::e;
  using boost::math::constants::two_pi;
  using std::cos;
  using std::exp;
  using std::sqrt;
  Real x = v[0];
  Real y = v[1];
  Real arg1 = -sqrt((x * x + y * y) / 2) / 5;
  Real arg2 = cos(two_pi<Real>() * x) + cos(two_pi<Real>() * y);
  return -20 * exp(arg1) - exp(arg2 / 2) + 20 + e<Real>();
}

template<typename Real>
auto
rosenbrock_saddle(std::array<Real, 2> const& v) -> Real
{
  Real x{ v[0] };
  Real y{ v[1] };
  return static_cast<Real>(100 * (x * x - y) * (x * x - y) + (1 - x) * (1 - x));
}

template<class Real>
Real
rastrigin(std::vector<Real> const& v)
{
  using boost::math::constants::two_pi;
  using std::cos;
  auto A = static_cast<Real>(10);
  auto y = static_cast<Real>(10 * v.size());
  for (auto x : v) {
    y += x * x - A * cos(two_pi<Real>() * x);
  }
  return y;
}

// Useful for testing return-type != scalar argument type,
// and robustness to NaNs:
double
sphere(std::vector<float> const& v)
{
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
Real
three_hump_camel(std::array<Real, 2> const& v)
{
  Real x = v[0];
  Real y = v[1];
  auto xsq = x * x;
  return 2 * xsq - (1 + Real(1) / Real(20)) * xsq * xsq + xsq * xsq * xsq / 6 +
         x * y + y * y;
}

// Minima occurs at (3, 1/2) with value 0:
template<typename Real>
Real
beale(std::array<Real, 2> const& v)
{
  Real x = v[0];
  Real y = v[1];
  Real t1 = Real(3) / Real(2) - x + x * y;
  Real t2 = Real(9) / Real(4) - x + x * y * y;
  Real t3 = Real(21) / Real(8) - x + x * y * y * y;
  return t1 * t1 + t2 * t2 + t3 * t3;
}

// Multiobjective problems, every objective to be minimized.
// Zitzler, Deb and Thiele (2000), Comparison of multiobjective evolutionary algorithms: empirical results.
// ZDT1, 2 and 3 are posed on [0, 1]^30 and ZDT6 on [0, 1]^10.
// The Pareto front is at g = 1, which is x[i] = 0 for i > 0, where f2 = h(f1).
template<class Container>
auto
zdt(Container const& x, int problem) -> std::array<typename Container::value_type, 2>
{
  using Real = typename Container::value_type;
  using boost::math::constants::pi;
  using std::exp;
  using std::pow;
  using std::sin;
  using std::sqrt;
  Real s = 0;
  for (size_t i = 1; i < x.size(); ++i) {
    s += x[i];
  }
  auto n = static_cast<Real>(x.size() - 1);
  Real f1 = x[0];
  Real g = 1 + 9 * s / n;
  if (problem == 6) {
    f1 = 1 - exp(-4 * x[0]) * pow(sin(6 * pi<Real>() * x[0]), 6);
    g = 1 + 9 * pow(s / n, Real(1) / 4);
  }
  Real h = 1 - sqrt(f1 / g);
  if (problem == 2 || problem == 6) {
    h = 1 - (f1 / g) * (f1 / g);
  } else if (problem == 3) {
    h -= (f1 / g) * sin(10 * pi<Real>() * f1);
  }
  return { f1, g * h };
}

// Deb, Thiele, Laumanns and Zitzler (2005), Scalable test problems for evolutionary multiobjective optimization.
// Three objectives on [0, 1]^n; the Pareto front is the part of the unit sphere in the positive octant.
template<class Container>
auto
dtlz2(Container const& x) -> std::vector<typename Container::value_type>
{
  using Real = typename Container::value_type;
  using boost::math::constants::half_pi;
  using std::sin;
  Real g = 1;
  for (size_t i = 2; i < x.size(); ++i) {
    g += (x[i] - Real(1) / 2) * (x[i] - Real(1) / 2);
  }
  // cos(t) as sin(pi/2 - t), so that no objective is negative at x = 1, where cos(half_pi<float>()) < 0:
  Real c0 = sin(half_pi<Real>() * (1 - x[0]));
  Real c1 = sin(half_pi<Real>() * (1 - x[1]));
  return { g * c0 * c1, g * c0 * sin(half_pi<Real>() * x[1]), g * sin(half_pi<Real>() * x[0]) };
}

// Fonseca and Fleming (1995), on [-4, 4]^n. The Pareto set is x[i] = u / sqrt(n) for -1 <= u <= 1, where
// f1 = 1 - exp(-(u - 1)^2) and f2 = 1 - exp(-(u + 1)^2).
template<class Container>
auto
fonseca_fleming(Container const& x) -> std::array<typename Container::value_type, 2>
{
  using Real = typename Container::value_type;
  using std::exp;
  using std::sqrt;
  Real s = 1 / sqrt(static_cast<Real>(x.size()));
  Real d1 = 0;
  Real d2 = 0;
  for (auto xi : x) {
    d1 += (xi - s) * (xi - s);
    d2 += (xi + s) * (xi + s);
  }
  return { 1 - exp(-d1), 1 - exp(-d2) };
}

// Kursawe (1991), on [-5, 5]^3. The Pareto front is in three pieces, and is not convex.
template<class Container>
auto
kursawe(Container const& x) -> std::array<typename Container::value_type, 2>
{
  using Real = typename Container::value_type;
  using std::abs;
  using std::exp;
  using std::pow;
  using std::sin;
  using std::sqrt;
  Real f1 = 0;
  Real f2 = 0;
  for (size_t i = 0; i < x.size(); ++i) {
    if (i + 1 < x.size()) {
      f1 -= 10 * exp(-sqrt(x[i] * x[i] + x[i + 1] * x[i + 1]) / 5);
    }
    f2 += pow(abs(x[i]), Real(4) / 5) + 5 * sin(x[i] * x[i] * x[i]);
  }
  return { f1, f2 };
}

#endif
