/*
 * Copyright Nick Thompson, 2019
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */

#include "math_unit_test.hpp"
#include <vector>
#include <random>
#include <type_traits>
#include <boost/math/constants/constants.hpp>
#include <boost/math/interpolators/cardinal_trigonometric.hpp>
#include <boost/math/tools/test_value.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/cpp_complex.hpp>
#ifdef BOOST_HAS_FLOAT128
#include <boost/multiprecision/float128.hpp>
#include <boost/multiprecision/complex128.hpp>
#endif

#if !defined(TEST1) && !defined(TEST2) && !defined(TEST3) && !defined(TEST4) && !defined(TEST5)
#  define TEST1
#  define TEST2
#  define TEST3
#  define TEST4
#  define TEST5
#endif


using std::sin;
using boost::math::constants::two_pi;
using boost::math::interpolators::cardinal_trigonometric;
using boost::multiprecision::cpp_bin_float_50;

// <random> distributions only accept the built-in floating point types.
template<class Real>
using dist_type = typename std::conditional<std::is_floating_point<Real>::value, Real, double>::type;

template<class Real>
void test_constant()
{
    Real t0 = 0;
    Real h = 1;
    for(size_t n = 1; n < 20; ++n)
    {
      Real c = 8;
      std::vector<Real> v(n, c);
      auto ct = cardinal_trigonometric<decltype(v)>(v, t0, h);
      // The FFT of constant data leaves rounding error in the coefficients that should vanish,
      // and the derivatives multiply the k-th of them by k and k^2.
      CHECK_ULP_CLOSE(c, ct(0.3), 6);
      CHECK_ULP_CLOSE(c*h*n, ct.integrate(), 3);
      CHECK_ULP_CLOSE(c*c*h*n, ct.squared_l2(), 6);
      CHECK_MOLLIFIED_CLOSE(Real(0), ct.prime(0.8), 75*std::numeric_limits<Real>::epsilon());
      CHECK_MOLLIFIED_CLOSE(Real(0), ct.double_prime(0.8), 150*std::numeric_limits<Real>::epsilon());
    }
}

template<class Real>
void test_interpolation_condition()
{
  std::mt19937 gen(1234);
  std::uniform_real_distribution<dist_type<Real>> dis(1, 10);

  for(size_t n = 1; n < 20; ++n) {
    Real t0 = static_cast<Real>(dis(gen));
    Real h = static_cast<Real>(dis(gen));
    std::vector<Real> v(n);
    for (size_t i = 0; i < n; ++i) {
      v[i] = static_cast<Real>(dis(gen));
    }
    auto ct = cardinal_trigonometric<decltype(v)>(v, t0, h);
    for (size_t i = 0; i < n; ++i) {
      Real arg = t0 + i*h;
      Real expected = v[i];
      Real computed = ct(arg);
      // The abscissa t0 + i*h is itself rounded, and the phase error it causes grows with n;
      // 5n ulps was too tight for float with clang-cl at n = 11 (60 ulps).
      if(!CHECK_ULP_CLOSE(expected, computed, 8*n))
      {
        std::cerr << "  Samples: " << n << "\n";
      }
    }
  }

}


template<class Real>
void test_sampled_sine()
{
    using std::sin;
    using std::cos;
    for (unsigned n = 15; n < 50; ++n)
    {
      Real t0 = 0;
      Real T = 1;
      Real h = T/n;
      std::vector<Real> v(n);
      auto s = [&](Real t) { return sin(two_pi<Real>()*(t-t0)/T);};
      auto s_prime = [&](Real t) { return two_pi<Real>()*cos(two_pi<Real>()*(t-t0)/T)/T;};
      auto s_double_prime = [&](Real t) { return -two_pi<Real>()*two_pi<Real>()*sin(two_pi<Real>()*(t-t0)/T)/(T*T);};
      for(size_t j = 0; j < v.size(); ++j)
      {
          Real t = t0 + j*h;
          v[j] = s(t);
      }
      auto ct = cardinal_trigonometric<decltype(v)>(v, t0, h);
      CHECK_ULP_CLOSE(T, ct.period(), 3);
      std::mt19937 gen(1234);
      std::uniform_real_distribution<dist_type<Real>> dist(0, 500);

      unsigned j = 0;
      while (j++ < 50) {
        Real arg = static_cast<Real>(dist(gen));
        Real expected = s(arg);
        Real computed = ct(arg);
        CHECK_MOLLIFIED_CLOSE(expected, computed, std::numeric_limits<Real>::epsilon()*4000);

        expected = s_prime(arg);
        computed = ct.prime(arg);
        CHECK_MOLLIFIED_CLOSE(expected, computed, 18000*std::numeric_limits<Real>::epsilon());

        expected = s_double_prime(arg);
        computed = ct.double_prime(arg);
        CHECK_MOLLIFIED_CLOSE(expected, computed, 100000*std::numeric_limits<Real>::epsilon());

      }
      CHECK_MOLLIFIED_CLOSE(Real(0), ct.integrate(), std::numeric_limits<Real>::epsilon());
    }
}

template<class Real>
void test_bump()
{
  using std::exp;
  using std::abs;
  using std::sqrt;
  using std::pow;
  auto bump = [](Real x)->Real { if (abs(x) >= 1) { return Real(0); } return exp(-Real(1)/(Real(1)-x*x)); };
  auto bump_prime = [](Real x)->Real {
      if (abs(x) >= 1) { return Real(0); }

      return -2*x*exp(-Real(1)/(Real(1)-x*x))/pow(1-x*x,2);
  };

  auto bump_double_prime = [](Real x)->Real {
      if (abs(x) >= 1) { return Real(0); }

      return (6*pow(x,4)-2)*exp(-Real(1)/(Real(1)-x*x))/pow(1-x*x,4);
  };


  Real t0 = -1;
  size_t n = 4096;
  Real h = Real(2)/Real(n);

  std::vector<Real> v(n);
  for(size_t i = 0; i < n; ++i)
  {
      Real t = t0 + i*h;
      v[i] = bump(t);
  }

  auto ct = cardinal_trigonometric<decltype(v)>(v, t0, h);
  // With 4096 samples the interpolant differs from the bump by about 1e-37 (1e-33 and 1e-30 for the derivatives),
  // which only matters for types with more than 113 bits.
  using std::max;
  Real eps = std::numeric_limits<Real>::epsilon();
  Real tol0 = (max)(Real(2*eps), Real(1e-35));
  Real tol1 = (max)(Real(4000*eps), Real(1e-31));
  Real tol2 = (max)(Real(4000*4000*eps), Real(1e-28));
  std::mt19937 gen(323723);
  std::uniform_real_distribution<long double> dis(-0.9, 0.9);

  size_t i = 0;
  while (i++ < 1000)
  {
      Real t = static_cast<Real>(dis(gen));
      Real expected = bump(t);
      Real computed = ct(t);
      if(!CHECK_MOLLIFIED_CLOSE(expected, computed, tol0)) {
          std::cerr << "  Problem occurred at abscissa " << t << "\n";
      }

      expected = bump_prime(t);
      computed = ct.prime(t);
      if(!CHECK_MOLLIFIED_CLOSE(expected, computed, tol1)) {
          std::cerr << "  Problem occurred at abscissa " << t << "\n";
      }

      expected = bump_double_prime(t);
      computed = ct.double_prime(t);
      if(!CHECK_MOLLIFIED_CLOSE(expected, computed, tol2)) {
          std::cerr << "  Problem occurred at abscissa " << t << "\n";
      }


  }

  // Wolfram Alpha:
  // NIntegrate[Exp[-1/(1-x*x)],{x,-1,1}]
  CHECK_ULP_CLOSE(BOOST_MATH_TEST_VALUE(Real, 0.443993816168079437823048921170552663761201789045697497307484553947041), ct.integrate(), 3);

  // NIntegrate[Exp[-2/(1-x*x)],{x,-1,1}]
  CHECK_ULP_CLOSE(BOOST_MATH_TEST_VALUE(Real, 0.133086120844994271556947327955328571362579155162813005534500258889539), ct.squared_l2(), 1);


}


int main()
{

#ifdef TEST1
    test_constant<float>();
    test_sampled_sine<float>();
    test_bump<float>();
    test_interpolation_condition<float>();
#endif


#ifdef TEST2
    test_constant<double>();
    test_sampled_sine<double>();
    test_bump<double>();
    test_interpolation_condition<double>();
#endif

#ifdef TEST3
    test_constant<long double>();
    test_sampled_sine<long double>();
    test_bump<long double>();
    test_interpolation_condition<long double>();
#endif

#ifdef TEST4
#ifdef BOOST_HAS_FLOAT128
    using boost::multiprecision::float128;
    test_constant<float128>();
    test_sampled_sine<float128>();
    test_bump<float128>();
    test_interpolation_condition<float128>();
#endif
#endif

#ifdef TEST5
    test_constant<cpp_bin_float_50>();
    test_sampled_sine<cpp_bin_float_50>();
    test_bump<cpp_bin_float_50>();
    test_interpolation_condition<cpp_bin_float_50>();
#endif

    return boost::math::test::report_errors();
}
