/*
 * Copyright Nick Thompson, 2017
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */
#include "math_unit_test.hpp"
#include <algorithm>
#include <cmath>
#include <random>
#include <vector>
#include <boost/type_index.hpp>
#include <boost/math/special_functions/chebyshev.hpp>
#include <boost/math/special_functions/chebyshev_transform.hpp>
#include <boost/math/special_functions/sinc.hpp>
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

using boost::math::chebyshev_t;
using boost::math::chebyshev_t_prime;
using boost::math::chebyshev_u;
using boost::math::chebyshev_transform;


template<class Real>
void test_sin_chebyshev_transform()
{
    using boost::math::chebyshev_transform;
    using boost::math::constants::half_pi;
    using std::sin;
    using std::cos;
    using std::abs;

    Real tol = std::numeric_limits<Real>::epsilon();
    auto f = [](Real x)->Real { return sin(x); };
    Real a = 0;
    Real b = 1;
    chebyshev_transform<Real> cheb(f, a, b, tol);
    // Higher precision needs a longer series, and the rounding error of the derivative grows like its length squared.
    Real slack = std::numeric_limits<Real>::digits > 64 ? 10 : 1;

    Real x = a;
    while (x < b)
    {
        Real s = sin(x);
        Real c = cos(x);
        CHECK_ABSOLUTE_ERROR(s, cheb(x), slack*tol);
        CHECK_ABSOLUTE_ERROR(c, cheb.prime(x), 200*slack*tol);
        x += static_cast<Real>(1)/static_cast<Real>(1 << 7);
    }

    Real Q = cheb.integrate();

    CHECK_ABSOLUTE_ERROR(1 - cos(static_cast<Real>(1)), Q, 100*tol);
}


template<class Real>
void test_sinc_chebyshev_transform()
{
    using std::cos;
    using std::sin;
    using std::abs;
    using boost::math::sinc_pi;
    using boost::math::chebyshev_transform;
    using boost::math::constants::half_pi;

    Real tol = 100*std::numeric_limits<Real>::epsilon();
    auto f = [](Real x) { return boost::math::sinc_pi(x); };
    Real a = 0;
    Real b = 1;
    chebyshev_transform<Real> cheb(f, a, b, tol/50);
    // See test_sin_chebyshev_transform.
    Real slack = std::numeric_limits<Real>::digits > 64 ? 10 : 1;

    Real x = a;
    while (x < b)
    {
        Real s = sinc_pi(x);
        Real ds = (cos(x)-sinc_pi(x))/x;
        if (x == 0) { ds = 0; }

        CHECK_ABSOLUTE_ERROR(s, cheb(x), tol);
        CHECK_ABSOLUTE_ERROR(ds, cheb.prime(x), 10*slack*tol);
        x += static_cast<Real>(1)/static_cast<Real>(1 << 7);
    }

    Real Q = cheb.integrate();
    //NIntegrate[Sinc[x], {x, 0, 1}, WorkingPrecision -> 200, AccuracyGoal -> 150, PrecisionGoal -> 150, MaxRecursion -> 150]
    Real Q_exp = BOOST_MATH_TEST_VALUE(Real, 0.94608307036718301494135331382317965781233795473811179047145477356668);
    CHECK_ABSOLUTE_ERROR(Q_exp, Q, tol);
}



//Examples taken from "Approximation Theory and Approximation Practice", by Trefethen
template<class Real>
void test_atap_examples()
{
    using std::sin;
    using std::exp;
    using std::sqrt;
    using boost::math::constants::half;
    using boost::math::sinc_pi;
    using boost::math::chebyshev_transform;
    using boost::math::constants::half_pi;

    Real tol = 10*std::numeric_limits<Real>::epsilon();
    auto f1 = [](Real x) { return ((0 < x) - (x < 0)) - x/2; };
    auto f2 = [](Real x) { Real t = sin(6*x); Real s = sin(x + exp(2*x));
                           Real u = (0 < s) - (s < 0);
                           return t + u; };

    //auto f3 = [](Real x) { return sin(6*x) + sin(60*exp(x)); };
    //auto f4 = [](Real x) { return 1/(1+1000*(x+half<Real>())*(x+half<Real>())) + 1/sqrt(1+1000*(x-Real(1)/Real(2))*(x-Real(1)/Real(2)));};
    Real a = -1;
    Real b = 1;
    chebyshev_transform<Real> cheb1(f1, a, b, tol);
    chebyshev_transform<Real> cheb2(f2, a, b, tol);
    //chebyshev_transform<Real> cheb3(f3, a, b, tol);

    Real x = a;
    while (x < b)
    {
        // f1 and f2 are not differentiable; standard convergence rate theorems don't apply.
        // Basically, the max refinements are always hit; so the error is not related to the precision of the type.
        Real acceptable_error = sqrt(tol);
        Real acceptable_error_2 = 9e-4;
        // The error is that of truncating the series of a discontinuous function (about 1/n at the
        // maximum length n = 2^24), so it does not shrink with the precision beyond that of long double.
        if (std::is_same<Real, long double>::value || std::numeric_limits<Real>::digits > 64)
        {
            acceptable_error = 1.6e-5;
        }
        if (std::is_same<Real, double>::value)
        {
            acceptable_error *= 500;
        }
        CHECK_ABSOLUTE_ERROR(f1(x), cheb1(x), acceptable_error);

        CHECK_ABSOLUTE_ERROR(f2(x), cheb2(x), acceptable_error_2);
        x += static_cast<Real>(1)/static_cast<Real>(1 << 7);
    }
}


//Validate that the Chebyshev polynomials are well approximated by the Chebyshev transform.
template<class Real>
void test_chebyshev_chebyshev_transform()
{
    Real tol = 500*std::numeric_limits<Real>::epsilon();
    // T_0 = 1:
    auto t0 = [](Real) { return 1; };
    chebyshev_transform<Real> cheb0(t0, -1, 1);
    CHECK_ABSOLUTE_ERROR(2, cheb0.coefficients()[0], tol);

    Real x = -1;
    while (x < 1)
    {
        CHECK_ABSOLUTE_ERROR(1, cheb0(x), tol);
        CHECK_ABSOLUTE_ERROR(Real(0), cheb0.prime(x), tol);
        x += static_cast<Real>(1)/static_cast<Real>(1 << 7);
    }

    // T_1 = x:
    auto t1 = [](Real x) { return x; };
    chebyshev_transform<Real> cheb1(t1, -1, 1);
    CHECK_ABSOLUTE_ERROR(Real(1), cheb1.coefficients()[1], tol);

    x = -1;
    while (x < 1)
    {
        CHECK_ABSOLUTE_ERROR(x, cheb1(x), tol);
        CHECK_ABSOLUTE_ERROR(Real(1), cheb1.prime(x), tol);
        x += static_cast<Real>(1)/static_cast<Real>(1 << 7);
    }


    auto t2 = [](Real x) { return 2*x*x-1; };
    chebyshev_transform<Real> cheb2(t2, -1, 1);
    CHECK_ABSOLUTE_ERROR(Real(1), cheb2.coefficients()[2], tol);

    x = -1;
    while (x < 1)
    {
        CHECK_ABSOLUTE_ERROR(t2(x), cheb2(x), tol);
        CHECK_ABSOLUTE_ERROR(4*x, cheb2.prime(x), tol);
        x += static_cast<Real>(1)/static_cast<Real>(1 << 7);
    }
}

// The DCT-II behind chebyshev_transform, for odd and even n, against a naive sum in 50 digits;
// the in-place transform agrees with the out-of-place one to 0 ulps.
template <class Real>
void test_dct_ii()
{
   using boost::multiprecision::cpp_bin_float_50;
   std::mt19937_64 gen(5);
   std::uniform_real_distribution<double> dis(-1, 1);
   for (std::size_t n : {1, 2, 3, 4, 5, 7, 8, 12, 17, 64, 100, 127, 1000})
   {
      std::vector<Real> x(n);
      for (auto& t : x)
      {
         t = Real(dis(gen));
      }
      boost::math::detail::dct_ii<Real> dct(n);
      std::vector<Real> y(n);
      dct.execute(x.data(), y.data());
      std::vector<Real> z = x;
      dct.execute(z.data());
      CHECK_EQUAL(z == y, true);
      // y[k] = 2 sum_j x[j] cos(pi k (2j+1)/(2n)); the error bound is normwise, like the FFT's.
      Real scale = 0;
      std::vector<cpp_bin_float_50> ref(n);
      for (std::size_t k = 0; k < n; ++k)
      {
         cpp_bin_float_50 sum = 0;
         for (std::size_t j = 0; j < n; ++j)
         {
            sum += cpp_bin_float_50(x[j]) * boost::math::cos_pi(cpp_bin_float_50((2 * j + 1) * k) / (2 * n));
         }
         ref[k] = 2 * sum;
         scale = (std::max)(scale, static_cast<Real>(abs(ref[k])));
      }
      Real L = static_cast<Real>(std::log2(static_cast<double>(n)) + 1);
      for (std::size_t k = 0; k < n; ++k)
      {
         CHECK_ABSOLUTE_ERROR(static_cast<Real>(ref[k]), y[k], 4 * std::numeric_limits<Real>::epsilon() * L * scale);
      }
   }
}

int main()
{
#ifdef TEST1
    test_dct_ii<float>();
    test_chebyshev_chebyshev_transform<float>();
    test_sin_chebyshev_transform<float>();
    test_atap_examples<float>();
    test_sinc_chebyshev_transform<float>();
#endif
#ifdef TEST2
    test_dct_ii<double>();
    test_chebyshev_chebyshev_transform<double>();
    test_sin_chebyshev_transform<double>();
    test_atap_examples<double>();
    test_sinc_chebyshev_transform<double>();
#endif
#ifdef TEST3
    test_chebyshev_chebyshev_transform<long double>();
    test_sin_chebyshev_transform<long double>();
    test_atap_examples<long double>();
    test_sinc_chebyshev_transform<long double>();
#endif
#ifdef TEST4
#ifdef BOOST_HAS_FLOAT128
    using boost::multiprecision::float128;
    test_chebyshev_chebyshev_transform<float128>();
    test_sin_chebyshev_transform<float128>();
    test_atap_examples<float128>();
    test_sinc_chebyshev_transform<float128>();
#endif
#endif
#ifdef TEST5
    // The atap examples are not smooth, so they run to the maximum length (2^23 points); too slow in software arithmetic.
    using boost::multiprecision::cpp_bin_float_50;
    test_chebyshev_chebyshev_transform<cpp_bin_float_50>();
    test_sin_chebyshev_transform<cpp_bin_float_50>();
    test_sinc_chebyshev_transform<cpp_bin_float_50>();
#endif

    return boost::math::test::report_errors();
}


