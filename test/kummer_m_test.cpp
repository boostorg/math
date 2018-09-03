// Copyright Nick Thompson 2017.
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#define BOOST_TEST_MODULE KummerMTest
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <iostream>
#include <iomanip>
#include <boost/math/special_functions/kummer_m.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>

using std::exp;
using boost::math::kummer_m;
using boost::multiprecision::cpp_bin_float_quad;

/*
 * Interesting identities to test:
 * M(1/2, 1, x) = e^(x/2)I_{0}(x/2), where I_{0} is the modified Bessel function
 * M(2,1,x) = e^x (1+x)
 * M(1,2,x) = (e^x-1)/x
 * M(a,a,x) = e^x
 * M(1, 1/2, x) = 1+e^x\sqrt(pi x)erf(\sqrt(x))
 * M(0, b, x) = 1
 * M(-1, b,x) = 1 - x/b
 * Kummer's transformation: M(a,b,x) = e^x M(b-a, b, -z)
 */

template<class Real>
void test_polynomial_reduction()
{
    //For a a negative integer, M(a,b,x) is a polynomial.
    // M(0, b, x) = 1:
    Real tol = 10*std::numeric_limits<Real>::epsilon();
    Real x = -0.5;
    Real b = 1;
    while (x < 0.5) {
      Real expected = 1;
      Real computed = kummer_m(0, b, x);
      BOOST_CHECK_CLOSE_FRACTION(expected, computed, tol);

      expected = 1 - x/b;
      computed = kummer_m(-1, b, x);
      BOOST_CHECK_CLOSE_FRACTION(expected, computed, tol);

      expected = 1 - 2*x/b + x*x/(b*(b+1));
      computed = kummer_m(-2, b, x);
      BOOST_CHECK_CLOSE_FRACTION(expected, computed, tol);
      x += 0.1;
    }
}

template<class Real>
void test_exp_identity()
{
    std::cout << std::setprecision(std::numeric_limits<Real>::digits10);

    Real tol = 10*std::numeric_limits<Real>::epsilon();
    Real x = -0.5;
    while (x < 1)
    {
        Real expected = exp(x);
        Real computed = kummer_m(1,1, x);
        BOOST_CHECK_CLOSE_FRACTION(expected, computed, tol);
        x += 0.1;
    }
}


BOOST_AUTO_TEST_CASE(KummerMTest)
{
    test_polynomial_reduction<float>();
    test_polynomial_reduction<double>();
    test_polynomial_reduction<long double>();

    test_exp_identity<float>();
    test_exp_identity<double>();
    test_exp_identity<long double>();
    test_exp_identity<cpp_bin_float_quad>();
}
