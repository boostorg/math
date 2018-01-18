/*
 *  (C) Copyright Nick Thompson 2018.
 *  Use, modification and distribution are subject to the
 *  Boost Software License, Version 1.0. (See accompanying file
 *  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 *
 */

#define BOOST_TEST_MODULE van_der_corput_test
#include <iostream>
#include <iomanip>
#include <boost/test/included/unit_test.hpp>
#include <boost/math/tools/van_der_corput.hpp>
#include <boost/math/special_functions/prime.hpp>

using boost::math::van_der_corput;
using boost::math::modified_van_der_corput;
using boost::math::prime;
using std::numeric_limits;


template<class Real>
void test_base_2()
{
    Real tol = numeric_limits<Real>::epsilon();
    Real z = van_der_corput<Real, int>(1, 2);
    BOOST_CHECK_CLOSE_FRACTION(z, (Real) 1/ (Real) 2, tol);
    z = van_der_corput<Real, int>(2, 2);
    BOOST_CHECK_CLOSE_FRACTION(z, (Real) 1/ (Real) 4, tol);
    z = van_der_corput<Real, int>(3, 2);
    BOOST_CHECK_CLOSE_FRACTION(z, (Real) 3/ (Real) 4, tol);
    z = van_der_corput<Real, int>(4, 2);
    BOOST_CHECK_CLOSE_FRACTION(z, (Real) 1/ (Real) 8, tol);
    z = van_der_corput<Real, int>(5, 2);
    BOOST_CHECK_CLOSE_FRACTION(z, (Real) 5/ (Real) 8, tol);
}

template<class Real>
void test_van_der_corput()
{
    Real tol = 5*numeric_limits<Real>::epsilon();
    for (int b = 3; b < 10000; ++b)
    {
        Real z = van_der_corput<Real, int>(0, b);
        BOOST_CHECK_EQUAL(z, 0);
        z = van_der_corput<Real, int>(1, b);
        BOOST_CHECK_CLOSE_FRACTION(z, (Real) 1/ (Real) b, tol);
        z = van_der_corput<Real, int>(2, b);
        BOOST_CHECK_CLOSE_FRACTION(z, (Real) 2/ (Real) b, tol);
        z = van_der_corput<Real, int>(b, b);
        Real expected = (Real) 1/ (Real) b;
        expected /= (Real) b;
        BOOST_CHECK_CLOSE_FRACTION(z, expected, tol);
        z = van_der_corput<Real, int>(b*b, b);
        expected /= (Real) b;
        BOOST_CHECK_CLOSE_FRACTION(z, expected, tol);
    }
}

template<class Real>
void test_modified_van_der_corput()
{
    Real tol = 5*numeric_limits<Real>::epsilon();
    for (int b = 3; b < 10000; ++b)
    {
        for (int x = 0; x < 100; ++x)
        {
          Real z1 = van_der_corput<Real, int>(x, b);
          Real z2 = modified_van_der_corput<Real, int>(x, b, 1);
          BOOST_CHECK_EQUAL(z1, z2);
        }
    }
}



BOOST_AUTO_TEST_CASE(van_der_corput_test)
{
    test_base_2<double>();
    test_van_der_corput<double>();
    test_modified_van_der_corput<double>();
}
