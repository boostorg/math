/*
 * Copyright Nick Thompson, 2017
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */
#define BOOST_TEST_MODULE chebyshev_test

#include <boost/type_index.hpp>
#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/math/special_functions/chebyshev.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>

using boost::multiprecision::cpp_bin_float_quad;
using boost::multiprecision::cpp_bin_float_50;
using boost::multiprecision::cpp_bin_float_100;
using boost::math::chebyshev_t;
using boost::math::chebyshev_t_prime;
using boost::math::chebyshev_u;

template<class Real>
void test_polynomials()
{
    std::cout << "Testing explicit polynomial representations of the Chebyshev polynomials on type " << boost::typeindex::type_id<Real>().pretty_name()  << "\n";

    Real x = -2;
    Real tol = 100*std::numeric_limits<Real>::epsilon();
    while (x < 2)
    {
        BOOST_CHECK_CLOSE_FRACTION(chebyshev_t(0, x), 1, tol);
        BOOST_CHECK_CLOSE_FRACTION(chebyshev_t(1, x), x, tol);
        BOOST_CHECK_CLOSE_FRACTION(chebyshev_t(2, x), 2*x*x - 1, tol);
        BOOST_CHECK_CLOSE_FRACTION(chebyshev_t(3, x), x*(4*x*x-3), tol);
        BOOST_CHECK_CLOSE_FRACTION(chebyshev_t(4, x), 8*x*x*x*x - 8*x*x + 1, tol);
        BOOST_CHECK_CLOSE_FRACTION(chebyshev_t(5, x), x*(16*x*x*x*x - 20*x*x + 5), tol);
        x += 1/static_cast<Real>(1<<7);
    }

    x = -2;
    tol = 10*tol;
    while (x < 2)
    {
        BOOST_CHECK_CLOSE_FRACTION(chebyshev_u(0, x), 1, tol);
        BOOST_CHECK_CLOSE_FRACTION(chebyshev_u(1, x), 2*x, tol);
        BOOST_CHECK_CLOSE_FRACTION(chebyshev_u(2, x), 4*x*x - 1, tol);
        BOOST_CHECK_CLOSE_FRACTION(chebyshev_u(3, x), 4*x*(2*x*x - 1), tol);
        x += 1/static_cast<Real>(1<<7);
    }
}


template<class Real>
void test_derivatives()
{
    std::cout << "Testing explicit polynomial representations of the Chebyshev polynomial derivatives on type " << boost::typeindex::type_id<Real>().pretty_name()  << "\n";

    Real x = -2;
    Real tol = 1000*std::numeric_limits<Real>::epsilon();
    while (x < 2)
    {
        BOOST_CHECK_CLOSE_FRACTION(chebyshev_t_prime(0, x), 0, tol);
        BOOST_CHECK_CLOSE_FRACTION(chebyshev_t_prime(1, x), 1, tol);
        BOOST_CHECK_CLOSE_FRACTION(chebyshev_t_prime(2, x), 4*x, tol);
        BOOST_CHECK_CLOSE_FRACTION(chebyshev_t_prime(3, x), 3*(4*x*x - 1), tol);
        BOOST_CHECK_CLOSE_FRACTION(chebyshev_t_prime(4, x), 16*x*(2*x*x - 1), tol);
        // This one makes the tolerance have to grow too large; the Chebyshev recurrence is more stable than naive polynomial evaluation anyway.
        //BOOST_CHECK_CLOSE_FRACTION(chebyshev_t_prime(5, x), 5*(4*x*x*(4*x*x - 3) + 1), tol);
        x += 1/static_cast<Real>(1<<7);
    }
}

BOOST_AUTO_TEST_CASE(chebyshev_test)
{
    test_polynomials<float>();
    test_polynomials<double>();
    test_polynomials<long double>();
    test_polynomials<cpp_bin_float_quad>();

    test_derivatives<float>();
    test_derivatives<double>();
    test_derivatives<long double>();
    test_derivatives<cpp_bin_float_quad>();
}
