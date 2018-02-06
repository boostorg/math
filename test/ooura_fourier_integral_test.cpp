// Copyright Nick Thompson, 2018
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)
#define BOOST_TEST_MODULE test_ooura_fourier_transform

#include <cmath>
#include <iostream>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/math_fwd.hpp>
#include <boost/type_index.hpp>
#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/math/quadrature/ooura_fourier_integrals.hpp>

using boost::multiprecision::cpp_bin_float_50;
using boost::multiprecision::cpp_bin_float_quad;
using boost::math::constants::third;
using boost::math::constants::half;
using boost::math::constants::pi;
using boost::math::quadrature::ooura_fourier_sin;
using boost::math::quadrature::ooura_fourier_cos;
using std::sqrt;
using std::log;

template<class Real>
void test_sinc()
{
    std::cout << "Testing sinc integral on type " << boost::typeindex::type_id<Real>().pretty_name()  << "\n";
    using std::numeric_limits;
    Real tol = 50*numeric_limits<Real>::epsilon();
    auto f = [](Real x)->Real { return 1/x; };
    Real omega = 1;
    while (omega < 5)
    {
        Real Is = ooura_fourier_sin<decltype(f), Real>(f, omega);
        BOOST_CHECK_CLOSE_FRACTION(Is, pi<Real>()/2, tol);

        Is = ooura_fourier_sin<decltype(f), Real>(f, -omega);
        BOOST_CHECK_CLOSE_FRACTION(Is, -pi<Real>()/2, tol);
        omega += 1;
    }
}

template<class Real>
void test_exp()
{
    std::cout << "Testing exponential integral on type " << boost::typeindex::type_id<Real>().pretty_name()  << "\n";
    using std::exp;
    using std::numeric_limits;
    Real tol = 50*numeric_limits<Real>::epsilon();
    auto f = [](Real x)->Real {return exp(-x);};
    Real omega = 1;
    while (omega < 5)
    {
        Real Is = ooura_fourier_sin<decltype(f), Real>(f, omega);
        Real exact = omega/(1+omega*omega);
        BOOST_CHECK_CLOSE_FRACTION(Is, exact, tol);
        omega += 1;
    }
}

template<class Real>
void test_root()
{
    std::cout << "Testing integral of sin(kx)/sqrt(x) on type " << boost::typeindex::type_id<Real>().pretty_name()  << "\n";
    using std::sqrt;
    using std::numeric_limits;
    Real tol = 1000*numeric_limits<Real>::epsilon();
    auto f = [](Real x)->Real { return 1/sqrt(x);};
    Real omega = 1;
    while (omega < 5)
    {
        Real Is = ooura_fourier_sin<decltype(f), Real>(f, omega);
        Real exact = sqrt(pi<Real>()/(2*omega));
        BOOST_CHECK_CLOSE_FRACTION(Is, exact, tol);
        omega += 1;
    }
}


// This works, but doesn't recover the precision you want in a unit test:
// template<class Real>
// void test_log()
// {
//     std::cout << "Testing integral of log(x)sin(x) on type " << boost::typeindex::type_id<Real>().pretty_name()  << "\n";
//     using std::log;
//     using std::exp;
//     using std::numeric_limits;
//     using boost::math::constants::euler;
//     Real tol = 1000*numeric_limits<Real>::epsilon();
//     auto f = [](Real x)->Real { return exp(-100*numeric_limits<Real>::epsilon()*x)*log(x);};
//     Real omega = 1;
//     Real Is = ooura_fourier_sin<decltype(f), Real>(f, omega, sqrt(numeric_limits<Real>::epsilon())/100);
//     BOOST_CHECK_CLOSE_FRACTION(Is, -euler<Real>(), tol);
// }


template<class Real>
void test_cos_integral1()
{
    std::cout << "Testing integral of cos(x)/(x*x+1) on type " << boost::typeindex::type_id<Real>().pretty_name()  << "\n";
    using std::exp;
    using boost::math::constants::half_pi;
    using boost::math::constants::e;
    using std::numeric_limits;
    Real tol = 10*numeric_limits<Real>::epsilon();
    auto f = [](Real x)->Real { return 1/(x*x+1);};
    Real omega = 1;
    Real Is = ooura_fourier_cos<decltype(f), Real>(f, omega);
    Real exact = half_pi<Real>()/e<Real>();
    BOOST_CHECK_CLOSE_FRACTION(Is, exact, tol);
}


BOOST_AUTO_TEST_CASE(fourier_integral_test)
{
    test_sinc<float>();
    test_sinc<double>();
    test_sinc<long double>();
    test_sinc<cpp_bin_float_quad>();

    test_exp<float>();
    test_exp<double>();
    test_exp<long double>();

    test_root<float>();
    test_root<double>();
    test_root<long double>();

    test_cos_integral1<float>();
    test_cos_integral1<double>();
    test_cos_integral1<long double>();
}
