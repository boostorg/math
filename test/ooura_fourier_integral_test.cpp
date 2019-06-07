// Copyright Nick Thompson, 2019
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)
#define BOOST_TEST_MODULE test_ooura_fourier_transform

#include <cmath>
#include <iostream>
#include <boost/type_index.hpp>
#include <boost/test/included/unit_test.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/math/quadrature/ooura_fourier_integrals.hpp>

using boost::math::quadrature::ooura_fourier_sin;
using boost::math::constants::pi;

template<class Real>
void test_ooura_eta()
{
    using boost::math::quadrature::detail::ooura_eta;
    std::cout << "Testing eta function on type " << boost::typeindex::type_id<Real>().pretty_name()  << "\n";
    {
        Real x = 0;
        Real alpha = 7;
        auto [eta, eta_prime] = ooura_eta(x, alpha);
        BOOST_CHECK_SMALL(eta, std::numeric_limits<Real>::min());
        BOOST_CHECK_CLOSE_FRACTION(eta_prime, 2 + alpha + Real(1)/Real(4), 10*std::numeric_limits<Real>::epsilon());
    }

    {
        Real alpha = 4;
        for (Real z = 0.125; z < 500; z += 0.125) {
            Real x = std::log(z);
            auto [eta, eta_prime] = ooura_eta(x, alpha);
            BOOST_CHECK_CLOSE_FRACTION(eta, 2*x + alpha*(1-1/z) + (z-1)/4, 10*std::numeric_limits<Real>::epsilon());
            BOOST_CHECK_CLOSE_FRACTION(eta_prime, 2 + alpha/z + z/4, 10*std::numeric_limits<Real>::epsilon());
        }
    }
}

template<class Real>
void test_ooura_sin_nodes_and_weights()
{
    using boost::math::quadrature::detail::ooura_sin_node_and_weight;
    using boost::math::quadrature::detail::ooura_eta;
    std::cout << "Testing nodes and weights on type " << boost::typeindex::type_id<Real>().pretty_name()  << "\n";
    {
        long n = 1;
        Real alpha = 1;
        Real h = 1;
        auto [node, weight] = ooura_sin_node_and_weight(n, h, alpha);
        Real expected_node = pi<Real>()/(1-exp(-ooura_eta(n*h, alpha).first));
        BOOST_CHECK_CLOSE_FRACTION(node,  expected_node,10*std::numeric_limits<Real>::epsilon());
    }
}

template<class Real>
void test_ooura_alpha() {
    std::cout << "Testing Ooura alpha on type " << boost::typeindex::type_id<Real>().pretty_name()  << "\n";
    using std::sqrt;
    using std::log1p;
    using boost::math::quadrature::detail::calculate_ooura_alpha;
    Real alpha = calculate_ooura_alpha(Real(1));
    Real expected = 1/sqrt(16 + 4*log1p(pi<Real>()));
    BOOST_CHECK_CLOSE_FRACTION(alpha, expected, 10*std::numeric_limits<Real>::epsilon());
}

template<class Real>
void test_sinc()
{
    std::cout << "Testing sinc integral on type " << boost::typeindex::type_id<Real>().pretty_name()  << "\n";
    using std::numeric_limits;
    Real tol = 50*numeric_limits<Real>::epsilon();
    ooura_fourier_sin<Real> integrator(tol);
    auto f = [](Real x)->Real { return 1/x; };
    Real omega = 1;
    while (omega < 10)
    {
        Real Is = integrator.integrate(f, omega);
        BOOST_CHECK_CLOSE_FRACTION(Is, pi<Real>()/2, tol);

        Is = integrator.integrate(f, -omega);
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
    ooura_fourier_sin<Real> integrator(tol);
    auto f = [](Real x)->Real {return exp(-x);};
    Real omega = 1;
    while (omega < 5)
    {
        Real Is = integrator.integrate(f, omega);
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
    Real tol = 10*numeric_limits<Real>::epsilon();
    ooura_fourier_sin<Real> integrator(tol);
    auto f = [](Real x)->Real { return 1/sqrt(x);};
    Real omega = 1;
    while (omega < 5)
    {

        Real Is = integrator.integrate(f, omega);
        Real exact = sqrt(pi<Real>()/(2*omega));
        BOOST_CHECK_CLOSE_FRACTION(Is, exact, tol);
        omega += 1;
    }
}

// See: https://scicomp.stackexchange.com/questions/32790/numerical-evaluation-of-highly-oscillatory-integral/32799#32799
template<class Real>
Real asymptotic(Real lambda) {
    using std::sin;
    using std::cos;
    using boost::math::constants::pi;
    Real I1 = cos(lambda - pi<Real>()/4)*sqrt(2*pi<Real>()/lambda);
    Real I2 = sin(lambda - pi<Real>()/4)*sqrt(2*pi<Real>()/(lambda*lambda*lambda))/8;
    return I1 + I2;
}

template<class Real>
void test_double_osc()
{
    std::cout << "Testing double oscillation on type " << boost::typeindex::type_id<Real>().pretty_name()  << "\n";
    using std::sqrt;
    using std::numeric_limits;
    Real tol = 10*numeric_limits<Real>::epsilon();
    ooura_fourier_sin<Real> integrator(tol);
    Real lambda = 7;
    auto f = [&lambda](Real x)->Real { return cos(lambda*cos(x))/x; };
    Real omega = 1;
    Real Is = 2*integrator.integrate(f, omega);
    Real exact = asymptotic(lambda);
    BOOST_CHECK_CLOSE_FRACTION(Is, exact, 0.01);
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

/*
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

*/

BOOST_AUTO_TEST_CASE(ooura_fourier_transform_test)
{
    test_ooura_eta<float>();
    test_ooura_eta<double>();
    test_ooura_eta<long double>();

    test_ooura_sin_nodes_and_weights<float>();
    test_ooura_sin_nodes_and_weights<double>();
    test_ooura_sin_nodes_and_weights<long double>();

    test_ooura_alpha<float>();
    test_ooura_alpha<double>();
    test_ooura_alpha<long double>();

    test_sinc<float>();
    test_sinc<double>();
    test_sinc<long double>();

    test_exp<float>();
    test_exp<double>();
    test_exp<long double>();

    test_root<float>();
    test_root<double>();
    test_root<long double>();

    test_double_osc<float>();
    test_double_osc<double>();
    test_double_osc<long double>();
}
