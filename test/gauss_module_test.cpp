// Copyright Nick Thompson, 2017
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Module build check for the Gauss-Legendre quadrature component. It consumes
// <boost/math/quadrature/gauss.hpp> through import boost.math and exercises the
// public gauss<Real, Points>::integrate entry point. Expected values are the
// mathematically-exact integrals and known-good results taken from
// gauss_quadrature_test.cpp (double/float paths only).

#ifndef BOOST_MATH_BUILD_MODULE
#include <boost/math/quadrature/gauss.hpp>
#else
import boost.math;
#endif

#ifndef BOOST_MATH_BUILD_MODULE
#include <cmath>
#endif
#include "math_unit_test.hpp"

using boost::math::quadrature::gauss;

template <class Real>
void test_spots(const char* type_name)
{
    std::cout << "Testing gauss quadrature for type " << type_name << std::endl;

    using std::exp;
    using std::sqrt;
    using std::cos;

    const Real pi {boost::math::constants::pi<Real>()};

    // A Gauss rule integrates polynomials of degree <= 2N-1 exactly, so a
    // 7-point rule reproduces linear and quadratic integrands to rounding.

    // Integral of 5x + 7 over [0, 1] is 9.5; the L1 norm matches on a positive
    // integrand.
    auto linear = [](const Real& x) { return 5 * x + 7; };
    Real L1 {};
    Real Q {gauss<Real, 7>::integrate(linear, static_cast<Real>(0), static_cast<Real>(1), &L1)};
    CHECK_ULP_CLOSE(static_cast<Real>(9.5), Q, 32);
    CHECK_ULP_CLOSE(static_cast<Real>(9.5), L1, 32);

    // Swapping the limits negates the result.
    Q = gauss<Real, 7>::integrate(linear, static_cast<Real>(1), static_cast<Real>(0));
    CHECK_ULP_CLOSE(static_cast<Real>(-9.5), Q, 32);

    // Integral of 5x^2 + 7x + 12 over [0, 1] is 103/6.
    auto quadratic = [](const Real& x) { return 5 * x * x + 7 * x + 12; };
    Q = gauss<Real, 7>::integrate(quadratic, static_cast<Real>(0), static_cast<Real>(1), &L1);
    CHECK_ULP_CLOSE(static_cast<Real>(103) / 6, Q, 32);
    CHECK_ULP_CLOSE(static_cast<Real>(103) / 6, L1, 32);

    // Smooth analytic integrand: exp(t) cos(t) over [0, pi/2] is (e^(pi/2) - 1)/2.
    auto smooth = [](const Real& t) { return exp(t) * cos(t); };
    Q = gauss<Real, 7>::integrate(smooth, static_cast<Real>(0), pi / 2);
    const Real smooth_expected {(exp(pi / 2) - 1) / 2};
    CHECK_ABSOLUTE_ERROR(smooth_expected, Q, static_cast<Real>(1e-5));

    // Area of the quarter unit circle: integral of sqrt(1 - t^2) over [0, 1] is pi/4.
    auto circle = [](const Real& t) { return sqrt(1 - t * t); };
    Q = gauss<Real, 7>::integrate(circle, static_cast<Real>(0), static_cast<Real>(1));
    CHECK_ABSOLUTE_ERROR(pi / 4, Q, static_cast<Real>(3.5e-3));
}

int main()
{
    test_spots<double>("double");
    test_spots<float>("float");

    return boost::math::test::report_errors();
}
