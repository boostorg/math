// Copyright Nick Thompson, 2019
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Module build check for the Ooura Fourier integrals. Exercises the public
// entry points ooura_fourier_sin<Real>::integrate and
// ooura_fourier_cos<Real>::integrate against known analytic values when the
// test is built against the boost.math C++20 module.

#ifndef BOOST_MATH_BUILD_MODULE
#include <boost/math/quadrature/ooura_fourier_integrals.hpp>
#include <boost/math/constants/constants.hpp>
#else
import boost.math;
#endif

#ifndef BOOST_MATH_BUILD_MODULE
#include <cmath>
#include <limits>
#include <utility>
#endif
#include "math_unit_test.hpp"

template <class Real>
void test_sin_integrals(const char* type_name)
{
    std::cout << "Testing Ooura sin integrals for type " << type_name << std::endl;

    using boost::math::quadrature::ooura_fourier_sin;

    const Real tol {10 * std::numeric_limits<Real>::epsilon()};
    ooura_fourier_sin<Real> integrator {tol};

    // Dirichlet integral: integral of sin(omega x)/x on (0, inf) is pi/2 for omega > 0.
    auto sinc = [](Real x) -> Real { return 1 / x; };
    {
        auto [Is, err] = integrator.integrate(sinc, Real(1));
        (void)err;
        CHECK_MOLLIFIED_CLOSE(boost::math::constants::pi<Real>() / 2, Is,
                              50 * std::numeric_limits<Real>::epsilon());
    }
    {
        // Odd in omega, so a negative frequency flips the sign.
        auto [Is, err] = integrator.integrate(sinc, Real(-1));
        (void)err;
        CHECK_MOLLIFIED_CLOSE(-boost::math::constants::pi<Real>() / 2, Is,
                              50 * std::numeric_limits<Real>::epsilon());
    }

    // integral of exp(-x) sin(omega x) on (0, inf) is omega/(1 + omega^2).
    auto expo = [](Real x) -> Real { return std::exp(-x); };
    for (Real omega = 1; omega < 3; ++omega)
    {
        auto [Is, err] = integrator.integrate(expo, omega);
        (void)err;
        const Real exact {omega / (1 + omega * omega)};
        CHECK_MOLLIFIED_CLOSE(exact, Is, 50 * std::numeric_limits<Real>::epsilon());
    }

    // integral of sin(omega x)/sqrt(x) on (0, inf) is sqrt(pi/(2 omega)).
    auto root = [](Real x) -> Real { return 1 / std::sqrt(x); };
    {
        const Real omega {1};
        auto [Is, err] = integrator.integrate(root, omega);
        (void)err;
        const Real exact {std::sqrt(boost::math::constants::pi<Real>() / (2 * omega))};
        CHECK_MOLLIFIED_CLOSE(exact, Is, 100 * std::numeric_limits<Real>::epsilon());
    }

    // A zero integrand must integrate to exactly zero.
    auto zero = [](Real) -> Real { return Real(0); };
    {
        auto [Is, err] = integrator.integrate(zero, Real(1));
        (void)err;
        CHECK_EQUAL(Is, Real(0));
    }
}

template <class Real>
void test_cos_integrals(const char* type_name)
{
    std::cout << "Testing Ooura cos integrals for type " << type_name << std::endl;

    using boost::math::quadrature::ooura_fourier_cos;

    const Real tol {10 * std::numeric_limits<Real>::epsilon()};
    ooura_fourier_cos<Real> integrator {tol};

    // integral of cos(x)/(x^2 + 1) on (0, inf) is (pi/2)/e (omega = 1).
    auto lorentz = [](Real x) -> Real { return 1 / (x * x + 1); };
    {
        auto [Is, err] = integrator.integrate(lorentz, Real(1));
        (void)err;
        const Real exact {boost::math::constants::half_pi<Real>() / boost::math::constants::e<Real>()};
        CHECK_MOLLIFIED_CLOSE(exact, Is, 10 * std::numeric_limits<Real>::epsilon());
    }

    // integral of exp(-a x) cos(omega x) on (0, inf) is a/(a^2 + omega^2).
    for (Real a = 1; a < 4; ++a)
    {
        auto expo = [&a](Real x) -> Real { return std::exp(-a * x); };
        for (Real omega = 1; omega < 3; ++omega)
        {
            auto [Is, err] = integrator.integrate(expo, omega);
            (void)err;
            const Real exact {a / (a * a + omega * omega)};
            CHECK_MOLLIFIED_CLOSE(exact, Is, 500 * std::numeric_limits<Real>::epsilon());
        }
    }
}

int main()
{
    test_sin_integrals<float>("float");
    test_sin_integrals<double>("double");

    test_cos_integrals<float>("float");
    test_cos_integrals<double>("double");

    return boost::math::test::report_errors();
}
