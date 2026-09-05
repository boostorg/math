/*
 * Copyright Nick Thompson, 2020
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */

#include "math_unit_test.hpp"
#include <numeric>
#include <utility>
#include <iomanip>
#include <iostream>
#include <random>
#include <cmath>
#include <boost/assert.hpp>
#include <boost/core/demangle.hpp>
#include <boost/hana/for_each.hpp>
#include <boost/hana/ext/std/integer_sequence.hpp>
#include <boost/math/tools/condition_numbers.hpp>
#include <boost/math/differentiation/finite_difference.hpp>
#include <boost/math/special_functions/daubechies_wavelet.hpp>
#include <boost/math/special_functions/next.hpp>
#include <boost/math/quadrature/trapezoidal.hpp>

#ifdef BOOST_HAS_FLOAT128
#include <boost/multiprecision/float128.hpp>
using boost::multiprecision::float128;
#endif


using boost::math::constants::pi;
using boost::math::constants::root_two;

template<typename Real>
void test_exact_value()
{
    // The global phase of the wavelet is not constrained by anything other than convention.
    // Make sure that our conventions match the rest of the world:
    auto psi = boost::math::daubechies_wavelet<Real, 2>(2);
    Real computed = psi(1);
    Real expected = -1.366025403784439;
    CHECK_MOLLIFIED_CLOSE(expected, computed, 0.0001);
}

template<typename Real, int p>
void test_quadratures()
{
    std::cout << "Testing quadratures of " << p << " vanishing moment Daubechies wavelet on type " << boost::core::demangle(typeid(Real).name()) << "\n";
    using boost::math::quadrature::trapezoidal;
    auto psi = boost::math::daubechies_wavelet<Real, p>();
    std::cout << "Wavelet functor size is " << psi.bytes() << " bytes" << std::endl;
        
    Real tol = std::numeric_limits<Real>::epsilon();
    Real error_estimate = std::numeric_limits<Real>::quiet_NaN();
    Real L1 = std::numeric_limits<Real>::quiet_NaN();
    auto [a, b] = psi.support();
    CHECK_ULP_CLOSE(Real(-p+1), a, 0);
    CHECK_ULP_CLOSE(Real(p), b, 0);
    // A wavelet is a function of zero average; ensure the quadrature over its support is zero.
    Real Q = trapezoidal(psi, a, b, tol, 15, &error_estimate, &L1);
    if (!CHECK_MOLLIFIED_CLOSE(Real(0), Q, Real(0.0001)))
    {
        std::cerr << "  Quadrature of " << p << " vanishing moment wavelet does not vanish.\n";
        std::cerr << "  Error estimate: " << error_estimate << ", L1 norm: " << L1 << "\n";
    }
    auto psi_sq = [psi](Real x) {
        Real t = psi(x);
        return t*t;
    };
    Q = trapezoidal(psi_sq, a, b, tol, 15, &error_estimate, &L1);
    Real quad_tol = 2000*std::sqrt(std::numeric_limits<Real>::epsilon())/(p*p*p);
    if (!CHECK_MOLLIFIED_CLOSE(Real(1), Q, quad_tol))
    {
        std::cerr << "  L2 norm of " << p << " vanishing moment wavelet does not vanish.\n";
        std::cerr << "  Error estimate: " << error_estimate << ", L1 norm: " << L1 << "\n";
    }
    // psi is orthogonal to its integer translates: \int \psi(x-k) \psi(x) \, \mathrm{d}x = 0
    // g_n = 1/sqrt(2) <psi(t/2), phi(t-n)> (Mallat, 7.55)

    // Now hit the boundary. Much can go wrong here; this just tests for segfaults:
    int samples = 500;
    Real xlo = a;
    Real xhi = b;
    for (int i = 0; i < samples; ++i)
    {
        CHECK_ULP_CLOSE(Real(0), psi(xlo), 0);
        CHECK_ULP_CLOSE(Real(0), psi(xhi), 0);
        if constexpr (p > 2)
        {
            CHECK_ULP_CLOSE(Real(0), psi.prime(xlo), 0);
            CHECK_ULP_CLOSE(Real(0), psi.prime(xhi), 0);
            if constexpr (p >= 6) {
                CHECK_ULP_CLOSE(Real(0), psi.double_prime(xlo), 0);
                CHECK_ULP_CLOSE(Real(0), psi.double_prime(xhi), 0);
            }
        }
        xlo = std::nextafter(xlo, std::numeric_limits<Real>::lowest());
        xhi = std::nextafter(xhi, (std::numeric_limits<Real>::max)());
    }

    xlo = a;
    xhi = b;
    for (int i = 0; i < samples; ++i) {
        std::cout << std::setprecision(std::numeric_limits<Real>::max_digits10);
        BOOST_ASSERT(abs(psi(xlo)) <= 5);
        BOOST_ASSERT(abs(psi(xhi)) <= 5);
        if constexpr (p > 2)
        {
            BOOST_ASSERT(abs(psi.prime(xlo)) <= 5);
            BOOST_ASSERT(abs(psi.prime(xhi)) <= 5);
            if constexpr (p >= 6)
            {
                BOOST_ASSERT(abs(psi.double_prime(xlo)) <= 5);
                BOOST_ASSERT(abs(psi.double_prime(xhi)) <= 5);
            }
        }
        xlo = std::nextafter(xlo, (std::numeric_limits<Real>::max)());
        xhi = std::nextafter(xhi, std::numeric_limits<Real>::lowest());
    }
}

// Exercise both derivative-capable interpolators at points between grid nodes.
// Differentiating prime() checks double_prime() independently of the grid generator.
// Use coarse grids so the finite-difference stencil fits inside an interpolation interval.
template<class Real, int p>
void test_refinement_derivatives(int refinements)
{
    using boost::math::differentiation::finite_difference_derivative;
    auto f = boost::math::daubechies_wavelet<Real, p>(refinements);
    auto value = [&](Real x) { return f(x); };
    auto prime = [&](Real x) { return f.prime(x); };
    const Real tolerance = std::is_same_v<Real, float> ? Real(0.005) : Real(0.000001);
    for (int i = 0; i < 8; ++i)
    {
        Real x = Real(0) + Real(i)/8 + Real(1)/32;
        auto d1 = finite_difference_derivative<decltype(value), Real, 2>(value, x);
        auto d2 = finite_difference_derivative<decltype(prime), Real, 2>(prime, x);
        if constexpr (std::is_same_v<Real, float>)
        {
            // The automatic float step is too large for the curvature of these
            // interpolants. This dyadic step stays well inside the grid cell.
            Real h = Real(1)/1024;
            d1 = (value(x + h) - value(x - h))/(2*h);
            d2 = (prime(x + h) - prime(x - h))/(2*h);
        }
        CHECK_MOLLIFIED_CLOSE(d1, f.prime(x), tolerance);
        CHECK_MOLLIFIED_CLOSE(d2, f.double_prime(x), tolerance);
    }
}

// At shared dyadic nodes, changing the refinement must preserve function values.
template<class Real, int p>
void test_explicit_refinement()
{
    auto coarse = boost::math::daubechies_wavelet<Real, p>(3);
    auto fine = boost::math::daubechies_wavelet<Real, p>(4);
    CHECK_EQUAL(true, fine.bytes() > coarse.bytes());
    auto [a, b] = coarse.support();
    for (Real x = a; x <= b; x += Real(1)/8)
    {
        CHECK_MOLLIFIED_CLOSE(coarse(x), fine(x), 32*std::numeric_limits<Real>::epsilon());
        if constexpr (p > 2)
        {
            CHECK_MOLLIFIED_CLOSE(coarse.prime(x), fine.prime(x), 32*std::numeric_limits<Real>::epsilon());
        }
        if constexpr (p >= 6)
        {
            CHECK_MOLLIFIED_CLOSE(coarse.double_prime(x), fine.double_prime(x), 32*std::numeric_limits<Real>::epsilon());
        }
    }
}

template<class Real>
void test_haar_wavelet()
{
    for (int refinements : {-2, -1, 3})
    {
        auto psi = boost::math::daubechies_wavelet<Real, 1>(refinements);
        CHECK_EQUAL(Real(0), psi.support().first);
        CHECK_EQUAL(Real(1), psi.support().second);
        for (Real x : {Real(-1), Real(0), Real(0.5), Real(1), Real(2)})
        {
            CHECK_EQUAL(Real(0), psi(x));
        }
        CHECK_EQUAL(Real(1), psi(Real(0.25)));
        CHECK_EQUAL(Real(-1), psi(Real(0.75)));
        CHECK_EQUAL(Real(1), psi(std::nextafter(Real(0.5), Real(0))));
        CHECK_EQUAL(Real(-1), psi(std::nextafter(Real(0.5), Real(1))));
    }
    CHECK_THROW((boost::math::daubechies_wavelet<Real, 1>(0)), std::domain_error);
}

int main()
{
    #ifndef __MINGW32__
    try
    {
      test_haar_wavelet<float>();
      test_haar_wavelet<double>();
      test_exact_value<double>();

      boost::hana::for_each(std::make_index_sequence<17>(), [&](auto i) {
         test_quadratures<float, i + 3>();
         test_quadratures<double, i + 3>();
         });
      test_refinement_derivatives<float, 6>(4);
      test_refinement_derivatives<float, 19>(-2);
      test_refinement_derivatives<double, 6>(4);
      test_refinement_derivatives<double, 19>(4);
      boost::hana::for_each(std::make_index_sequence<18>(), [&](auto i) {
         test_explicit_refinement<float, i + 2>();
         test_explicit_refinement<double, i + 2>();
      });
    }
    catch (std::bad_alloc)
    {
       // not much we can do about this, this test uses lots of memory!
    }
    #endif
    return boost::math::test::report_errors();
}
