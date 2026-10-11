//  Copyright John Maddock 2006.
//  Copyright Paul A. Bristow 2007.
//  Copyright Matt Borland 2026.

//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
//  Module build check for boost::math::tools::brent_find_minima. It exercises
//  both exported overloads of the Brent minimiser when the test is built
//  against the boost.math module. The parabola case is mathematically exact;
//  the gamma-function case mirrors the expected minimum abscissa used in
//  test_minima.cpp.

#ifndef BOOST_MATH_BUILD_MODULE
#include <boost/math/tools/minima.hpp>
#include <boost/math/special_functions/gamma.hpp>
#else
import boost.math;
#endif

#include "math_unit_test.hpp"

// f(v) = 3 (v - 3)^2 + 4 has a unique minimum at v = 3 with value 4.
// Brent's parabolic interpolation is exact for a quadratic, so both overloads
// recover the vertex to well within the requested tolerances.
template <class Real>
void test_parabola(Real loc_tol, Real val_tol)
{
    using boost::math::tools::brent_find_minima;

    auto f = [](Real v) -> Real
    {
        const Real a {v - Real(3)};
        return Real(3) * a * a + Real(4);
    };

    // Overload reporting only (location, value).
    auto m = brent_find_minima(f, Real(-10), Real(10), 50);
    CHECK_ABSOLUTE_ERROR(Real(3), m.first, loc_tol);
    CHECK_ABSOLUTE_ERROR(Real(4), m.second, val_tol);
    // The reported value can never dip below the true minimum of 4.
    CHECK_GE(m.second, Real(4));

    // Overload with an explicit iteration budget: max_iter is written back in
    // place as the number of iterations actually consumed.
    boost::math::uintmax_t max_iter {1000};
    auto m2 = brent_find_minima(f, Real(-10), Real(10), 50, max_iter);
    CHECK_ABSOLUTE_ERROR(Real(3), m2.first, loc_tol);
    CHECK_ABSOLUTE_ERROR(Real(4), m2.second, val_tol);
    // Convergence is rapid, and the budget must have been decremented in place.
    CHECK_LE(max_iter, static_cast<boost::math::uintmax_t>(100));
}

// tgamma and lgamma both attain their minimum on the positive real axis at the
// same abscissa x0 = 1.4616321449683623..., where tgamma(x0) = 0.8856031944....
template <class Real>
void test_gamma_minimum(Real loc_tol, Real val_tol)
{
    using boost::math::tools::brent_find_minima;

    const Real x0 {static_cast<Real>(1.4616321449683623)};

    auto lg = [](Real x) -> Real { return boost::math::lgamma(x); };
    auto ml = brent_find_minima(lg, Real(0.5), Real(10), 50);
    CHECK_ABSOLUTE_ERROR(x0, ml.first, loc_tol);

    auto tg = [](Real x) -> Real { return boost::math::tgamma(x); };
    auto mt = brent_find_minima(tg, Real(0.5), Real(10), 50);
    CHECK_ABSOLUTE_ERROR(x0, mt.first, loc_tol);
    CHECK_ABSOLUTE_ERROR(static_cast<Real>(0.8856031944108887), mt.second, val_tol);
}

int main()
{
    test_parabola<double>(1e-5, 1e-6);
    test_gamma_minimum<double>(1e-6, 1e-6);

    test_parabola<float>(1e-3f, 1e-4f);
    test_gamma_minimum<float>(1e-3f, 1e-4f);

    return boost::math::test::report_errors();
}
