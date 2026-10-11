//  (C) Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Module build check for the boost::math::tools root finders. It exercises the
// exported entry points of <boost/math/tools/roots.hpp> (quadratic_roots, bisect,
// newton_raphson_iterate, halley_iterate, schroder_iterate) when the test is built
// against the boost.math module. The quadratic_roots expected values mirror the
// assertions in test_roots.cpp; the iterative solvers locate sqrt(2) as the
// mathematically certain root of x^2 - 2.

#ifndef BOOST_MATH_BUILD_MODULE
#include <boost/math/tools/roots.hpp>
#else
import boost.math;
#endif

#include "math_unit_test.hpp"

#ifndef BOOST_MATH_BUILD_MODULE
#include <cmath>
#include <limits>
#include <utility>
#include <tuple>
#endif

// quadratic_roots: expected values copied from test_solve_real_quadratic in test_roots.cpp.
template <class Real>
void test_quadratic_roots()
{
    using boost::math::tools::quadratic_roots;

    // quadratic_roots promotes float inputs to double, so cast the computed
    // roots back to Real before comparing against the Real expected values.

    // x^2 - 1 = 0 has roots -1 and 1.
    auto p = quadratic_roots<Real>(1, 0, -1);
    CHECK_ULP_CLOSE(static_cast<Real>(-1), static_cast<Real>(p.first), 2);
    CHECK_ULP_CLOSE(static_cast<Real>(1), static_cast<Real>(p.second), 2);

    // 7 x^2 = 0 has a double root at 0.
    p = quadratic_roots<Real>(7, 0, 0);
    CHECK_EQUAL(static_cast<Real>(0), static_cast<Real>(p.first));
    CHECK_EQUAL(static_cast<Real>(0), static_cast<Real>(p.second));

    // (x - 7)^2 = x^2 - 14 x + 49 has a double root at 7.
    p = quadratic_roots<Real>(1, -14, 49);
    CHECK_ULP_CLOSE(static_cast<Real>(7), static_cast<Real>(p.first), 4);
    CHECK_ULP_CLOSE(static_cast<Real>(7), static_cast<Real>(p.second), 4);
}

// Kahan's example demonstrates the necessity of an fma; double-only known-good
// value taken directly from test_roots.cpp.
void test_quadratic_roots_kahan()
{
    using boost::math::tools::quadratic_roots;

    auto p = quadratic_roots<double>(94906265.625, -189812534.0, 94906268.375);
    CHECK_ULP_CLOSE(1.0, p.first, 8);
    CHECK_ULP_CLOSE(1.000000028975958, p.second, 8);
}

// bisect / newton / halley / schroder all locate sqrt(2), the root of x^2 - 2.
template <class Real>
void test_iterative_solvers()
{
    using boost::math::tools::bisect;
    using boost::math::tools::newton_raphson_iterate;
    using boost::math::tools::halley_iterate;
    using boost::math::tools::schroder_iterate;

    const Real expected = std::sqrt(static_cast<Real>(2));
    const int digits = std::numeric_limits<Real>::digits;

    // bisect: the bracket [0, 2] contains a sign change at sqrt(2).
    auto f0 = [](Real x) { return x * x - static_cast<Real>(2); };
    boost::math::tools::eps_tolerance<Real> tol(digits);
    auto bracket = bisect(f0, static_cast<Real>(0), static_cast<Real>(2), tol);
    Real root = bracket.first + (bracket.second - bracket.first) / 2;
    CHECK_ULP_CLOSE(expected, root, 8);

    // Newton-Raphson: functor returns (f, f').
    auto f1 = [](Real x) { return std::make_pair(x * x - static_cast<Real>(2), static_cast<Real>(2) * x); };
    root = newton_raphson_iterate(f1, static_cast<Real>(1.5), static_cast<Real>(1), static_cast<Real>(2), digits);
    CHECK_ULP_CLOSE(expected, root, 4);

    // Halley and Schroder: functor returns (f, f', f'').
    auto f2 = [](Real x) { return std::make_tuple(x * x - static_cast<Real>(2), static_cast<Real>(2) * x, static_cast<Real>(2)); };
    root = halley_iterate(f2, static_cast<Real>(1.5), static_cast<Real>(1), static_cast<Real>(2), digits);
    CHECK_ULP_CLOSE(expected, root, 4);

    root = schroder_iterate(f2, static_cast<Real>(1.5), static_cast<Real>(1), static_cast<Real>(2), digits);
    CHECK_ULP_CLOSE(expected, root, 4);
}

int main()
{
    test_quadratic_roots<float>();
    test_quadratic_roots<double>();
    test_quadratic_roots_kahan();

    test_iterative_solvers<float>();
    test_iterative_solvers<double>();

    return boost::math::test::report_errors();
}
