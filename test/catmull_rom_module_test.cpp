/*
 * Copyright Nick Thompson, 2017
 * Copyright Matt Borland, 2026
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */
//
// Module build check for catmull_rom. Consuming the interpolator through
// `import boost.math;` confirms the exported class template and its
// max_parameter/parameter_at_point/operator()/prime members are usable from a
// module consumer, along with both the move-vector and initializer_list
// constructors. The assertions exercise only mathematically-certain properties:
// a Catmull-Rom spline passes through every control point (the interpolation
// condition), and on collinear, equally spaced points it reproduces the
// straight line through them exactly, so both the value and the tangent are
// exact. Expected values are taken from the test_linear case of
// catmull_rom_test.cpp.

#ifndef BOOST_MATH_BUILD_MODULE
#include <boost/math/interpolators/catmull_rom.hpp>
#else
import boost.math;
#endif

#include "math_unit_test.hpp"

#ifndef BOOST_MATH_BUILD_MODULE
#include <array>
#include <vector>
#include <utility>
#include <limits>
#endif

using boost::math::catmull_rom;

// Four equally spaced, collinear points on the x-axis. The spline passes through
// every node, is the straight line s -> (s, 0, 0) on the interior segments, and
// has the constant unit tangent (1, 0, 0) there.
template <class Real>
void test_linear(const char* type_name)
{
    std::cout << "Testing catmull_rom linear reproduction for type " << type_name << std::endl;

    const Real eps {std::numeric_limits<Real>::epsilon()};

    std::vector<std::array<Real, 3>> v(4);
    v[0] = {Real(0), Real(0), Real(0)};
    v[1] = {Real(1), Real(0), Real(0)};
    v[2] = {Real(2), Real(0), Real(0)};
    v[3] = {Real(3), Real(0), Real(0)};

    catmull_rom<std::array<Real, 3>> cr(std::move(v));

    // One unit of parameter between consecutive nodes, total length 3.
    CHECK_ULP_CLOSE(Real(3), cr.max_parameter(), 4);
    CHECK_ABSOLUTE_ERROR(Real(0), cr.parameter_at_point(0), eps);
    CHECK_ULP_CLOSE(Real(1), cr.parameter_at_point(1), 4);
    CHECK_ULP_CLOSE(Real(2), cr.parameter_at_point(2), 4);
    CHECK_ULP_CLOSE(Real(3), cr.parameter_at_point(3), 4);

    // Interpolation condition: the curve passes through each control point.
    const auto p0 = cr(Real(0));
    CHECK_ABSOLUTE_ERROR(Real(0), p0[0], eps);
    CHECK_ABSOLUTE_ERROR(Real(0), p0[1], eps);
    CHECK_ABSOLUTE_ERROR(Real(0), p0[2], eps);

    const auto p1 = cr(Real(1));
    CHECK_ULP_CLOSE(Real(1), p1[0], 4);
    CHECK_ABSOLUTE_ERROR(Real(0), p1[1], eps);
    CHECK_ABSOLUTE_ERROR(Real(0), p1[2], eps);

    const auto p2 = cr(Real(2));
    CHECK_ULP_CLOSE(Real(2), p2[0], 4);
    CHECK_ABSOLUTE_ERROR(Real(0), p2[1], eps);
    CHECK_ABSOLUTE_ERROR(Real(0), p2[2], eps);

    const auto p3 = cr(Real(3));
    CHECK_ULP_CLOSE(Real(3), p3[0], 4);
    CHECK_ABSOLUTE_ERROR(Real(0), p3[1], eps);
    CHECK_ABSOLUTE_ERROR(Real(0), p3[2], eps);

    // On the interior segment [1, 2] the curve is the line s -> (s, 0, 0) and
    // its tangent is the constant (1, 0, 0).
    const Real interior[3] {Real(1.25), Real(1.5), Real(1.75)};
    for (const Real s : interior)
    {
        const auto p = cr(s);
        CHECK_ULP_CLOSE(s, p[0], 4);
        CHECK_ABSOLUTE_ERROR(Real(0), p[1], eps);
        CHECK_ABSOLUTE_ERROR(Real(0), p[2], eps);

        const auto tangent = cr.prime(s);
        CHECK_ULP_CLOSE(Real(1), tangent[0], 4);
        CHECK_ABSOLUTE_ERROR(Real(0), tangent[1], eps);
        CHECK_ABSOLUTE_ERROR(Real(0), tangent[2], eps);
    }
}

// The initializer_list constructor produces the same interpolant, and the curve
// reproduces its node values at the parameter reported by parameter_at_point.
template <class Real>
void test_initializer_list(const char* type_name)
{
    std::cout << "Testing catmull_rom initializer_list construction for type " << type_name << std::endl;

    const Real eps {std::numeric_limits<Real>::epsilon()};

    const std::array<Real, 3> a0 {Real(0), Real(0), Real(0)};
    const std::array<Real, 3> a1 {Real(1), Real(0), Real(0)};
    const std::array<Real, 3> a2 {Real(2), Real(0), Real(0)};
    const std::array<Real, 3> a3 {Real(3), Real(0), Real(0)};

    catmull_rom<std::array<Real, 3>> cr({a0, a1, a2, a3});

    const auto q0 = cr(cr.parameter_at_point(0));
    CHECK_ABSOLUTE_ERROR(Real(0), q0[0], eps);
    CHECK_ABSOLUTE_ERROR(Real(0), q0[1], eps);

    const auto q1 = cr(cr.parameter_at_point(1));
    CHECK_ULP_CLOSE(Real(1), q1[0], 4);
    CHECK_ABSOLUTE_ERROR(Real(0), q1[1], eps);

    const auto q2 = cr(cr.parameter_at_point(2));
    CHECK_ULP_CLOSE(Real(2), q2[0], 4);
    CHECK_ABSOLUTE_ERROR(Real(0), q2[1], eps);
}

int main()
{
    test_linear<float>("float");
    test_linear<double>("double");

    test_initializer_list<float>("float");
    test_initializer_list<double>("double");

    return boost::math::test::report_errors();
}
