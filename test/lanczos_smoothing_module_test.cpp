//  (C) Copyright Nick Thompson 2019.
//  (C) Copyright Matt Borland 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Module build check for boost::math::differentiation::discrete_lanczos_derivative.
// The exported class template builds its interior and boundary filters by calling
// the namespace-scope templates in differentiation::detail (interior_velocity_filter,
// boundary_velocity_filter, acceleration_filter). Constructing and evaluating the
// derivative here forces those internal templates to be instantiated and resolved
// from a module consumer. Data and expected values are copied from
// lanczos_smoothing_test.cpp; only mathematically-certain results are checked.

#ifndef BOOST_MATH_BUILD_MODULE
#include <boost/math/differentiation/lanczos_smoothing.hpp>
#else
import boost.math;
#endif

#ifndef BOOST_MATH_BUILD_MODULE
#include <vector>
#include <cstddef>
#endif
#include "math_unit_test.hpp"

// First derivative (order 1) on unit-spacing samples, double precision.
void test_velocity()
{
    using boost::math::differentiation::discrete_lanczos_derivative;

    // Default filter: n = 18, approximation_order = 3.
    auto d = discrete_lanczos_derivative<double>(1.0);
    CHECK_EQUAL(1.0, d.get_spacing());

    // Constant signal: the first derivative is identically zero.
    std::vector<double> constant(500, 7.0);
    auto dc = d(constant);
    for (std::size_t i = 0; i < dc.size(); ++i)
    {
        CHECK_ABSOLUTE_ERROR(0.0, dc[i], 1e-9);
    }

    // Linear signal v = 9 i + 6: the derivative is the slope 9 everywhere
    // (data and expected value from test_data_representations). A generous
    // absolute tolerance keeps the check robust where the internal long double
    // filter build is only 64-bit wide.
    std::vector<double> linear(500);
    for (std::size_t i = 0; i < linear.size(); ++i)
    {
        linear[i] = 9.0 * static_cast<double>(i) + 6.0;
    }
    auto dl = d(linear);
    for (std::size_t i = 0; i < dl.size(); ++i)
    {
        CHECK_ABSOLUTE_ERROR(9.0, dl[i], 1e-6);
    }

    // Quadratic signal v = 15 i^2 + 7 i + 8: the derivative is 30 i + 7 (test_boundary_lanczos).
    // Checked on the well-conditioned interior where the symmetric filter applies.
    std::vector<double> quad(500);
    for (std::size_t i = 0; i < quad.size(); ++i)
    {
        const double x = static_cast<double>(i);
        quad[i] = 15.0 * x * x + 7.0 * x + 8.0;
    }
    auto dq = d(quad);
    const std::size_t interior[] = {std::size_t(50), std::size_t(100), std::size_t(250), std::size_t(400)};
    for (const std::size_t i : interior)
    {
        const double expected = 30.0 * static_cast<double>(i) + 7.0;
        CHECK_MOLLIFIED_CLOSE(expected, dq[i], 1e-9);
        // The single-index operator() must agree with the whole-vector transform.
        CHECK_MOLLIFIED_CLOSE(expected, d(quad, i), 1e-9);
    }
}

// Rescaling: halving the sample rate doubles the reported derivative (test_rescaling).
void test_rescaling()
{
    using boost::math::differentiation::discrete_lanczos_derivative;

    std::vector<double> v(500);
    for (std::size_t i = 0; i < v.size(); ++i)
    {
        const double x = static_cast<double>(i);
        v[i] = 7.0 * x * x + 9.0 * x + 6.0;
    }

    auto lanczos1 = discrete_lanczos_derivative<double>(1.0);
    auto lanczos2 = discrete_lanczos_derivative<double>(2.0);
    auto d1 = lanczos1(v);
    auto d2 = lanczos2(v);
    for (std::size_t i = 0; i < v.size(); ++i)
    {
        CHECK_MOLLIFIED_CLOSE(d1[i], 2.0 * d2[i], 1e-12);
    }
}

// Second derivative (order 2, acceleration filter), double precision (test_lanczos_acceleration).
void test_acceleration()
{
    using boost::math::differentiation::discrete_lanczos_derivative;

    auto lanczos = discrete_lanczos_derivative<double, 2>(1.0, 4, 3);
    CHECK_EQUAL(1.0, lanczos.get_spacing());

    // Constant signal: the second derivative is zero.
    std::vector<double> constant(100, 7.0);
    for (std::size_t i = 0; i < constant.size(); ++i)
    {
        CHECK_ABSOLUTE_ERROR(0.0, lanczos(constant, i), 1e-9);
    }

    // Quadratic signal v = 7 i^2 + 9 i + 6: the second derivative is exactly 14.
    std::vector<double> quad(100);
    for (std::size_t i = 0; i < quad.size(); ++i)
    {
        const double x = static_cast<double>(i);
        quad[i] = 7.0 * x * x + 9.0 * x + 6.0;
    }
    for (std::size_t i = 0; i < quad.size(); ++i)
    {
        CHECK_MOLLIFIED_CLOSE(14.0, lanczos(quad, i), 1e-8);
    }
}

// Exercise the float instantiation of the exported template through the module.
// The interior first derivative of a linear ramp reproduces the slope.
void test_velocity_float()
{
    using boost::math::differentiation::discrete_lanczos_derivative;

    auto d = discrete_lanczos_derivative<float>(1.0f);

    std::vector<float> linear(100);
    for (std::size_t i = 0; i < linear.size(); ++i)
    {
        linear[i] = 9.0f * static_cast<float>(i) + 6.0f;
    }

    // Interior indices for the default n = 18 filter satisfy 18 <= i <= 81.
    const std::size_t interior[] = {std::size_t(20), std::size_t(40), std::size_t(60), std::size_t(80)};
    for (const std::size_t i : interior)
    {
        CHECK_MOLLIFIED_CLOSE(9.0f, d(linear, i), 1e-3f);
    }
}

int main()
{
    test_velocity();
    test_rescaling();
    test_acceleration();
    test_velocity_float();

    return boost::math::test::report_errors();
}
