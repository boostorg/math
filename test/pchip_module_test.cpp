/*
 * Copyright Nick Thompson, 2020
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */

// Module build check for the pchip interpolator. Consuming it through
// `import boost.math;` confirms the exported class template and its
// operator()/prime/domain members are usable from a module consumer.
// The assertions below exercise only mathematically-certain pchip
// properties (exact reproduction of constants and affine functions and the
// interpolation condition at the nodes), with data and tolerances taken
// from pchip_test.cpp.

#ifndef BOOST_MATH_BUILD_MODULE
#include <boost/math/interpolators/pchip.hpp>
#else
import boost.math;
#endif

#include "math_unit_test.hpp"

#ifndef BOOST_MATH_BUILD_MODULE
#include <vector>
#include <utility>
#include <cstddef>
#include <iostream>
#endif

using boost::math::interpolators::pchip;

// A pchip built from a constant sample reproduces that constant exactly and
// its derivative vanishes; the domain spans the first and last abscissa.
template <class Real>
void test_constant(const char* type_name)
{
    std::cout << "Testing pchip constant reproduction for type " << type_name << std::endl;

    std::vector<Real> x {0, 1, 2, 3, 9, 22, 81};
    std::vector<Real> y(x.size(), static_cast<Real>(7));

    auto spline = pchip(std::move(x), std::move(y));

    const auto dom = spline.domain();
    CHECK_ULP_CLOSE(static_cast<Real>(0), dom.first, 1);
    CHECK_ULP_CLOSE(static_cast<Real>(81), dom.second, 1);

    const std::vector<Real> sample {static_cast<Real>(0), static_cast<Real>(0.25),
                                    static_cast<Real>(1.5), static_cast<Real>(5.5),
                                    static_cast<Real>(15), static_cast<Real>(50),
                                    static_cast<Real>(81)};
    for (std::size_t i {0}; i < sample.size(); ++i)
    {
        CHECK_ULP_CLOSE(static_cast<Real>(7), spline(sample[i]), 2);
        CHECK_ULP_CLOSE(static_cast<Real>(0), spline.prime(sample[i]), 2);
    }
}

// pchip reproduces an affine function exactly: nodes and interior points land
// on the line and the derivative equals the constant slope.
template <class Real>
void test_linear(const char* type_name)
{
    std::cout << "Testing pchip affine reproduction for type " << type_name << std::endl;

    // Four-node line y = x: node values and interior samples (tolerances
    // copied from pchip_test.cpp).
    {
        const std::vector<Real> x {0, 1, 2, 3};
        const std::vector<Real> y {0, 1, 2, 3};

        auto xc = x;
        auto yc = y;
        auto spline = pchip(std::move(xc), std::move(yc));

        CHECK_ULP_CLOSE(y[0], spline(x[0]), 0);
        CHECK_ULP_CLOSE(static_cast<Real>(0.5), spline(static_cast<Real>(0.5)), 10);
        CHECK_ULP_CLOSE(y[1], spline(x[1]), 0);
        CHECK_ULP_CLOSE(static_cast<Real>(1.5), spline(static_cast<Real>(1.5)), 10);
        CHECK_ULP_CLOSE(y[2], spline(x[2]), 0);
        CHECK_ULP_CLOSE(static_cast<Real>(2.5), spline(static_cast<Real>(2.5)), 10);
        CHECK_ULP_CLOSE(y[3], spline(x[3]), 0);
    }

    // A finely sampled line y = x has slope exactly one everywhere; pchip
    // reproduces both the value and the derivative exactly.
    {
        std::vector<Real> x(45);
        std::vector<Real> y(45);
        for (std::size_t i {0}; i < x.size(); ++i)
        {
            x[i] = static_cast<Real>(i);
            y[i] = static_cast<Real>(i);
        }

        auto spline = pchip(std::move(x), std::move(y));
        for (Real t {0}; t < static_cast<Real>(44); t += static_cast<Real>(0.5))
        {
            CHECK_ULP_CLOSE(t, spline(t), 0);
            CHECK_ULP_CLOSE(static_cast<Real>(1), spline.prime(t), 0);
        }
    }
}

// The interpolation condition s(x_j) == y_j holds at every node regardless of
// the (here non-linear) sample or the endpoint-derivative choice.
template <class Real>
void test_interpolation_condition(const char* type_name)
{
    std::cout << "Testing pchip interpolation condition for type " << type_name << std::endl;

    // y = x*x, exact at integer nodes.
    const std::vector<Real> x {0, 1, 2, 3, 4};
    const std::vector<Real> y {0, 1, 4, 9, 16};

    auto xc = x;
    auto yc = y;
    auto spline = pchip(std::move(xc), std::move(yc));
    for (std::size_t i {0}; i < x.size(); ++i)
    {
        CHECK_ULP_CLOSE(y[i], spline(x[i]), 2);
    }

    // Endpoint derivatives do not affect the interpolation condition.
    xc = x;
    yc = y;
    spline = pchip(std::move(xc), std::move(yc), static_cast<Real>(0), static_cast<Real>(0));
    for (std::size_t i {0}; i < x.size(); ++i)
    {
        CHECK_ULP_CLOSE(y[i], spline(x[i]), 2);
    }
}

int main()
{
    test_constant<float>("float");
    test_constant<double>("double");

    test_linear<float>("float");
    test_linear<double>("double");

    test_interpolation_condition<float>("float");
    test_interpolation_condition<double>("double");

    return boost::math::test::report_errors();
}
