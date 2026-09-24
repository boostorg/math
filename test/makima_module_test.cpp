// Copyright Nick Thompson, 2020
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Module build check for makima (modified Akima) interpolation. Consuming the
// interpolator through `import boost.math;` confirms the exported class
// template and its operator()/prime members are usable from a module consumer.
// The assertions exercise only mathematically-certain makima properties (exact
// reproduction of constants and straight lines, node reproduction), with
// expected values and tolerances taken from makima_test.cpp. The circular
// buffer, long double, and multiprecision paths of that test are dropped: a
// module consumer may not include Boost.CircularBuffer or Boost.Multiprecision
// textually.

#ifndef BOOST_MATH_BUILD_MODULE
#include <boost/math/interpolators/makima.hpp>
#else
import boost.math;
#endif

#include "math_unit_test.hpp"

#ifndef BOOST_MATH_BUILD_MODULE
#include <vector>
#include <utility>
#include <cstddef>
#endif

using boost::math::interpolators::makima;

// A makima spline through constant data reproduces the constant exactly, and
// its derivative vanishes everywhere on the interval.
template <class Real>
void test_constant(const char* type_name)
{
    std::cout << "Testing makima constant reproduction for type " << type_name << std::endl;

    std::vector<Real> x {0, 1, 2, 3, 9, 22, 81};
    std::vector<Real> y(x.size(), static_cast<Real>(7));

    auto spline = makima(std::move(x), std::move(y));

    for (Real t {static_cast<Real>(0)}; t <= static_cast<Real>(81); t += static_cast<Real>(0.5))
    {
        CHECK_ULP_CLOSE(static_cast<Real>(7), spline(t), 2);
        CHECK_ULP_CLOSE(static_cast<Real>(0), spline.prime(t), 2);
    }
}

// makima reproduces straight-line data exactly: values match at and between the
// nodes, and the derivative equals the constant slope of the line.
template <class Real>
void test_linear(const char* type_name)
{
    std::cout << "Testing makima straight-line reproduction for type " << type_name << std::endl;

    std::vector<Real> x {0, 1, 2, 3};
    std::vector<Real> y {0, 1, 2, 3};

    auto xc = x;
    auto yc = y;
    auto spline = makima(std::move(xc), std::move(yc));

    // Node values are reproduced to the bit; interior points lie on the line.
    CHECK_ULP_CLOSE(static_cast<Real>(0), spline(static_cast<Real>(0)), 0);
    CHECK_ULP_CLOSE(static_cast<Real>(0.5), spline(static_cast<Real>(0.5)), 10);
    CHECK_ULP_CLOSE(static_cast<Real>(1), spline(static_cast<Real>(1)), 0);
    CHECK_ULP_CLOSE(static_cast<Real>(1.5), spline(static_cast<Real>(1.5)), 10);
    CHECK_ULP_CLOSE(static_cast<Real>(2), spline(static_cast<Real>(2)), 0);
    CHECK_ULP_CLOSE(static_cast<Real>(2.5), spline(static_cast<Real>(2.5)), 10);
    CHECK_ULP_CLOSE(static_cast<Real>(3), spline(static_cast<Real>(3)), 0);

    // A longer ramp x[i] = y[i] = i; the interpolant is the identity, slope 1.
    std::vector<Real> xr(45);
    std::vector<Real> yr(45);
    for (std::size_t i {0}; i < xr.size(); ++i)
    {
        xr[i] = static_cast<Real>(i);
        yr[i] = static_cast<Real>(i);
    }

    auto ramp = makima(std::move(xr), std::move(yr));
    for (Real t {static_cast<Real>(0)}; t < static_cast<Real>(44); t += static_cast<Real>(0.5))
    {
        CHECK_ULP_CLOSE(t, ramp(t), 0);
        CHECK_ULP_CLOSE(static_cast<Real>(1), ramp.prime(t), 0);
    }

    // The endpoint-derivative overload, given the true slope, is also exact.
    std::vector<Real> xe(45);
    std::vector<Real> ye(45);
    for (std::size_t i {0}; i < xe.size(); ++i)
    {
        xe[i] = static_cast<Real>(i);
        ye[i] = static_cast<Real>(i);
    }

    auto ends = makima(std::move(xe), std::move(ye), static_cast<Real>(1), static_cast<Real>(1));
    for (Real t {static_cast<Real>(0)}; t < static_cast<Real>(44); t += static_cast<Real>(0.5))
    {
        CHECK_ULP_CLOSE(t, ends(t), 0);
        CHECK_ULP_CLOSE(static_cast<Real>(1), ends.prime(t), 0);
    }
}

// The defining property of the interpolant: it reproduces the sample values at
// the nodes, independent of the endpoint-derivative choice.
template <class Real>
void test_interpolation_condition(const char* type_name)
{
    std::cout << "Testing makima interpolation condition for type " << type_name << std::endl;

    std::vector<Real> x {static_cast<Real>(0.1), static_cast<Real>(0.7), static_cast<Real>(1.3),
                         static_cast<Real>(2.0), static_cast<Real>(3.6), static_cast<Real>(5.1)};
    std::vector<Real> y {static_cast<Real>(1.2), static_cast<Real>(-0.3), static_cast<Real>(2.5),
                         static_cast<Real>(0.8), static_cast<Real>(1.9), static_cast<Real>(-1.1)};

    auto xc = x;
    auto yc = y;
    auto spline = makima(std::move(xc), std::move(yc));
    for (std::size_t i {0}; i < x.size(); ++i)
    {
        CHECK_ULP_CLOSE(y[i], spline(x[i]), 2);
    }

    // Prescribing endpoint derivatives must not disturb the interpolation.
    auto xc2 = x;
    auto yc2 = y;
    auto spline2 = makima(std::move(xc2), std::move(yc2), static_cast<Real>(0), static_cast<Real>(0));
    for (std::size_t i {0}; i < x.size(); ++i)
    {
        CHECK_ULP_CLOSE(y[i], spline2(x[i]), 2);
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
