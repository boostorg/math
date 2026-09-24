// Copyright Nick Thompson, 2020
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Module build check for cubic_hermite. Consuming the interpolators through
// `import boost.math;` confirms the exported class templates cubic_hermite,
// cardinal_cubic_hermite and cardinal_cubic_hermite_aos, along with their
// operator()/prime/domain members, are usable from a module consumer.
// The assertions exercise only mathematically-certain spline properties:
// exact reproduction of constant, affine and quadratic functions (a cubic
// Hermite segment is the unique cubic matching value and derivative at both
// endpoints, so it reproduces any polynomial of degree <= 3), plus node
// reproduction. Expected values and tolerances are taken from
// cubic_hermite_test.cpp.

#ifndef BOOST_MATH_BUILD_MODULE
#include <boost/math/interpolators/cubic_hermite.hpp>
#else
import boost.math;
#endif

#include "math_unit_test.hpp"

#ifndef BOOST_MATH_BUILD_MODULE
#include <vector>
#include <array>
#include <utility>
#include <cstddef>
#endif

using boost::math::interpolators::cubic_hermite;
using boost::math::interpolators::cardinal_cubic_hermite;
using boost::math::interpolators::cardinal_cubic_hermite_aos;

// A constant is reproduced exactly and its derivative vanishes everywhere.
template <class Real>
void test_constant()
{
    std::vector<Real> x {0, 1, 2, 3, 9, 22, 81};
    std::vector<Real> y(x.size(), static_cast<Real>(7));
    std::vector<Real> dydx(x.size(), static_cast<Real>(0));

    auto s = cubic_hermite(std::move(x), std::move(y), std::move(dydx));

    const Real nodes[] {0, 1, 2, 3, 9, 22, 81};
    for (auto node : nodes)
    {
        CHECK_ULP_CLOSE(static_cast<Real>(7), s(node), 2);
        CHECK_ULP_CLOSE(static_cast<Real>(0), s.prime(node), 2);
    }

    // Interior samples: the constant and its zero slope are reproduced there too.
    CHECK_ULP_CLOSE(static_cast<Real>(7), s(static_cast<Real>(1.5)), 2);
    CHECK_ULP_CLOSE(static_cast<Real>(7), s(static_cast<Real>(15)), 2);
    CHECK_ULP_CLOSE(static_cast<Real>(0), s.prime(static_cast<Real>(1.5)), 2);
    CHECK_ULP_CLOSE(static_cast<Real>(0), s.prime(static_cast<Real>(15)), 2);

    auto [lo, hi] = s.domain();
    CHECK_ULP_CLOSE(static_cast<Real>(0), lo, 0);
    CHECK_ULP_CLOSE(static_cast<Real>(81), hi, 0);
}

// f(x) = x with unit slope: nodes are reproduced exactly and the straight line
// is reproduced at interior points.
template <class Real>
void test_linear()
{
    std::vector<Real> x {0, 1, 2, 3};
    std::vector<Real> y {0, 1, 2, 3};
    std::vector<Real> dydx {1, 1, 1, 1};

    auto s = cubic_hermite(std::move(x), std::move(y), std::move(dydx));

    CHECK_ULP_CLOSE(static_cast<Real>(0), s(static_cast<Real>(0)), 0);
    CHECK_ULP_CLOSE(static_cast<Real>(1) / static_cast<Real>(2), s(static_cast<Real>(1) / static_cast<Real>(2)), 10);
    CHECK_ULP_CLOSE(static_cast<Real>(1), s(static_cast<Real>(1)), 0);
    CHECK_ULP_CLOSE(static_cast<Real>(3) / static_cast<Real>(2), s(static_cast<Real>(3) / static_cast<Real>(2)), 10);
    CHECK_ULP_CLOSE(static_cast<Real>(2), s(static_cast<Real>(2)), 0);
    CHECK_ULP_CLOSE(static_cast<Real>(5) / static_cast<Real>(2), s(static_cast<Real>(5) / static_cast<Real>(2)), 10);
    CHECK_ULP_CLOSE(static_cast<Real>(3), s(static_cast<Real>(3)), 0);

    // The derivative of an affine function is the constant slope.
    CHECK_ULP_CLOSE(static_cast<Real>(1), s.prime(static_cast<Real>(1) / static_cast<Real>(2)), 2);
    CHECK_ULP_CLOSE(static_cast<Real>(1), s.prime(static_cast<Real>(5) / static_cast<Real>(2)), 2);
}

// A cubic Hermite interpolant reproduces quadratics exactly, so with
// y = x*x/2 and dydx = x we recover s(t) = t*t/2 and s'(t) = t between nodes.
template <class Real>
void test_quadratic()
{
    std::vector<Real> x {0, 1, 2, 3, 4};
    std::vector<Real> y {0, static_cast<Real>(0.5), 2, static_cast<Real>(4.5), 8};
    std::vector<Real> dydx {0, 1, 2, 3, 4};

    auto s = cubic_hermite(std::move(x), std::move(y), std::move(dydx));

    const Real samples[] {static_cast<Real>(0.5), static_cast<Real>(1.5), static_cast<Real>(2.5), static_cast<Real>(3.5)};
    for (auto t : samples)
    {
        CHECK_ULP_CLOSE(t * t / static_cast<Real>(2), s(t), 10);
        CHECK_ULP_CLOSE(t, s.prime(t), 150);
    }
}

// The evenly spaced (cardinal) forms reproduce the same affine function, both
// the split y/dydx interface and the array-of-structs interface.
template <class Real>
void test_cardinal_linear()
{
    const Real x0 {static_cast<Real>(0)};
    const Real dx {static_cast<Real>(1)};

    std::vector<Real> y {0, 1, 2, 3};
    std::vector<Real> dydx {1, 1, 1, 1};

    auto s = cardinal_cubic_hermite(std::move(y), std::move(dydx), x0, dx);

    CHECK_ULP_CLOSE(static_cast<Real>(0), s(static_cast<Real>(0)), 0);
    CHECK_ULP_CLOSE(static_cast<Real>(1) / static_cast<Real>(2), s(static_cast<Real>(1) / static_cast<Real>(2)), 10);
    CHECK_ULP_CLOSE(static_cast<Real>(1), s(static_cast<Real>(1)), 0);
    CHECK_ULP_CLOSE(static_cast<Real>(3) / static_cast<Real>(2), s(static_cast<Real>(3) / static_cast<Real>(2)), 10);
    CHECK_ULP_CLOSE(static_cast<Real>(2), s(static_cast<Real>(2)), 0);
    CHECK_ULP_CLOSE(static_cast<Real>(5) / static_cast<Real>(2), s(static_cast<Real>(5) / static_cast<Real>(2)), 10);
    CHECK_ULP_CLOSE(static_cast<Real>(3), s(static_cast<Real>(3)), 0);
    CHECK_ULP_CLOSE(static_cast<Real>(1), s.prime(static_cast<Real>(3) / static_cast<Real>(2)), 2);

    // Array-of-structs interface: data[i] = {y_i, dydx_i}.
    std::vector<std::array<Real, 2>> data(4);
    for (std::size_t i {0}; i < data.size(); ++i)
    {
        data[i][0] = static_cast<Real>(i);
        data[i][1] = static_cast<Real>(1);
    }

    auto saos = cardinal_cubic_hermite_aos(std::move(data), x0, dx);

    CHECK_ULP_CLOSE(static_cast<Real>(0), saos(static_cast<Real>(0)), 0);
    CHECK_ULP_CLOSE(static_cast<Real>(1) / static_cast<Real>(2), saos(static_cast<Real>(1) / static_cast<Real>(2)), 10);
    CHECK_ULP_CLOSE(static_cast<Real>(1), saos(static_cast<Real>(1)), 0);
    CHECK_ULP_CLOSE(static_cast<Real>(2), saos(static_cast<Real>(2)), 0);
    CHECK_ULP_CLOSE(static_cast<Real>(3), saos(static_cast<Real>(3)), 0);
    CHECK_ULP_CLOSE(static_cast<Real>(1), saos.prime(static_cast<Real>(3) / static_cast<Real>(2)), 2);
}

int main()
{
    test_constant<float>();
    test_constant<double>();

    test_linear<float>();
    test_linear<double>();

    test_quadratic<float>();
    test_quadratic<double>();

    test_cardinal_linear<float>();
    test_cardinal_linear<double>();

    return boost::math::test::report_errors();
}
