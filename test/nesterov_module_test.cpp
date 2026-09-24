//           Copyright Maksym Zhelyenzyakov 2025-2026.
// Distributed under the Boost Software License, Version 1.0.
//      (See accompanying file LICENSE_1_0.txt or copy at
//           https://www.boost.org/LICENSE_1_0.txt)
//
// Module build check for the Nesterov accelerated gradient component. The
// make_nag factory drives a reverse-mode autodiff objective whose rvar argument
// type the module does not export, so the exercised public entry point here is
// nesterov_update_policy: its scalar update rule is the mathematical core of the
// algorithm and crosses the module boundary. The values below are exact in
// binary floating point.

#ifndef BOOST_MATH_BUILD_MODULE
#include <boost/math/optimization/nesterov.hpp>
#else
import boost.math;
#endif

#include "math_unit_test.hpp"

template <typename T>
void test_spots()
{
    using boost::math::optimization::nesterov_update_policy;

    // lr = mu = 1/2 (exact). The scalar update, starting from velocity v = 0, is:
    //   v_prev = v;  v = mu * v - lr * g;  x += -mu * v_prev + (1 + mu) * v
    nesterov_update_policy<T> up {static_cast<T>(0.5), static_cast<T>(0.5)};
    CHECK_EQUAL(up.lr(), static_cast<T>(0.5));
    CHECK_EQUAL(up.mu(), static_cast<T>(0.5));

    T x {static_cast<T>(8)};
    T g {static_cast<T>(4)};
    T v {static_cast<T>(0)};

    // Step 1: v = -0.5 * 4 = -2, x = 8 + 1.5 * (-2) = 5.
    up(x, g, v);
    CHECK_EQUAL(v, static_cast<T>(-2));
    CHECK_EQUAL(x, static_cast<T>(5));

    // Step 2: v = 0.5 * (-2) - 0.5 * 4 = -3, x = 5 + (1 - 4.5) = 1.5.
    up(x, g, v);
    CHECK_EQUAL(v, static_cast<T>(-3));
    CHECK_EQUAL(x, static_cast<T>(1.5));
}

int main()
{
    test_spots<float>();
    test_spots<double>();

    return boost::math::test::report_errors();
}
