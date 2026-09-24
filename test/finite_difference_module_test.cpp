//  (C) Copyright Nick Thompson 2018.
//  (C) Copyright Matt Borland 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Module build check for the finite-difference numerical differentiation
// component. Consuming boost.math as a named module, this exercises both public
// entry points, finite_difference_derivative across every supported order and
// complex_step_derivative, using functions whose derivatives are exact and
// mathematically certain (polynomials) so no library helpers are needed to
// produce the expected values.

#ifndef BOOST_MATH_BUILD_MODULE
#include <boost/math/differentiation/finite_difference.hpp>
#else
import boost.math;
#endif

#ifndef BOOST_MATH_BUILD_MODULE
#include <complex>
#include <iostream>
#endif
#include "math_unit_test.hpp"

using boost::math::differentiation::finite_difference_derivative;
using boost::math::differentiation::complex_step_derivative;

// tol_fd bounds the finite-difference approximation error (order 1 is the
// loosest at ~2*sqrt(eps)); tol_cs bounds the essentially exact complex step.
template <class T>
void test_spots(const char* type_name, T tol_fd, T tol_cs)
{
    std::cout << "Testing finite_difference for type " << type_name << std::endl;

    // f(x) = x*x has f'(x) = 2x, so f'(3) = 6. Exercise every supported order.
    auto square = [](T t) { return t * t; };

    CHECK_ABSOLUTE_ERROR(static_cast<T>(6), (finite_difference_derivative<decltype(square), T, 1>(square, static_cast<T>(3))), tol_fd);
    CHECK_ABSOLUTE_ERROR(static_cast<T>(6), (finite_difference_derivative<decltype(square), T, 2>(square, static_cast<T>(3))), tol_fd);
    CHECK_ABSOLUTE_ERROR(static_cast<T>(6), (finite_difference_derivative<decltype(square), T, 4>(square, static_cast<T>(3))), tol_fd);
    CHECK_ABSOLUTE_ERROR(static_cast<T>(6), (finite_difference_derivative<decltype(square), T, 6>(square, static_cast<T>(3))), tol_fd);
    CHECK_ABSOLUTE_ERROR(static_cast<T>(6), (finite_difference_derivative<decltype(square), T, 8>(square, static_cast<T>(3))), tol_fd);

    // The error-estimate output argument is populated on the default (order 6) path.
    T error_estimate {};
    const T d {finite_difference_derivative<decltype(square), T, 6>(square, static_cast<T>(3), &error_estimate)};
    CHECK_ABSOLUTE_ERROR(static_cast<T>(6), d, tol_fd);
    CHECK_GE(error_estimate, static_cast<T>(0));

    // f(x) = x*x*x has f'(x) = 3x*x, so f'(2) = 12 (via the default order).
    auto cube = [](T t) { return t * t * t; };
    CHECK_ABSOLUTE_ERROR(static_cast<T>(12), (finite_difference_derivative<decltype(cube), T>(cube, static_cast<T>(2))), tol_fd);

    // A straight line f(x) = 2x + 1 has f'(x) = 2 everywhere.
    auto line = [](T t) { return static_cast<T>(2) * t + static_cast<T>(1); };
    CHECK_ABSOLUTE_ERROR(static_cast<T>(2), (finite_difference_derivative<decltype(line), T, 4>(line, static_cast<T>(5))), tol_fd);

    // Complex-step differentiation is well conditioned and exact for polynomials.
    auto complex_square = [](std::complex<T> z) { return z * z; };
    CHECK_ABSOLUTE_ERROR(static_cast<T>(6), complex_step_derivative(complex_square, static_cast<T>(3)), tol_cs);

    auto complex_cube = [](std::complex<T> z) { return z * z * z; };
    CHECK_ABSOLUTE_ERROR(static_cast<T>(12), complex_step_derivative(complex_cube, static_cast<T>(2)), tol_cs);
}

int main()
{
    test_spots<double>("double", 1e-5, 1e-11);
    test_spots<float>("float", 1e-2f, 1e-5f);

    return boost::math::test::report_errors();
}
