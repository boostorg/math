// Copyright Nick Thompson, 2017
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Module build check for cardinal_cubic_b_spline. Consuming the interpolator
// through `import boost.math;` confirms the exported class template and its
// operator()/prime/double_prime members are usable from a module consumer.
// The assertions below exercise only mathematically-certain spline properties
// (node reproduction, exact reproduction of constants and affine functions),
// with expected values taken from cardinal_cubic_b_spline_test.cpp.

#ifndef BOOST_MATH_BUILD_MODULE
#include <boost/math/interpolators/cardinal_cubic_b_spline.hpp>
#else
import boost.math;
#endif

#include "math_unit_test.hpp"

#ifndef BOOST_MATH_BUILD_MODULE
#include <vector>
#include <limits>
#include <cmath>
#include <cstddef>
#endif

// s(x_j) == f(x_j) at every node: the defining interpolation condition, which
// holds regardless of the endpoint-derivative choice.
template <class Real>
void test_interpolation_condition(const char* type_name)
{
    std::cout << "Testing cardinal_cubic_b_spline interpolation condition for type " << type_name << std::endl;

    using std::sin;

    const Real x0 {static_cast<Real>(1)};
    const Real step {static_cast<Real>(0.125)};

    std::vector<Real> v(20);
    for (std::size_t i {0}; i < v.size(); ++i)
    {
        v[i] = sin(x0 + step * static_cast<Real>(i));
    }

    boost::math::interpolators::cardinal_cubic_b_spline<Real> spline(v.data(), v.size(), x0, step);

    for (std::size_t i {0}; i < v.size(); ++i)
    {
        const Real y {spline(x0 + step * static_cast<Real>(i))};
        CHECK_ULP_CLOSE(v[i], y, 64);
    }
}

// Constants are reproduced exactly and their derivatives vanish. Specifying
// zero endpoint derivatives makes the reproduction exact.
template <class Real>
void test_constant_function(const char* type_name)
{
    std::cout << "Testing cardinal_cubic_b_spline constant reproduction for type " << type_name << std::endl;

    const Real constant {static_cast<Real>(50.2)};
    const Real a {static_cast<Real>(5)};
    const Real step {static_cast<Real>(0.02)};

    std::vector<Real> v(20, constant);

    boost::math::interpolators::cardinal_cubic_b_spline<Real> spline(v.data(), v.size(), a, step, static_cast<Real>(0), static_cast<Real>(0));

    const Real eps {std::numeric_limits<Real>::epsilon()};
    for (std::size_t i {0}; i < v.size(); ++i)
    {
        const Real arg {a + step * static_cast<Real>(i) + static_cast<Real>(0.002)};
        CHECK_ULP_CLOSE(constant, spline(arg), 16);
        CHECK_ABSOLUTE_ERROR(static_cast<Real>(0), spline.prime(arg), 16 * eps);
        // double_prime carries a 1/h^2 factor, so the residual of the (exactly
        // zero) second derivative is amplified; keep this tolerance loose.
        CHECK_ABSOLUTE_ERROR(static_cast<Real>(0), spline.double_prime(arg), 5000 * eps);
    }
}

// Cubic b splines reproduce affine functions exactly, and the derivative of
// f(x) = m*x + b is the constant slope m.
template <class Real>
void test_affine_function(const char* type_name)
{
    std::cout << "Testing cardinal_cubic_b_spline affine reproduction for type " << type_name << std::endl;

    using std::sqrt;

    const Real m {static_cast<Real>(3)};
    const Real b {static_cast<Real>(2)};
    const Real step {static_cast<Real>(0.05)};

    std::vector<Real> v(20);
    for (std::size_t i {0}; i < v.size(); ++i)
    {
        v[i] = m * (step * static_cast<Real>(i)) + b;
    }

    // Specify the exact endpoint derivatives so the affine fit is exact.
    boost::math::interpolators::cardinal_cubic_b_spline<Real> spline(v.data(), v.size(), static_cast<Real>(0), step, m, m);

    const Real tol {100 * sqrt(std::numeric_limits<Real>::epsilon())};
    for (std::size_t i {0}; i + 1 < v.size(); ++i)
    {
        const Real arg {step * static_cast<Real>(i) + static_cast<Real>(0.001)};
        CHECK_ABSOLUTE_ERROR(m * arg + b, spline(arg), tol);
        CHECK_ABSOLUTE_ERROR(m, spline.prime(arg), tol);
    }
}

int main()
{
    test_interpolation_condition<float>("float");
    test_interpolation_condition<double>("double");

    test_constant_function<float>("float");
    test_constant_function<double>("double");

    test_affine_function<float>("float");
    test_affine_function<double>("double");

    return boost::math::test::report_errors();
}
