// Copyright Nick Thompson, 2021
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Module build check for bilinear_uniform. Exercises the exported interpolator
// through `import boost.math;` so that its templates (and the std::shared_ptr /
// std::vector members they rely on) resolve when the test is compiled against
// the boost.math named module. The assertions use only mathematically certain
// properties of bilinear interpolation: a constant field is reproduced
// everywhere, the interpolant reproduces its node values exactly, and bilinear
// interpolation is exact for functions of the form a + b*x + c*y + d*x*y. The
// unit-square expected values reuse the formula from bilinear_uniform_test.cpp.

#ifndef BOOST_MATH_BUILD_MODULE
#include <boost/math/interpolators/bilinear_uniform.hpp>
#else
import boost.math;
#endif

#ifndef BOOST_MATH_BUILD_MODULE
#include <vector>
#include <utility>
#include <iostream>
#endif
#include "math_unit_test.hpp"

using boost::math::interpolators::bilinear_uniform;

template <class Real>
void test_bilinear(const char* type_name)
{
    std::cout << "Testing bilinear_uniform on type " << type_name << std::endl;

    // A constant field is reproduced exactly at every query point.
    {
        const Real value = static_cast<Real>(1.5);
        std::vector<Real> v(2 * 2, value);
        auto ub = bilinear_uniform<std::vector<Real>>(std::move(v), 2, 2);
        CHECK_ULP_CLOSE(value, ub(static_cast<Real>(0), static_cast<Real>(0)), 0);
        CHECK_ULP_CLOSE(value, ub(static_cast<Real>(0.5), static_cast<Real>(0.25)), 1);
        CHECK_ULP_CLOSE(value, ub(static_cast<Real>(1), static_cast<Real>(1)), 0);
    }

    // Unit square with corner values {1, 2, 3, 4} stored row major:
    // f(x0 + i*dx, y0 + j*dy) lives at f[j*cols + i], so this is
    //   f[0]=1 at (0,0), f[1]=2 at (1,0), f[2]=3 at (0,1), f[3]=4 at (1,1).
    // Expected values follow the unit-square bilinear formula used in
    // bilinear_uniform_test.cpp:
    //   f(x,y) = v0*(1-x)*(1-y) + v1*x*(1-y) + v2*(1-x)*y + v3*x*y
    {
        std::vector<Real> v { static_cast<Real>(1), static_cast<Real>(2),
                              static_cast<Real>(3), static_cast<Real>(4) };
        auto ub = bilinear_uniform<std::vector<Real>>(std::move(v), 2, 2);

        // Node reproduction at the four corners.
        CHECK_ULP_CLOSE(static_cast<Real>(1), ub(static_cast<Real>(0), static_cast<Real>(0)), 0);
        CHECK_ULP_CLOSE(static_cast<Real>(2), ub(static_cast<Real>(1), static_cast<Real>(0)), 0);
        CHECK_ULP_CLOSE(static_cast<Real>(3), ub(static_cast<Real>(0), static_cast<Real>(1)), 0);
        CHECK_ULP_CLOSE(static_cast<Real>(4), ub(static_cast<Real>(1), static_cast<Real>(1)), 0);

        // Center and edge midpoints (all dyadic, so bit-exact in practice).
        CHECK_ULP_CLOSE(static_cast<Real>(2.5), ub(static_cast<Real>(0.5), static_cast<Real>(0.5)), 1);
        CHECK_ULP_CLOSE(static_cast<Real>(1.5), ub(static_cast<Real>(0.5), static_cast<Real>(0)), 1);
        CHECK_ULP_CLOSE(static_cast<Real>(3.5), ub(static_cast<Real>(0.5), static_cast<Real>(1)), 1);
        CHECK_ULP_CLOSE(static_cast<Real>(2), ub(static_cast<Real>(0), static_cast<Real>(0.5)), 1);
        CHECK_ULP_CLOSE(static_cast<Real>(3), ub(static_cast<Real>(1), static_cast<Real>(0.5)), 1);

        // Interior dyadic point (0.25, 0.75) -> 2.75.
        CHECK_ULP_CLOSE(static_cast<Real>(2.75), ub(static_cast<Real>(0.25), static_cast<Real>(0.75)), 2);
    }

    // Bilinear interpolation is exact for f(x,y) = 1 + 2x + 3y + 4xy.
    // Sample it on a 3x3 unit grid (x,y in {0,1,2}); interior queries must
    // reproduce f exactly. Node layout f[j*cols + i]:
    //   j=0: 1, 3, 5   j=1: 4, 10, 16   j=2: 7, 17, 27
    {
        std::vector<Real> v { static_cast<Real>(1),  static_cast<Real>(3),  static_cast<Real>(5),
                              static_cast<Real>(4),  static_cast<Real>(10), static_cast<Real>(16),
                              static_cast<Real>(7),  static_cast<Real>(17), static_cast<Real>(27) };
        auto ub = bilinear_uniform<std::vector<Real>>(std::move(v), 3, 3);

        // Node reproduction.
        CHECK_ULP_CLOSE(static_cast<Real>(10), ub(static_cast<Real>(1), static_cast<Real>(1)), 0);
        CHECK_ULP_CLOSE(static_cast<Real>(27), ub(static_cast<Real>(2), static_cast<Real>(2)), 0);

        // Interior points where bilinear interpolation reproduces f exactly.
        CHECK_ULP_CLOSE(static_cast<Real>(4.5),  ub(static_cast<Real>(0.5), static_cast<Real>(0.5)), 2);
        CHECK_ULP_CLOSE(static_cast<Real>(8.5),  ub(static_cast<Real>(1.5), static_cast<Real>(0.5)), 2);
        CHECK_ULP_CLOSE(static_cast<Real>(9.5),  ub(static_cast<Real>(0.5), static_cast<Real>(1.5)), 2);
        CHECK_ULP_CLOSE(static_cast<Real>(17.5), ub(static_cast<Real>(1.5), static_cast<Real>(1.5)), 2);
    }
}

int main()
{
    test_bilinear<float>("float");
    test_bilinear<double>("double");

    return boost::math::test::report_errors();
}
