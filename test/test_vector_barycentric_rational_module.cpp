// Copyright Nick Thompson, 2019
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Module build check for vector_barycentric_rational. Exercises the exported
// interpolator through `import boost.math;` so that its templates (and the
// std::shared_ptr / std::vector members they rely on) resolve when the test is
// compiled against the boost.math named module. The point type is a plain
// std::array<Real, 2> so the check stays free of Eigen and uBLAS, neither of
// which the module exports. The assertions use only mathematically certain
// properties: the interpolant reproduces its node values exactly, and a
// constant sample set is interpolated by that constant with a zero derivative.

#ifndef BOOST_MATH_BUILD_MODULE
#include <boost/math/interpolators/vector_barycentric_rational.hpp>
#else
import boost.math;
#endif

#ifndef BOOST_MATH_BUILD_MODULE
#include <array>
#include <vector>
#include <utility>
#include <limits>
#include <cmath>
#include <cstddef>
#endif
#include "math_unit_test.hpp"

template <class Real>
void test_vector_barycentric(const char* type_name)
{
    std::cout << "Testing vector_barycentric_rational on type " << type_name << std::endl;

    using std::sqrt;
    using Point = std::array<Real, 2>;
    using Interpolator = boost::math::interpolators::vector_barycentric_rational<std::vector<Real>, std::vector<Point>>;

    // Node abscissas, strictly increasing as the constructor requires.
    std::vector<Real> t { static_cast<Real>(0), static_cast<Real>(1), static_cast<Real>(2),
                          static_cast<Real>(3), static_cast<Real>(4), static_cast<Real>(5),
                          static_cast<Real>(6), static_cast<Real>(7) };
    // Two independent components per node.
    std::vector<Real> a { static_cast<Real>(2.5), static_cast<Real>(-1.0), static_cast<Real>(3.25),
                          static_cast<Real>(0.5), static_cast<Real>(-4.0), static_cast<Real>(6.75),
                          static_cast<Real>(1.0), static_cast<Real>(-2.0) };
    std::vector<Real> b { static_cast<Real>(1.0), static_cast<Real>(2.0), static_cast<Real>(4.0),
                          static_cast<Real>(8.0), static_cast<Real>(16.0), static_cast<Real>(32.0),
                          static_cast<Real>(64.0), static_cast<Real>(128.0) };

    std::vector<Point> y(t.size());
    for (std::size_t i = 0; i < t.size(); ++i)
    {
        y[i][0] = a[i];
        y[i][1] = b[i];
    }

    // Interpolation condition: evaluating at a node returns the stored point
    // exactly (operator() short-circuits when t == t_[i]), so this is bit-exact.
    {
        std::vector<Real> tc = t;
        std::vector<Point> yc = y;
        Interpolator interp(std::move(tc), std::move(yc));

        for (std::size_t i = 0; i < t.size(); ++i)
        {
            // void out-parameter overload.
            Point z;
            interp(z, t[i]);
            CHECK_ULP_CLOSE(y[i][0], z[0], 0);
            CHECK_ULP_CLOSE(y[i][1], z[1], 0);

            // value-returning overload.
            Point z2 = interp(t[i]);
            CHECK_ULP_CLOSE(y[i][0], z2[0], 0);
            CHECK_ULP_CLOSE(y[i][1], z2[1], 0);

            // eval_with_prime reproduces the node value in its output argument.
            Point v, dvdt;
            interp.eval_with_prime(v, dvdt, t[i]);
            CHECK_ULP_CLOSE(y[i][0], v[0], 0);
            CHECK_ULP_CLOSE(y[i][1], v[1], 0);

            // pair-returning eval_with_prime reproduces the node value too.
            std::pair<Point, Point> pr = interp.eval_with_prime(t[i]);
            CHECK_ULP_CLOSE(y[i][0], pr.first[0], 0);
            CHECK_ULP_CLOSE(y[i][1], pr.first[1], 0);
        }
    }

    // Same interpolation condition with an explicit approximation order of 5,
    // which must be strictly less than the node count of 8.
    {
        std::vector<Real> tc = t;
        std::vector<Point> yc = y;
        Interpolator interp5(std::move(tc), std::move(yc), 5);

        for (std::size_t i = 0; i < t.size(); ++i)
        {
            Point z;
            interp5(z, t[i]);
            CHECK_ULP_CLOSE(y[i][0], z[0], 0);
            CHECK_ULP_CLOSE(y[i][1], z[1], 0);
        }
    }

    // A constant sample set is interpolated by that constant, with a zero
    // derivative, when evaluated away from the nodes.
    {
        const Real c0 = static_cast<Real>(-8);
        const Real c1 = static_cast<Real>(3.5);
        std::vector<Real> tc { static_cast<Real>(0), static_cast<Real>(1), static_cast<Real>(2),
                               static_cast<Real>(3), static_cast<Real>(4), static_cast<Real>(5) };
        std::vector<Point> yc(tc.size());
        for (std::size_t i = 0; i < tc.size(); ++i)
        {
            yc[i][0] = c0;
            yc[i][1] = c1;
        }
        Interpolator interpc(std::move(tc), std::move(yc));

        const Real prime_tol = sqrt(std::numeric_limits<Real>::epsilon());
        for (const Real s : { static_cast<Real>(0.5), static_cast<Real>(2.5), static_cast<Real>(4.5) })
        {
            Point z = interpc(s);
            CHECK_ULP_CLOSE(c0, z[0], 4);
            CHECK_ULP_CLOSE(c1, z[1], 4);

            Point dzdt = interpc.prime(s);
            CHECK_ABSOLUTE_ERROR(static_cast<Real>(0), dzdt[0], prime_tol);
            CHECK_ABSOLUTE_ERROR(static_cast<Real>(0), dzdt[1], prime_tol);
        }
    }
}

int main()
{
    test_vector_barycentric<float>("float");
    test_vector_barycentric<double>("double");

    return boost::math::test::report_errors();
}
