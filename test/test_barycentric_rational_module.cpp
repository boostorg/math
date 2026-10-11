// Copyright Nick Thompson, 2017
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Module build check for barycentric_rational. Exercises the exported
// interpolator through `import boost.math;` so that its templates (and the
// std::shared_ptr / std::vector members they rely on) resolve when the test is
// compiled against the boost.math named module. The assertions use only
// mathematically certain properties: an interpolant reproduces its node values
// exactly, and a constant sample set is interpolated by that constant with a
// zero derivative.

#ifndef BOOST_MATH_BUILD_MODULE
#include <boost/math/interpolators/barycentric_rational.hpp>
#else
import boost.math;
#endif

#ifndef BOOST_MATH_BUILD_MODULE
#include <vector>
#include <utility>
#include <limits>
#include <cmath>
#include <cstddef>
#endif
#include "math_unit_test.hpp"

template <class Real>
void test_barycentric(const char* type_name)
{
    std::cout << "Testing barycentric_rational on type " << type_name << std::endl;

    using std::sqrt;

    // Interpolation condition: the interpolant reproduces its node values.
    // operator()(x_i) returns the stored y_i exactly, so this is bit-exact.
    std::vector<Real> x { static_cast<Real>(0), static_cast<Real>(1), static_cast<Real>(2),
                          static_cast<Real>(3), static_cast<Real>(4), static_cast<Real>(5),
                          static_cast<Real>(6) };
    std::vector<Real> y { static_cast<Real>(2.5), static_cast<Real>(-1.0), static_cast<Real>(3.25),
                          static_cast<Real>(0.5), static_cast<Real>(-4.0), static_cast<Real>(6.75),
                          static_cast<Real>(1.0) };

    boost::math::interpolators::barycentric_rational<Real> interp(x.data(), y.data(), y.size());
    for (std::size_t i = 0; i < x.size(); ++i)
    {
        CHECK_ULP_CLOSE(y[i], interp(x[i]), 0);
    }

    // Same condition with an explicit approximation order of 5 (needs order < node count).
    std::vector<Real> xh { static_cast<Real>(0), static_cast<Real>(1), static_cast<Real>(2),
                           static_cast<Real>(3), static_cast<Real>(4), static_cast<Real>(5),
                           static_cast<Real>(6), static_cast<Real>(7) };
    std::vector<Real> yh { static_cast<Real>(1), static_cast<Real>(2), static_cast<Real>(4),
                           static_cast<Real>(8), static_cast<Real>(16), static_cast<Real>(32),
                           static_cast<Real>(64), static_cast<Real>(128) };
    boost::math::interpolators::barycentric_rational<Real> interp5(xh.data(), yh.data(), yh.size(), 5);
    for (std::size_t i = 0; i < xh.size(); ++i)
    {
        CHECK_ULP_CLOSE(yh[i], interp5(xh[i]), 0);
    }

    // The move constructor produces an equivalent interpolant.
    std::vector<Real> xm = x;
    std::vector<Real> ym = y;
    boost::math::interpolators::barycentric_rational<Real> minterp(std::move(xm), std::move(ym));
    for (std::size_t i = 0; i < x.size(); ++i)
    {
        CHECK_ULP_CLOSE(y[i], minterp(x[i]), 0);
    }

    // A constant sample set is interpolated by that constant, with zero derivative.
    const Real c = static_cast<Real>(-8);
    std::vector<Real> xc { static_cast<Real>(0), static_cast<Real>(1), static_cast<Real>(2),
                           static_cast<Real>(3), static_cast<Real>(4), static_cast<Real>(5) };
    std::vector<Real> yc(xc.size(), c);
    boost::math::interpolators::barycentric_rational<Real> interpc(xc.data(), yc.data(), yc.size());

    const Real prime_tol = sqrt(std::numeric_limits<Real>::epsilon());
    for (const Real t : { static_cast<Real>(0.5), static_cast<Real>(2.5), static_cast<Real>(4.5) })
    {
        CHECK_ULP_CLOSE(c, interpc(t), 4);
        CHECK_ABSOLUTE_ERROR(static_cast<Real>(0), interpc.prime(t), prime_tol);
    }
}

int main()
{
    test_barycentric<float>("float");
    test_barycentric<double>("double");

    return boost::math::test::report_errors();
}
