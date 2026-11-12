//  Copyright Nick Thompson, 2026
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
//  See: https://github.com/boostorg/math/issues/1637
//
//  For small v the quadrature used for the tail of the non-central t cdf lost the mass in the
//  s^(v-1) singularity at the origin: the complemented cdf was out by a factor of 2 for v = 0.001,
//  which in turn broke find_non_centrality.

#ifndef BOOST_MATH_BUILD_MODULE
#include <boost/math/distributions/non_central_t.hpp>
#include <boost/math/distributions/complement.hpp>
#include <boost/math/tools/precision.hpp>
#include <algorithm>
#else
import boost.math;
#endif

#include "math_unit_test.hpp"

template <class Real>
void test()
{
    using dist = boost::math::non_central_t_distribution<Real>;
    using boost::math::complement;

    // Rounding v = 0.001 to float already moves the result by 3 epsilon:
    const Real tol = (std::max)(20 * boost::math::tools::epsilon<Real>(), static_cast<Real>(1e-12));

    const Real v = static_cast<Real>(0.001L);
    const Real nc = static_cast<Real>(-3.086951096691082L);
    const Real t = 300;
    // Reference values are from mpmath at 50 digits, with d = -nc and H(s) = exp(-v s^2 / 2) ncdf(-(t s + d)):
    //   2 (v/2)^(v/2) / (v gamma(v/2)) * quad(lambda u: H(u**(1/v)), [0, ..., 1])
    const Real q = static_cast<Real>(1.000000448365218891339442e-3L);

    CHECK_RELATIVE_ERROR(q, cdf(complement(dist(v, nc), t)), tol);
    // The same tail by reflection:
    CHECK_RELATIVE_ERROR(q, cdf(dist(v, -nc), -t), tol);

    // The inverse is limited by how well v and q are represented in type Real:
    const Real find_tol = (std::max)(1000 * boost::math::tools::epsilon<Real>(), static_cast<Real>(1e-10));
    CHECK_RELATIVE_ERROR(nc, dist::find_non_centrality(complement(v, t, q)), find_tol);

    if (boost::math::tools::digits<Real>() >= 53)
    {
        // As reported: 1 - p loses 3 digits here.
        const Real p = static_cast<Real>(0.99899999955163478L);
        CHECK_RELATIVE_ERROR(nc, dist::find_non_centrality(v, t, p), static_cast<Real>(1e-9));
    }

    // The quadrature is used whenever the tail is small:
    const Real loose = (std::max)(20 * boost::math::tools::epsilon<Real>(), static_cast<Real>(1e-11));
    CHECK_RELATIVE_ERROR(static_cast<Real>(9.159270549984962703551346e-4L), cdf(complement(dist(static_cast<Real>(0.01L), nc), t)), loose);
    CHECK_RELATIVE_ERROR(static_cast<Real>(4.242090926480852647116002e-4L), cdf(complement(dist(static_cast<Real>(0.1L), nc), t)), loose);
    // Either side of v = 1, where the singularity goes, and of v = 2, where the quadrature changes:
    CHECK_RELATIVE_ERROR(static_cast<Real>(7.498258893378480277293159e-7L), cdf(complement(dist(static_cast<Real>(0.9990234375L), nc), t)), loose);
    CHECK_RELATIVE_ERROR(static_cast<Real>(7.451097257679060976066953e-7L), cdf(complement(dist(static_cast<Real>(1), nc), t)), loose);
    CHECK_RELATIVE_ERROR(static_cast<Real>(1.634578958865346598692263e-9L), cdf(complement(dist(static_cast<Real>(1.9990234375L), nc), t)), loose);
    CHECK_RELATIVE_ERROR(static_cast<Real>(1.625247361717577498979435e-9L), cdf(complement(dist(static_cast<Real>(2), nc), t)), loose);
}

int main()
{
    test<float>();
    test<double>();

    return boost::math::test::report_errors();
}
