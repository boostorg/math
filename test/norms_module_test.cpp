/*
 *  (C) Copyright Nick Thompson 2018.
 *  Copyright Matt Borland 2026.
 *  Use, modification and distribution are subject to the
 *  Boost Software License, Version 1.0. (See accompanying file
 *  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */
//
// Module build check for the vector norms in <boost/math/tools/norms.hpp>.
// Consuming the exported free functions through the boost.math module exercises
// sup_norm, l1_norm, l2_norm, lp_norm, l0_pseudo_norm, hamming_distance,
// total_variation and the l1/l2/lp/sup distances. Every expected value here is
// mathematically certain (sums, a perfect square root, telescoping differences)
// and mirrors the double/float assertions in norms_test.cpp. The
// multiprecision, complex, integer and Boost.Test parts of that test are dropped.

#ifndef BOOST_MATH_BUILD_MODULE
#include <boost/math/tools/norms.hpp>
#else
import boost.math;
#endif

#ifndef BOOST_MATH_BUILD_MODULE
#include <vector>
#include <array>
#include <cstddef>
#endif
#include "math_unit_test.hpp"

template <class Real>
void test_norms()
{
    // sup_norm: largest magnitude. |-2| beats |1| and |0|.
    std::vector<Real> v {-2, 1, 0};
    CHECK_ULP_CLOSE(Real(2), boost::math::tools::sup_norm(v), 1);
    std::array<Real, 3> va {-2, 1, 0};
    CHECK_ULP_CLOSE(Real(2), boost::math::tools::sup_norm(va), 1);

    // l1_norm: sum of magnitudes.
    std::vector<Real> ones {1, 1, 1};
    CHECK_ULP_CLOSE(Real(3), boost::math::tools::l1_norm(ones), 1);

    // l2_norm: sqrt of the sum of squares; four ones give sqrt(4) = 2 exactly.
    std::vector<Real> four_ones {1, 1, 1, 1};
    CHECK_ULP_CLOSE(Real(2), boost::math::tools::l2_norm(four_ones), 1);

    // lp_norm: a single nonzero entry reproduces its magnitude for any p.
    std::array<Real, 3> u {1, 0, 0};
    CHECK_ULP_CLOSE(Real(1), boost::math::tools::lp_norm(u.begin(), u.end(), 3), 2);
    u[0] = -8;
    CHECK_ULP_CLOSE(Real(8), boost::math::tools::lp_norm(u.cbegin(), u.cend(), 3), 16);

    // l0_pseudo_norm: count of nonzero entries.
    std::vector<Real> sparse {0, 0, 1};
    CHECK_EQUAL(std::size_t(1), boost::math::tools::l0_pseudo_norm(sparse));

    // hamming_distance: number of differing positions.
    std::vector<Real> h1 {1, 2, 3};
    std::vector<Real> h2 {1, 2, 4};
    CHECK_EQUAL(std::size_t(1), boost::math::tools::hamming_distance(h1, h2));
    CHECK_EQUAL(std::size_t(0), boost::math::tools::hamming_distance(h1, h1));

    // total_variation: sum of consecutive absolute differences.
    std::vector<Real> flat {1, 1};
    CHECK_ULP_CLOSE(Real(0), boost::math::tools::total_variation(flat), 0);
    std::vector<Real> step {1, 2};
    CHECK_ULP_CLOSE(Real(1), boost::math::tools::total_variation(step), 1);
    std::vector<Real> ramp {0, 1, 2, 3, 4, 5};
    CHECK_ULP_CLOSE(Real(5), boost::math::tools::total_variation(ramp), 1);

    // l1_distance: distance from a vector to itself is zero; {1,2,3} to {1,1,1}
    // is 0 + 1 + 2 = 3.
    std::vector<Real> a {1, 2, 3};
    std::vector<Real> b {1, 1, 1};
    CHECK_ULP_CLOSE(Real(0), boost::math::tools::l1_distance(a, a), 0);
    CHECK_ULP_CLOSE(Real(3), boost::math::tools::l1_distance(a, b), 1);

    // l2_distance: a vector is at zero distance from itself.
    CHECK_ULP_CLOSE(Real(0), boost::math::tools::l2_distance(a, a), 0);

    // sup_distance: {1,1,1,1} to the zero vector is 1.
    std::vector<Real> f1 {1, 1, 1, 1};
    std::vector<Real> f0 {0, 0, 0, 0};
    CHECK_ULP_CLOSE(Real(0), boost::math::tools::sup_distance(f1, f1), 0);
    CHECK_ULP_CLOSE(Real(1), boost::math::tools::sup_distance(f1, f0), 1);

    // lp_distance: {1,0,0} to the zero vector is 1 for any p.
    std::vector<Real> p1 {1, 0, 0};
    std::vector<Real> p0 {0, 0, 0};
    CHECK_ULP_CLOSE(Real(0), boost::math::tools::lp_distance(p1, p1, 3), 0);
    CHECK_ULP_CLOSE(Real(1), boost::math::tools::lp_distance(p1, p0, 3), 2);
}

int main()
{
    test_norms<float>();
    test_norms<double>();

    return boost::math::test::report_errors();
}
