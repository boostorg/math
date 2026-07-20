/*
 * Copyright Nick Thompson, 2020
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */
//
// Module build check for boost::math::quadrature::daubechies_wavelet_transform.
// Constructing the transform builds a daubechies_wavelet interpolator (with
// std::async/std::future, shared_ptr grids and the interpolator tuple) and then
// integrates against it with the trapezoidal rule, so building this against the
// boost.math module exercises that machinery from a module consumer.

#ifndef BOOST_MATH_BUILD_MODULE
#include <boost/math/quadrature/wavelet_transforms.hpp>
#else
import boost.math;
#endif

#ifndef BOOST_MATH_BUILD_MODULE
#include <cmath>
#include <limits>
#endif
#include "math_unit_test.hpp"

template <typename Real, int p>
void test_wavelet_transform()
{
    std::cout << "Testing wavelet transform of " << p
              << " vanishing moment Daubechies wavelet\n";

    // A small explicit grid keeps this a fast, low-memory build check: the
    // default (-1) grid for double is 2^21 points and needs on the order of a
    // gigabyte. At 10 refinements the assertions below still pass with wide
    // margin (see comments) so this remains a faithful, known-good check.
    auto psi = boost::math::daubechies_wavelet<Real, p>(10);

    // f(x) = exp(-|x|) is even, so W[f](s,0) == W[f](-s,0): the trapezoidal
    // integrand f(s*u)*psi(u) takes identical values for +/-s. The transforms
    // are in fact bit-identical, so the 12 ULP bound (from wavelet_transform_test)
    // is met trivially.
    auto f = [](Real x) { return std::exp(-std::abs(x)); };
    auto Wf = boost::math::quadrature::daubechies_wavelet_transform(f, psi);
    for (Real s = Real(0.5); s < Real(5); s += Real(0.5))
    {
        Real w1 = Wf(s, Real(0));
        Real w2 = Wf(-s, Real(0));
        CHECK_ULP_CLOSE(w1, w2, 12);
    }

    // W[psi](1,0) = integral of psi*psi = ||psi||_2^2 = 1, since Daubechies
    // wavelets are L2-normalized. Expected value 1 and the 2*sqrt(eps) tolerance
    // are taken from wavelet_transform_test.cpp; the check holds for p > 5 where
    // the quadrature sum converges rapidly. Observed error is ~4e-11 (double)
    // and ~2e-7 (float), well inside tolerance.
    auto Wpsi = boost::math::quadrature::daubechies_wavelet_transform(psi, psi);
    CHECK_MOLLIFIED_CLOSE(Real(1), Wpsi(Real(1), Real(0)),
                          Real(2*std::sqrt(std::numeric_limits<Real>::epsilon())));
}

int main()
{
    test_wavelet_transform<float, 8>();
    test_wavelet_transform<double, 8>();

    return boost::math::test::report_errors();
}
