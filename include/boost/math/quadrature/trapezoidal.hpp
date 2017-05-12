/*
 * Copyright Nick Thompson, 2017
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 *
 * Use the adaptive trapezoidal rule to estimate the integral of periodic functions over a period,
 * or to integrate a function whose derivative vanishes at the endpoints.
 *
 * If your function does not satisfy these conditions, and instead is simply continuous and bounded
 * over the whole interval, then this routine will still converge, albeit slowly. However, there
 * are much more efficient methods in this case, including Romberg, Simpson, and double exponential quadrature.
 */

#ifndef BOOST_MATH_QUADRATURE_TRAPEZOIDAL_HPP
#define BOOST_MATH_QUADRATURE_TRAPEZOIDAL_HPP

#include <cmath>
#include <limits>
#include <stdexcept>
#include <boost/math/constants/constants.hpp>

namespace boost{ namespace math{ namespace quadrature {

template<class F, class Real>
Real trapezoidal(F f, Real a, Real b, Real tol = sqrt(std::numeric_limits<Real>::epsilon()), size_t max_refinements = 10, Real* error_estimate = nullptr, Real* L1 = nullptr)
{
    using std::abs;
    using std::isfinite;
    using boost::math::constants::half;
    if(a >= b)
    {
        throw std::domain_error("a < b for integration over the region [a, b] is required.\n");
    }
    if (!isfinite(a))
    {
        throw std::domain_error("Left endpoint of integration must be finite for adaptive trapezoidal integration.\n");
    }
    if (!isfinite(b))
    {
        throw std::domain_error("Right endpoint of integration must be finite for adaptive trapedzoidal integration.\n");
    }

    Real ya = f(a);
    Real yb = f(b);
    Real h = (b - a)*half<Real>();
    Real I0 = (ya + yb)*h;
    Real IL0 = (abs(ya) + abs(yb))*h;

    Real yh = f(a + h);
    Real I1 = half<Real>()*I0 + yh*h;
    Real IL1 = half<Real>()*IL0 + abs(yh)*h;

    // The recursion is:
    // I_k = 1/2 I_{k-1} + 1/2^k \sum_{j=1; j odd, j < 2^k} f(a + j(b-a)/2^k)
    size_t k = 2;
    // We want to go through at least 4 levels so we have sampled the function at least 10 times.
    // Otherwise, we could terminate prematurely and miss essential features.
    // This is of course possible anyway, but 10 samples seems to be a reasonable compromise.
    Real error = abs(I0 - I1);
    while (k < 4 || (k < max_refinements && error > tol*IL1) )
    {
        I0 = I1;
        IL0 = IL1;

        I1 = half<Real>()*I0;
        IL1 = half<Real>()*IL0;
        size_t p = 1 << k;
        h *= half<Real>();
        Real sum = 0;
        Real absum = 0;

        for(size_t j = 1; j < p; j += 2)
        {
            Real y = f(a + j*h);
            sum += y;
            absum += abs(y);
        }

        I1 += sum*h;
        IL1 += absum*h;
        ++k;
        error = abs(I0 - I1);
    }

    if (error_estimate)
    {
        *error_estimate = error;
    }

    if (L1)
    {
        *L1 = IL1;
    }

    return I1;
}

}}}
#endif
