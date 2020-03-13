/*
 * Copyright Nick Thompson, 2020
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef BOOST_MATH_QUADRATURE_WAVELET_TRANSFORMS_HPP
#define BOOST_MATH_QUADRATURE_WAVELET_TRANSFORMS_HPP
#include <boost/math/special_functions/daubechies_wavelet.hpp>
#include <boost/math/quadrature/trapezoidal.hpp>

namespace boost::math::quadrature {

template<class F, typename Real, int p>
class daubechies_wavelet_transform
{
public:
    daubechies_wavelet_transform(F f, int grid_refinements = -1, Real tol = boost::math::tools::root_epsilon<Real>(),
    int max_refinements = 12) : f_{f}, psi_(grid_refinements), tol_{tol}, max_refinements_{max_refinements}
    {}

    Real operator()(Real s, Real t) const
    {
        using std::sqrt;
        using boost::math::quadrature::trapezoidal;
        if (s == 0)
        {
            return std::numeric_limits<Real>::quiet_NaN();
        }

        // -p + 1 < (x-t)/s < p so if s > 0 then s(1-p) + t < x < sp + t
        // if s < 0, then
        // -sp + s >= x-t >= sp <=> sp + t < x < s(1-p) + t
        Real a = -s*p + s + t;
        Real b = s*p + t;
        if (s < 0)
        {
            std::swap(a, b);
        }
        Real Q = trapezoidal(f_, a, b, tol_, max_refinements_);
        return Q/sqrt(s);
    }

private:
    F f_;
    boost::math::daubechies_wavelet<Real, p> psi_;
    Real tol_;
    int max_refinements_;
};


}
#endif