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

    auto operator()(Real s, Real t)->decltype(std::declval<F>()(std::declval<Real>())) const
    {
        using std::sqrt;
        using std::abs;
        using boost::math::quadrature::trapezoidal;

        // -p + 1 < (x-t)/s < p so if s > 0 then s(1-p) + t < x < sp + t
        // if s < 0, then
        // -sp + s >= x-t >= sp <=> sp + t < x < s(1-p) + t
        Real a = -s*p + s + t;
        Real b = s*p + t;
        if (s < 0)
        {
            std::swap(a, b);
        }
        if (s == 0)
        {
            return std::numeric_limits<Real>::quiet_NaN();
        }

        return trapezoidal(f_, a, b, tol_, max_refinements_)/sqrt(abs(s));
    }

private:
    F f_;
    boost::math::daubechies_wavelet<Real, p> psi_;
    Real tol_;
    int max_refinements_;
};


}
#endif