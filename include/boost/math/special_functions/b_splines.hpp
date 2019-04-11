//  (C) Copyright Nick Thompson 2019.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_SPECIAL_B_SPLINES_HPP
#define BOOST_MATH_SPECIAL_B_SPLINES_HPP

#include <boost/math/special_functions/factorials.hpp>

namespace boost { namespace math {


template<class Real, class Z = int>
class b_spline {

public:
    b_spline(Z  n) : n_{n}
    {
        static_assert(std::is_integral_v<Z>, "Z must be an integer type.");

        if (n < 1)
        {
            throw std::domain_error("B-splines are only defined for n >= 0.");
        }
        coeffs_.resize(n_ + 1);
        // There's faster ways to do this, but n is generally not large.
        for (Z v = 0; v < n_ + 1; ++v)
        {
            coeffs_[v] = Real(n_)/(boost::math::factorial<Real>(v)*boost::math::factorial<Real>(n-v));
            if (v & 1)
            {
                coeffs_[v] = -coeffs_[v];
            }
        }
    }
    // Schoenberg, Cardinal Spline Interpolation, Equation 1.4:
    Real operator()(Real x)
    {
        if (x <= 0 || x >= n_)
        {
            return Real(0);
        }

        if (x > Real(n_)/Real(2))
        {
            return this->operator()(Real(n_) - x);
        }

        Real sum = 0;
        using std::pow;
        for (Z v = 0; v < n_; ++v)
        {
            Real t = x - v;
            if (t <= 0)
            {
                return sum;
            }
            sum += coeffs_[v]*pow(t, n_ - 1);

        }
        return sum;
    }

    // Schoenberg, Cardinal Spline Interpolation, Equation 1.5:
    Real centered(Real x)
    {
        return this->operator()(x + Real(n_)/Real(2));
    }

private:
    std::vector<Real> coeffs_;
    Z n_;
};

}}
#endif
