// Copyright Nick Thompson, 2017
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

// This computes the Catmull-Rom spline from a list of points.

#ifndef BOOST_MATH_INTERPOLATORS_WHITAKKER_SHANNON_HPP
#define BOOST_MATH_INTERPOLATORS_WHITAKKER_SHANNON_HPP
#include <boost/math/special_functions/sinc.hpp>
#include <boost/math/constants/constants.hpp>

namespace boost::math::interpolators {

template<class RandomAccessContainer>
class whittaker_shannon {
public:

    using Real = typename RandomAccessContainer::value_type;
    whittaker_shannon(RandomAccessContainer&& y, Real t0, Real h) : m_y{std::move(y)}, m_t0{t0}, m_inv_h{1/h}
    {
    }

    Real operator()(Real t) const {
        using boost::math::constants::pi;
        Real y = 0;
        Real x = (t - m_t0)*m_inv_h;
        for (size_t i = 0; i < m_y.size(); ++i)
        {
            Real arg = pi<Real>()*(x - i);
            // We can probably speed this up using the recurrence (5.4.6) in Numerical Recipes, Third Edition.
            y += m_y[i]*boost::math::sinc_pi(arg);
        }
        return y;
    }


    Real operator[](size_t i) const {
        return m_y[i];
    }

private:
    RandomAccessContainer m_y;
    Real m_t0;
    Real m_inv_h;
};
}
#endif
