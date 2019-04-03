// Copyright Nick Thompson, 2019
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)
#ifndef BOOST_MATH_INTERPOLATORS_WHITAKKER_SHANNON_HPP
#define BOOST_MATH_INTERPOLATORS_WHITAKKER_SHANNON_HPP
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/sin_pi.hpp>

namespace boost { namespace math { namespace interpolators {

template<class RandomAccessContainer>
class whittaker_shannon {
public:

    using Real = typename RandomAccessContainer::value_type;
    whittaker_shannon(RandomAccessContainer&& y, Real const & t0, Real const & h) : m_y{std::move(y)}, m_t0{t0}, m_h{h}
    {
    }

    Real operator()(Real t) const {
        using boost::math::constants::pi;
        using std::sin;
        Real y = 0;
        Real x = (t - m_t0)/m_h;
        
        for (size_t i = 0; i < m_y.size(); ++i)
        {
            Real denom = (x - i);
            if (denom == 0) {
                return m_y[i];
            }
            if (i & 1) {
                y -= m_y[i]/denom;
            }
            else {
                y += m_y[i]/denom;
            }
        }
        return y*sin_pi(x)/pi<Real>();
    }


    Real operator[](size_t i) const {
        return m_y[i];
    }

    RandomAccessContainer&& return_data() {
        return std::move(m_y);
    }


private:
    RandomAccessContainer m_y;
    Real m_t0;
    Real m_h;
};
}}}
#endif
