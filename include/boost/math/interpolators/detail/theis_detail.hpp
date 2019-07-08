// Copyright Nick Thompson, 2019
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)
#ifndef BOOST_MATH_INTERPOLATORS_THEIS_DETAIL_HPP
#define BOOST_MATH_INTERPOLATORS_THEIS_DETAIL_HPP
#include <cmath>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/sin_pi.hpp>

namespace boost { namespace math { namespace interpolators { namespace detail {

template<class RandomAccessContainer>
class theis_detail {
public:

    using Real = typename RandomAccessContainer::value_type;
    theis_detail(RandomAccessContainer&& y, Real const & t0, Real const & h) : m_y{std::move(y)}, m_t0{t0}, m_h{h}
    {}

    inline Real operator()(Real t) const {
      using boost::math::constants::pi;
      using std::isfinite;
      using std::floor;
      Real y = 0;
      Real x = (t - m_t0)/m_h;
      if (x == floor(x)) {
        size_t i = static_cast<size_t>(floor(x));
        return m_y[i];
      }

      Real z = x;
      auto it = m_y.begin();
      auto end = m_y.end();
      while(it != end)
      {
          y += *it++/(z*z);
          z -= 1;
      }
      Real coeff = boost::math::sin_pi(x)/pi<Real>();
      return y*coeff*coeff;
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
}}}}
#endif
