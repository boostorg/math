//  Copyright Nick Thompson 2019.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_DISTRIBUTIONS_EMPIRICAL_CUMULATIVE_DISTRIBUTION_FUNCTION_HPP
#define BOOST_MATH_DISTRIBUTIONS_EMPIRICAL_CUMULATIVE_DISTRIBUTION_FUNCTION_HPP
#include <algorithm>

namespace boost { namespace math{

template<class RandomAccessContainer>
class empirical_cumulative_distribution_function {
    using Real = typename RandomAccessContainer::value_type;
public:
    empirical_cumulative_distribution_function(RandomAccessContainer && v)
    {
        static_assert(!std::is_integral_v<Real>, "Integer data not implemented.");
        m_v = std::move(v);
        if (!std::is_sorted(m_v.begin(), m_v.end())) {
            std::sort(m_v.begin(), m_v.end());
        }
    }

    auto operator()(Real x) const {
       if (x < m_v[0]) {
         return Real(0);
       }
       if (x >= m_v[m_v.size()-1]) {
         return Real(1);
       }
       auto it = std::upper_bound(m_v.begin(), m_v.end(), x);
       return static_cast<Real>(std::distance(m_v.begin(), it))/static_cast<Real>(m_v.size());
    }


private:
    RandomAccessContainer m_v;
};

}}
#endif
