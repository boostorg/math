//  (C) Copyright Matt Borland 2022.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_STATISTICS_CHATTERJEE_CORRELATION_HPP
#define BOOST_MATH_STATISTICS_CHATTERJEE_CORRELATION_HPP

#include <cstdint>
#include <cmath>
#include <algorithm>
#include <iterator>
#include <vector>
#include <limits>
#include <type_traits>
#include <boost/math/tools/assert.hpp>
#include <boost/math/tools/config.hpp>
#include <boost/math/statistics/detail/rank.hpp>

namespace boost { namespace math { namespace statistics {

namespace detail {

template <typename ReturnType, typename ForwardIterator>
ReturnType chatterjee_correlation_seq_impl(ForwardIterator u_begin, ForwardIterator u_end, ForwardIterator v_begin, ForwardIterator v_end)
{
    using std::abs;
    
    BOOST_MATH_ASSERT_MSG(std::is_sorted(u_begin, u_end), "The x values must be sorted in order to use this funtionality");

    const std::vector<std::size_t> rank_vector = rank(v_begin, v_end);

    std::size_t sum = 0;
    for (std::size_t i = 1; i < rank_vector.size(); ++i)
    {
        // Equivalent to abs(rank_vector[i] - rank_vector[i-1]) but avoids unsigned underflow in the intermediate step
        if (rank_vector[i] > rank_vector[i-1])
        {
            sum += rank_vector[i] - rank_vector[i-1];
        }
        else
        {
            sum += rank_vector[i-1] - rank_vector[i];
        }
    }

    ReturnType result = static_cast<ReturnType>(1) - (static_cast<ReturnType>(3 * sum) / static_cast<ReturnType>(rank_vector.size() * rank_vector.size() - 1));

    // If the result is 1 then Y is constant and all of the elements must be ties
    if (abs(result - static_cast<ReturnType>(1)) < std::numeric_limits<ReturnType>::epsilon())
    {
        return std::numeric_limits<ReturnType>::quiet_NaN();
    }

    return result;
}

} // Namespace detail

template <typename Container, typename Real = typename Container::value_type, 
          typename ReturnType = typename std::conditional<std::is_integral<Real>::value, double, Real>::type>
inline ReturnType chatterjee_correlation(const Container& u, const Container& v)
{
    return detail::chatterjee_correlation_seq_impl<ReturnType>(std::begin(u), std::end(u), std::begin(v), std::end(v));
}

}}} // Namespace boost::math::statistics

#endif // BOOST_MATH_STATISTICS_CHATTERJEE_CORRELATION_HPP
