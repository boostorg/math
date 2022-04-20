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
#include <utility>
#include <type_traits>
#include <boost/math/tools/assert.hpp>
#include <boost/math/tools/config.hpp>
#include <boost/math/statistics/detail/rank.hpp>

#ifdef BOOST_MATH_EXEC_COMPATIBLE
#include <execution>
#include <future>
#include <thread>
#endif

namespace boost { namespace math { namespace statistics {

namespace detail {

template <typename ReturnType, typename ForwardIterator>
ReturnType chatterjee_correlation_seq_impl(ForwardIterator u_begin, ForwardIterator u_end, ForwardIterator v_begin, ForwardIterator v_end)
{
    using std::abs;
    
    BOOST_MATH_ASSERT_MSG(std::is_sorted(u_begin, u_end), "The x values must be sorted in order to use this functionality");

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

    // If the result is 1 then Y is constant and all the elements must be ties
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

#ifdef BOOST_MATH_EXEC_COMPATIBLE

namespace boost::math::statistics {

namespace detail {

template <typename ReturnType, typename ExecutionPolicy, typename ForwardIterator>
ReturnType chatterjee_correlation_par_impl(ExecutionPolicy&& exec, ForwardIterator u_begin, ForwardIterator u_end,
                                                                   ForwardIterator v_begin, ForwardIterator v_end)
{
    using std::abs;
    BOOST_MATH_ASSERT_MSG(std::is_sorted(std::forward<ExecutionPolicy>(exec), u_begin, u_end), "The x values must be sorted in order to use this functionality");

    const auto rank_vector = rank(std::forward<ExecutionPolicy>(exec), v_begin, v_end);

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

    // If the result is 1 then Y is constant and all the elements must be ties
    if (abs(result - static_cast<ReturnType>(1)) < std::numeric_limits<ReturnType>::epsilon())
    {
        return std::numeric_limits<ReturnType>::quiet_NaN();
    }

    return result;
}

} // Namespace detail

template <typename ExecutionPolicy, typename Container, typename Real = typename Container::value_type,
          typename ReturnType = std::conditional_t<std::is_integral_v<Real>, double, Real>>
inline ReturnType chatterjee_correlation(ExecutionPolicy&& exec, const Container& u, const Container& v)
{
    if constexpr (std::is_same_v<std::remove_reference_t<decltype(exec)>, decltype(std::execution::seq)>)
    {
        return detail::chatterjee_correlation_seq_impl<ReturnType>(std::cbegin(u), std::cend(u),
                                                                   std::cbegin(v), std::cend(v));
    }
    else
    {
        return detail::chatterjee_correlation_par_impl<ReturnType>(std::forward<ExecutionPolicy>(exec),
                                                                   std::cbegin(u), std::cend(u),
                                                                   std::cbegin(v), std::cend(v));
    }
}

} // Namespace boost::math::statistics

#endif

#endif // BOOST_MATH_STATISTICS_CHATTERJEE_CORRELATION_HPP
