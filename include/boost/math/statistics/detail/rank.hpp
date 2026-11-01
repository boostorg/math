//  (C) Copyright Matt Borland 2022
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_STATISTICS_DETAIL_RANK_HPP
#define BOOST_MATH_STATISTICS_DETAIL_RANK_HPP

#ifndef BOOST_MATH_BUILD_MODULE
#include <cstdint>
#include <vector>
#include <numeric>
#include <utility>
#include <iterator>
#include <algorithm>
#endif
#include <boost/math/tools/config.hpp>

#ifdef BOOST_MATH_EXEC_COMPATIBLE
#ifndef BOOST_MATH_BUILD_MODULE
#include <execution>
#endif
#endif

#ifndef BOOST_MATH_EXEC_COMPATIBLE

BOOST_MATH_NAMESPACE_BEGIN namespace statistics { namespace detail {

template <typename ForwardIterator, typename T = typename std::iterator_traits<ForwardIterator>::value_type>
auto rank(ForwardIterator first, ForwardIterator last) -> std::vector<std::size_t>
{
    std::size_t elements = std::distance(first, last);

    std::vector<std::pair<T, std::size_t>> rank_vector(elements);
    std::size_t i = 0;
    while (first != last)
    {
        rank_vector[i] = std::make_pair(*first, i);
        ++i;
        ++first;
    }

    std::sort(rank_vector.begin(), rank_vector.end());

    // Tied values all get the largest rank of their run, i.e. result[i] + 1 = #{j : y_j <= y_i}.
    std::size_t rank = 0;
    std::vector<std::size_t> result(elements);
    for (i = elements; i-- > 0;)
    {
        if (i + 1 == elements || rank_vector[i].first != rank_vector[i + 1].first)
        {
            rank = i;
        }
        result[rank_vector[i].second] = rank;
    }

    return result;
}

template <typename Container>
inline auto rank(const Container& c) -> std::vector<std::size_t>
{
    return rank(std::begin(c), std::end(c));
}

}} BOOST_MATH_NAMESPACE_END

#else

BOOST_MATH_NAMESPACE_BEGIN namespace statistics::detail {

template <typename ExecutionPolicy, typename ForwardIterator, typename T = typename std::iterator_traits<ForwardIterator>::value_type>
auto rank(ExecutionPolicy&& exec, ForwardIterator first, ForwardIterator last)
{
    std::size_t elements = std::distance(first, last);

    std::vector<std::pair<T, std::size_t>> rank_vector(elements);
    std::size_t i = 0;
    while (first != last)
    {
        rank_vector[i] = std::make_pair(*first, i);
        ++i;
        ++first;
    }

    std::sort(exec, rank_vector.begin(), rank_vector.end());

    // Tied values all get the largest rank of their run, i.e. result[i] + 1 = #{j : y_j <= y_i}.
    std::size_t rank = 0;
    std::vector<std::size_t> result(elements);
    for (i = elements; i-- > 0;)
    {
        if (i + 1 == elements || rank_vector[i].first != rank_vector[i + 1].first)
        {
            rank = i;
        }
        result[rank_vector[i].second] = rank;
    }

    return result;
}

template <typename ExecutionPolicy, typename Container>
inline auto rank(ExecutionPolicy&& exec, const Container& c)
{
    return rank(exec, std::cbegin(c), std::cend(c));
}

template <typename ForwardIterator, typename T = typename std::iterator_traits<ForwardIterator>::value_type>
inline auto rank(ForwardIterator first, ForwardIterator last)
{
    return rank(std::execution::seq, first, last);
}

template <typename Container>
inline auto rank(const Container& c)
{
    return rank(std::execution::seq, std::cbegin(c), std::cend(c));
}

} BOOST_MATH_NAMESPACE_END

#endif // BOOST_MATH_EXEC_COMPATIBLE

#endif // BOOST_MATH_STATISTICS_DETAIL_RANK_HPP
