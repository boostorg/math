//  (C) Copyright Matt Borland 2022
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_STATISTICS_DETAIL_RANK_HPP
#define BOOST_MATH_STATISTICS_DETAIL_RANK_HPP

#include <cstdint>
#include <vector>
#include <numeric>
#include <utility>
#include <iterator>
#include <boost/math/tools/config.hpp>

#ifndef BOOST_MATH_EXEC_COMPATIBLE

#include <execution>

namespace boost { namespace math { namespace statistics { namespace detail {

template <typename ForwardIterator, typename T = typename std::iterator_traits<ForwardIterator>::value_type>
auto rank(ForwardIterator first, ForwardIterator last) -> std::vector<std::size_t>
{
    const auto elements = std::distance(first, last);

    std::vector<std::pair<T, std::size_t>> rank_vector(elements);
    std::size_t i = 0;
    while (first != last)
    {
        rank_vector[i] = std::make_pair(*first, i);
        ++i;
        ++first;
    }

    std::sort(rank_vector.begin(), rank_vector.end());

    std::pair<T, std::size_t> rank;
    std::vector<std::size_t> result(elements);
    for (i = 0; i < elements; ++i)
    {
        if (rank_vector[i].first != rank.first)
        {
            rank = std::make_pair(rank_vector[i].first, i);
        }
        result[rank_vector[i].second] = rank.second;
    }

    return result;
}

template <typename Container>
inline auto rank(const Container& c) -> std::vector<std::size_t>
{
    return rank(std::begin(c), std::end(c));
}

}}}} // Namespaces

#else

namespace boost::math::statistics::detail {

template <typename ExecutionPolicy, typename ForwardIterator, typename T = typename std::iterator_traits<ForwardIterator>::value_type>
auto rank(ExecutionPolicy&& exec, ForwardIterator first, ForwardIterator last)
{
    const auto elements = std::distance(first, last);

    std::vector<std::pair<T, std::size_t>> rank_vector(elements);
    std::size_t i = 0;
    while (first != last)
    {
        rank_vector[i] = std::make_pair(*first, i);
        ++i;
        ++first;
    }

    std::sort(exec, rank_vector.begin(), rank_vector.end());

    std::pair<T, std::size_t> rank;
    std::vector<std::size_t> result(elements);
    for (i = 0; i < elements; ++i)
    {
        if (rank_vector[i].first != rank.first)
        {
            rank = std::make_pair(rank_vector[i].first, i);
        }
        result[rank_vector[i].second] = rank.second;
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

} // Namespaces

#endif // BOOST_MATH_EXEC_COMPATIBLE

#endif // BOOST_MATH_STATISTICS_DETAIL_RANK_HPP
