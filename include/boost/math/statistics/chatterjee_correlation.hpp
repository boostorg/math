//  (C) Copyright Matt Borland 2022.
//  (C) Copyright Oleksandr Kornijcuk 2026.
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

template <typename BDIter>
std::size_t chatterjee_transform(BDIter begin, BDIter end)
{
    std::size_t sum = 0;

    while(++begin != end)
    {
        if(*begin > *std::prev(begin))
        {
            sum += *begin - *std::prev(begin);
        }
        else
        {
            sum += *std::prev(begin) - *begin;
        }
    }

    return sum;
}

template <typename ReturnType, typename ForwardIterator>
ReturnType chatterjee_correlation_seq_impl(ForwardIterator u_begin, ForwardIterator u_end, ForwardIterator v_begin, ForwardIterator v_end)
{
    using std::abs;

    BOOST_MATH_ASSERT_MSG(std::is_sorted(u_begin, u_end), "The x values must be sorted in order to use this functionality");

    const std::vector<std::size_t> rank_vector = rank(v_begin, v_end);

    std::size_t sum = chatterjee_transform(rank_vector.begin(), rank_vector.end());

    ReturnType result = static_cast<ReturnType>(1) - (static_cast<ReturnType>(3 * sum) / static_cast<ReturnType>(rank_vector.size() * rank_vector.size() - 1));

    // If the result is 1 then Y is constant and all the elements must be ties
    if (abs(result - static_cast<ReturnType>(1)) < std::numeric_limits<ReturnType>::epsilon())
    {
        return std::numeric_limits<ReturnType>::quiet_NaN();
    }

    return result;
}

// Generalised (M right-nearest-neighbour) Chatterjee correlation of Lin and Han (2021),
// "On boosting the power of Chatterjee's rank correlation", https://arxiv.org/abs/2108.06828.
//
// For a sample sorted by X, the inner sum of equation (2.4) is
//
//     S = sum_{i} sum_{m=1}^{M} min(R_i, R_{j_m(i)}),
//
// where R_i is the (1-based) rank of Y_i and j_m(i) is the m-th right nearest neighbour of X_i.
// After sorting by X the m-th right neighbour of element i is element i + m when it exists
// (i + m < n); otherwise j_m(i) = i (the sentinel of the paper), contributing min(R_i, R_i) = R_i.
//
// Note: rank() returns 0-based ranks, whereas the paper's formula uses 1-based ranks. The +1 offset
// cancels in the M = 1 statistic above (which uses |R_{i+1} - R_i|) but does NOT cancel under
// min(.,.), so it is applied explicitly here.
inline std::uint64_t chatterjee_mnn_min_sum(const std::vector<std::size_t>& rank_vector,
                                            std::size_t lo, std::size_t hi, std::size_t M)
{
    const std::size_t n = rank_vector.size();
    std::uint64_t sum = 0;

    for (std::size_t i = lo; i < hi; ++i)
    {
        const std::size_t r_i = rank_vector[i] + 1;

        // Number of right neighbours of i that actually exist within [0, n).
        const std::size_t existing = (i + M < n) ? M : (n - 1 - i);

        for (std::size_t m = 1; m <= existing; ++m)
        {
            const std::size_t r_j = rank_vector[i + m] + 1;
            sum += static_cast<std::uint64_t>((std::min)(r_i, r_j));
        }

        // The remaining (M - existing) neighbours are sentinels j_m(i) = i, each contributing r_i.
        sum += static_cast<std::uint64_t>(M - existing) * static_cast<std::uint64_t>(r_i);
    }

    return sum;
}

template <typename ReturnType>
ReturnType chatterjee_mnn_result(std::uint64_t sum, std::size_t n, std::size_t M)
{
    const ReturnType n_r = static_cast<ReturnType>(n);
    const ReturnType m_r = static_cast<ReturnType>(M);

    // Denominator (n + 1) * [nM + M(M + 1) / 4] of equation (2.4); chosen so that E[xi_{n,M}] = 0 under H0.
    const ReturnType denominator = (n_r + static_cast<ReturnType>(1)) *
        (n_r * m_r + m_r * (m_r + static_cast<ReturnType>(1)) / static_cast<ReturnType>(4));

    return static_cast<ReturnType>(-2) + static_cast<ReturnType>(6) * static_cast<ReturnType>(sum) / denominator;
}

template <typename ReturnType, typename ForwardIterator>
ReturnType chatterjee_correlation_mnn_seq_impl(ForwardIterator u_begin, ForwardIterator u_end,
                                               ForwardIterator v_begin, ForwardIterator v_end, std::size_t M)
{
    BOOST_MATH_ASSERT_MSG(std::is_sorted(u_begin, u_end), "The x values must be sorted in order to use this functionality");

    const std::size_t n = static_cast<std::size_t>(std::distance(v_begin, v_end));

    BOOST_MATH_ASSERT_MSG(M >= 1, "The number of right nearest neighbours M must be at least 1");
    BOOST_MATH_ASSERT_MSG(M <= n, "The number of right nearest neighbours M must not exceed the number of samples");

    // If Y is constant the statistic is undefined (Y must be non-constant); return a quiet NaN.
    // This is checked on the input directly: rank() collapses tied values, so a constant Y would
    // otherwise reduce the effective sample size to one.
    if (v_begin != v_end)
    {
        using value_type = typename std::iterator_traits<ForwardIterator>::value_type;
        const value_type first_value = *v_begin;
        if (std::all_of(std::next(v_begin), v_end,
                        [&first_value](const value_type& value) { return value == first_value; }))
        {
            return std::numeric_limits<ReturnType>::quiet_NaN();
        }
    }

    const std::vector<std::size_t> rank_vector = rank(v_begin, v_end);

    const std::uint64_t sum = chatterjee_mnn_min_sum(rank_vector, 0, n, M);

    return chatterjee_mnn_result<ReturnType>(sum, n, M);
}

} // Namespace detail

template <typename Container, typename Real = typename Container::value_type,
          typename ReturnType = typename std::conditional<std::is_integral<Real>::value, double, Real>::type>
inline ReturnType chatterjee_correlation(const Container& u, const Container& v)
{
    return detail::chatterjee_correlation_seq_impl<ReturnType>(std::begin(u), std::end(u), std::begin(v), std::end(v));
}

template <typename Container, typename Real = typename Container::value_type,
          typename ReturnType = typename std::conditional<std::is_integral<Real>::value, double, Real>::type>
inline ReturnType chatterjee_correlation_mnn(const Container& u, const Container& v, std::size_t M)
{
    return detail::chatterjee_correlation_mnn_seq_impl<ReturnType>(std::begin(u), std::end(u), std::begin(v), std::end(v), M);
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

    auto rank_vector = rank(std::forward<ExecutionPolicy>(exec), v_begin, v_end);

    const auto num_threads = std::thread::hardware_concurrency() == 0 ? 2u : std::thread::hardware_concurrency();
    std::vector<std::future<std::size_t>> future_manager {};
    const auto elements_per_thread = std::ceil(static_cast<double>(rank_vector.size()) / num_threads);

    auto it = rank_vector.begin();
    auto end = rank_vector.end();
    for(std::size_t i {}; i < num_threads - 1; ++i)
    {
        future_manager.emplace_back(std::async(std::launch::async | std::launch::deferred, [it, elements_per_thread]() -> std::size_t
        {
            return chatterjee_transform(it, std::next(it, elements_per_thread));
        }));
        it = std::next(it, elements_per_thread - 1);
    }

    future_manager.emplace_back(std::async(std::launch::async | std::launch::deferred, [it, end]() -> std::size_t
    {
        return chatterjee_transform(it, end);
    }));

    std::size_t sum {};
    for(std::size_t i {}; i < future_manager.size(); ++i)
    {
        sum += future_manager[i].get();
    }

    ReturnType result = static_cast<ReturnType>(1) - (static_cast<ReturnType>(3 * sum) / static_cast<ReturnType>(rank_vector.size() * rank_vector.size() - 1));

    // If the result is 1 then Y is constant and all the elements must be ties
    if (abs(result - static_cast<ReturnType>(1)) < std::numeric_limits<ReturnType>::epsilon())
    {
        return std::numeric_limits<ReturnType>::quiet_NaN();
    }

    return result;
}

template <typename ReturnType, typename ExecutionPolicy, typename ForwardIterator>
ReturnType chatterjee_correlation_mnn_par_impl(ExecutionPolicy&& exec, ForwardIterator u_begin, ForwardIterator u_end,
                                                                       ForwardIterator v_begin, ForwardIterator v_end, std::size_t M)
{
    BOOST_MATH_ASSERT_MSG(std::is_sorted(std::forward<ExecutionPolicy>(exec), u_begin, u_end), "The x values must be sorted in order to use this functionality");

    const std::size_t n = static_cast<std::size_t>(std::distance(v_begin, v_end));

    BOOST_MATH_ASSERT_MSG(M >= 1, "The number of right nearest neighbours M must be at least 1");
    BOOST_MATH_ASSERT_MSG(M <= n, "The number of right nearest neighbours M must not exceed the number of samples");

    // If Y is constant the statistic is undefined (Y must be non-constant); return a quiet NaN.
    // This is checked on the input directly: rank() collapses tied values, so a constant Y would
    // otherwise reduce the effective sample size to one.
    if (v_begin != v_end)
    {
        using value_type = typename std::iterator_traits<ForwardIterator>::value_type;
        const value_type first_value = *v_begin;
        if (std::all_of(std::next(v_begin), v_end,
                        [&first_value](const value_type& value) { return value == first_value; }))
        {
            return std::numeric_limits<ReturnType>::quiet_NaN();
        }
    }

    const auto rank_vector = rank(std::forward<ExecutionPolicy>(exec), v_begin, v_end);

    const auto num_threads = std::thread::hardware_concurrency() == 0 ? 2u : std::thread::hardware_concurrency();

    // The i-loop is embarrassingly parallel: each task sums over a disjoint contiguous range of
    // indices [lo, hi) and reads (read-only) rank_vector entries up to i + M, which may fall in a
    // neighbouring task's range. Because there are no writes, the partition needs no overlap; this
    // differs from the M = 1 parallel path above, which splits the data array (and overlaps by one
    // element) for the difference-based transform.
    std::vector<std::future<std::uint64_t>> future_manager {};
    future_manager.reserve(num_threads);

    const std::size_t chunk = static_cast<std::size_t>(std::ceil(static_cast<double>(n) / num_threads));

    std::size_t lo = 0;
    while (lo < n)
    {
        const std::size_t hi = (std::min)(lo + chunk, n);
        future_manager.emplace_back(std::async(std::launch::async | std::launch::deferred,
            [&rank_vector, lo, hi, M]() -> std::uint64_t
            {
                return chatterjee_mnn_min_sum(rank_vector, lo, hi, M);
            }));
        lo = hi;
    }

    std::uint64_t sum = 0;
    for (std::size_t i {}; i < future_manager.size(); ++i)
    {
        sum += future_manager[i].get();
    }

    return chatterjee_mnn_result<ReturnType>(sum, n, M);
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

template <typename ExecutionPolicy, typename Container, typename Real = typename Container::value_type,
          typename ReturnType = std::conditional_t<std::is_integral_v<Real>, double, Real>>
inline ReturnType chatterjee_correlation_mnn(ExecutionPolicy&& exec, const Container& u, const Container& v, std::size_t M)
{
    if constexpr (std::is_same_v<std::remove_reference_t<decltype(exec)>, decltype(std::execution::seq)>)
    {
        return detail::chatterjee_correlation_mnn_seq_impl<ReturnType>(std::cbegin(u), std::cend(u),
                                                                       std::cbegin(v), std::cend(v), M);
    }
    else
    {
        return detail::chatterjee_correlation_mnn_par_impl<ReturnType>(std::forward<ExecutionPolicy>(exec),
                                                                       std::cbegin(u), std::cend(u),
                                                                       std::cbegin(v), std::cend(v), M);
    }
}

} // Namespace boost::math::statistics

#endif

#endif // BOOST_MATH_STATISTICS_CHATTERJEE_CORRELATION_HPP
