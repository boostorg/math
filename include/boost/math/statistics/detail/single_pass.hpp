//  (C) Copyright Nick Thompson 2018
//  (C) Copyright Matt Borland 2020
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_STATISTICS_UNIVARIATE_STATISTICS_DETAIL_SINGLE_PASS_HPP
#define BOOST_MATH_STATISTICS_UNIVARIATE_STATISTICS_DETAIL_SINGLE_PASS_HPP

#include <tuple>
#include <iterator>
#include <atomic>
#include <thread>
#include <type_traits>
#include <future>

namespace boost::math::statistics::detail
{
// Higham, Accuracy and Stability, equation 1.6a and 1.6b:
// Calculates Mean, M2, and variance
template<typename ForwardIterator>
auto variance_integeral_impl(ForwardIterator first, ForwardIterator last) -> std::tuple<double, double, double>
{
    double M = *first;
    double Q = 0;
    double k = 2;
    double M2 = 0;
    for (auto it = std::next(first); it != last; ++it)
    {
        double delta_1 = *it - M;
        Q = Q + ((k-1)*delta_1*delta_1)/k;
        M = M + delta_1/k;
        k += 1;
        double delta_2 = *it - M;
        M2 = M2 + (delta_1 * delta_2);
    }
    return std::make_tuple(M, M2, Q/(k-1));
}

template<typename ForwardIterator>
auto variance_real_impl(ForwardIterator first, ForwardIterator last) ->
    std::tuple<typename std::iterator_traits<ForwardIterator>::value_type,
               typename std::iterator_traits<ForwardIterator>::value_type,
               typename std::iterator_traits<ForwardIterator>::value_type>
{
    using Real = typename std::iterator_traits<ForwardIterator>::value_type;

    Real M = *first;
    Real Q = 0;
    Real k = 2;
    Real M2 = 0;
    for (auto it = std::next(first); it != last; ++it)
    {
        Real tmp = (*it - M)/k;
        Real delta_1 = *it - M;
        Q += k*(k-1)*tmp*tmp;
        M += tmp;
        k += 1;
        Real delta_2 = *it - M;
        M2 = M2 + (delta_1 * delta_2);
    }
    return std::make_tuple(M, M2, Q/(k-1));
}

// Global thread counter required for recursive calls. Reset before each parallel method is called
static std::atomic<unsigned> thread_counter {1};

// http://i.stanford.edu/pub/cstr/reports/cs/tr/79/773/CS-TR-79-773.pdf
// http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.214.8508&rep=rep1&type=pdf
template<typename ReturnType, typename ForwardIterator>
ReturnType parallel_variance_impl(ForwardIterator first, ForwardIterator last)
{
    using Real = typename std::iterator_traits<ForwardIterator>::value_type;
    static unsigned num_threads {std::thread::hardware_concurrency()};

    const auto elements {std::distance(first, last)};
    const auto range_a {std::floor(elements / 2)};
    const auto range_b {elements - range_a};

    thread_counter.fetch_add(2);
    auto future_a {std::async(std::launch::async, [first, range_a, &thread_counter, num_threads]() -> ReturnType
    {
        if constexpr (std::is_integral_v<Real>)
        {
            if(thread_counter + 2 <= num_threads && range_a > 10)
                return parallel_variance_impl<ReturnType>(first, std::next(first, range_a));
            else
                return variance_integeral_impl(first, std::next(first, range_a));
        }
        else
        {
            if(thread_counter + 2 <= num_threads && range_a > 10)
                return parallel_variance_impl<ReturnType>(first, std::next(first, range_a));
            else
                return variance_real_impl(first, std::next(first, range_a));
        }
    })};
    auto future_b {std::async(std::launch::async, [first, last, range_a, &thread_counter, num_threads]() -> ReturnType
    {
        if constexpr (std::is_integral_v<Real>)
        {
            if(thread_counter + 2 <= num_threads && range_a > 10)
                return parallel_variance_impl<ReturnType>(std::next(first, range_a), last);
            else
                return variance_integeral_impl(std::next(first, range_a), last);
        }
        else
        {
            if(thread_counter + 2 <= num_threads && range_a > 10)
                return parallel_variance_impl<ReturnType>(std::next(first, range_a), last);
            else
                return variance_real_impl(std::next(first, range_a), last);
        }
    })};

    const auto results_a {future_a.get()};
    const auto results_b {future_b.get()};

    const auto mean_a = std::get<0>(results_a);
    const auto M2_a = std::get<1>(results_a);
    const auto mean_b = std::get<0>(results_b);
    const auto M2_b = std::get<1>(results_b);

    const auto n_ab = elements;
    const auto delta = mean_b - mean_a;
    const auto mean_ab = (range_a * mean_a + range_b * mean_b) / n_ab;
    const auto M2_ab = M2_a + M2_b + delta * delta * (range_a * range_b / n_ab);

    return std::make_tuple(mean_ab, M2_ab, n_ab);
}
}

#endif // BOOST_MATH_STATISTICS_UNIVARIATE_STATISTICS_DETAIL_SINGLE_PASS_HPP
