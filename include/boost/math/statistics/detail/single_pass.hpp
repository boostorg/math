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
#include <memory>
#include <thread>
#include <type_traits>
#include <future>
#include <cmath>
#include <algorithm>
#include <iostream>

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
    auto future_a {std::async(std::launch::async, [first, range_a]() -> ReturnType
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
    auto future_b {std::async(std::launch::async, [first, last, range_a]() -> ReturnType
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

// https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Higher-order_statistics
// Use first four moments of disjoint sets for policies with skewness and kurtosis
template<typename ForwardIterator>
auto first_four_moments_integral_impl(ForwardIterator first, ForwardIterator last)
{
    double M1 = *first;
    double M2 = 0;
    double M3 = 0;
    double M4 = 0;
    double n = 2;
    for (auto it = std::next(first); it != last; ++it)
    {
        double delta21 = *it - M1;
        double tmp = delta21/n;
        M4 = M4 + tmp*(tmp*tmp*delta21*((n-1)*(n*n-3*n+3)) + 6*tmp*M2 - 4*M3);
        M3 = M3 + tmp*((n-1)*(n-2)*delta21*tmp - 3*M2);
        M2 = M2 + tmp*(n-1)*delta21;
        M1 = M1 + tmp;
        n += 1;
    }

    return std::make_tuple(M1, M2, M3, M4, n-1);
}

template<typename ForwardIterator>
auto first_four_moments_real_impl(ForwardIterator first, ForwardIterator last)
{
    using Real = typename std::iterator_traits<ForwardIterator>::value_type;

    Real M1 = *first;
    Real M2 = 0;
    Real M3 = 0;
    Real M4 = 0;
    Real n = 2;
    for (auto it = std::next(first); it != last; ++it)
    {
        Real delta21 = *it - M1;
        Real tmp = delta21/n;
        M4 = M4 + tmp*(tmp*tmp*delta21*((n-1)*(n*n-3*n+3)) + 6*tmp*M2 - 4*M3);
        M3 = M3 + tmp*((n-1)*(n-2)*delta21*tmp - 3*M2);
        M2 = M2 + tmp*(n-1)*delta21;
        M1 = M1 + tmp;
        n += 1;
    }

    return std::make_tuple(M1, M2, M3, M4, n-1);
}

template<typename ReturnType, typename ForwardIterator>
ReturnType parallel_first_four_moments_impl(ForwardIterator first, ForwardIterator last)
{
    using Real = typename std::iterator_traits<ForwardIterator>::value_type;
    static unsigned num_threads {std::thread::hardware_concurrency()};

    const auto elements {std::distance(first, last)};
    const auto range_a {std::floor(elements / 2)};
    const auto range_b {elements - range_a};

    thread_counter.fetch_add(2);
    auto future_a {std::async(std::launch::async, [first, range_a]() -> ReturnType
    {
        if constexpr (std::is_integral_v<Real>)
        {
            if(thread_counter + 2 <= num_threads && range_a > 10)
                return parallel_first_four_moments_impl<ReturnType>(first, std::next(first, range_a));
            else
                return first_four_moments_integral_impl(first, std::next(first, range_a));
        }
        else
        {
            if(thread_counter + 2 <= num_threads && range_a > 10)
                return parallel_first_four_moments_impl<ReturnType>(first, std::next(first, range_a));
            else
                return first_four_moments_real_impl(first, std::next(first, range_a));
        }
    })};
    auto future_b {std::async(std::launch::async, [first, last, range_a]() -> ReturnType
    {
        if constexpr (std::is_integral_v<Real>)
        {
            if(thread_counter + 2 <= num_threads && range_a > 10)
                return parallel_first_four_moments_impl<ReturnType>(std::next(first, range_a), last);
            else
                return first_four_moments_integral_impl(std::next(first, range_a), last);
        }
        else
        {
            if(thread_counter + 2 <= num_threads && range_a > 10)
                return parallel_first_four_moments_impl<ReturnType>(std::next(first, range_a), last);
            else
                return first_four_moments_real_impl(std::next(first, range_a), last);
        }
    })};

    const auto results_a {future_a.get()};
    const auto results_b {future_b.get()};

    const auto M1_a = std::get<0>(results_a);
    const auto M2_a = std::get<1>(results_a);
    const auto M3_a = std::get<2>(results_a);
    const auto M4_a = std::get<3>(results_a);

    const auto M1_b = std::get<0>(results_b);
    const auto M2_b = std::get<1>(results_b);
    const auto M3_b = std::get<2>(results_b);
    const auto M4_b = std::get<3>(results_b);

    const auto n_ab = elements;
    const auto delta = M1_b - M1_a;
    
    const auto M1_ab = (range_a * M1_a + range_b * M1_b) / n_ab;
    const auto M2_ab = M2_a + M2_b + delta * delta * (range_a * range_b / n_ab);
    const auto M3_ab = M3_a + M3_b + (delta * delta * delta) * range_a * range_b * (range_a - range_b) / (n_ab * n_ab)    
                       + 3 * delta * (range_a * M2_b - range_b * M2_a) / n_ab;
    const auto M4_ab = M4_a + M4_b + (delta * delta * delta * delta) * range_a * range_b * (range_a * range_a - range_a * range_b + range_b * range_b) / (n_ab * n_ab * n_ab)
                       + 6 * delta * delta * (range_a * range_a * M2_b + range_b * range_b * M2_a) / (n_ab * n_ab) 
                       + 4 * delta * (range_a * M3_b - range_b * M3_a) / n_ab;

    return std::make_tuple(M1_ab, M2_ab, M3_ab, M4_ab, n_ab);
}

// Follows equation 1.5 of:
// https://prod.sandia.gov/techlib-noauth/access-control.cgi/2008/086212.pdf
template<class ForwardIterator>
auto skewness_sequential_impl(ForwardIterator first, ForwardIterator last)
{
    using Real = typename std::iterator_traits<ForwardIterator>::value_type;
    using std::sqrt;
    BOOST_ASSERT_MSG(first != last, "At least one sample is required to compute skewness.");
    if constexpr (std::is_integral<Real>::value)
    {
        double M1 = *first;
        double M2 = 0;
        double M3 = 0;
        double n = 2;
        for (auto it = std::next(first); it != last; ++it)
        {
            double delta21 = *it - M1;
            double tmp = delta21/n;
            M3 = M3 + tmp*((n-1)*(n-2)*delta21*tmp - 3*M2);
            M2 = M2 + tmp*(n-1)*delta21;
            M1 = M1 + tmp;
            n += 1;
        }

        double var = M2/(n-1);
        if (var == 0)
        {
            // The limit is technically undefined, but the interpretation here is clear:
            // A constant dataset has no skewness.
            return double(0);
        }
        double skew = M3/(M2*sqrt(var));
        return skew;
    }
    else
    {
        Real M1 = *first;
        Real M2 = 0;
        Real M3 = 0;
        Real n = 2;
        for (auto it = std::next(first); it != last; ++it)
        {
            Real delta21 = *it - M1;
            Real tmp = delta21/n;
            M3 += tmp*((n-1)*(n-2)*delta21*tmp - 3*M2);
            M2 += tmp*(n-1)*delta21;
            M1 += tmp;
            n += 1;
        }

        Real var = M2/(n-1);
        if (var == 0)
        {
            // The limit is technically undefined, but the interpretation here is clear:
            // A constant dataset has no skewness.
            return Real(0);
        }
        Real skew = M3/(M2*sqrt(var));
        return skew;
    }
}

// https://htor.inf.ethz.ch/publications/img/atomic-bench.pdf
template<typename ReturnType, typename ExecutionPolicy, typename RandomAccessIterator>
ReturnType gini_coefficient_parallel_impl(ExecutionPolicy&& exec, RandomAccessIterator first, RandomAccessIterator last)
{
    ReturnType i = 1;
    ReturnType num = 0;
    ReturnType denom = 0;
    
    std::for_each(exec, first, last, [&i, &num, &denom](auto val)
    {
        num = num + val * i;
        denom = denom + val;
        i = i + 1;
    });

    if(denom == 0)
    {
        return ReturnType(0);
    }

    return ((2*num)/denom - i)/(i-1);
}

template<class RandomAccessIterator>
auto gini_coefficient_sequential_impl(RandomAccessIterator first, RandomAccessIterator last)
{
    using Real = typename std::iterator_traits<RandomAccessIterator>::value_type;
    BOOST_ASSERT_MSG(first != last && std::next(first) != last, "Computation of the Gini coefficient requires at least two samples.");

    std::sort(first, last);
    if constexpr (std::is_integral<Real>::value)
    {
        double i = 1;
        double num = 0;
        double denom = 0;
        for (auto it = first; it != last; ++it)
        {
            num += *it*i;
            denom += *it;
            ++i;
        }

        // If the l1 norm is zero, all elements are zero, so every element is the same.
        if (denom == 0)
        {
            return double(0);
        }

        return ((2*num)/denom - i)/(i-1);
    }
    else
    {
        Real i = 1;
        Real num = 0;
        Real denom = 0;
        for (auto it = first; it != last; ++it)
        {
            num += *it*i;
            denom += *it;
            ++i;
        }

        // If the l1 norm is zero, all elements are zero, so every element is the same.
        if (denom == 0)
        {
            return Real(0);
        }

        return ((2*num)/denom - i)/(i-1);
    }
}
}

#endif // BOOST_MATH_STATISTICS_UNIVARIATE_STATISTICS_DETAIL_SINGLE_PASS_HPP
