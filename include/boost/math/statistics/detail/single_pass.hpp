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
#include <cmath>
#include <algorithm>
#include <valarray>

namespace boost { namespace math { namespace statistics { namespace detail {

template<typename ReturnType, typename ForwardIterator>
ReturnType mean_sequential_impl(ForwardIterator first, ForwardIterator last)
{
    const std::size_t elements {static_cast<std::size_t>(std::distance(first, last))};
    std::valarray<ReturnType> mu {0, 0, 0, 0};
    std::valarray<ReturnType> temp {0, 0, 0, 0};
    ReturnType i {1};
    const ForwardIterator end {std::next(first, elements - (elements % 4))};
    ForwardIterator it {first};

    while(it != end)
    {
        const ReturnType inv {ReturnType(1) / i};
        temp = {static_cast<ReturnType>(*it++), static_cast<ReturnType>(*it++), static_cast<ReturnType>(*it++), static_cast<ReturnType>(*it++)};
        temp -= mu;
        mu += (temp *= inv);
        i += 1;
    }

    const ReturnType num1 {ReturnType(elements - (elements % 4))/ReturnType(4)};
    const ReturnType num2 {num1 + ReturnType(elements % 4)};

    while(it != last)
    {
        mu[3] += (*it-mu[3])/i;
        i += 1;
        ++it;
    }

    return (num1 * std::valarray<ReturnType>(mu[std::slice(0,3,1)]).sum() + num2 * mu[3]) / ReturnType(elements);
}

// Higham, Accuracy and Stability, equation 1.6a and 1.6b:
// Calculates Mean, M2, and variance
template<typename ReturnType, typename ForwardIterator>
ReturnType variance_sequential_impl(ForwardIterator first, ForwardIterator last)
{
    using Real = typename std::tuple_element<0, ReturnType>::type;

    Real M = *first;
    Real Q = 0;
    Real k = 2;
    Real M2 = 0;

    for(auto it = std::next(first); it != last; ++it)
    {
        Real tmp = (*it - M) / k;
        Real delta_1 = *it - M;
        Q += k*(k-1)*tmp*tmp;
        M += tmp;
        k += 1;
        Real delta_2 = *it - M;
        M2 += delta_1 * delta_2;
    }

    return std::make_tuple(M, M2, Q/(k-1));
}

// https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Higher-order_statistics
template<typename ReturnType, typename ForwardIterator>
ReturnType first_four_moments_sequential_impl(ForwardIterator first, ForwardIterator last)
{
    using Real = typename std::tuple_element<0, ReturnType>::type;
    using Size = typename std::tuple_element<4, ReturnType>::type;

    Real M1 = *first;
    Real M2 = 0;
    Real M3 = 0;
    Real M4 = 0;
    Size n = 2;
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

// TODO(mborland): replace recursive decomposition with linear
// TODO(mborland): refactor name to match the rest e.g. first_four_moments_parallel_impl
template<typename ReturnType, typename ForwardIterator>
ReturnType parallel_first_four_moments_impl(ForwardIterator first, ForwardIterator last)
{
    using Real = typename std::tuple_element<0, ReturnType>::type;
    
    static unsigned num_threads {std::thread::hardware_concurrency()};

    const auto elements {std::distance(first, last)};
    const auto range_a {std::floor(elements / 2)};
    const auto range_b {elements - range_a};
    static std::atomic<unsigned> thread_counter {};
    thread_counter.fetch_add(2);
    auto future_a {std::async(std::launch::async | std::launch::deferred, [first, range_a]() -> ReturnType
    {
        if(thread_counter + 2 <= num_threads && range_a > 10)
        {
            return parallel_first_four_moments_impl<ReturnType>(first, std::next(first, range_a));
        }
        else
        {
            return first_four_moments_sequential_impl<ReturnType>(first, std::next(first, range_a));
        }
    })};
    auto future_b {std::async(std::launch::async | std::launch::deferred, [first, last, range_a]() -> ReturnType
    {
        if(thread_counter + 2 <= num_threads && range_a > 10)
        {
            return parallel_first_four_moments_impl<ReturnType>(std::next(first, range_a), last);
        }
        else
        {
            return first_four_moments_sequential_impl<ReturnType>(std::next(first, range_a), last);
        }
    })};

    const auto results_a {future_a.get()};
    const auto results_b {future_b.get()};
    thread_counter.fetch_sub(2);

    const Real M1_a = std::get<0>(results_a);
    const Real M2_a = std::get<1>(results_a);
    const Real M3_a = std::get<2>(results_a);
    const Real M4_a = std::get<3>(results_a);

    const Real M1_b = std::get<0>(results_b);
    const Real M2_b = std::get<1>(results_b);
    const Real M3_b = std::get<2>(results_b);
    const Real M4_b = std::get<3>(results_b);

    const Real n_ab = elements;
    const Real delta = M1_b - M1_a;
    
    const Real M1_ab = (range_a * M1_a + range_b * M1_b) / n_ab;
    const Real M2_ab = M2_a + M2_b + delta * delta * (range_a * range_b / n_ab);
    const Real M3_ab = M3_a + M3_b + (delta * delta * delta) * range_a * range_b * (range_a - range_b) / (n_ab * n_ab)    
                       + 3 * delta * (range_a * M2_b - range_b * M2_a) / n_ab;
    const Real M4_ab = M4_a + M4_b + (delta * delta * delta * delta) * range_a * range_b * (range_a * range_a - range_a * range_b + range_b * range_b) / (n_ab * n_ab * n_ab)
                       + 6 * delta * delta * (range_a * range_a * M2_b + range_b * range_b * M2_a) / (n_ab * n_ab) 
                       + 4 * delta * (range_a * M3_b - range_b * M3_a) / n_ab;

    return std::make_tuple(M1_ab, M2_ab, M3_ab, M4_ab, n_ab);
}

// Follows equation 1.5 of:
// https://prod.sandia.gov/techlib-noauth/access-control.cgi/2008/086212.pdf
template<typename ReturnType, typename ForwardIterator>
ReturnType skewness_sequential_impl(ForwardIterator first, ForwardIterator last)
{
    using std::sqrt;
    BOOST_ASSERT_MSG(first != last, "At least one sample is required to compute skewness.");
    
    ReturnType M1 = *first;
    ReturnType M2 = 0;
    ReturnType M3 = 0;
    ReturnType n = 2;
        
    for (auto it = std::next(first); it != last; ++it)    
    {
        ReturnType delta21 = *it - M1;
        ReturnType tmp = delta21/n;
        M3 += tmp*((n-1)*(n-2)*delta21*tmp - 3*M2);
        M2 += tmp*(n-1)*delta21;
        M1 += tmp;
        n += 1;
    }
   
    ReturnType var = M2/(n-1);
    
    if (var == 0)
    {
        // The limit is technically undefined, but the interpretation here is clear:
        // A constant dataset has no skewness.
        return ReturnType(0);
    }
    
    ReturnType skew = M3/(M2*sqrt(var));
    return skew;
}

template<typename ReturnType, typename ExecutionPolicy, typename RandomAccessIterator>
ReturnType gini_coefficient_parallel_impl(ExecutionPolicy&& exec, RandomAccessIterator first, RandomAccessIterator last)
{
    using Real = typename std::iterator_traits<RandomAccessIterator>::value_type;
    
    ReturnType i = 1;
    ReturnType num = 0;
    ReturnType denom = 0;
    
    std::for_each(exec, first, last, [&i, &num, &denom](const Real& val)
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

template<typename ReturnType, typename ForwardIterator>
ReturnType gini_coefficient_sequential_impl(ForwardIterator first, ForwardIterator last)
{
    ReturnType i = 1;
    ReturnType num = 0;
    ReturnType denom = 0;

    for(auto it = first; it != last; ++it)
    {
        num += *it*i;
        denom += *it;
        ++i;
    }

    // If the l1 norm is zero, all elements are zero, so every element is the same.
    if(denom == 0)
    {
        return ReturnType(0);
    }
    else
    {
        return ((2*num)/denom - i)/(i-1);
    }
}

template<typename ForwardIterator, typename OutputIterator>
OutputIterator mode_impl(ForwardIterator first, ForwardIterator last, OutputIterator output)
{
    using Z = typename std::iterator_traits<ForwardIterator>::value_type;
    using Size = typename std::iterator_traits<ForwardIterator>::difference_type;

    std::vector<Z> modes {};
    modes.reserve(16);
    Size max_counter {0};

    while(first != last)
    {
        Size current_count {0};
        ForwardIterator end_it {first};
        while(end_it != last && *end_it == *first)
        {
            ++current_count;
            ++end_it;
        }

        if(current_count > max_counter)
        {
            modes.resize(1);
            modes[0] = *first;
            max_counter = current_count;
        }

        else if(current_count == max_counter)
        {
            modes.emplace_back(*first);
        }

        first = end_it;
    }

    return std::move(modes.begin(), modes.end(), output);
}
}}}}

#endif // BOOST_MATH_STATISTICS_UNIVARIATE_STATISTICS_DETAIL_SINGLE_PASS_HPP
