//  (C) Copyright Nick Thompson 2018.
//  (C) Copyright Matt Borland 2020.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_STATISTICS_UNIVARIATE_STATISTICS_HPP
#define BOOST_MATH_STATISTICS_UNIVARIATE_STATISTICS_HPP

#include <boost/math/statistics/detail/single_pass.hpp>
#include <boost/assert.hpp>
#include <algorithm>
#include <iterator>
#include <tuple>
#include <cmath>
#include <vector>
#include <type_traits>
#include <execution>
#include <numeric>
#include <list>

namespace boost::math::statistics {

template<class ExecutionPolicy, class ForwardIterator>
inline auto mean(ExecutionPolicy&& exec, ForwardIterator first, ForwardIterator last)
{
    using Real = typename std::iterator_traits<ForwardIterator>::value_type;
    BOOST_ASSERT_MSG(first != last, "At least one sample is required to compute the mean.");
    
    if constexpr (std::is_integral_v<Real>)
    {
        if constexpr (std::is_same_v<std::remove_reference_t<decltype(exec)>, decltype(std::execution::seq)>)
        {
            return detail::mean_sequential_impl<double>(first, last);
        }
        else
        {
            return std::reduce(exec, first, last, 0.0) / std::distance(first, last);
        }
    }
    else
    {
        if constexpr (std::is_same_v<std::remove_reference_t<decltype(exec)>, decltype(std::execution::seq)>)
        {
            return detail::mean_sequential_impl<Real>(first, last);
        }
        else
        {
            return std::reduce(exec, first, last, Real(0.0)) / Real(std::distance(first, last));
        }
    }
}

template<class ExecutionPolicy, class Container>
inline auto mean(ExecutionPolicy&& exec, Container const & v)
{
    return mean(exec, std::cbegin(v), std::cend(v));
}

template<class ForwardIterator>
inline auto mean(ForwardIterator first, ForwardIterator last)
{
    return mean(std::execution::seq, first, last);
}

template<class Container>
inline auto mean(Container const & v)
{
    return mean(std::execution::seq, std::cbegin(v), std::cend(v));
}

template<class ExecutionPolicy, class ForwardIterator>
inline auto variance(ExecutionPolicy&& exec, ForwardIterator first, ForwardIterator last)
{
    using Real = typename std::iterator_traits<ForwardIterator>::value_type;
    
    if constexpr (std::is_integral_v<Real>)
    {
        if constexpr (std::is_same_v<std::remove_reference_t<decltype(exec)>, decltype(std::execution::seq)>)
        {
           return std::get<2>(detail::variance_sequential_impl<std::tuple<double, double, double>>(first, last));
        }
        else
        {
            const auto results = detail::parallel_variance_impl<std::tuple<double, double, double>>(first, last);
            return std::get<1>(results) / std::get<2>(results);
        }
    }
    else
    {
        if constexpr (std::is_same_v<std::remove_reference_t<decltype(exec)>, decltype(std::execution::seq)>)
        {
            return std::get<2>(detail::variance_sequential_impl<std::tuple<Real, Real, Real>>(first, last));
        }
        else
        {
            const auto results = detail::parallel_variance_impl<std::tuple<Real, Real, Real>>(first, last);
            return std::get<1>(results) / std::get<2>(results);
        }
    }
}

template<class ExecutionPolicy, class Container>
inline auto variance(ExecutionPolicy&& exec, Container const & v)
{
    return variance(exec, std::cbegin(v), std::cend(v));
}

template<class ForwardIterator>
inline auto variance(ForwardIterator first, ForwardIterator last)
{
    return variance(std::execution::seq, first, last);
}

template<class Container>
inline auto variance(Container const & v)
{
    return variance(std::execution::seq, std::cbegin(v), std::cend(v));
}

template<class ExecutionPolicy, class ForwardIterator>
inline auto sample_variance(ExecutionPolicy&& exec, ForwardIterator first, ForwardIterator last)
{
    const auto n = std::distance(first, last);
    BOOST_ASSERT_MSG(n > 1, "At least two samples are required to compute the sample variance.");
    return n*variance(exec, first, last)/(n-1);
}

template<class ExecutionPolicy, class Container>
inline auto sample_variance(ExecutionPolicy&& exec, Container const & v)
{
    return sample_variance(exec, std::cbegin(v), std::cend(v));
}

template<class ForwardIterator>
inline auto sample_variance(ForwardIterator first, ForwardIterator last)
{
    return sample_variance(std::execution::seq, first, last);
}

template<class Container>
inline auto sample_variance(Container const & v)
{
    return sample_variance(std::execution::seq, std::cbegin(v), std::cend(v));
}

template<class ExecutionPolicy, class ForwardIterator>
inline auto mean_and_sample_variance(ExecutionPolicy&& exec, ForwardIterator first, ForwardIterator last)
{
    using Real = typename std::iterator_traits<ForwardIterator>::value_type;

    if constexpr (std::is_integral_v<Real>)
    {
        if constexpr (std::is_same_v<std::remove_reference_t<decltype(exec)>, decltype(std::execution::seq)>)
        {
            const auto results = detail::variance_sequential_impl<std::tuple<double, double, double>>(first, last);
            return std::make_pair(std::get<0>(results), std::get<2>(results));
        }
        else
        {
            const auto results = detail::parallel_variance_impl<std::tuple<double, double, double>>(first, last);
            return std::make_pair(std::get<0>(results), std::get<1>(results) / std::get<2>(results));
        }
    }
    else
    {
        if constexpr (std::is_same_v<std::remove_reference_t<decltype(exec)>, decltype(std::execution::seq)>)
        {
            const auto results = detail::variance_sequential_impl<std::tuple<Real, Real, Real>>(first, last);
            return std::make_pair(std::get<0>(results), std::get<2>(results));
        }
        else
        {
            const auto results = detail::parallel_variance_impl<std::tuple<Real, Real, Real>>(first, last);
            return std::make_pair(std::get<0>(results), std::get<1>(results) / std::get<2>(results));
        }
    }
}

template<class ExecutionPolicy, class Container>
inline auto mean_and_sample_variance(ExecutionPolicy&& exec, Container const & v)
{
    return mean_and_sample_variance(exec, std::cbegin(v), std::cend(v));
}

template<class ForwardIterator>
inline auto mean_and_sample_variance(ForwardIterator first, ForwardIterator last)
{
    return mean_and_sample_variance(std::execution::seq, first, last);
}

template<class Container>
inline auto mean_and_sample_variance(Container const & v)
{
    return mean_and_sample_variance(std::execution::seq, std::cbegin(v), std::cend(v));
}

template<class ExecutionPolicy, class ForwardIterator>
inline auto first_four_moments(ExecutionPolicy&& exec, ForwardIterator first, ForwardIterator last)
{
    using Real = typename std::iterator_traits<ForwardIterator>::value_type;

    if constexpr (std::is_integral_v<Real>)
    {
        if constexpr (std::is_same_v<std::remove_reference_t<decltype(exec)>, decltype(std::execution::seq)>)
        {
            const auto results = detail::first_four_moments_sequential_impl<std::tuple<double, double, double, double, double>>(first, last); 
            return std::make_tuple(std::get<0>(results), std::get<1>(results) / std::get<4>(results), std::get<2>(results) / std::get<4>(results), 
                                std::get<3>(results) / std::get<4>(results));
        }
        else
        {
            const auto results = detail::parallel_first_four_moments_impl<std::tuple<double, double, double, double, double>>(first, last);
            return std::make_tuple(std::get<0>(results), std::get<1>(results) / std::get<4>(results), std::get<2>(results) / std::get<4>(results), 
                                   std::get<3>(results) / std::get<4>(results));
        }
    }
    else
    {
        if constexpr (std::is_same_v<std::remove_reference_t<decltype(exec)>, decltype(std::execution::seq)>)
        {
            const auto results = detail::first_four_moments_sequential_impl<std::tuple<Real, Real, Real, Real, Real>>(first, last);
            return std::make_tuple(std::get<0>(results), std::get<1>(results) / std::get<4>(results), std::get<2>(results) / std::get<4>(results), 
                                   std::get<3>(results) / std::get<4>(results));
        }
        else
        {
            const auto results = detail::parallel_first_four_moments_impl<std::tuple<Real, Real, Real, Real, Real>>(first, last);
            return std::make_tuple(std::get<0>(results), std::get<1>(results) / std::get<4>(results), std::get<2>(results) / std::get<4>(results), 
                                   std::get<3>(results) / std::get<4>(results));
        }
    }
}

template<class ExecutionPolicy, class Container>
inline auto first_four_moments(ExecutionPolicy&& exec, Container const & v)
{
    return first_four_moments(exec, std::cbegin(v), std::cend(v));
}

template<class ForwardIterator>
inline auto first_four_moments(ForwardIterator first, ForwardIterator last)
{
    return first_four_moments(std::execution::seq, first, last);
}

template<class Container>
inline auto first_four_moments(Container const & v)
{
    return first_four_moments(std::execution::seq, std::cbegin(v), std::cend(v));
}

// https://prod.sandia.gov/techlib-noauth/access-control.cgi/2008/086212.pdf
template<class ExecutionPolicy, class ForwardIterator>
inline auto skewness(ExecutionPolicy&& exec, ForwardIterator first, ForwardIterator last)
{
    using Real = typename std::iterator_traits<ForwardIterator>::value_type;

    if constexpr (std::is_same_v<std::remove_reference_t<decltype(exec)>, decltype(std::execution::seq)>)
    {
        if constexpr (std::is_integral_v<Real>)
        {
            return detail::skewness_sequential_impl<double>(first, last);
        }
        else
        {
            return detail::skewness_sequential_impl<Real>(first, last);
        }
    }
    else 
    {
        const auto [M1, M2, M3, M4] = first_four_moments(exec, first, last);
        const auto n = std::distance(first, last);
        const auto var = M2/(n-1);

        if (M2 == 0)
        {
            // The limit is technically undefined, but the interpretation here is clear:
            // A constant dataset has no skewness.
            if constexpr (std::is_integral_v<Real>)
            {
                return double(0);
            }
            else
            {
                return Real(0);
            }
        }
        else
        {
            return M3/(M2*sqrt(var)) / Real(2);
        }
    }
}

template<class ExecutionPolicy, class Container>
inline auto skewness(ExecutionPolicy&& exec, Container & v)
{
    return skewness(exec, std::cbegin(v), std::cend(v));
}

template<class ForwardIterator>
inline auto skewness(ForwardIterator first, ForwardIterator last)
{
    return skewness(std::execution::seq, first, last);
}

template<class Container>
inline auto skewness(Container const & v)
{
    return skewness(std::execution::seq, std::cbegin(v), std::cend(v));
}

// Follows equation 1.6 of:
// https://prod.sandia.gov/techlib-noauth/access-control.cgi/2008/086212.pdf
template<class ExecutionPolicy, class ForwardIterator>
inline auto kurtosis(ExecutionPolicy&& exec, ForwardIterator first, ForwardIterator last)
{
    const auto [M1, M2, M3, M4] = first_four_moments(exec, first, last);
    if (M2 == 0)
    {
        return M2;
    }
    return M4/(M2*M2);
}

template<class ExecutionPolicy, class Container>
inline auto kurtosis(ExecutionPolicy&& exec, Container const & v)
{
    return kurtosis(exec, std::cbegin(v), std::cend(v));
}

template<class ForwardIterator>
inline auto kurtosis(ForwardIterator first, ForwardIterator last)
{
    return kurtosis(std::execution::seq, first, last);
}

template<class Container>
inline auto kurtosis(Container const & v)
{
    return kurtosis(std::execution::seq, std::cbegin(v), std::cend(v));
}

template<class ExecutionPolicy, class ForwardIterator>
inline auto excess_kurtosis(ExecutionPolicy&& exec, ForwardIterator first, ForwardIterator last)
{
    return kurtosis(exec, first, last) - 3;
}

template<class ExecutionPolicy, class Container>
inline auto excess_kurtosis(ExecutionPolicy&& exec, Container const & v)
{
    return excess_kurtosis(exec, std::cbegin(v), std::cend(v));
}

template<class ForwardIterator>
inline auto excess_kurtosis(ForwardIterator first, ForwardIterator last)
{
    return excess_kurtosis(std::execution::seq, first, last);
}

template<class Container>
inline auto excess_kurtosis(Container const & v)
{
    return excess_kurtosis(std::execution::seq, std::cbegin(v), std::cend(v));
}


template<class ExecutionPolicy, class RandomAccessIterator>
auto median(ExecutionPolicy&& exec, RandomAccessIterator first, RandomAccessIterator last)
{
    const auto num_elems = std::distance(first, last);
    BOOST_ASSERT_MSG(num_elems > 0, "The median of a zero length vector is undefined.");
    if (num_elems & 1)
    {
        auto middle = first + (num_elems - 1)/2;
        std::nth_element(exec, first, middle, last);
        return *middle;
    }
    else
    {
        auto middle = first + num_elems/2 - 1;
        std::nth_element(exec, first, middle, last);
        std::nth_element(exec, middle, middle+1, last);
        return (*middle + *(middle+1))/2;
    }
}


template<class ExecutionPolicy, class RandomAccessContainer>
inline auto median(ExecutionPolicy&& exec, RandomAccessContainer & v)
{
    return median(exec, std::begin(v), std::end(v));
}

template<class RandomAccessIterator>
inline auto median(RandomAccessIterator first, RandomAccessIterator last)
{
    return median(std::execution::seq, first, last);
}

template<class RandomAccessContainer>
inline auto median(RandomAccessContainer & v)
{
    return median(std::execution::seq, std::begin(v), std::end(v));
}

template<class ExecutionPolicy, class RandomAccessIterator>
inline auto gini_coefficient(ExecutionPolicy&& exec, RandomAccessIterator first, RandomAccessIterator last)
{
    using Real = typename std::iterator_traits<RandomAccessIterator>::value_type;

    // C++20 enables constexpr std::is_sorted and compile time std::sort
    #if __cplusplus > 201900 || _MSVC_LANG > 201900
    if constexpr (!std::is_sorted(first, last))
    {
        std::sort(first, last);
    }
    #else
    if(!std::is_sorted(exec, first, last))
    {
        std::sort(exec, first, last);
    }
    #endif

    if constexpr (std::is_same_v<std::remove_reference_t<decltype(exec)>, decltype(std::execution::seq)>)
    {
        if constexpr (std::is_integral_v<Real>)
        {
            return detail::gini_coefficient_sequential_impl<double>(first, last);
        }
        else
        {
            return detail::gini_coefficient_sequential_impl<Real>(first, last);
        }   
    }
    
    else if constexpr (std::is_integral_v<Real>)
    {
        return detail::gini_coefficient_parallel_impl<double>(exec, first, last);
    }

    else
    {
        return detail::gini_coefficient_parallel_impl<Real>(exec, first, last);
    }
}

template<class ExecutionPolicy, class RandomAccessContainer>
inline auto gini_coefficient(ExecutionPolicy&& exec, RandomAccessContainer & v)
{
    return gini_coefficient(exec, std::begin(v), std::end(v));
}

template<class RandomAccessIterator>
inline auto gini_coefficient(RandomAccessIterator first, RandomAccessIterator last)
{
    return gini_coefficient(std::execution::seq, first, last);
}

template<class RandomAccessContainer>
inline auto gini_coefficient(RandomAccessContainer & v)
{
    return gini_coefficient(std::execution::seq, std::begin(v), std::end(v));
}

template<class ExecutionPolicy, class RandomAccessIterator>
inline auto sample_gini_coefficient(ExecutionPolicy&& exec, RandomAccessIterator first, RandomAccessIterator last)
{
    const auto n = std::distance(first, last);
    return n*gini_coefficient(exec, first, last)/(n-1);
}

template<class ExecutionPolicy, class RandomAccessContainer>
inline auto sample_gini_coefficient(ExecutionPolicy&& exec, RandomAccessContainer & v)
{
    return sample_gini_coefficient(exec, std::begin(v), std::end(v));
}

template<class RandomAccessIterator>
inline auto sample_gini_coefficient(RandomAccessIterator first, RandomAccessIterator last)
{
    return sample_gini_coefficient(std::execution::seq, first, last);
}

template<class RandomAccessContainer>
inline auto sample_gini_coefficient(RandomAccessContainer & v)
{
    return sample_gini_coefficient(std::execution::seq, std::begin(v), std::end(v));
}

template<class ExecutionPolicy, class RandomAccessIterator>
auto median_absolute_deviation(ExecutionPolicy&& exec, RandomAccessIterator first, RandomAccessIterator last, 
    typename std::iterator_traits<RandomAccessIterator>::value_type center=std::numeric_limits<typename std::iterator_traits<RandomAccessIterator>::value_type>::quiet_NaN())
{
    using std::abs;
    using Real = typename std::iterator_traits<RandomAccessIterator>::value_type;
    using std::isnan;
    if (isnan(center))
    {
        center = boost::math::statistics::median(exec, first, last);
    }
    const auto num_elems = std::distance(first, last);
    BOOST_ASSERT_MSG(num_elems > 0, "The median of a zero-length vector is undefined.");
    auto comparator = [&center](Real a, Real b) { return abs(a-center) < abs(b-center);};
    if (num_elems & 1)
    {
        auto middle = first + (num_elems - 1)/2;
        std::nth_element(exec, first, middle, last, comparator);
        return abs(*middle);
    }
    else
    {
        auto middle = first + num_elems/2 - 1;
        std::nth_element(exec, first, middle, last, comparator);
        std::nth_element(exec, middle, middle+1, last, comparator);
        return (abs(*middle) + abs(*(middle+1)))/abs(static_cast<Real>(2));
    }
}

template<class ExecutionPolicy, class RandomAccessContainer>
inline auto median_absolute_deviation(ExecutionPolicy&& exec, RandomAccessContainer & v, 
    typename RandomAccessContainer::value_type center=std::numeric_limits<typename RandomAccessContainer::value_type>::quiet_NaN())
{
    return median_absolute_deviation(exec, std::begin(v), std::end(v), center);
}

template<class RandomAccessIterator>
inline auto median_absolute_deviation(RandomAccessIterator first, RandomAccessIterator last, 
    typename RandomAccessIterator::value_type center=std::numeric_limits<typename RandomAccessIterator::value_type>::quiet_NaN())
{
    return median_absolute_deviation(std::execution::seq, first, last, center);
}

template<class RandomAccessContainer>
inline auto median_absolute_deviation(RandomAccessContainer & v, 
    typename RandomAccessContainer::value_type center=std::numeric_limits<typename RandomAccessContainer::value_type>::quiet_NaN())
{
    return median_absolute_deviation(std::execution::seq, std::begin(v), std::end(v), center);
}

template<class ExecutionPolicy, class ForwardIterator>
auto interquartile_range(ExecutionPolicy&& exec, ForwardIterator first, ForwardIterator last)
{
    using Real = typename std::iterator_traits<ForwardIterator>::value_type;
    static_assert(!std::is_integral_v<Real>, "Integer values have not yet been implemented.");
    auto m = std::distance(first,last);
    BOOST_ASSERT_MSG(m >= 3, "At least 3 samples are required to compute the interquartile range.");
    auto k = m/4;
    auto j = m - (4*k);
    // m = 4k+j.
    // If j = 0 or j = 1, then there are an even number of samples below the median, and an even number above the median.
    //    Then we must average adjacent elements to get the quartiles.
    // If j = 2 or j = 3, there are an odd number of samples above and below the median, these elements may be directly extracted to get the quartiles.

    if (j==2 || j==3)
    {
        auto q1 = first + k;
        auto q3 = first + 3*k + j - 1;
        std::nth_element(exec, first, q1, last);
        Real Q1 = *q1;
        std::nth_element(exec, q1, q3, last);
        Real Q3 = *q3;
        return Q3 - Q1;
    } else {
        // j == 0 or j==1:
        auto q1 = first + k - 1;
        auto q3 = first + 3*k - 1 + j;
        std::nth_element(exec, first, q1, last);
        Real a = *q1;
        std::nth_element(exec, q1, q1 + 1, last);
        Real b = *(q1 + 1);
        Real Q1 = (a+b)/2;
        std::nth_element(exec, q1, q3, last);
        a = *q3;
        std::nth_element(exec, q3, q3 + 1, last);
        b = *(q3 + 1);
        Real Q3 = (a+b)/2;
        return Q3 - Q1;
    }
}

template<class ExecutionPolicy, class RandomAccessContainer>
inline auto interquartile_range(ExecutionPolicy&& exec, RandomAccessContainer & v)
{
    return interquartile_range(exec, std::begin(v), std::end(v));
}

template<class RandomAccessIterator>
inline auto interquartile_range(RandomAccessIterator first, RandomAccessIterator last)
{
    return interquartile_range(std::execution::seq, first, last);
}

template<class RandomAccessContainer>
inline auto interquartile_range(RandomAccessContainer & v)
{
    return interquartile_range(std::execution::seq, std::begin(v), std::end(v));
}

template<class ExecutionPolicy, class ForwardIterator, class OutputIterator>
inline OutputIterator mode(ExecutionPolicy&& exec, ForwardIterator first, ForwardIterator last, OutputIterator output)
{   
    #if __cplusplus > 201900 || _MSVC_LANG > 201900
    if constexpr (!std::is_sorted(first, last))
    {
        if constexpr (std::is_same_v<typename std::iterator_traits<ForwardIterator>::iterator_category(), std::random_access_iterator_tag>)
        {
            std::sort(first, last);
        }
        else
        {
            BOOST_ASSERT("Data must be sorted for sequential mode calculation");
        }
    }
    #else
    if(!std::is_sorted(exec, first, last))
    {
        if constexpr (std::is_same_v<typename std::iterator_traits<ForwardIterator>::iterator_category(), std::random_access_iterator_tag>)
        {
            std::sort(exec, first, last);
        }
        else
        {
            BOOST_ASSERT("Data must be sorted for sequential mode calculation");
        }
    }
    #endif

    return detail::mode_impl(first, last, output);
}

template<class ExecutionPolicy, class Container, class OutputIterator>
inline OutputIterator mode(ExecutionPolicy&& exec, Container & v, OutputIterator output)
{
    return mode(exec, std::begin(v), std::end(v), output);
}

template<class ForwardIterator, class OutputIterator>
inline OutputIterator mode(ForwardIterator first, ForwardIterator last, OutputIterator output)
{
    return mode(std::execution::seq, first, last, output);
}

// Requires enable_if_t to not clash with impl that returns std::list
template<class Container, class OutputIterator, std::enable_if_t<std::is_execution_policy_v<Container>, bool> = true>
inline OutputIterator mode(Container & v, OutputIterator output)
{
    return mode(std::execution::seq, std::begin(v), std::end(v), output);
}

// std::list is the return type for the proposed STL stats library

template<class ExecutionPolicy, class ForwardIterator, class Real = typename std::iterator_traits<ForwardIterator>::value_type>
inline auto mode(ExecutionPolicy&& exec, ForwardIterator first, ForwardIterator last)
{
    std::list<Real> modes;
    mode(exec, first, last, modes.begin());
    return modes;
}

template<class ExecutionPolicy, class Container>
inline auto mode(ExecutionPolicy&& exec, Container & v)
{
    return mode(exec, std::begin(v), std::end(v));
}

template<class ForwardIterator>
inline auto mode(ForwardIterator first, ForwardIterator last)
{
    return mode(std::execution::seq, first, last);
}

template<class Container>
inline auto mode(Container & v)
{
    return mode(std::execution::seq, std::begin(v), std::end(v));
}
}
#endif
