//  (C) Copyright Nick Thompson 2018.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_STATISTICS_UNIVARIATE_STATISTICS_HPP
#define BOOST_MATH_STATISTICS_UNIVARIATE_STATISTICS_HPP

#include <boost/math/statistics/detail/single_pass.hpp>
#include <algorithm>
#include <iterator>
#include <tuple>
#include <cmath>
#include <vector>
#include <boost/assert.hpp>
#include <atomic>
#include <type_traits>
#include <execution>
#include <numeric>
#include <valarray>
#include <thread>
#include <future>
#include <iostream>

namespace boost::math::statistics {

template<class ExecutionPolicy, class ForwardIterator>
auto mean(ExecutionPolicy&& exec, ForwardIterator first, ForwardIterator last)
{
    using Real = typename std::iterator_traits<ForwardIterator>::value_type;
    BOOST_ASSERT_MSG(first != last, "At least one sample is required to compute the mean.");
    if constexpr (std::is_integral_v<Real>)
    {
        if constexpr (std::is_same_v<std::remove_reference_t<decltype(exec)>, decltype(std::execution::seq)>)
        {
            double mu = 0;
            double i = 1;
            for(auto it = first; it != last; ++it) 
            {
                mu = mu + (*it - mu)/i;
                i += 1;
            }
            return mu;
        }
        else
        {
            return std::reduce(exec, first, last, 0.0) / std::distance(first, last);
        }
    }
    else if constexpr (std::is_same_v<typename std::iterator_traits<ForwardIterator>::iterator_category, std::random_access_iterator_tag>)
    {
        std::size_t elements {std::distance(first, last)};
        std::valarray<Real> mu {0, 0, 0, 0};
        std::valarray<Real> temp {0, 0, 0, 0};
        Real i = 1;
        auto end {last - (elements % 4)};

        for(auto it {first}; it != end; it += 4)
        {
            const Real inv {Real(1) / i};
            temp = {*it, *(it+1), *(it+2), *(it+3)};
            temp -= mu;
            mu += (temp *= inv);
            i += 1;
        }

        Real num1 = Real(elements - (elements %4))/Real(4);
        Real num2 = num1 + Real(elements % 4);

        for (auto it = end; it != last; ++it)
        {
            mu[3] += (*it-mu[3])/i;
            i += 1;
        }

        return (num1 * std::valarray<Real>(mu[std::slice(0,3,1)]).sum() + num2 * mu[3]) / Real(elements);
    }
    else
    {
        if constexpr (std::is_same_v<std::remove_reference_t<decltype(exec)>, decltype(std::execution::seq)>)
        {
            auto it = first;
            Real mu = *it;
            Real i = 2;
            while(++it != last)
            {
                mu += (*it - mu)/i;
                i += 1;
            }
            return mu;
        }
        else
        {
            return std::reduce(exec, first, last, static_cast<Real>(0.0)) / std::distance(first, last);
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
            return std::get<2>(detail::variance_integeral_impl(first, last));
        }
        else
        {
            detail::thread_counter = 1;
            const auto results = detail::parallel_variance_impl<std::tuple<double, double, double>>(first, last);
            return std::get<1>(results) / std::get<2>(results);
        }
    }
    else
    {
        if constexpr (std::is_same_v<std::remove_reference_t<decltype(exec)>, decltype(std::execution::seq)>)
        {
            return std::get<2>(detail::variance_real_impl(first, last));
        }
        else
        {
            detail::thread_counter = 1;
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
            const auto results = detail::variance_integeral_impl(first, last);
            return std::make_pair(std::get<0>(results), std::get<2>(results));
        }
        else
        {
            detail::thread_counter = 1;
            const auto results = detail::parallel_variance_impl<std::tuple<double, double, double>>(first, last);
            return std::make_pair(std::get<0>(results), std::get<1>(results) / std::get<2>(results));
        }
    }
    else
    {
        if constexpr (std::is_same_v<std::remove_reference_t<decltype(exec)>, decltype(std::execution::seq)>)
        {
            const auto results = detail::variance_real_impl(first, last);
            return std::make_pair(std::get<0>(results), std::get<2>(results));
        }
        else
        {
            detail::thread_counter = 1;
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

// TODO(mborland): find/try decompistion algo with impl split by type similar to variance
//                 Should satisfy this, first four momenets, kurtosis, and excess kurtosis
//
// TODO(mborland): Reorder skewness below first four moments and call it like kurtosis does
// Follows equation 1.5 of:
// https://prod.sandia.gov/techlib-noauth/access-control.cgi/2008/086212.pdf
template<class ForwardIterator>
auto skewness(ForwardIterator first, ForwardIterator last)
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

template<class Container>
inline auto skewness(Container const & v)
{
    return skewness(v.cbegin(), v.cend());
}

template<class ExecutionPolicy, class ForwardIterator>
inline auto first_four_moments(ExecutionPolicy&& exec, ForwardIterator first, ForwardIterator last)
{
    using Real = typename std::iterator_traits<ForwardIterator>::value_type;

    if constexpr (std::is_integral_v<Real>)
    {
        if constexpr (std::is_same_v<std::remove_reference_t<decltype(exec)>, decltype(std::execution::seq)>)
        {
            const auto results = detail::first_four_moments_integral_impl(first, last);
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
            const auto results = detail::first_four_moments_real_impl(first, last);
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


// Follows equation 1.6 of:
// https://prod.sandia.gov/techlib-noauth/access-control.cgi/2008/086212.pdf
template<class ForwardIterator>
auto kurtosis(ForwardIterator first, ForwardIterator last)
{
    auto [M1, M2, M3, M4] = first_four_moments(first, last);
    if (M2 == 0)
    {
        return M2;
    }
    return M4/(M2*M2);
}

template<class Container>
inline auto kurtosis(Container const & v)
{
    return kurtosis(v.cbegin(), v.cend());
}

template<class ForwardIterator>
auto excess_kurtosis(ForwardIterator first, ForwardIterator last)
{
    return kurtosis(first, last) - 3;
}

template<class Container>
inline auto excess_kurtosis(Container const & v)
{
    return excess_kurtosis(v.cbegin(), v.cend());
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

// https://htor.inf.ethz.ch/publications/img/atomic-bench.pdf
template<class RandomAccessIterator>
auto gini_coefficient(RandomAccessIterator first, RandomAccessIterator last)
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

template<class RandomAccessContainer>
inline auto gini_coefficient(RandomAccessContainer & v)
{
    return gini_coefficient(v.begin(), v.end());
}

template<class RandomAccessIterator>
inline auto sample_gini_coefficient(RandomAccessIterator first, RandomAccessIterator last)
{
    size_t n = std::distance(first, last);
    return n*gini_coefficient(first, last)/(n-1);
}

template<class RandomAccessContainer>
inline auto sample_gini_coefficient(RandomAccessContainer & v)
{
    return sample_gini_coefficient(v.begin(), v.end());
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

template<class ForwardIterator, class OutputIterator>
auto sorted_mode(ForwardIterator first, ForwardIterator last, OutputIterator output) -> decltype(output)
{
    using Z = typename std::iterator_traits<ForwardIterator>::value_type;
    static_assert(std::is_integral<Z>::value, "Floating point values have not yet been implemented.");
    using Size = typename std::iterator_traits<ForwardIterator>::difference_type;

    std::vector<Z> modes {};
    modes.reserve(16);
    Size max_counter {0};

    while(first != last)
    {
        Size current_count {0};
        auto end_it {first};
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

template<class Container, class OutputIterator>
inline auto sorted_mode(Container & v, OutputIterator output) -> decltype(output)
{
    return sorted_mode(v.begin(), v.end(), output);
}

template<class RandomAccessIterator, class OutputIterator>
auto mode(RandomAccessIterator first, RandomAccessIterator last, OutputIterator output) -> decltype(output)
{
    std::sort(first, last);
    return sorted_mode(first, last, output);
}

template<class RandomAccessContainer, class OutputIterator>
inline auto mode(RandomAccessContainer & v, OutputIterator output) -> decltype(output)
{
    return mode(v.begin(), v.end(), output);
}

}
#endif
