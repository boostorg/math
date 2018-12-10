//  (C) Copyright Nick Thompson 2018.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_TOOLS_UNIVARIATE_STATISTICS_HPP
#define BOOST_MATH_TOOLS_UNIVARIATE_STATISTICS_HPP

#include <algorithm>
#include <iterator>
#include <boost/type_traits.hpp>
#include <boost/assert.hpp>
#include <boost/multiprecision/detail/number_base.hpp>


namespace boost{ namespace math{ namespace tools {

template<class ForwardIterator>
auto
mean(ForwardIterator first, ForwardIterator last)
{
    using Real = typename std::iterator_traits<ForwardIterator>::value_type;
    BOOST_ASSERT_MSG(first != last, "At least one sample is required to compute the mean.");
    if constexpr (std::is_integral<Real>::value)
    {
        double mu = 0;
        double i = 1;
        for(auto it = first; it != last; ++it) {
            mu = mu + (*it - mu)/i;
            i += 1;
        }
        return mu;
    }
    else
    {
        Real mu = 0;
        Real i = 1;
        for(auto it = first; it != last; ++it) {
            mu = mu + (*it - mu)/i;
            i += 1;
        }
        return mu;
    }
}

template<class Container>
inline auto mean(Container const & v)
{
    return mean(v.cbegin(), v.cend());
}

template<class ForwardIterator>
auto
mean_and_population_variance(ForwardIterator first, ForwardIterator last)
{
    using Real = typename std::iterator_traits<ForwardIterator>::value_type;
    BOOST_ASSERT_MSG(first != last, "At least one sample is required to compute mean and variance.");
    // Higham, Accuracy and Stability, equation 1.6a and 1.6b:
    if constexpr (std::is_integral<Real>::value)
    {
        double M = *first;
        double Q = 0;
        double k = 2;
        for (auto it = first + 1; it != last; ++it)
        {
            double tmp = *it - M;
            Q = Q + ((k-1)*tmp*tmp)/k;
            M = M + tmp/k;
            k += 1;
        }
        return std::make_pair(M, Q/(k-1));
    }
    else
    {
        Real M = *first;
        Real Q = 0;
        Real k = 2;
        for (auto it = first + 1; it != last; ++it)
        {
            Real tmp = *it - M;
            Q = Q + ((k-1)*tmp*tmp)/k;
            M = M + tmp/k;
            k += 1;
        }

        return std::make_pair(M, Q/(k-1));
    }
}

template<class Container>
inline auto mean_and_population_variance(Container const & v)
{
    return mean_and_population_variance(v.cbegin(), v.cend());
}

template<class RandomAccessIterator>
auto median(RandomAccessIterator first, RandomAccessIterator last)
{
    size_t num_elems = std::distance(first, last);
    BOOST_ASSERT_MSG(num_elems > 0, "The median of a zero length vector is undefined.");
    if (num_elems & 1)
    {
        auto middle = first + (num_elems - 1)/2;
        std::nth_element(first, middle, last);
        return *middle;
    }
    else
    {
        auto middle = first + num_elems/2 - 1;
        std::nth_element(first, middle, last);
        std::nth_element(middle, middle+1, last);
        return (*middle + *(middle+1))/2;
    }
}


template<class RandomAccessContainer>
inline auto median(RandomAccessContainer & v)
{
    return median(v.begin(), v.end());
}


template<class ForwardIterator>
auto gini_coefficient(ForwardIterator first, ForwardIterator last)
{
    using Real = typename std::iterator_traits<ForwardIterator>::value_type;
    BOOST_ASSERT_MSG(first != last && std::next(first) != last, "Computation of the Gini coefficient requires at least two samples.");

    std::sort(first, last);

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

    return ((2*num)/denom - i)/(i-2);
}

template<class RandomAccessContainer>
inline auto gini_coefficient(RandomAccessContainer & v)
{
    return gini_coefficient(v.begin(), v.end());
}

}}}
#endif
