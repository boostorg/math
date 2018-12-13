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

// Follows equation 1.5 of:
// https://prod.sandia.gov/techlib-noauth/access-control.cgi/2008/086212.pdf
template<class ForwardIterator>
auto
population_skewness(ForwardIterator first, ForwardIterator last)
{
    using Real = typename std::iterator_traits<ForwardIterator>::value_type;
    BOOST_ASSERT_MSG(first != last, "At least one sample is required to compute skewness.");
    if constexpr (std::is_integral<Real>::value)
    {
        double M1 = *first;
        double M2 = 0;
        double M3 = 0;
        double n = 2;
        for (auto it = first + 1; it != last; ++it)
        {
            double delta21 = *it - M1;
            double tmp = delta21/n;
            M3 = M3 + tmp*((n-1)*(n-2)*delta21*tmp - 3*M2);
            M2 = M2 + tmp*(n-1)*delta21;
            M1 = M1 + tmp;
            n += 1;
        }

        double variance = M2/(n-1);
        if (variance == 0)
        {
            // The limit is technically undefined, but the interpretation here is clear:
            // A constant dataset has no skewness.
            return double(0);
        }
        double skewness = M3/((n-1)*variance*sqrt(variance));
        return skewness;
    }
    else
    {
        Real M1 = *first;
        Real M2 = 0;
        Real M3 = 0;
        Real n = 2;
        for (auto it = first + 1; it != last; ++it)
        {
            Real delta21 = *it - M1;
            Real tmp = delta21/n;
            M3 = M3 + tmp*((n-1)*(n-2)*delta21*tmp - 3*M2);
            M2 = M2 + tmp*(n-1)*delta21;
            M1 = M1 + tmp;
            n += 1;
        }

        Real variance = M2/(n-1);
        if (variance == 0)
        {
            // The limit is technically undefined, but the interpretation here is clear:
            // A constant dataset has no skewness.
            return Real(0);
        }
        Real skewness = M3/((n-1)*variance*sqrt(variance));
        return skewness;
    }
}

template<class Container>
inline auto population_skewness(Container const & v)
{
    return population_skewness(v.cbegin(), v.cend());
}

// Follows equation 1.6 of:
// https://prod.sandia.gov/techlib-noauth/access-control.cgi/2008/086212.pdf
template<class ForwardIterator>
auto
population_kurtosis(ForwardIterator first, ForwardIterator last)
{
    using Real = typename std::iterator_traits<ForwardIterator>::value_type;
    BOOST_ASSERT_MSG(first != last, "At least one sample is required to compute kurtosis.");
    if constexpr (std::is_integral<Real>::value)
    {
        double M1 = *first;
        double M2 = 0;
        double M3 = 0;
        double M4 = 0;
        double n = 2;
        for (auto it = first + 1; it != last; ++it)
        {
            double delta21 = *it - M1;
            double tmp = delta21/n;
            M4 = M4 + tmp*(tmp*tmp*delta21*((n-1)*(n*n-3*n+3)) + 6*tmp*M2 - 4*M3);
            M3 = M3 + tmp*((n-1)*(n-2)*delta21*tmp - 3*M2);
            M2 = M2 + tmp*(n-1)*delta21;
            M1 = M1 + tmp;
            n += 1;
        }

        double variance = M2/(n-1);
        if (variance == 0)
        {
            return double(0);
        }
        double kurtosis = M4/((n-1)*variance*variance);
        return kurtosis;
    }
    else
    {
        Real M1 = *first;
        Real M2 = 0;
        Real M3 = 0;
        Real M4 = 0;
        Real n = 2;
        for (auto it = first + 1; it != last; ++it)
        {
            Real delta21 = *it - M1;
            Real tmp = delta21/n;
            M4 = M4 + tmp*(tmp*tmp*delta21*((n-1)*(n*n-3*n+3)) + 6*tmp*M2 - 4*M3);
            M3 = M3 + tmp*((n-1)*(n-2)*delta21*tmp - 3*M2);
            M2 = M2 + tmp*(n-1)*delta21;
            M1 = M1 + tmp;
            n += 1;
        }

        Real variance = M2/(n-1);
        if (variance == 0)
        {
            // Again, the limit is technically undefined, but the interpretation here is clear:
            // A constant dataset has no kurtosis.
            return Real(0);
        }
        Real kurtosis = M4/((n-1)*variance*variance);
        return kurtosis;
    }
}

template<class Container>
inline auto population_kurtosis(Container const & v)
{
    return population_kurtosis(v.cbegin(), v.cend());
}


// Follows equation 1.5/1.6 of:
// https://prod.sandia.gov/techlib-noauth/access-control.cgi/2008/086212.pdf
template<class ForwardIterator>
auto
first_four_moments(ForwardIterator first, ForwardIterator last)
{
    using Real = typename std::iterator_traits<ForwardIterator>::value_type;
    BOOST_ASSERT_MSG(first != last, "At least one sample is required to compute the first four moments.");
    if constexpr (std::is_integral<Real>::value)
    {
        double M1 = *first;
        double M2 = 0;
        double M3 = 0;
        double M4 = 0;
        double n = 2;
        for (auto it = first + 1; it != last; ++it)
        {
            double delta21 = *it - M1;
            double tmp = delta21/n;
            M4 = M4 + tmp*(tmp*tmp*delta21*((n-1)*(n*n-3*n+3)) + 6*tmp*M2 - 4*M3);
            M3 = M3 + tmp*((n-1)*(n-2)*delta21*tmp - 3*M2);
            M2 = M2 + tmp*(n-1)*delta21;
            M1 = M1 + tmp;
            n += 1;
        }

        return std::make_tuple(M1, M2/(n-1), M3/(n-1), M4/(n-1));
    }
    else
    {
        Real M1 = *first;
        Real M2 = 0;
        Real M3 = 0;
        Real M4 = 0;
        Real n = 2;
        for (auto it = first + 1; it != last; ++it)
        {
            Real delta21 = *it - M1;
            Real tmp = delta21/n;
            M4 = M4 + tmp*(tmp*tmp*delta21*((n-1)*(n*n-3*n+3)) + 6*tmp*M2 - 4*M3);
            M3 = M3 + tmp*((n-1)*(n-2)*delta21*tmp - 3*M2);
            M2 = M2 + tmp*(n-1)*delta21;
            M1 = M1 + tmp;
            n += 1;
        }

        return std::make_tuple(M1, M2/(n-1), M3/(n-1), M4/(n-1));
    }
}

template<class Container>
inline auto first_four_moments(Container const & v)
{
    return first_four_moments(v.cbegin(), v.cend());
}


}}}
#endif
