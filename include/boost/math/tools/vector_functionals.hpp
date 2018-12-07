//  (C) Copyright Nick Thompson 2018.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_TOOLS_VECTOR_FUNCTIONALS_HPP
#define BOOST_MATH_TOOLS_VECTOR_FUNCTIONALS_HPP

#include <algorithm>
#include <iterator>
#include <iostream>
#include <boost/type_traits.hpp>
#include <boost/assert.hpp>
/*
 * A set of tools for computing scalar quantities associated with lists of numbers.
 */


namespace boost{ namespace math{ namespace tools {

template<class ForwardIterator>
auto
mean(ForwardIterator first, ForwardIterator last)
{
    typedef typename std::remove_const<typename std::remove_reference<decltype(*std::declval<ForwardIterator>())>::type>::type Real;
    BOOST_ASSERT_MSG(first != last, "At least one sample is required to compute the mean.");
    Real mu = 0;
    Real i = 1;
    for(auto it = first; it != last; ++it) {
        mu = mu + (*it - mu)/i;
        i += 1;
    }
    return mu;
}

template<class ForwardIterator>
auto
mean_and_population_variance(ForwardIterator first, ForwardIterator last)
{
    typedef typename std::remove_const<typename std::remove_reference<decltype(*std::declval<ForwardIterator>())>::type>::type Real;
    BOOST_ASSERT_MSG(first != last, "At least one sample is required to compute mean and variance.");
    // Higham, Accuracy and Stability, equation 1.6a and 1.6b:
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

template<class ForwardIterator>
auto median(ForwardIterator first, ForwardIterator last)
{
    typedef typename std::remove_const<typename std::remove_reference<decltype(*std::declval<ForwardIterator>())>::type>::type Real;
    Real m = 0;
    return m;
}

// Mallat, "A Wavelet Tour of Signal Processing", equation 2.60:
template<class ForwardIterator>
auto total_variation(ForwardIterator first, ForwardIterator last)
{
    typedef typename std::remove_const<typename std::remove_reference<decltype(*std::declval<ForwardIterator>())>::type>::type Real;
    BOOST_ASSERT_MSG(first != last && std::next(first) != last, "At least two samples are required to compute the total variation.");
    Real tot = 0;
    auto it = first;
    Real tmp = *it;
    while (++it != last)
    {
        tot += abs(*it - tmp);
        tmp = *it;
    }
    return tot;
}

// Mallat, equation 10.4 uses the base-2 logarithm.
template<class ForwardIterator>
auto shannon_entropy(ForwardIterator first, ForwardIterator last)
{
    typedef typename std::remove_const<typename std::remove_reference<decltype(*std::declval<ForwardIterator>())>::type>::type Real;
    using std::log2;
    Real entropy = 0;
    for (auto it = first; it != last; ++it)
    {
        Real tmp = *it;
        if (tmp != 0)
        {
            entropy += tmp*log2(tmp);
        }
    }
    return -entropy;
}

template<class ForwardIterator>
auto sup_norm(ForwardIterator first, ForwardIterator last)
{
    BOOST_ASSERT_MSG(first != last, "At least one value is required to compute the sup norm.");
    typedef typename std::remove_const<typename std::remove_reference<decltype(*std::declval<ForwardIterator>())>::type>::type RealOrComplex;
    using std::abs;
    if constexpr (boost::is_complex<RealOrComplex>::value)
    {
        auto it = max_element(first, last, [](RealOrComplex a, RealOrComplex b) { return abs(b) > abs(a); });
        return abs(*it);
    }
    else
    {
        auto pair = minmax_element(first, last);
        if (abs(*pair.first) > abs(*pair.second))
        {
            return abs(*pair.first);
        }
        else
        {
            return abs(*pair.second);
        }
    }
}

template<class ForwardIterator>
size_t l0_norm(ForwardIterator first, ForwardIterator last)
{
    typedef typename std::remove_const<typename std::remove_reference<decltype(*std::declval<ForwardIterator>())>::type>::type RealOrComplex;
    size_t count = 0;
    for (auto it = first; it != last; ++it)
    {
        if (*it != RealOrComplex(0))
        {
            ++count;
        }
    }
    return count;
}

template<class ForwardIterator>
auto lp_norm(ForwardIterator first, ForwardIterator last, typename std::remove_const<typename std::remove_reference<decltype(*std::declval<ForwardIterator>())>::type>::type p)
{
    using std::pow;
    using std::is_floating_point;
    using std::isfinite;
    typedef typename std::remove_const<typename std::remove_reference<decltype(*std::declval<ForwardIterator>())>::type>::type RealOrComplex;
    if constexpr (boost::is_complex<RealOrComplex>::value)
    {
        BOOST_ASSERT_MSG(p.real() >= 0, "For p < 0, the lp norm is not a norm.");
        BOOST_ASSERT_MSG(p.imag() == 0, "For imaginary p, the lp norm is not a norm.");
        using std::norm;
        decltype(p.real()) lp = 0;
        for (auto it = first; it != last; ++it)
        {
            lp += pow(norm(*it), p.real()/2);
        }

        auto result = pow(lp, 1/p.real());
        if (!isfinite(result))
        {
            auto a = boost::math::tools::sup_norm(first, last);
            decltype(p.real()) lp = 0;
            for (auto it = first; it != last; ++it)
            {
                lp += pow(abs(*it)/a, p.real());
            }
            result = a*pow(lp, 1/p.real());
        }
        return result;
    }
    else if constexpr (is_floating_point<RealOrComplex>::value)
    {
        BOOST_ASSERT_MSG(p >= 0, "For p < 0, the lp norm is not a norm");
        RealOrComplex lp = 0;
        for (auto it = first; it != last; ++it)
        {
            lp += pow(abs(*it), p);
        }
        RealOrComplex result = pow(lp, 1/p);
        if (!isfinite(result))
        {
            RealOrComplex a = boost::math::tools::sup_norm(first, last);
            lp = 0;
            for (auto it = first; it != last; ++it)
            {
                lp += pow(abs(*it)/a, p);
            }
            result = a*pow(lp, 1/p);
        }
        return result;
    }
    else
    {
        BOOST_ASSERT_MSG(false, "Unable to determine if the input type is real or complex.");
    }
}

template<class ForwardIterator>
auto gini_coefficient(ForwardIterator first, ForwardIterator last)
{
    typedef typename std::remove_const<typename std::remove_reference<decltype(*std::declval<ForwardIterator>())>::type>::type Real;
    BOOST_ASSERT_MSG(first != last && std::next(first) != last, "Computation of the Gini coefficient requires at least two samples.");

    std::sort(first, last);

    Real i = 1;
    Real num = 0;
    Real denom = 0;
    for (auto it = first; it != last; ++it) {
        num += *it*i;
        denom += *it;
        ++i;
    }
    return ((2*num)/denom - i)/(i-2);
}


template<class ForwardIterator>
auto absolute_gini_coefficient(ForwardIterator first, ForwardIterator last)
{
    typedef typename std::remove_const<typename std::remove_reference<decltype(*std::declval<ForwardIterator>())>::type>::type RealOrComplex;
    BOOST_ASSERT_MSG(first != last && std::next(first) != last, "Computation of the Gini coefficient requires at least two samples.");

    std::sort(first, last,  [](RealOrComplex a, RealOrComplex b) { return abs(b) > abs(a); });


    decltype(abs(*first)) i = 1;
    decltype(abs(*first)) num = 0;
    decltype(abs(*first)) denom = 0;
    std::cout <<  "{";
    for (auto it = first; it != last; ++it) {
        std::cout << abs(*it) << ", ";
        decltype(abs(*first)) tmp = abs(*it);
        num += tmp*i;
        denom += tmp;
        ++i;
    }
    std::cout << "}\n";
    return ((2*num)/denom - i)/(i-2);
}



}}}
#endif
