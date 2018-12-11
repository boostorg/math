//  (C) Copyright Nick Thompson 2018.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_TOOLS_SIGNAL_STATISTICS_HPP
#define BOOST_MATH_TOOLS_SIGNAL_STATISTICS_HPP

#include <algorithm>
#include <iterator>
#include <boost/type_traits.hpp>
#include <boost/assert.hpp>
#include <boost/multiprecision/detail/number_base.hpp>


namespace boost{ namespace math{ namespace tools {

template<class RandomAccessIterator>
auto absolute_median(RandomAccessIterator first, RandomAccessIterator last)
{
    using std::abs;
    using RealOrComplex = typename std::iterator_traits<RandomAccessIterator>::value_type;
    size_t num_elems = std::distance(first, last);
    BOOST_ASSERT_MSG(num_elems > 0, "The median of a zero-length vector is undefined.");
    auto comparator = [](RealOrComplex a, RealOrComplex b) { return abs(a) < abs(b);};
    if (num_elems & 1)
    {
        auto middle = first + (num_elems - 1)/2;
        std::nth_element(first, middle, last, comparator);
        return abs(*middle);
    }
    else
    {
        auto middle = first + num_elems/2 - 1;
        std::nth_element(first, middle, last, comparator);
        std::nth_element(middle, middle+1, last, comparator);
        return (abs(*middle) + abs(*(middle+1)))/abs(static_cast<RealOrComplex>(2));
    }
}

template<class RandomAccessContainer>
inline auto absolute_median(RandomAccessContainer & v)
{
    return absolute_median(v.begin(), v.end());
}


template<class ForwardIterator>
auto shannon_entropy(ForwardIterator first, ForwardIterator last)
{
    using Real = typename std::iterator_traits<ForwardIterator>::value_type;
    using std::log;
    Real entropy = 0;
    for (auto it = first; it != last; ++it)
    {
        if (*it != 0)
        {
            entropy += (*it)*log(*it);
        }
    }
    return -entropy;
}

template<class Container>
inline auto shannon_entropy(Container const & v)
{
    return shannon_entropy(v.cbegin(), v.cend());
}


template<class ForwardIterator>
auto shannon_cost(ForwardIterator first, ForwardIterator last)
{
    using Real = typename std::iterator_traits<ForwardIterator>::value_type;
    using std::log;
    Real cost = 0;
    for (auto it = first; it != last; ++it)
    {
        if (*it != 0)
        {
            Real tmp = abs(*it);
            cost += tmp*tmp*log(tmp*tmp);
        }
    }
    return -cost;
}

template<class Container>
inline auto shannon_cost(Container const & v)
{
    return shannon_cost(v.cbegin(), v.cend());
}


template<class ForwardIterator>
auto absolute_gini_coefficient(ForwardIterator first, ForwardIterator last)
{
    using RealOrComplex = typename std::iterator_traits<ForwardIterator>::value_type;
    BOOST_ASSERT_MSG(first != last && std::next(first) != last, "Computation of the Gini coefficient requires at least two samples.");

    std::sort(first, last,  [](RealOrComplex a, RealOrComplex b) { return abs(b) > abs(a); });


    decltype(abs(*first)) i = 1;
    decltype(abs(*first)) num = 0;
    decltype(abs(*first)) denom = 0;
    for (auto it = first; it != last; ++it)
    {
        decltype(abs(*first)) tmp = abs(*it);
        num += tmp*i;
        denom += tmp;
        ++i;
    }

    // If the l1 norm is zero, all elements are zero, so every element is the same.
    if (denom == 0)
    {
        decltype(abs(*first)) zero = 0;
        return zero;
    }
    return ((2*num)/denom - i)/(i-2);
}

template<class RandomAccessContainer>
inline auto absolute_gini_coefficient(RandomAccessContainer & v)
{
    return absolute_gini_coefficient(v.begin(), v.end());
}

// The Hoyer sparsity measure is defined in:
// https://arxiv.org/pdf/0811.4706.pdf
template<class ForwardIterator>
auto hoyer_sparsity(const ForwardIterator first, const ForwardIterator last)
{
    using RealIntOrComplex = typename std::iterator_traits<ForwardIterator>::value_type;
    using std::abs;
    using std::sqrt;
    BOOST_ASSERT_MSG(first != last, "Computation of the Hoyer sparsity requires at least one sample.");

    decltype(abs(*first)) l1 = 0;
    decltype(abs(*first)) l2 = 0;
    decltype(abs(*first)) n = 0;
    for (auto it = first; it != last; ++it)
    {
        decltype(abs(*first)) tmp = abs(*it);
        l1 += tmp;
        l2 += tmp*tmp;
        n += 1;
    }
    if constexpr (std::is_integral<RealIntOrComplex>::value)
    {
        double rootn = sqrt(n);
        return (rootn - l1/sqrt(l2) )/ (rootn - 1);
    }
    else
    {
        decltype(abs(*first)) rootn = sqrt(n);
        return (rootn - l1/sqrt(l2) )/ (rootn - 1);
    }
}

template<class Container>
inline auto hoyer_sparsity(Container const & v)
{
    return hoyer_sparsity(v.cbegin(), v.cend());
}


template<class Container>
auto oracle_snr(Container const & signal, Container const & noise)
{
    using Real = typename Container::value_type;
    BOOST_ASSERT_MSG(signal.size() == noise.size(), "Signal and noise must be have the same number of elements.");
    if constexpr (std::is_integral<Real>::value)
    {
        double numerator = 0;
        double denominator = 0;
        for (size_t i = 0; i < signal.size(); ++i)
        {
            numerator += signal[i]*signal[i];
            denominator += noise[i]*noise[i];
        }
        if (numerator == 0 && denominator == 0)
        {
            return std::numeric_limits<double>::quiet_NaN();
        }
        if (denominator == 0)
        {
            return std::numeric_limits<double>::infinity();
        }
        return numerator/denominator;
    }
    else if constexpr (boost::is_complex<Real>::value ||
                       boost::multiprecision::number_category<Real>::value == boost::multiprecision::number_kind_complex)

    {
        using std::norm;
        typename Real::value_type numerator = 0;
        typename Real::value_type denominator = 0;
        for (size_t i = 0; i < signal.size(); ++i)
        {
            numerator += norm(signal[i]);
            denominator += norm(noise[i]);
        }
        if (numerator == 0 && denominator == 0)
        {
            return std::numeric_limits<typename Real::value_type>::quiet_NaN();
        }
        if (denominator == 0)
        {
            return std::numeric_limits<typename Real::value_type>::infinity();
        }

        return numerator/denominator;
    }
    else
    {
        Real numerator = 0;
        Real denominator = 0;
        for (size_t i = 0; i < signal.size(); ++i)
        {
            numerator += signal[i]*signal[i];
            denominator += noise[i]*noise[i];
        }
        if (numerator == 0 && denominator == 0)
        {
            return std::numeric_limits<Real>::quiet_NaN();
        }
        if (denominator == 0)
        {
            return std::numeric_limits<Real>::infinity();
        }

        return numerator/denominator;
    }
}

// Follows the definition of SNR given in Mallat, A Wavelet Tour of Signal Processing, equation 11.16.
template<class Container>
auto oracle_snr_db(Container const & signal, Container const & noise)
{
    using std::log10;
    return 10*log10(oracle_snr(signal, noise));
}

// Of course since we have an oracle snr estimator, we should have an snr estimator not requiring oracle data.
// The M2M4 estimator is reputed to be quite good, as is the SVR measure.
// A good reference is:
// D. R. Pauluzzi and N. C. Beaulieu, "A comparison of SNR estimation techniques for the AWGN channel," IEEE Trans. Communications, Vol. 48, No. 10, pp. 1681-1691, 2000.
// A nice python implementation:
// https://github.com/gnuradio/gnuradio/blob/master/gr-digital/examples/snr_estimators.py
// However, we have not implemented kurtosis and kurtosis, which is required of the method.


}}}
#endif
