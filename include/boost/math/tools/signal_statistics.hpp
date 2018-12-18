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
#include <boost/math/tools/univariate_statistics.hpp>


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
auto absolute_gini_coefficient(ForwardIterator first, ForwardIterator last)
{
    using std::abs;
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
    return ((2*num)/denom - i)/(i-1);
}

template<class RandomAccessContainer>
inline auto absolute_gini_coefficient(RandomAccessContainer & v)
{
    return boost::math::tools::absolute_gini_coefficient(v.begin(), v.end());
}

template<class ForwardIterator>
auto sample_absolute_gini_coefficient(ForwardIterator first, ForwardIterator last)
{
    size_t n = std::distance(first, last);
    return n*boost::math::tools::absolute_gini_coefficient(first, last)/(n-1);
}

template<class RandomAccessContainer>
inline auto sample_absolute_gini_coefficient(RandomAccessContainer & v)
{
    return boost::math::tools::sample_absolute_gini_coefficient(v.begin(), v.end());
}


// The Hoyer sparsity measure is defined in:
// https://arxiv.org/pdf/0811.4706.pdf
template<class ForwardIterator>
auto hoyer_sparsity(const ForwardIterator first, const ForwardIterator last)
{
    using RealIntOrComplex = typename std::iterator_traits<ForwardIterator>::value_type;
    using std::abs;
    using std::sqrt;
    BOOST_ASSERT_MSG(first != last && std::next(first) != last, "Computation of the Hoyer sparsity requires at least two samples.");

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
    return boost::math::tools::hoyer_sparsity(v.cbegin(), v.cend());
}


template<class Container>
auto oracle_snr(Container const & signal, Container const & noisy_signal)
{
    using Real = typename Container::value_type;
    BOOST_ASSERT_MSG(signal.size() == noisy_signal.size(),
                     "Signal and noisy_signal must be have the same number of elements.");
    if constexpr (std::is_integral<Real>::value)
    {
        double numerator = 0;
        double denominator = 0;
        for (size_t i = 0; i < signal.size(); ++i)
        {
            numerator += signal[i]*signal[i];
            denominator += (noisy_signal[i] - signal[i])*(noisy_signal[i] - signal[i]);
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
            denominator += norm(noisy_signal[i] - signal[i]);
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
            denominator += (signal[i] - noisy_signal[i])*(signal[i] - noisy_signal[i]);
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

template<class Container>
auto mean_invariant_oracle_snr(Container const & signal, Container const & noisy_signal)
{
    using Real = typename Container::value_type;
    BOOST_ASSERT_MSG(signal.size() == noisy_signal.size(), "Signal and noisy signal must be have the same number of elements.");

    Real mu = boost::math::tools::mean(signal);
    Real numerator = 0;
    Real denominator = 0;
    for (size_t i = 0; i < signal.size(); ++i)
    {
        Real tmp = signal[i] - mu;
        numerator += tmp*tmp;
        denominator += (signal[i] - noisy_signal[i])*(signal[i] - noisy_signal[i]);
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

template<class Container>
auto mean_invariant_oracle_snr_db(Container const & signal, Container const & noisy_signal)
{
    using std::log10;
    return 10*log10(boost::math::tools::mean_invariant_oracle_snr(signal, noisy_signal));
}


// Follows the definition of SNR given in Mallat, A Wavelet Tour of Signal Processing, equation 11.16.
template<class Container>
auto oracle_snr_db(Container const & signal, Container const & noisy_signal)
{
    using std::log10;
    return 10*log10(boost::math::tools::oracle_snr(signal, noisy_signal));
}

// A good reference on the M2M4 estimator:
// D. R. Pauluzzi and N. C. Beaulieu, "A comparison of SNR estimation techniques for the AWGN channel," IEEE Trans. Communications, Vol. 48, No. 10, pp. 1681-1691, 2000.
// A nice python implementation:
// https://github.com/gnuradio/gnuradio/blob/master/gr-digital/examples/snr_estimators.py

template<class Container>
auto m2m4_snr_estimator(Container const & noisy_signal,  typename Container::value_type estimated_signal_kurtosis=1, typename Container::value_type estimated_noise_kurtosis=3)
{
    BOOST_ASSERT_MSG(estimated_signal_kurtosis > 0, "The estimated signal kurtosis must be positive");
    BOOST_ASSERT_MSG(estimated_noise_kurtosis > 0, "The estimated noise kurtosis must be positive.");
    using Real = typename Container::value_type;
    using std::sqrt;
    if constexpr (std::is_floating_point<Real>::value ||
       boost::multiprecision::number_category<Real>::value == boost::multiprecision::number_kind_floating_point)
    {
        // If we first eliminate N, we obtain the quadratic equation:
        // (ka+kw-6)S^2 + 2M2(3-kw)S + kw*M2^2 - M4 = 0 =: a*S^2 + bs*N + cs = 0
        // If we first eliminate S, we obtain the quadratic equation:
        // (ka+kw-6)N^2 + 2M2(3-ka)N + ka*M2^2 - M4 = 0 =: a*N^2 + bn*N + cn = 0
        // I believe these equations are totally independent quadratics;
        // if one has a complex solution it is not necessarily the case that the other must also.
        // However, I can't prove that, so there is a chance that this does unnecessary work.
        // Future improvements: There are algorithms which can solve quadratics much more effectively than the naive implementation found here.
        // See: https://stackoverflow.com/questions/48979861/numerically-stable-method-for-solving-quadratic-equations/50065711#50065711
        auto [M1, M2, M3, M4] = boost::math::tools::first_four_moments(noisy_signal);
        if (M4 == 0)
        {
            // The signal is constant. There is no noise:
            return std::numeric_limits<Real>::infinity();
        }
        // Change to notation in Pauluzzi, equation 41:
        auto kw = estimated_noise_kurtosis;
        auto ka = estimated_signal_kurtosis;
        // A common case, since it's the default:
        Real a = (ka+kw-6);
        Real bs = 2*M2*(3-kw);
        Real cs = kw*M2*M2 - M4;
        Real bn = 2*M2*(3-ka);
        Real cn = ka*M2*M2 - M4;
        Real N, S;
        if(kw == 3)
        {
            if (ka == 3)
            {
                // When ka = kw = 3, then either the system is inconsistent, or the system does not have a unique solution:
                return std::numeric_limits<Real>::quiet_NaN();
            }
            Real Ssq = -cs/a;
            if (Ssq < 0)
            {
                Real radicand = bn*bn - 4*a*cn;
                if (radicand < 0)
                {
                    return std::numeric_limits<Real>::quiet_NaN();
                }
                N = (-bn + sqrt(radicand))/(2*a);
                if (N < 0)
                {
                    N =  (-bn - sqrt(radicand))/(2*a);
                    if (N < 0)
                    {
                        return std::numeric_limits<Real>::quiet_NaN();
                    }
                    S = M2 - N;
                    if (S < 0)
                    {
                        return std::numeric_limits<Real>::quiet_NaN();
                    }
                    return S/N;
                }
                S = M2 - N;
                if (S < 0)
                {
                    return std::numeric_limits<Real>::quiet_NaN();
                }
                return S/N;
            }
            S = sqrt(Ssq);
            N = M2 - S;
            if (N < 0)
            {
                return std::numeric_limits<Real>::quiet_NaN();
            }
            return S/N;
        }

        // Maybe I should look for some very small distance from 6, but . . .
        if (ka+kw == 6)
        {
            // In this case we don't need to solve a quadratic equation:
            S = -cs/bs;
            N = -cn/bn;
            if (S/N < 0)
            {
                return std::numeric_limits<Real>::quiet_NaN();
            }
            return S/N;
        }

        // The special cases have been taken care of.
        // Now we must resort to solving a full quadratic.
        Real radicand = bs*bs - 4*a*cs;
        if (radicand < 0)
        {
            // See if we have a solution for N:
            radicand = bn*bn - 4*a*cn;
            if (radicand < 0)
            {
                // Both S and N are complex:
                return std::numeric_limits<Real>::quiet_NaN();
            }
            // N is real. Can it be made positive?
            N = (-bn + sqrt(radicand))/(2*a);
            if (N < 0)
            {
                N = (-bn - sqrt(radicand))/(2*a);
                if (N < 0)
                {
                    return std::numeric_limits<Real>::quiet_NaN();
                }
            }
            S = M2 - N;
            if (S < 0)
            {
                return std::numeric_limits<Real>::quiet_NaN();
            }
            return S/N;
        }

        S = (-bs + sqrt(radicand))/(2*a);
        if (S < 0)
        {
            S = (-bs - sqrt(radicand))/(2*a);
            if (S < 0)
            {
                return std::numeric_limits<Real>::quiet_NaN();
            }
        }
        N = M2 - S;
        if (N < 0)
        {
            return std::numeric_limits<Real>::quiet_NaN();
        }
        return S/N;
    }
    else
    {
        BOOST_ASSERT_MSG(false, "The M2M4 estimator has not been implemented for this type.");
    }
}

template<class Container>
auto m2m4_snr_estimator_db(Container const & noisy_signal,  typename Container::value_type estimated_signal_kurtosis=1, typename Container::value_type estimated_noise_kurtosis=3)
{
    using std::log10;
    return 10*log10(m2m4_snr_estimator(noisy_signal, estimated_signal_kurtosis, estimated_noise_kurtosis));
}

}}}
#endif
