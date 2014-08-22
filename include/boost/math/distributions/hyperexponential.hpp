//  Copyright 2014 Marco Guazzone (marco.guazzone@gmail.com)
//
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// This module implements the Hyper-Exponential distribution.
//
// References:
// - "Queueing Theory in Manufacturing Systems Analysis and Design" by H.T. Papadopolous, C. Heavey and J. Browne (Chapman & Hall/CRC, 1993)
// - http://reference.wolfram.com/language/ref/HyperexponentialDistribution.html
// - http://en.wikipedia.org/wiki/Hyperexponential_distribution
//

#ifndef BOOST_MATH_DISTRIBUTIONS_HYPEREXPONENTIAL_HPP
#define BOOST_MATH_DISTRIBUTIONS_HYPEREXPONENTIAL_HPP


#include <numeric>
#include <boost/math/distributions/complement.hpp>
#include <boost/math/distributions/detail/common_error_handling.hpp>
#include <boost/math/distributions/exponential.hpp>
#include <boost/math/policies/policy.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/math/tools/precision.hpp>
#include <boost/math/tools/roots.hpp>
//#include <boost/math/tools/tuple.hpp>
#include <boost/range/begin.hpp>
#include <boost/range/end.hpp>
#include <cstddef>
//#include <iostream>
#include <limits>
#include <utility>
#include <vector>


namespace boost { namespace math {

namespace detail {

template <typename Dist>
typename Dist::value_type generic_quantile(const Dist& dist, const typename Dist::value_type& p, const typename Dist::value_type& guess, bool comp, const char* function);

} // Namespace detail


template <typename RealT, typename PolicyT>
class hyperexponential_distribution;


namespace /*<unnamed>*/ { namespace hyperexp_detail {

template <typename T>
void normalize(std::vector<T>& v)
{
    const T sum = std::accumulate(v.begin(), v.end(), static_cast<T>(0));

    const typename std::vector<T>::iterator end = v.end();
    for (typename std::vector<T>::iterator it = v.begin();
         it != end;
         ++it)
    {
        *it /= sum;
    }
}

template <typename T>
bool iszero(T x)
{
#ifdef FP_ZERO
    return boost::math::fpclassify(x) == FP_ZERO;
    // Alternatively, we could use std::fpclassify (but this is available from ISO C99)
    //return std::fpclassify(x) == FP_ZERO;
#else
    return = ((x < 0) ? bool(-x < (std::numeric_limits<T>::min)())
                      : bool(+x < (std::numeric_limits<T>::min)()));
#endif // FP_ZERO
}

/*
template <typename RealT, typename PolicyT>
class quantile_functor
{
    public: quantile_functor(hyperexponential_distribution<RealT,PolicyT> const& dist, RealT const& p, bool complement)
    : dist_(dist),
      p_(p),
      comp_(complement)
    {
    }

    public: boost::math::tuple<RealT,RealT> operator()(RealT const& x)
    {
        const RealT c = comp_ ? boost::math::cdf(boost::math::complement(dist_, x)) : boost::math::cdf(dist_, x);
        const RealT fx = c - p_;  // Difference (cdf - value) to minimize.
        const RealT dx = comp_ ? -boost::math::pdf(dist_, x) : boost::math::pdf(dist_, x); // pdf is 1st derivative.
        // return both function evaluation difference f(x) and 1st derivative f'(x).
        return boost::math::make_tuple(fx, dx);
    }

    private: const hyperexponential_distribution<RealT,PolicyT> dist_;
    private: const RealT p_;
    private: const bool comp_;
}; // quantile_functor
*/

template <typename RealT, typename PolicyT>
bool check_probabilities(char const* function, std::vector<RealT> const& probabilities, RealT* presult, PolicyT const& pol)
{
    const std::size_t n = probabilities.size();
    RealT sum = 0;
    for (std::size_t i = 0; i < n; ++i)
    {
        if (probabilities[i] < 0
            || probabilities[i] > 1
            || !boost::math::isfinite(probabilities[i]))
        {
            *presult = policies::raise_domain_error<RealT>(function,
                                                           "The elements of parameter \"probabilities\" must be >= 0 and <= 1, but at least one of them was: %1%.",
                                                           probabilities[i],
                                                           pol);
            return false;
        }
        sum += probabilities[i];
    }

    if (!iszero(sum-RealT(1)))
    {
        *presult = policies::raise_domain_error<RealT>(function,
                                                       "The elements of parameter \"probabilities\" must sum to 1, but their sum is: %1%.",
                                                       sum,
                                                       pol);
        return false;
    }

    return true;
}

template <typename RealT, typename PolicyT>
bool check_rates(char const* function, std::vector<RealT> const& rates, RealT* presult, PolicyT const& pol)
{
    const std::size_t n = rates.size();
    for (std::size_t i = 0; i < n; ++i)
    {
        if (rates[i] <= 0
            || !boost::math::isfinite(rates[i]))
        {
            *presult = policies::raise_domain_error<RealT>(function,
                                                           "The elements of parameter \"rates\" must be > 0, but at least one of them is: %1%.",
                                                           rates[i],
                                                           pol);
            return false;
        }
    }
    return true;
}

template <typename RealT, typename PolicyT>
bool check_dist(char const* function, std::vector<RealT> const& probabilities, std::vector<RealT> const& rates, RealT* presult, PolicyT const& pol)
{
    if (probabilities.size() != rates.size())
    {
        *presult = policies::raise_domain_error<RealT>(function,
                                                       "The parameters \"probabilities\" and \"rates\" must have the same length, but their size differ by: %1%.",
                                                       std::abs(static_cast<RealT>(probabilities.size())-static_cast<RealT>(rates.size())),
                                                       pol);
        return false;
    }

    return check_probabilities(function, probabilities, presult, pol)
           && check_rates(function, rates, presult, pol);
}

template <typename RealT, typename PolicyT>
bool check_x(char const* function, RealT x, RealT* presult, PolicyT const& pol)
{
    if (x < 0 || boost::math::isnan(x))
    {
        *presult = policies::raise_domain_error<RealT>(function, "The random variable must be >= 0, but is: %1%.", x, pol);
        return false;
    }
    return true;
}

template <typename RealT, typename PolicyT>
bool check_probability(char const* function, RealT p, RealT* presult, PolicyT const& pol)
{
    if (p < 0 || p > 1 || boost::math::isnan(p))
    {
        *presult = policies::raise_domain_error<RealT>(function, "The probability be >= 0 and <= 1, but is: %1%.", p, pol);
        return false;
    }
    return true;
}

template <typename RealT, typename PolicyT>
RealT quantile_impl(hyperexponential_distribution<RealT, PolicyT> const& dist, RealT const& p, bool comp)
{
    // Don't have a closed form so try to numerically solve the inverse CDF...

    typedef typename policies::evaluation<RealT, PolicyT>::type value_type;
    typedef typename policies::normalise<PolicyT,
                                         policies::promote_float<false>,
                                         policies::promote_double<false>,
                                         policies::discrete_quantile<>,
                                         policies::assert_undefined<> >::type forwarding_policy;

    static const char* function = comp ? "boost::math::quantile(const boost::math::complemented2_type<boost::math::hyperexponential_distribution<%1%>, %1%>&)"
                                       : "boost::math::quantile(const boost::math::hyperexponential_distribution<%1%>&, %1%)";

    RealT result = 0;

    if (!check_probability(function, p, &result, PolicyT()))
    {
        return result;
    }

    const std::size_t n = dist.num_phases();
    const std::vector<RealT> probs = dist.probabilities();
    const std::vector<RealT> rates = dist.rates();

    // A possible (but inaccurate) approximation is given below, where the
    // quantile is given by the weighted sum of exponential quantiles:
    RealT guess = 0;
    if (comp)
    {
        for (std::size_t i = 0; i < n; ++i)
        {
            const exponential_distribution<RealT,PolicyT> exp(rates[i]);

            guess += probs[i]*quantile(complement(exp, p));
        }
    }
    else
    {
        for (std::size_t i = 0; i < n; ++i)
        {
            const exponential_distribution<RealT,PolicyT> exp(rates[i]);

            guess += probs[i]*quantile(exp, p);
        }
    }

    // Fast return in case the Hyper-Exponential is essentially and Exponential
    if (n == 1)
    {
        return guess;
    }

    value_type q;
    q = detail::generic_quantile(hyperexponential_distribution<RealT,forwarding_policy>(probs, rates),
                                 p,
                                 guess,
                                 comp,
                                 function);

    result = policies::checked_narrowing_cast<RealT,forwarding_policy>(q, function);

    return result;
}

}} // Namespace <unnamed>::hyperexp_detail


template <typename RealT = double, typename PolicyT = policies::policy<> >
class hyperexponential_distribution
{
    public: typedef RealT value_type;
    public: typedef PolicyT policy_type;


    public: template <typename ProbIterT, typename RateIterT>
            hyperexponential_distribution(ProbIterT prob_first, ProbIterT prob_last,
                                          RateIterT rate_first, RateIterT rate_last)
    : probs_(prob_first, prob_last),
      rates_(rate_first, rate_last)
    {
        RealT err;
        hyperexp_detail::check_dist("boost::math::hyperexponential_distribution<%1%>::hyperexponential_distribution",
                                    probs_,
                                    rates_,
                                    &err,
                                    PolicyT());
    }

    public: template <typename ProbRangeT, typename RateRangeT>
            hyperexponential_distribution(ProbRangeT prob_range, RateRangeT rate_range)
    : probs_(boost::begin(prob_range), boost::end(prob_range)),
      rates_(boost::begin(rate_range), boost::end(rate_range))
    {
        assert(probs_.size() == rates_.size());

        hyperexp_detail::normalize(probs_);

        RealT err;
        hyperexp_detail::check_dist("boost::math::hyperexponential_distribution<%1%>::hyperexponential_distribution",
                                    probs_,
                                    rates_,
                                    &err,
                                    PolicyT());
    }

    public: std::vector<RealT> probabilities() const
    {
        return probs_;
    }

    public: std::vector<RealT> rates() const
    {
        return rates_;
    }

    public: std::size_t num_phases() const
    {
        return rates_.size();
    }


    private: std::vector<RealT> probs_;
    private: std::vector<RealT> rates_;
}; // class hyperexponential_distribution


// Convenient type synonym for double.
typedef hyperexponential_distribution<double> hyperexponential;


// Range of permissible values for random variable x
template <typename RealT, typename PolicyT>
std::pair<RealT,RealT> range(hyperexponential_distribution<RealT,PolicyT> const&)
{
    if (std::numeric_limits<RealT>::has_infinity)
    {
        return std::make_pair(static_cast<RealT>(0), std::numeric_limits<RealT>::infinity()); // 0 to +inf.
    }

    return std::make_pair(static_cast<RealT>(0), tools::max_value<RealT>()); // 0 to +<max value>
}

// Range of supported values for random variable x.
// This is range where cdf rises from 0 to 1, and outside it, the pdf is zero.
template <typename RealT, typename PolicyT>
std::pair<RealT,RealT> support(hyperexponential_distribution<RealT,PolicyT> const&)
{
    return std::make_pair(tools::min_value<RealT>(), tools::max_value<RealT>()); // <min value> to +<max value>.
}

template <typename RealT, typename PolicyT>
RealT pdf(hyperexponential_distribution<RealT, PolicyT> const& dist, RealT const& x)
{
    RealT result = 0;

    if (!hyperexp_detail::check_x("boost::math::pdf(const boost::math::hyperexponential_distribution<%1%>&, %1%)", x, &result, PolicyT()))
    {
        return result;
    }

    const std::size_t n = dist.num_phases();
    const std::vector<RealT> probs = dist.probabilities();
    const std::vector<RealT> rates = dist.rates();

    for (std::size_t i = 0; i < n; ++i)
    {
        //const exponential_distribution<RealT,PolicyT> exp(rates[i]);

        //result += probs[i]*pdf(exp, x);
        result += probs[i]*rates[i]*std::exp(-rates[i]*x);
    }

    return result;
}

template <typename RealT, typename PolicyT>
RealT cdf(hyperexponential_distribution<RealT, PolicyT> const& dist, RealT const& x)
{
    RealT result = 0;

    if (!hyperexp_detail::check_x("boost::math::cdf(const boost::math::hyperexponential_distribution<%1%>&, %1%)", x, &result, PolicyT()))
    {
        return result;
    }

    const std::size_t n = dist.num_phases();
    const std::vector<RealT> probs = dist.probabilities();
    const std::vector<RealT> rates = dist.rates();

    for (std::size_t i = 0; i < n; ++i)
    {
        const exponential_distribution<RealT,PolicyT> exp(rates[i]);

        result += probs[i]*cdf(exp, x);
    }

    return result;
}

template <typename RealT, typename PolicyT>
RealT quantile(hyperexponential_distribution<RealT, PolicyT> const& dist, RealT const& p)
{
    return hyperexp_detail::quantile_impl(dist, p , false);
}

template <typename RealT, typename PolicyT>
RealT cdf(complemented2_type<hyperexponential_distribution<RealT,PolicyT>, RealT> const& c)
{
    RealT const& x = c.param;
    hyperexponential_distribution<RealT,PolicyT> const& dist = c.dist;

    RealT result = 0;

    if (!hyperexp_detail::check_x("boost::math::cdf(boost::math::complemented2_type<const boost::math::hyperexponential_distribution<%1%>&, %1%>)", x, &result, PolicyT()))
    {
        return result;
    }

    const std::size_t n = dist.num_phases();
    const std::vector<RealT> probs = dist.probabilities();
    const std::vector<RealT> rates = dist.rates();

    for (std::size_t i = 0; i < n; ++i)
    {
        const exponential_distribution<RealT,PolicyT> exp(rates[i]);

        result += probs[i]*cdf(complement(exp, x));
    }

    return result;
}


template <typename RealT, typename PolicyT>
RealT quantile(complemented2_type<hyperexponential_distribution<RealT, PolicyT>, RealT> const& c)
{
    RealT const& p = c.param;
    hyperexponential_distribution<RealT,PolicyT> const& dist = c.dist;

    return hyperexp_detail::quantile_impl(dist, p , true);
}

template <typename RealT, typename PolicyT>
RealT mean(hyperexponential_distribution<RealT, PolicyT> const& dist)
{
    RealT result = 0;

    const std::size_t n = dist.num_phases();
    const std::vector<RealT> probs = dist.probabilities();
    const std::vector<RealT> rates = dist.rates();

    for (std::size_t i = 0; i < n; ++i)
    {
        const exponential_distribution<RealT,PolicyT> exp(rates[i]);

        result += probs[i]*mean(exp);
    }

    return result;
}

template <typename RealT, typename PolicyT>
RealT variance(hyperexponential_distribution<RealT, PolicyT> const& dist)
{
    RealT result = 0;

    const std::size_t n = dist.num_phases();
    const std::vector<RealT> probs = dist.probabilities();
    const std::vector<RealT> rates = dist.rates();

    // 2(\sum_{i=1}^n \frac{p_i}{\lambda_i^2})-(\sum_{i=1}^n \frac{p_i}{\lambda_i})^2
    for (std::size_t i = 0; i < n; ++i)
    {
        result += probs[i]/(rates[i]*rates[i]);
    }

    const RealT mean = mean(dist);

    result = 2.0*result-mean*mean;

    return result;
}

template <typename RealT, typename PolicyT>
RealT skewness(hyperexponential_distribution<RealT,PolicyT> const& dist)
{
    // (6(\sum_{i=1}^n\frac{p_i}{\lambda_i^3}) - (3(2(\sum_{i=1}^n\frac{p_i}{\lambda_i^2) - (\sum_{i=1}^n\frac{p_i}{\lambda_i})^2) + (\sum_{i=1}^n\frac{p_i}{\lambda_i})^2)(\sum_{i=1}^n\frac{p_i}{\lambda_i}))/((\sum_{i=1}^n\frac{p_i}{\lambda_i^2})-(\sum_{i=1}^n\frac{p_i}{\lambda_i})^2)^(3/2)

    const std::size_t n = dist.num_phases();
    const std::vector<RealT> probs = dist.probabilities();
    const std::vector<RealT> rates = dist.rates();

    RealT s1 = 0; // \sum_{i=1}^n \frac{p_i}{\lambda_i}
    RealT s2 = 0; // \sum_{i=1}^n \frac{p_i}{\lambda_i^2}
    RealT s3 = 0; // \sum_{i=1}^n \frac{p_i}{\lambda_i^3}
    for (std::size_t i = 0; i < n; ++i)
    {
        const RealT p = probs[i];
        const RealT r = rates[i];
        const RealT r2 = r*r;
        const RealT r3 = r2*r;

        s1 += p/r;
        s2 += p/r2;
        s3 += p/r3;
    }

    const RealT s1s1 = s1*s1;

    const RealT num = (6.0*s3 - (3.0*(2.0*s2 - s1s1) + s1s1)*s1);
    const RealT den = (2.0*s2 - s1s1);

    return num/std::pow(den, 1.5);
}

template <typename RealT, typename PolicyT>
RealT kurtosis(hyperexponential_distribution<RealT,PolicyT> const& dist)
{
    // (24(\sum_{i=1}^n\frac{p_i}{\lambda_i^4})-24(\sum_{i=1}^n\frac{p_i}{\lambda_i^3})(\sum_{i=1}^n\frac{p_i}{\lambda_i})+3(2(2(\sum_{i=1}^n\frac{p_i}{\lambda_i^2})-(\sum_{i=1}^n\frac{p_i}{\lambda_i})^2) + (\sum_{i=1}^n\frac{p_i}{\lambda_i})^2)(\sum_{i=1}^n\frac{p_i}{\lambda_i})^2)/((\sum_{i=1}^n\frac{p_i}{\lambda_i^2})-(\sum_{i=1}^n\frac{p_i}{\lambda_i})^2)^2

    const std::size_t n = dist.num_phases();
    const std::vector<RealT> probs = dist.probabilities();
    const std::vector<RealT> rates = dist.rates();

    RealT s1 = 0; // \sum_{i=1}^n \frac{p_i}{\lambda_i}
    RealT s2 = 0; // \sum_{i=1}^n \frac{p_i}{\lambda_i^2}
    RealT s3 = 0; // \sum_{i=1}^n \frac{p_i}{\lambda_i^3}
    RealT s4 = 0; // \sum_{i=1}^n \frac{p_i}{\lambda_i^4}
    for (std::size_t i = 0; i < n; ++i)
    {
        const RealT p = probs[i];
        const RealT r = rates[i];
        const RealT r2 = r*r;
        const RealT r3 = r2*r;
        const RealT r4 = r3*r;

        s1 += p/r;
        s2 += p/r2;
        s3 += p/r3;
        s4 += p/r4;
    }

    const RealT s1s1 = s1*s1;

    const RealT num = (24.0*s4 - 24.0*s3*s1 + 3.0*(2.0*(2.0*s2 - s1s1) + s1s1)*s1s1);
    const RealT den = (2.0*s2 - s1s1);

    return num/(den*den);
}

template <typename RealT, typename PolicyT>
RealT kurtosis_excess(hyperexponential_distribution<RealT,PolicyT> const& dist)
{
    return kurtosis(dist) - 3;
}

template <typename RealT, typename PolicyT>
RealT mode(hyperexponential_distribution<RealT,PolicyT> const& dist)
{
    return 0;
}

}} // namespace boost::math

// This include must be at the end, *after* the accessors
// for this distribution have been defined, in order to
// keep compilers that support two-phase lookup happy.
#include <boost/math/distributions/detail/derived_accessors.hpp>
#include <boost/math/distributions/detail/generic_quantile.hpp>

#endif // BOOST_MATH_DISTRIBUTIONS_HYPEREXPONENTIAL
