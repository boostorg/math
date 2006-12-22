// boost\math\distributions\bernoulli.hpp

// Copyright John Maddock 2006.
// Copyright Paul A. Bristow 2006.

// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

// http://en.wikipedia.org/wiki/bernoulli_distribution

// bernoulli distribution is the discrete probability distribution of
// the number (k) of successes, in a single Bernoulli trials.
// It is a version of the binomial distribution when n = 1.

// But note that the bernoulli distribution
// (like others including the poisson, binomial & negative binomial)
// is strictly defined as a discrete function: only integral values of k are envisaged.
// However because of the method of calculation using a continuous gamma function,
// it is convenient to treat it as if a continous function,
// and permit non-integral values of k.
// To enforce the strict mathematical model, users should use floor or ceil functions
// on k outside this function to ensure that k is integral.

#ifndef BOOST_MATH_SPECIAL_BERNOULLI_HPP
#define BOOST_MATH_SPECIAL_BERNOULLI_HPP

#pragma warning (disable: 4180) // temporary.

#include <boost/math/special_functions/beta.hpp> // for incomplete beta.
#include <boost/math/distributions/complement.hpp> // complements
#include <boost/math/distributions/detail/common_error_handling.hpp> // error checks
#include <boost/math/special_functions/fpclassify.hpp> // isnan.
#include <boost/math/tools/roots.hpp> // for root finding.

#include <utility>

#if defined (BOOST_MSVC) && defined(BOOST_MATH_THROW_ON_DOMAIN_ERROR)
#pragma warning (disable: 4180)
#  pragma warning(push)
//#  pragma warning(disable: 4702) // unreachable code
// in domain_error_imp in error_handling
#endif

namespace boost
{
  namespace math
  {
    namespace bernoulli_detail
    {
      // Common error checking routines for bernoulli distribution functions:
      template <class RealType>
      inline bool check_success_fraction(const char* function, const RealType& p, RealType* result)
      {
        if(!(boost::math::isfinite)(p) || (p < 0) || (p > 1))
        {
          *result = tools::domain_error<RealType>(
            function,
            "Success fraction argument is %1%, but must be >= 0 and <= 1 !", p);
          return false;
        }
        return true;
      }
      template <class RealType>
      inline bool check_dist(const char* function, const RealType& p, RealType* result)
      {
        return check_success_fraction(function, p, result);

      }
      template <class RealType>
      bool check_dist_and_k(const char* function, const RealType& p, RealType k, RealType* result)
      {
        if(check_dist(function, p, result) == false)
        {
          return false;
        }
        if((k < 0) || !(boost::math::isfinite)(k))
        {
          *result = tools::domain_error<RealType>(
            function,
            "Number of successes argument is %1%, but must be 0 or 1 !", k);
          return false;
        }
        if(k > 1)
        {
          *result = tools::domain_error<RealType>(
            function,
            "Number of successes argument is %1%, but must be 0 or 1 !", k);
          return false;
        }
        return true;
      }
      template <class RealType>
      inline bool check_dist_and_prob(const char* function, RealType p, RealType prob, RealType* result)
      {
        if(check_dist(function, p, result) && detail::check_probability(function, prob, result) == false)
        {
          return false;
        }
        return true;
      }
    } // namespace bernoulli_detail


    template <class RealType = double>
    class bernoulli_distribution
    {
    public:
      typedef RealType value_type;

      bernoulli_distribution(RealType p) : m_p(p)
      {
        RealType result; // of checks.
        bernoulli_detail::check_dist(
          BOOST_CURRENT_FUNCTION,
          m_p,
          &result);
      } // bernoulli_distribution constructor.

      RealType success_fraction() const
      { // Probability.
        return m_p;
      }

      enum interval_type
      { //
        clopper_pearson_exact_interval,
        jeffreys_prior_interval
      }; // class bernoulli_distribution

      // Estimation of the success fraction parameter.
      // The best estimate is actually simply successes/trials,
      // these functions are used
      // to obtain confidence intervals for the success fraction.
      //
      static RealType estimate_lower_bound_on_p(
        RealType trials,
        RealType successes,
        RealType probability,
        interval_type t = clopper_pearson_exact_interval)
      {
        RealType result; // of checks.
        if(false == bernoulli_detail::check_dist_and_k(
          BOOST_CURRENT_FUNCTION, trials, RealType(0), successes, &result)
          &&
          bernoulli_detail::check_dist_and_prob(
          BOOST_CURRENT_FUNCTION, trials, RealType(0), probability, &result))
        { return result; }

        if(successes == 0)
        {
          return 0;
        }

        // NOTE!!! The Clopper Pearson formula uses "successes" not
        // "successes+1" as usual to get the lower bound,
        // see http://www.itl.nist.gov/div898/handbook/prc/section2/prc241.htm
        return (t == clopper_pearson_exact_interval) ? ibeta_inv(successes, trials - successes + 1, probability)
          : ibeta_inv(successes + 0.5f, trials - successes + 0.5f, probability);
      }
      static RealType estimate_upper_bound_on_p(
        RealType trials,
        RealType successes,
        RealType probability,
        interval_type t = clopper_pearson_exact_interval)
      {
        // Error checks:
        RealType result;
        if(false == bernoulli_detail::check_dist_and_k(
          BOOST_CURRENT_FUNCTION, trials, RealType(0), successes, &result)
          &&
          bernoulli_detail::check_dist_and_prob(
          BOOST_CURRENT_FUNCTION, trials, RealType(0), probability, &result))
        { return result; }

        if(trials == successes)
          return 1;

        return (t == clopper_pearson_exact_interval) ? ibetac_inv(successes + 1, trials - successes, probability)
          : ibetac_inv(successes + 0.5f, trials - successes + 0.5f, probability);
      }
      // Estimate number of trials parameter:
      //
      // "How many trials do I need to be P% sure of seeing k events?"
      //    or
      // "How many trials can I have to be P% sure of seeing fewer than k events?"
      //
      static RealType estimate_minimum_number_of_trials(
        RealType k,     // number of events
        RealType p,     // success fraction
        RealType alpha) // risk level
      {
        // Error checks:
        RealType result;
        if(false == bernoulli_detail::check_dist_and_k(
          BOOST_CURRENT_FUNCTION, k, p, k, &result)
          &&
          bernoulli_detail::check_dist_and_prob(
          BOOST_CURRENT_FUNCTION, k, p, alpha, &result))
        { return result; }

        result = ibetac_invb(k + 1, p, alpha);  // returns n - k
        return result + k;
      }

      static RealType estimate_maximum_number_of_trials(
        RealType k,     // number of events
        RealType p,     // success fraction
        RealType alpha) // risk level
      {
        // Error checks:
        RealType result;
        if(false == bernoulli_detail::check_dist_and_k(
          BOOST_CURRENT_FUNCTION, k, p, k, &result)
          &&
          bernoulli_detail::check_dist_and_prob(
          BOOST_CURRENT_FUNCTION, k, p, alpha, &result))
        { return result; }

        result = ibeta_invb(k + 1, p, alpha);  // returns n - k
        return result + k;
      }

    private:
      RealType m_p; // success_fraction
    }; // template <class RealType> class bernoulli_distribution

    typedef bernoulli_distribution<double> bernoulli;

    template <class RealType>
    const std::pair<RealType, RealType> range(const bernoulli_distribution<RealType>& dist)
    { // Range of permissible values for random variable k = {0, 1}.
      using boost::math::tools::max_value;
      return std::pair<RealType, RealType>(0, 1);
    }

    template <class RealType>
    const std::pair<RealType, RealType> support(const bernoulli_distribution<RealType>& dist)
    { // Range of supported values for random variable k = {0, 1}.
      // This is range where cdf rises from 0 to 1, and outside it, the pdf is zero.
      return std::pair<RealType, RealType>(0, 1);
    }

    template <class RealType>
    inline RealType mean(const bernoulli_distribution<RealType>& dist)
    { // Mean of bernoulli distribution = p (n = 1).
      return  dist.success_fraction();
    } // mean

    template <class RealType>
    inline RealType median(const bernoulli_distribution<RealType>& dist)
    { // Median of bernoulli distribution is not defined.
      return tools::domain_error<RealType>(BOOST_CURRENT_FUNCTION, "Median is not implemented, result is %1%!", std::numeric_limits<RealType>::quiet_NaN());
    } // median

    template <class RealType>
    inline RealType variance(const bernoulli_distribution<RealType>& dist)
    { // Variance of bernoulli distribution =p * q.
      return  dist.success_fraction() * (1 - dist.success_fraction());
    } // variance

    template <class RealType>
    RealType pdf(const bernoulli_distribution<RealType>& dist, const RealType k)
    { // Probability Density/Mass Function.
      BOOST_FPU_EXCEPTION_GUARD
      // Error check:
      RealType result; // of checks.
      if(false == bernoulli_detail::check_dist_and_k(
        BOOST_CURRENT_FUNCTION,
        dist.success_fraction(),
        k,
        &result))
      {
        return result;
      }
      // Assume k is integral.
      if (k == 1)
      {
        return dist.success_fraction(); // p
      }
      else if (k == 0)
      {
        return 1 - dist.success_fraction(); // q = 1-p
      }
      else
      {
        return 0;
      }
    } // pdf

    template <class RealType>
    RealType cdf(const bernoulli_distribution<RealType>& dist, const RealType k)
    { // Cumulative Distribution Function Bernoulli.
      RealType p = dist.success_fraction();
      // Error check:
      RealType result;
      if(false == bernoulli_detail::check_dist_and_k(
        BOOST_CURRENT_FUNCTION,
        p,
        k,
        &result))
      {
        return result;
      }
      if (k < 0)
      {
        return 0;
      }
      else if (k < 1)
      { // includes k == 0.
        1 - p;
      }
      else
      { // k >= 1
        return 1;
      }
    } // bernoulli cdf

    template <class RealType>
    RealType cdf(const complemented2_type<bernoulli_distribution<RealType>, RealType>& c)
    { // Complemented Cumulative Distribution Function bernoulli.
      RealType const& k = c.param;
      bernoulli_distribution<RealType> const& dist = c.dist;
      RealType p = dist.success_fraction();
      // Error checks:
      RealType result;
      if(false == bernoulli_detail::check_dist_and_k(
        BOOST_CURRENT_FUNCTION,
        p,
        k,
        &result))
      {
        return result;
      }
      if (k < 0)
      {
        return 1;
      }
      else if (k < 1)
      { // includes k == 0.
        p;
      }
      else
      { // k >= 1
        return 0;
      }
    } // bernoulli cdf complement

    namespace detail
    {
      template <class RealType>
      struct bernoulli_functor
      {
        bernoulli_functor(const bernoulli_distribution<RealType>& d, const RealType& target, bool c = false)
          : dist(d), t(target), complement(c) {}

        RealType operator()(const RealType k)
        {
          if(k >= dist.trials())
            return 1; // any positive value will do.
          return complement ? t - cdf(boost::math::complement(dist, k)) : cdf(dist, k) - t;
        }
      private:
        const bernoulli_distribution<RealType>& dist;
        RealType t;
        bool complement;
      }; // struct bernoulli_functor
    } // namespace detail

    template <class RealType>
    RealType quantile(const bernoulli_distribution<RealType>& dist, const RealType& p)
    { // Quantile or Percent Point bernoulli function.
      // Return the number of expected successes k,
      // for a given probability p.
      //
      // Error checks:
      using namespace std;  // ADL of std names
      RealType result;
      if(false == bernoulli_detail::check_dist_and_prob(
        BOOST_CURRENT_FUNCTION,
        dist.trials(),
        dist.success_fraction(),
        p,
        &result))
      {
        return result;
      }

      // Special cases:
      //
      if(p == 0)
      {  // There may actually be no answer to this question,
        // since the probability of zero successes may be non-zero,
        // but zero is the best we can do:
        return 0;
      }
      if(p == 1)
      {  // Probability of n or fewer successes is always one,
        // so n is the most sensible answer here:
        return dist.trials();
      }
      if (p <= pow(1 - dist.success_fraction(), dist.trials()))
      { // p <= pdf(dist, 0) == cdf(dist, 0)
        return 0; // So the only reasonable result is zero.
      } // And root finder would fail otherwise.

      // Solve for quantile numerically:
      //
      detail::bernoulli_functor<RealType> f(dist, p);
      tools::eps_tolerance<RealType> tol(tools::digits<RealType>());
      boost::uintmax_t max_iter = 1000;
      std::pair<RealType, RealType> r = tools::bracket_and_solve_root(
        f,
        dist.trials() / 2,
        static_cast<RealType>(8),
        true,
        tol,
        max_iter);
      if(max_iter >= 1000)
        tools::logic_error<RealType>(BOOST_CURRENT_FUNCTION, "Unable to locate the root within a reasonable number of iterations, closest approximation so far was %1%", r.first);
      // return centre point of range found:
      return r.first + (r.second - r.first) / 2;
    } // quantile

    template <class RealType>
    RealType quantile(const complemented2_type<bernoulli_distribution<RealType>, RealType>& c)
    { // Quantile or Percent Point bernoulli function.
      // Return the number of expected successes k for a given
      // complement of the probability q.
      //
      // Error checks:
      RealType q = c.param;
      const bernoulli_distribution<RealType>& dist = c.dist;
      RealType result;
      if(false == bernoulli_detail::check_dist_and_prob(
        BOOST_CURRENT_FUNCTION,
        dist.trials(),
        dist.success_fraction(),
        q,
        &result))
      {
        return result;
      }

      // Special cases:
      //
      if(q == 1)
      {  // There may actually be no answer to this question,
        // since the probability of zero successes may be non-zero,
        // but zero is the best we can do:
        return 0;
      }
      if(q == 0)
      {  // Probability of greater than n successes is always zero,
        // so n is the most sensible answer here:
        return dist.trials();
      }

      if (-q <= powm1(1 - dist.success_fraction(), dist.trials()))
      { // // q <= cdf(complement(dist, 0)) == pdf(dist, 0)
        return 0; // So the only reasonable result is zero.
      } // And root finder would fail otherwise.

      // Need to consider the case where
      //
      // Solve for quantile numerically:
      //
      detail::bernoulli_functor<RealType> f(dist, q, true);
      tools::eps_tolerance<RealType> tol(tools::digits<RealType>());
      boost::uintmax_t max_iter = 1000;
      std::pair<RealType, RealType> r = tools::bracket_and_solve_root(
        f,
        dist.trials() / 2,
        static_cast<RealType>(8),
        true,
        tol,
        max_iter);
      if(max_iter >= 1000)
        tools::logic_error<RealType>(BOOST_CURRENT_FUNCTION, "Unable to locate the root within a reasonable number of iterations, closest approximation so far was %1%", r.first);
      // return centre point of range found:
      return r.first + (r.second - r.first) / 2;
    } // quantile

    template <class RealType>
    inline RealType mode(const bernoulli_distribution<RealType>& dist)
    {
      return (dist.success_fraction() <= 0.5) ? 0 : 1; // p = 0.5 can be 0 or 1
    }

    template <class RealType>
    inline RealType skewness(const bernoulli_distribution<RealType>& dist)
    {
      using namespace std;; // Aid ADL for sqrt.
      RealType p = dist.success_fraction();
      RealType q = 1 - dist.success_fraction();
      RealType result = (q - p) / sqrt(p * q);
      return (1 - 2 * p) / sqrt(p * (1 - p));
    }

    template <class RealType>
    inline RealType kurtosis_excess(const bernoulli_distribution<RealType>& dist)
    {
      RealType p = dist.success_fraction();
      return (6 * p * p - 6 * p + 1) / (p * (1 - p)) - 3;
    }

    template <class RealType>
    inline RealType kurtosis(const bernoulli_distribution<RealType>& dist)
    {
      RealType p = dist.success_fraction();
      return (6 * p * p - 6 * p + 1) / (p * (1 - p));
    }

  } // namespace math
} // namespace boost

// This include must be at the end, *after* the accessors
// for this distribution have been defined, in order to
// keep compilers that support two-phase lookup happy.
#include <boost/math/distributions/detail/derived_accessors.hpp>

#endif // BOOST_MATH_SPECIAL_BERNOULLI_HPP

