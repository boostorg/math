// boost\math\special_functions\negative_binomial.hpp

// Copyright Paul A. Bristow 2006.

// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

// http://en.wikipedia.org/wiki/negative_binomial_distribution
// http://mathworld.wolfram.com/NegativeBinomialDistribution.html
// http://documents.wolfram.com/teachersedition/Teacher/Statistics/DiscreteDistributions.html
// The negative binomial distribution NegativeBinomialDistribution[n, p]
// is the distribution of the number (k) of failures that occur in a sequence of trials before
// r successes have occurred, where the probability of success in each trial is p.

// In a sequence of Bernoulli trials or events
// (independent, yes or no, succeed or fail) with success_fraction probability p,
// negative_binomial is the probability that k or fewer failures
// preceed the r th trial's success.

// Negative_binomial distribution is a discrete probability distribution.
// But note that the negative binomial distribution
// (like others including the binomial, Poisson & Bernoulli)
// is strictly defined as a discrete function: only integral values of k are envisaged.
// However because of the method of calculation using a continuous gamma function,
// it is convenient to treat it as if a continous function,
// and permit non-integral values of k.
// To enforce the strict mathematical model, users should use floor or ceil functions
// on k outside this function to ensure that k is integral.

// MATHCAD cumulative negative binomial pnbinom(k, n, p)

// random variable x is the number of success NOT the probability.

#ifndef BOOST_MATH_SPECIAL_NEGATIVE_BINOMIAL_HPP
#define BOOST_MATH_SPECIAL_NEGATIVE_BINOMIAL_HPP

#include <boost/math/special_functions/beta.hpp> // for ibeta(a, b, x) == Ix(a, b).
#include <boost/math/distributions/complement.hpp> // complement.
#include <boost/math/distributions/detail/common_error_handling.hpp> // error checks domain_error & logic_error.
#include <boost/math/special_functions/fpclassify.hpp> // isnan.
#include <boost/math/tools/roots.hpp> // for root finding.

#include <boost/type_traits/is_floating_point.hpp>
#include <boost/type_traits/is_integral.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/mpl/if.hpp>

#include <limits> // using std::numeric_limits;

#if defined (BOOST_MSVC) && defined(BOOST_MATH_THROW_ON_DOMAIN_ERROR)
#  pragma warning(push)
#  pragma warning(disable: 4702) // unreachable code
// in domain_error_imp in error_handling
#endif

namespace boost
{
  namespace math
  {
    namespace negative_binomial_detail
    { 
      // Common error checking routines for negative binomial distribution functions:
      template <class RealType>
      inline bool check_successes(const char* function, const RealType& r, RealType* result)
      {
        if( !(boost::math::isfinite)(r) || (r <= 0) )
        {
          *result = tools::domain_error<RealType>(
            function, 
            "Number of successes argument is %1%, but must be > 0 !", r);
          return false;
        }
        return true;
      }
      template <class RealType>
      inline bool check_success_fraction(const char* function, const RealType& p, RealType* result)
      {
        if( !(boost::math::isfinite)(p) || (p < 0) || (p > 1) )
        {
          *result = tools::domain_error<RealType>(
            function, 
            "Success fraction argument is %1%, but must be >= 0 and <= 1 !", p);
          return false;
        }
        return true;
      }
      template <class RealType>
      inline bool check_dist(const char* function, const RealType& r, const RealType& p, RealType* result)
      {
        return check_success_fraction(function, p, result) 
          && check_successes(function, r, result);
      }
      template <class RealType>
      bool check_dist_and_k(const char* function, const RealType& r, const RealType& p, RealType k, RealType* result)
      {
        if(check_dist(function, r, p, result) == false)
        {
          return false;
        }
        if( !(boost::math::isfinite)(k) || (k < 0) )
        { // Check k failures.
          *result = tools::domain_error<RealType>(
            function, 
            "Number of failures argument is %1%, but must be >= 0 !", k);
          return false;
        }
        return true;
      } // Check_dist_and_k

      template <class RealType>
      inline bool check_dist_and_prob(const char* function, const RealType& r, RealType p, RealType prob, RealType* result)
      {
        if(check_dist(function, r, p, result) && detail::check_probability(function, prob, result) == false)
        {
          return false;
        }
        return true;
      } // check_dist_and_prob
    } //  namespace negative_binomial_detail

    template <class RealType = double>
    class negative_binomial_distribution
    {
    public:
      typedef RealType value_type;

      negative_binomial_distribution(RealType r, RealType p) : m_r(r), m_p(p)
      { // Constructor.
        RealType result;
        negative_binomial_detail::check_dist(
          BOOST_CURRENT_FUNCTION,
          m_r, // check successes r > 0.
          m_p, // Check success_fraction 0 <= p <= 1.
          &result);
      } // negative_binomial_distribution constructor.

      // Private data getter class member functions.
      RealType success_fraction() const
      { // Probability of success as fraction in range 0 to 1.
        return m_p;
      }
      RealType successes() const
      { // Total number of successes r.
        return m_r;
      }

      // Estimation of the success_fraction parameter in a negative binomial distribution.
      // The best estimate is actually simply successes/trials:
      // these functions are used to obtain 
      // confidence intervals for the success fraction.
      // Use ibeta_inc function for lower bound but complement ibetac_inv for upper bound.

      static RealType estimate_lower_bound_on_p(
        RealType failures, 
        RealType successes,
        RealType probability)
      {
        // Error checks:
        RealType result;
        if(false == negative_binomial_detail::check_dist_and_k(
          BOOST_CURRENT_FUNCTION, failures, RealType(0), successes, &result)
          && 
          negative_binomial_detail::check_dist_and_prob(
          BOOST_CURRENT_FUNCTION, failures, RealType(0), probability, &result))
        {
          return result;
        }

        return ibeta_inv(successes + 1, failures - successes, probability); // TODO negative
      } // estimate_lower_bound_on_p

      static RealType estimate_upper_bound_on_p(
        RealType failures, 
        RealType successes,
        RealType probability)
      {
        // Error checks:
        RealType result;
        if(false == negative_binomial_detail::check_dist_and_k(
          BOOST_CURRENT_FUNCTION, failures, RealType(0), successes, &result)
          && 
          negative_binomial_detail::check_dist_and_prob(
          BOOST_CURRENT_FUNCTION, failures, RealType(0), probability, &result))
        {
          return result;
        }
        // Use complement ibetac_inv function for upper bound.
        return ibetac_inv(successes + 1, failures - successes, probability); // TODO negative
      } // estimate_upper_bound_on_p

      // Estimate number of failures parameter:  
      //
      // "How many failures do I need to be P% sure of seeing k successes?"  
      //    or
      // "How many failures can I have to be P% sure of seeing fewer than k events?"

      static RealType estimate_number_of_trials(
        RealType k,     // number of events.
        RealType p,     // success fraction.
        RealType probability) // probability threshold.
      {
        // Error checks:
        RealType result;
        if(false == negative_binomial_detail::check_dist_and_k(
          BOOST_CURRENT_FUNCTION, k, p, k, &result)
          && 
          negative_binomial_detail::check_dist_and_prob(
          BOOST_CURRENT_FUNCTION, k, p, probability, &result))
        { return result; }

        result = ibetac_invb(k + 1, p, probability);  // returns n - k
        return result + k;
      } // estimate_number_of_failures

      template <class P1, class P2, class P3>
      static RealType estimate_number_of_trials(
        const complemented3_type<P1, P2, P3>& c) 
      {
        // extract args:
        const RealType k = c.dist;     // number of events.
        const RealType p = c.param1;   // success fraction.
        const RealType Q = c.param2;   // probability threshold.

        // Error checks:
        RealType result;
        if(false == negative_binomial_detail::check_dist_and_k(
          BOOST_CURRENT_FUNCTION, k, p, k, &result)
          && 
          negative_binomial_detail::check_dist_and_prob(
          BOOST_CURRENT_FUNCTION, k, p, Q, &result))
        { return result; }

        result = ibeta_invb(k + 1, p, Q);  // returns n - k
        return result + k;
      } // RealType estimate_number_of_trials complemented

    private:
      RealType m_r; // successes.
      RealType m_p; // success_fraction
    }; // template <class RealType> class negative_binomial_distribution

    typedef negative_binomial_distribution<double> negative_binomial; // Reserved name of type double.

    template <class RealType>
    inline RealType mean(const negative_binomial_distribution<RealType>& dist)
    { // Mean of Negative Binomial distribution = r(1-p)/p.
      return dist.successes() * (1 - dist.success_fraction() ) / dist.success_fraction();
    } // mean

    template <class RealType>
    inline RealType mode(const negative_binomial_distribution<RealType>& dist)
    { // Mode of Negative Binomial distribution = floor[(r-1) * (1 - p)/p]
      return floor((dist.successes() -1) * (1 - dist.success_fraction()) / dist.success_fraction());
    } // mode

    template <class RealType>
    inline RealType median(const negative_binomial_distribution<RealType>& dist)
    { // Median of Negative Binomial distribution undefined
      return numeric_limits<RealType>::quiet_NaN();
    } // median

    template <class RealType>
    inline RealType hazard(const negative_binomial_distribution<RealType>& dist)
    { // hazard of Negative Binomial distribution ???  TODO
      return numeric_limits<RealType>::quiet_NaN();
    } // hazard

    template <class RealType>
    inline RealType chf(const negative_binomial_distribution<RealType>& dist)
    { // chf of Negative Binomial distribution ???  TODO
      return numeric_limits<RealType>::quiet_NaN();
    } // chf
    template <class RealType>
    inline RealType skewness(const negative_binomial_distribution<RealType>& dist)
    { // skewness of Negative Binomial distribution = 2-p / (sqrt(r(1-p))
      RealType p = dist.success_fraction();
      RealType r = dist.successes();

      return (2 - p) /
        sqrt(r * (1 - p));
    } // skewness

    template <class RealType>
    inline RealType kurtosis(const negative_binomial_distribution<RealType>& dist)
    { // kurtosis of Negative Binomial distribution = r(1-p)/p.
      RealType p = dist.success_fraction();
      RealType r = dist.successes();
      return 6 / r +
       (p * p) / r * (1 - p );
    } // kurtosis

    template <class RealType>
    inline RealType variance(const negative_binomial_distribution<RealType>& dist)
    { // Variance of Binomial distribution = r (1-p) / p^2.
      return  dist.successes() * (1 - dist.success_fraction())
        / (dist.success_fraction() * dist.success_fraction());
    } // variance

    template <class RealType>
    inline RealType standard_deviation(const negative_binomial_distribution<RealType>& dist)
    { // standard_deviation = sqrt (Variance of Binomial distribution = r (1-p) / p^2).
      return sqrt(variance(dist));
    } // variance


    template <class RealType>
    RealType pdf(const negative_binomial_distribution<RealType>& dist, const RealType k)
    { // Probability Density/Mass Function.
      BOOST_FPU_EXCEPTION_GUARD

      RealType r = dist.successes();
      RealType p = dist.success_fraction();

      RealType result = (p/(r + k)) * ibeta_derivative(r, static_cast<RealType>(k+1), p);
      // Numerical errors might cause probability to be slightly outside the range < 0 or > 1.
      // constrain_probability here?
      // This might cause trouble downstream, so warn, possibly throw exception, but constrain to the limits.
      // equivalent to
      // return exp(lgamma(r + k) - lgamma(r) - lgamma(k+1)) * pow(p, r) * pow((1-p), k);
      return result;
    } // negative_binomial_pdf

    template <class RealType>
    RealType cdf(const negative_binomial_distribution<RealType>& dist, const RealType k) 
    { // Cumulative Distribution Function of Negative Binomial.
      using boost::math::ibeta; // Regularized incomplete beta function.
      // k argument may be integral, signed, or unsigned, or floating point.
      // If necessary, it has already been promoted from an integral type.
      RealType p = dist.success_fraction();
      RealType r = dist.successes();
      // Error check:
      RealType result;
      if(false == negative_binomial_detail::check_dist_and_k(
        BOOST_CURRENT_FUNCTION,
        r,
        dist.success_fraction(),
        k,
        &result))
      {
        return result;
      }

      RealType probability = ibeta(r, static_cast<RealType>(k+1), p);
      // Ip(r, k+1) = ibeta(r, k+1, p)
      // constrain_probability here?
      // Numerical errors might cause probability to be slightly outside the range < 0 or > 1.
      // This might cause trouble downstream, so warn, possibly throw exception, but constrain to the limits.
      return probability;
    } // cdf Cumulative Distribution Function Negative Binomial.

      template <class RealType>
      RealType cdf(const complemented2_type<negative_binomial_distribution<RealType>, RealType>& c)
      { // Complemented Cumulative Distribution Function Negative Binomial.

      using boost::math::ibetac; // Regularized incomplete beta function complement.
      // k argument may be integral, signed, or unsigned, or floating point.
      // If necessary, it has already been promoted from an integral type.
      RealType const& k = c.param;
      negative_binomial_distribution<RealType> const& dist = c.dist;
      RealType p = dist.success_fraction();
      RealType r = dist.successes();
      // Error check:
      RealType result;
      if(false == negative_binomial_detail::check_dist_and_k(
        BOOST_CURRENT_FUNCTION,
        r,
        p,
        k,
        &result))
      {
        return result;
      }
      // Calculate cdf negative binomial using the incomplete beta function.
      // Use of ibeta here prevents cancellation errors in calculating
      // 1-p if p is very small, perhaps smaller than machine epsilon.
      // Ip(k+1, r) = ibetac(r, k+1, p)
      // constrain_probability here?
     RealType probability = ibetac(r, static_cast<RealType>(k+1), p);
      // Numerical errors might cause probability to be slightly outside the range < 0 or > 1.
      // This might cause trouble downstream, so warn, possibly throw exception, but constrain to the limits.
      return probability;
    } // cdf Cumulative Distribution Function Negative Binomial.

    template <class RealType>
    RealType quantile(const negative_binomial_distribution<RealType>& dist, const RealType& P)
    { // Quantile, percentile/100 or Percent Point Negative Binomial function.
      // Return the number of expected failures k for a given probability p.

      // Inverse cumulative Distribution Function or Quantile (percentile / 100) of negative_binomial Probability.
      // MAthCAD pnbinom return smallest k such that negative_binomial(k, n, p) >= probability.
      // k argument may be integral, signed, or unsigned, or floating point.
      // BUT Cephes/CodeCogs says: finds argument p (0 to 1) such that cdf(k, n, p) = y

      RealType p = dist.success_fraction();
      RealType r = dist.successes();

      //check dist and p 
      RealType result;
      if(false == negative_binomial_detail::check_dist_and_prob
        (BOOST_CURRENT_FUNCTION, r, p, P, &result))
      {
        return result;
      }

      // Special cases.
      if (P == 1)
      {  // Would need +infinity failures to have no probability of success.
       using std::numeric_limits;
       //return +numeric_limits<RealType>::infinity(); 
       return +numeric_limits<RealType>::max(); 
      }
      if (P == 0)
      { // No failures are allowed if P = 0.
        return 0;
      }
      using boost::math::ibeta_invb;
      // Calculate quantile of negative_binomial using the inverse incomplete beta function.
      return ibeta_invb(r, p, P) - 1; // 
    } // RealType quantile(const negative_binomial_distribution dist, p)

      template <class RealType>
      RealType quantile(const complemented2_type<negative_binomial_distribution<RealType>, RealType>& c)
      { // Quantile or Percent Point Binomial function.
        // Return the number of expected successes k for a given 
        // complement of the probability Q.
        //
      RealType p = dist.success_fraction();
      RealType r = dist.successes();

        // Error checks:
        RealType Q = c.param;
        const negative_binomial_distribution<RealType>& dist = c.dist;
        RealType result;
        if(false == negative_binomial_detail::check_dist_and_prob(
           BOOST_CURRENT_FUNCTION,
           r,
           p,
           q,
           &result))
        {
           return result;
        }
        // Special cases:
        //
        if(q == 1)
        {  // There may actually be no answer to this question,
           // since the probability of zero failures may be non-zero,
           // but zero is the best we can do:
           return 0;
        }
        if(q == 0)
        {  // Probability of greater than n failures is always zero,
           // so r is the most sensible answer here:  TODO check these special cases.
           return r;
        }

        return ibetac_invb(r, p, Q);
      } // quantile complement

 } // namespace math
} // namespace boost

// This include must be at the end, *after* the accessors
// for this distribution have been defined, in order to
// keep compilers that support two-phase lookup happy.
#include <boost/math/distributions/detail/derived_accessors.hpp>

#endif // BOOST_MATH_SPECIAL_NEGATIVE_BINOMIAL_HPP
