// boost\math\distributions\binomial.hpp

// Copyright John Maddock 2006.
// Copyright Paul A. Bristow 2006.

// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

// http://en.wikipedia.org/wiki/binomial_distribution

// Binomial distribution is the discrete probability distribution of
// the number (k) of successes in a sequence of
// n independent (yes or no, success or failure) Bernoulli trials.

// It expresses the probability of a number of events occurring in a fixed time
// if these events occur with a known average rate (probability of success),
// and are independent of the time since the last event.

// The binomial distribution was discovered by Siméon-Denis (1781–1840).

// The number of cars that pass through a certain point on a road during a given period of time.
// The number of spelling mistakes a secretary makes while typing a single page.
// The number of phone calls at a call center per minute.
// The number of times a web server is accessed per minute.
// The number of light bulbs that burn out in a certain amount of time.
// The number of roadkill found per unit length of road

// http:/en.wikipedia.org/wiki/binomial_distribution

// Given a sample of N measured values k[i],
// we wish to estimate the value of the parameter x (mean)
// of the binomial population from which the sample was drawn.
// To calculate the maximum likelihood value = 1/N sum i = 1 to N of k[i]

// Also may want a function for EXACTLY k.

// And probability that there are EXACTLY k occurrences is
// exp(-x) * pow(x, k) / factorial(k)
// where x is expected occurrences (mean) during the given interval.
// For example, if events occur, on average, every 4 min,
// and we are interested in number of events occurring in 10 min,
// then x = 10/4 = 2.5

// http://www.itl.nist.gov/div898/handbook/eda/section3/eda366i.htm

// The binomial distribution is used when there are
// exactly two mutually exclusive outcomes of a trial.
// These outcomes are appropriately labeled "success" and "failure".
// The binomial distribution is used to obtain
// the probability of observing x successes in N trials,
// with the probability of success on a single trial denoted by p.
// The binomial distribution assumes that p is fixed for all trials.

// P(x, p, n) = n!/(x! * (n-x)!) * p^x * (1-p)^(n-x)

// http://mathworld.wolfram.com/BinomialCoefficient.html

// The binomial coefficient (n; k) is the number of ways of picking
// k unordered outcomes from n possibilities,
// also known as a combination or combinatorial number.
// The symbols _nC_k and (n; k) are used to denote a binomial coefficient,
// and are sometimes read as "n choose k."
// (n; k) therefore gives the number of k-subsets  possible out of a set of n distinct items.

// For example:
//  The 2-subsets of {1,2,3,4} are the six pairs {1,2}, {1,3}, {1,4}, {2,3}, {2,4}, and {3,4}, so (4; 2)==6.

// http://functions.wolfram.com/GammaBetaErf/Binomial/ for evaluation.

#ifndef BOOST_MATH_SPECIAL_BINOMIAL_HPP
#define BOOST_MATH_SPECIAL_BINOMIAL_HPP

#include <boost/math/special_functions/beta.hpp> // for incomplete beta.
#include <boost/math/distributions/complement.hpp> // complements
#include <boost/math/distributions/detail/common_error_handling.hpp> // error checks
#include <boost/math/special_functions/fpclassify.hpp> // isnan.
#include <boost/math/tools/roots.hpp> // for root finding.

#if defined (BOOST_MSVC) && defined(BOOST_MATH_THROW_ON_DOMAIN_ERROR)
#  pragma warning(push)
#  pragma warning(disable: 4702) // unreachable code
// in domain_error_imp in error_handling
#endif

namespace boost
{
  namespace math
  {
     namespace binomial_detail{
        // common error checking routines for binomial distribution functions:
        template <class RealType>
        inline bool check_N(const char* function, const RealType& N, RealType* result)
        {
           if((N < 0) || !(boost::math::isfinite)(N))
           {
               *result = tools::domain_error<RealType>(
                  function, 
                  "Number of Trials argument is %1%, but must be >= 0 !", N);
               return false;
           }
           return true;
        }
        template <class RealType>
        inline bool check_success_fraction(const char* function, const RealType& p, RealType* result)
        {
           if((p < 0) || (p > 1) || !(boost::math::isfinite)(p))
           {
               *result = tools::domain_error<RealType>(
                  function, 
                  "Success fraction argument is %1%, but must be >= 0 and <= 1 !", p);
               return false;
           }
           return true;
        }
        template <class RealType>
        inline bool check_dist(const char* function, const RealType& N, const RealType& p, RealType* result)
        {
           return check_success_fraction(
              function, p, result) 
              && check_N(
               function, N, result);
        }
        template <class RealType>
        bool check_dist_and_k(const char* function, const RealType& N, const RealType& p, RealType k, RealType* result)
        {
           if(check_dist(function, N, p, result) == false)
              return false;
           if((k < 0) || !(boost::math::isfinite)(k))
           {
               *result = tools::domain_error<RealType>(
                  function, 
                  "Number of Successes argument is %1%, but must be >= 0 !", k);
               return false;
           }
           if(k > N)
           {
               *result = tools::domain_error<RealType>(
                  function, 
                  "Number of Successes argument is %1%, but must be <= Number of Trials !", k);
               return false;
           }
           return true;
        }
        template <class RealType>
        inline bool check_dist_and_prob(const char* function, const RealType& N, RealType p, RealType prob, RealType* result)
        {
           if(check_dist(function, N, p, result) && detail::check_probability(function, prob, result) == false)
              return false;
           return true;
        }
     }

    template <class RealType>
    class binomial_distribution
    {
    public:
      typedef RealType value_type;

      binomial_distribution(RealType n, RealType p) : m_n(n), m_p(p)
      {
         RealType r;
         binomial_detail::check_dist(
            BOOST_CURRENT_FUNCTION,
            m_n,
            m_p,
            &r);
      } // binomial_distribution constructor.

      RealType success_fraction() const
      { // Probability.
        return m_p;
      }
      RealType trials() const
      { // Total number of trials.
        return m_n;
      }

      //
      // Estimation of the success parameter.
      // The best estimate is actually simply
      // successes/trials, these functions are used
      // to obtain confidence intervals for the success
      // fraction:
      //
      static RealType estimate_lower_bound_on_p(
         RealType trials, 
         RealType successes,
         RealType probability)
      {
        // Error checks:
        RealType result;
        if(false == binomial_detail::check_dist_and_k(
           BOOST_CURRENT_FUNCTION, trials, RealType(0), successes, &result)
            && 
           binomial_detail::check_dist_and_prob(
           BOOST_CURRENT_FUNCTION, trials, RealType(0), probability, &result))
        { return result; }

        return ibeta_inv(successes + 1, trials - successes, probability);
      }
      static RealType estimate_upper_bound_on_p(
         RealType trials, 
         RealType successes,
         RealType probability)
      {
        // Error checks:
        RealType result;
        if(false == binomial_detail::check_dist_and_k(
           BOOST_CURRENT_FUNCTION, trials, RealType(0), successes, &result)
            && 
           binomial_detail::check_dist_and_prob(
           BOOST_CURRENT_FUNCTION, trials, RealType(0), probability, &result))
        { return result; }

        return ibetac_inv(successes + 1, trials - successes, probability);
      }
      //
      // Estimate number of trials parameter:
      //
      // "How many trials do I need to be P% sure of seeing k events?"
      //    or
      // "How many trials can I have to be P% sure of seeing fewer than k events?"
      //
      static RealType estimate_number_of_trials(
         RealType k,     // number of events
         RealType p,     // success fraction
         RealType probability) // probability threshold
      {
        // Error checks:
        RealType result;
        if(false == binomial_detail::check_dist_and_k(
           BOOST_CURRENT_FUNCTION, k, p, k, &result)
            && 
           binomial_detail::check_dist_and_prob(
           BOOST_CURRENT_FUNCTION, k, p, probability, &result))
        { return result; }

        result = ibetac_invb(k + 1, p, probability);  // returns n - k
        return result + k;
      }

      template <class P1, class P2, class P3>
      static RealType estimate_number_of_trials(
         const complemented3_type<P1, P2, P3>& c) 
      {
        // extract args:
        const RealType k = c.dist;     // number of events
        const RealType p = c.param1;   // success fraction
        const RealType Q = c.param2;   // probability threshold

        // Error checks:
        RealType result;
        if(false == binomial_detail::check_dist_and_k(
           BOOST_CURRENT_FUNCTION, k, p, k, &result)
            && 
           binomial_detail::check_dist_and_prob(
           BOOST_CURRENT_FUNCTION, k, p, Q, &result))
        { return result; }

        result = ibeta_invb(k + 1, p, Q);  // returns n - k
        return result + k;
      }

    private:
        RealType m_n; // Not sure if this shouldn't be an int?
        RealType m_p; // success_fraction
      }; // template <class RealType> class binomial_distribution

      typedef binomial_distribution<double> binomial; // Reserved name of type double.

      template <class RealType>
      inline RealType mean(const binomial_distribution<RealType>& dist)
      { // Mean of Binomial distribution = np.
        return  dist.trials() * dist.success_fraction();
      } // mean

      template <class RealType>
      inline RealType variance(const binomial_distribution<RealType>& dist)
      { // Variance of Binomial distribution = np(1-p).
        return  dist.trials() * dist.success_fraction() * (1 - dist.success_fraction());
      } // variance

      template <class RealType>
      RealType pdf(const binomial_distribution<RealType>& dist, const RealType k)
      { // Probability Density/Mass Function.
        using boost::math::tools::domain_error;
        using namespace std; // for ADL of std functions

        RealType n = dist.trials();

        // Error check:
        RealType result;
        if(false == binomial_detail::check_dist_and_k(
           BOOST_CURRENT_FUNCTION,
           n,
           dist.success_fraction(),
           k,
           &result))
        {
           return result;
        }

        // Special cases of success_fraction, regardless of k successes and regardless of n trials.
        if (dist.success_fraction() == 0)
        {
           // probability of zero successes is 1:
           return static_cast<RealType>(k == 0 ? 1 : 0);
        }
        if (dist.success_fraction() == 1)
        {
           // probability of n successes is 1:
           return static_cast<RealType>(k == n ? 1 : 0);
        }
        // k argument may be integral, signed, or unsigned, or floating point.
        // If necessary, it has already been promoted from an integral type.
        if (n == 0)
        {
          return 1;
        }
        if (k == 0)
        { // binomial coeffic (n 0) = 1,
          // n ^ 0 = 1
          return pow(1 - dist.success_fraction(), n);
        }
        if (k == n)
        { // binomial coeffic (n n) = 1,
          // n ^ 0 = 1
          return pow(dist.success_fraction(), k);  // * pow((1 - dist.success_fraction()), (n - k)) = 1
        }

        // Probability of getting exactly k successes 
        // if C(n, k) is the binomial coefficient then:
        //
        // f(k; n,p) = C(n, k) * p^k * (1-p)^(n-k) 
        //           = (n!/(k!(n-k)!)) * p^k * (1-p)^(n-k)
        //           = (tgamma(n+1) / (tgamma(k+1)*tgamma(n-k+1))) * p^k * (1-p)^(n-k)
        //           = p^k (1-p)^(n-k) / (beta(k+1, n-k+1) * (n+1))
        //           = ibeta_derivative(k+1, n-k+1, p) / (n+1)
        //
        using boost::math::ibeta_derivative; // a, b, x
        return ibeta_derivative(k+1, n-k+1, dist.success_fraction()) / (n+1);

      } // pdf

      template <class RealType>
      RealType cdf(const binomial_distribution<RealType>& dist, const RealType k)
      { // Cumulative Distribution Function Binomial.
        // The random variate k is the number of successes in n trials.
        // k argument may be integral, signed, or unsigned, or floating point.
        // If necessary, it has already been promoted from an integral type.

        // Returns the sum of the terms 0 through k of the Binomial Probability Density/Mass:
        //
        //   i=k
        //   --  ( n )   i      n-i
        //   >   |   |  p  (1-p)
        //   --  ( i )
        //   i=0

        // The terms are not summed directly (at least for larger k)
        // instead the incomplete beta integral is employed,
        // according to the formula:
        // P = I[1-p]( n-k, k+1).
        //   = 1 - I[p](k + 1, n - k)

        using boost::math::tools::domain_error;
        using namespace std; // for ADL of std functions

        RealType n = dist.trials();
        RealType p = dist.success_fraction();

        // Error check:
        RealType result;
        if(false == binomial_detail::check_dist_and_k(
           BOOST_CURRENT_FUNCTION,
           n,
           p,
           k,
           &result))
        {
           return result;
        }
        if (k == n)
        {
          return 1;
        }

        // Special cases, regardless of k.
        if (p == 0)
        { 
           // This need explanation: the pdf is zero for all
           // cases except when k == 0.  For zero p the probability
           // of zero successes is one.  Therefore the cdf is always
           // 1: the probability of k or *fewer* successes is always 1
           // if there are never any successes!
           return 1;
        }
        if (p == 1)
        {
          // This is correct but needs explanation, when k = 1
          // all the cdf and pdf values are zero *except* when
          // k == n, and that case has been handled above already.
          return 0;
        }
        if((k < 20) && (floor(k) == k))
        {
          // For small k use a finite sum, it's cheaper
          // than the incomplete beta:
          RealType result = 0;
          for(unsigned i = 0; i <= k; ++i)
             result += pdf(dist, static_cast<RealType>(i));
          return result;
        }
        // Calculate cdf binomial using the incomplete beta function.
        // P = I[1-p](n - k, k + 1)
        //   = 1 - I[p](k + 1, n - k)
        // Use of ibetac here prevents cancellation errors in calculating
        // 1-p if p is very small, perhaps smaller than machine epsilon.
        return ibetac(k + 1, n - k, p);
      } // binomial cdf

      template <class RealType>
      RealType cdf(const complemented2_type<binomial_distribution<RealType>, RealType>& c)
      { // Complemented Cumulative Distribution Function Binomial.
        // The random variate k is the number of successes in n trials.
        // k argument may be integral, signed, or unsigned, or floating point.
        // If necessary, it has already been promoted from an integral type.

        // Returns the sum of the terms k+1 through n of the Binomial Probability Density/Mass:
        //
        //   i=n
        //   --  ( n )   i      n-i
        //   >   |   |  p  (1-p)
        //   --  ( i )
        //   i=k+1

        // The terms are not summed directly (at least for larger k)
        // instead the incomplete beta integral is employed,
        // according to the formula:
        // Q = 1 -I[1-p]( n-k, k+1).
        //   = I[p](k + 1, n - k)

        using boost::math::tools::domain_error;
        using namespace std; // for ADL of std functions

        RealType const& k = c.param;
        binomial_distribution<RealType> const& dist = c.dist;
        RealType n = dist.trials();
        RealType p = dist.success_fraction();

        // Error checks:
        RealType result;
        if(false == binomial_detail::check_dist_and_k(
           BOOST_CURRENT_FUNCTION,
           n,
           p,
           k,
           &result))
        {
           return result;
        }

        if (k == n)
        {
          // Probability of greater than n successes is necessarily zero:
          return 0;
        }

        // Special cases, regardless of k.
        if (p == 0)
        { 
           // This need explanation: the pdf is zero for all
           // cases except when k == 0.  For zero p the probability
           // of zero successes is one.  Therefore the cdf is always
           // 1: the probability of *more than* k successes is always 0
           // if there are never any successes!
           return 0;
        }
        if (p == 1)
        {
          // This needs explanation, when p = 1
          // we always have n successes, so the probability
          // of more than k successes is 1 as long as k < n.
          // The k == n case has already been handled above.
          return 1;
        }
        if((n - k < 20) && (floor(k) == k) && (floor(n) == n))
        {
          // For small n-k use a finite sum, it's cheaper
          // than the incomplete beta:
          RealType result = 0;
          for(RealType i = n; i > k; i -= 1)
             result += pdf(dist, i);
          return result;
        }
        // Calculate cdf binomial using the incomplete beta function.
        // Q = 1 -I[1-p](n - k, k + 1)
        //   = I[p](k + 1, n - k)
        // Use of ibeta here prevents cancellation errors in calculating
        // 1-p if p is very small, perhaps smaller than machine epsilon.
        return ibeta(k + 1, n - k, p);
      } // binomial cdf

      namespace detail{

         template <class RealType>
         struct binomial_functor
         {
            binomial_functor(const binomial_distribution<RealType>& d, const RealType& target, bool c = false)
               : dist(d), t(target), complement(c) {}

            RealType operator()(const RealType k)
            {
               if(k >= dist.trials())
                  return 1; // any positive value will do.
               return complement ? t - cdf(boost::math::complement(dist, k)) : cdf(dist, k) - t;
            }
         private:
            const binomial_distribution<RealType>& dist;
            RealType t;
            bool complement;
         };

      }

      template <class RealType>
      RealType quantile(const binomial_distribution<RealType>& dist, const RealType& p)
      { // Quantile or Percent Point Binomial function.
        // Return the number of expected successes k for a given 
        // probability p.
        //
        // Error checks:
        RealType result;
        if(false == binomial_detail::check_dist_and_prob(
           BOOST_CURRENT_FUNCTION,
           dist.trials(),
           dist.success_fraction(),
           p,
           &result))
        {
           return result;
        }
        //
        // Special cases:
        //
        if(p == 0)
        {
           // There may actually be no answer to this question,
           // since the probability of zero successes may be non-zero,
           // but zero is the best we can do:
           return 0;
        }
        if(p == 1)
        {
           // probability of n or fewer successes is always one,
           // so n is the most sensible answer here:
           return dist.trials();
        }

        //
        // Solve for quantile numerically:
        //
        detail::binomial_functor<RealType> f(dist, p);
        tools::eps_tolerance<RealType> tol(tools::digits<RealType>());
        boost::uintmax_t max_iter = 200;
        std::pair<RealType, RealType> r = tools::bracket_and_solve_root(
           f, 
           dist.trials() / 2, 
           static_cast<RealType>(2),
           true,
           tol,
           max_iter);
        if(max_iter >= 200)
           tools::logic_error<RealType>(BOOST_CURRENT_FUNCTION, "Unable to locate the root within a reasonable number of iterations, closest approximation so far was %1%", r.first);
        // return centre point of range found:
        return r.first + (r.second - r.first) / 2;
      } // quantile

      template <class RealType>
      RealType quantile(const complemented2_type<binomial_distribution<RealType>, RealType>& c)
      { // Quantile or Percent Point Binomial function.
        // Return the number of expected successes k for a given 
        // complement of the probability q.
        //
        // Error checks:
        RealType q = c.param;
        const binomial_distribution<RealType>& dist = c.dist;
        RealType result;
        if(false == binomial_detail::check_dist_and_prob(
           BOOST_CURRENT_FUNCTION,
           dist.trials(),
           dist.success_fraction(),
           q,
           &result))
        {
           return result;
        }
        //
        // Special cases:
        //
        if(q == 1)
        {
           // There may actually be no answer to this question,
           // since the probability of zero successes may be non-zero,
           // but zero is the best we can do:
           return 0;
        }
        if(q == 0)
        {
           // probability of greater than n successes is always zero,
           // so n is the most sensible answer here:
           return dist.trials();
        }

        //
        // Solve for quantile numerically:
        //
        detail::binomial_functor<RealType> f(dist, q, true);
        tools::eps_tolerance<RealType> tol(tools::digits<RealType>());
        boost::uintmax_t max_iter = 200;
        std::pair<RealType, RealType> r = tools::bracket_and_solve_root(
           f, 
           dist.trials() / 2, 
           static_cast<RealType>(2),
           true,
           tol,
           max_iter);
        if(max_iter >= 200)
           tools::logic_error<RealType>(BOOST_CURRENT_FUNCTION, "Unable to locate the root within a reasonable number of iterations, closest approximation so far was %1%", r.first);
        // return centre point of range found:
        return r.first + (r.second - r.first) / 2;
      } // quantile

    } // namespace math
  } // namespace boost

// This include must be at the end, *after* the accessors
// for this distribution have been defined, in order to
// keep compilers that support two-phase lookup happy.
#include <boost/math/distributions/detail/derived_accessors.hpp>

#endif // BOOST_MATH_SPECIAL_BINOMIAL_HPP

