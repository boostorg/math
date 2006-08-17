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
#include <boost/math/distributions/complement.hpp>
#include <boost/math/special_functions/factorials.hpp> // for factorials.
#include <boost/math/special_functions/log1p.hpp> // for log1p
#include <boost/math/tools/roots.hpp> // for ibeta_derivative.

//#include <boost/math/tools/roots.hpp> // for domain_error & logic_error.
//#include <boost/math/tools/promotion.hpp> // for argument promotion.
//
//#include <boost/type_traits/is_floating_point.hpp>
//#include <boost/type_traits/is_integral.hpp>
//#include <boost/type_traits/is_same.hpp>
//#include <boost/mpl/if.hpp>
//
#if defined (BOOST_MSVC) && defined(BOOST_MATH_THROW_ON_DOMAIN_ERROR)
#  pragma warning(push)
#  pragma warning(disable: 4702) // unreachable code
// in domain_error_imp in error_handling
#endif

namespace boost
{
  namespace math
  {
    template <class RealType>
    class binomial_distribution
    {
    public:
      binomial_distribution(RealType n, RealType p) : m_n(n), m_p(p)
      {
        if(m_n < 0)
        { // n must be >= 0!
          m_n = tools::domain_error<RealType>(BOOST_CURRENT_FUNCTION, "n argument is %1%, but must be >= 0 !", n);
          // If domain_error does NOT throw, it will return NaN and m_n = NaN.
        }

        if ((m_p < 0) || (m_p > 1)) // success_fraction or probability of success 
        { // Check 0 >= success fraction <= 1.
          m_p = tools::domain_error<RealType>(BOOST_CURRENT_FUNCTION, "success fraction is %1%, but must be >= 0 and <= 1 !", m_p);
           // If domain_error does NOT throw, it will return NaN and m_p = NaN.
        }
      } // binomial_distribution constructor.

      RealType success_fraction() const
      { // Probability.
        return m_p;
      }
      RealType trials() const
      { // Total number of trials.
        return m_n;
      }
    private:
        RealType m_n; // Not sure if this shouldn't be an int?
        RealType m_p; // success_fraction
      }; // template <class RealType> class binomial_distribution

      typedef binomial_distribution<double> binomial; // Reserved name of type double.

      template <class RealType>
      RealType mean(const binomial_distribution<RealType>& dist)
      { // Mean of Binomial distribution = np.
        return  dist.trials() * dist.success_fraction();
      } // mean

      template <class RealType>
      RealType variance(const binomial_distribution<RealType>& dist)
      { // Mean of Binomial distribution = np.
        return  dist.trials() * dist.success_fraction() * (1 - dist.success_fraction());
      } // mean

      template <class RealType>
      RealType pdf(const binomial_distribution<RealType>& dist, const RealType k)
      { // Probability Density/Mass Function.
        using boost::math::tools::domain_error;
        // Special cases of success_fraction, regardless of k successes and regardless of n trials.
        if (dist.success_fraction() == 0)
        {
          return 0;
        }
        if (dist.success_fraction() == 1)
        {
          return 0;
        }
        RealType n = dist.trials();
        if(n < 0)
        { // k must be <= n!
          return tools::domain_error<RealType>(BOOST_CURRENT_FUNCTION, "n argument is %1%, but must be >= 0 !", n);
        }
        // k argument may be integral, signed, or unsigned, or floating point.
        // If necessary, it has already been promoted from an integral type.
        if (n == 0)
        {
          return 1;
        }
        if(k < 0)
        { // k must be >= 0!
          return tools::domain_error<RealType>(BOOST_CURRENT_FUNCTION, "k argument is %1%, but must be >= 0 !", k);
        }
        if(k > n)
        { // k must be <= n!
          return tools::domain_error<RealType>(BOOST_CURRENT_FUNCTION, "n argument is %1%, but must be >= k !", n);
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

        // Probability of getting exactly k successes f(k; n,p) = binomial coefficient C(n, k) = * p ^ k (1-p)^(n-k) 
        //  n!/(k! * !(n-k)) * p ^ k * (1-p) ^ (n-k)
        // Use binomial_coefficient function here?
        //return exp(lgamma(n+1) - (lgamma(k+1) + lgamma(n-k+1))) // binomial coefficient C(n, k)
        //  * pow(dist.success_fraction(), k) * pow((1 - dist.success_fraction()), (n - k));
        // // Might use a table for binomial coefficients for small n?
        // But get overflow for large n because lgamma(n) overflows, so instead use:
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
        // y = bdtr( k, n, p ) = incbet( n-k, k+1, 1-p).

        using boost::math::tools::domain_error;
        using std::numeric_limits;

        // k argument may be integral, signed, or unsigned, or floating-point.
        if(k < 0)
        { // k must be >= 0!
          return tools::domain_error<RealType>(BOOST_CURRENT_FUNCTION, "k argument is %1%, but must be >= 0 !", k);
          //  warning C4702: unreachable code ???
        }
        RealType n = dist.trials();
        if(n < k)
        { // n must be <= k!
          return tools::domain_error<RealType>(BOOST_CURRENT_FUNCTION, "n argument is %1%, but must be >= k !", n);
          //  warning C4702: unreachable code ???
        }
        if(k > n)
        { // k should be <= n.  TODO is this the best - or return 1?
          return tools::domain_error<RealType>(BOOST_CURRENT_FUNCTION, "k argument is %1%, but should be <= n !", k);
          //  warning C4702: unreachable code ???
        }
        if (k == n)
        {
          return 1;
        }
        RealType p = dist.success_fraction();

        // Special cases, regardless of k.
        if (p == 0)
        { // TODO should this be very close to zero?
          return 0;
        }
        if (p == 1)
        {
          return 0;
        }
        if (k == 0)
        {
          return pow(1 - p, n - k);
        }
        // Calculate cdf binomial using the incomplete beta function.
        // y = bdtr( k, n, p ) = incbet( n-k, k+1, 1-p). ibeta (a, b, x)
        return ibeta(n - k, k + 1, 1 - p);
      } // binomial cdf

      template <class RealType>
      RealType binomial_coefficient(RealType n, RealType k)
      { // Binomial coefficient C(n, k)
        using ::boost::math::lgamma;

        if(k < 0)
        { // k must be >= 0!
          return tools::domain_error<RealType>(BOOST_CURRENT_FUNCTION, "k argument is %1%, but must be >= 0 !", k);
          //  warning C4702: unreachable code ???
        }
        if(n < 0)
        { // n must be >= 0!
          return tools::domain_error<RealType>(BOOST_CURRENT_FUNCTION, "n argument is %1%, but must be >= 0 !", n);
          //  warning C4702: unreachable code ???
        }
        // Fixed n
        // http://functions.wolfram.com/GammaBetaErf/Binomial/03/01/01/ 
        if (n == 0)
        {
          return 1;
        }
        if (k >= n)
        {
          return 1;
        }
        if (k == 0)
        {
          return 1;
        }
        if (k == 1)
        {
          return n;
        }
        if (k == 2)
        {
          return (n-1) * n/2;
        }
        if (k == 3)
        {
          return (n-2) * (n - 1) * n / 6;
        }
        if (k == 4)
        {
          return (n-3) * (n-2) * (n-1) * n / 24;
        }
        // more? divisor is factorial(k) MPL?
        return exp(lgamma(n+1) - (lgamma(k+1) + lgamma(n-k+1)));
      } //       template <class RealType> RealType binomial_coefficient(RealType n, RealType k)

      template <class RealType>
      RealType quantile(const binomial_distribution<RealType>& dist, const RealType& k, const RealType& y)
      { // Quantile or Percent Point Binomial function.
        // Returns the cdf y that would give k successes in n trials.

        // Moshier's bdtri(k, n, y)
        // Finds the event probability q such that the sum of the
        // terms 0 through k of the Binomial probability density (cdf)
        // is equal to the given cumulative probability y.

        // This is accomplished using the inverse beta integral function:
        // 1 - q = incbetinv(n-k, k+1, y)

        using boost::math::tools::domain_error;
        using boost::math::expm1;
        using boost::math::log1p;
        using boost::math::ibeta;
        using boost::math::ibeta_inv;

        if ((y < 0) || (y > 1))
        { // Check 0 <= cdf <= 1. 
          return tools::domain_error<RealType>(BOOST_CURRENT_FUNCTION, "probability is %1%, but must be >= 0 and <= 1 !", y);
        }
        if (k < 0)
        { // k must be >= 0!
          return tools::domain_error<RealType>(BOOST_CURRENT_FUNCTION, "k argument is %1%, but must be >= 0 !", k);
        }
        RealType n = dist.trials();

        if (k > n)
        { // n must be <= k!
          return tools::domain_error<RealType>(BOOST_CURRENT_FUNCTION, "k argument is %1%, but must be <= n !", k);
        }

        RealType p = dist.success_fraction();
        // Special cases, regardless of k. 
        if (p == 0)
        {
          return 0;
        }
        if (p == 1)
        {
          return 0;
        }
        RealType dn = n - k; // Leave this optimisation to the compiler?
        if (k == 0)
        {
          if (y > 0.8)
          {
            return -expm1(log1p(y - 1) / dn);
          }
          else
          {
            return 1 - pow(y, 1 / dn);
          }
        }
        else
        {
          RealType q = ibeta(dn, k+1, 0.5); // quantile result.
          return ((q > 0.5) ? ibeta_inv(k+1, dn, 1 - y) : 1 - ibeta_inv(dn, k+1, y));
        }
      } // quantile

    } // namespace math
  } // namespace boost

#endif // BOOST_MATH_SPECIAL_BINOMIAL_HPP

  /*

  TODO remove this

  Old stuff parked

  // k argument may be integral, signed, or unsigned, or floating point.
  // If necessary, it has already been promoted from an integral type.
  //if(k < 0)
  //{ // k must be >= 0!
  //  return tools::domain_error<RealType>(BOOST_CURRENT_FUNCTION, "k argument is %1%, but must be >= 0 !", k);
  //}

  //if(k > n)
  //{ // n must be < k!
  //  return tools::domain_error<RealType>(BOOST_CURRENT_FUNCTION, "n argument is %1%, but n must be >= k !", n);
  //}
  // TODO Could use a 4 arg here for
  //if(n >= k)
  //{ // n must be < k!
  //  return tools::domain_error<RealType>(BOOST_CURRENT_FUNCTION, "n argument is %1%, but must be < %2% !", n, k);
  //}

  using boost::math::ibeta; // Regularized incomplete beta function.
  using boost::math::tools::domain_error;
  using boost::math::tools::logic_error;
  using std::numeric_limits;


  RealType probability = ibeta(static_cast<RealType>(n) - static_cast<RealType>(k), static_cast<RealType>(k)+1, 1 - success_fraction);
  // Cephes incbet(n-k, k+1, 1-p);
  // Numerical errors might cause probability to be slightly outside the range < 0 or > 1.
  // This might cause trouble downstream, so warn, possibly throw exception, but constrain to the limits.
  if (probability < 0)
  {
  logic_error<RealType>(BOOST_CURRENT_FUNCTION, "probability %1% is < 0, so has been constrained to zero !", probability);
  return static_cast<RealType>(0.); // Constrain to zero if logic_error does not throw.
  }
  if(probability > 1)
  {
  logic_error<RealType>(BOOST_CURRENT_FUNCTION, "probability %1% is > 1, so has been constrained to unity!", probability);
  return static_cast<RealType>(1.); // Constrain to unity if logic_error does not throw.
  }

  RealType probability = ibeta(static_cast<RealType>(k)+1, static_cast<RealType>(n) - static_cast<RealType>(k), success_fraction);
  // Cephes incbet(k+1, n-k, p);
  // Numerical errors might cause probability to be slightly outside the range < 0 or > 1.
  // This might cause trouble downstream, so warn, possibly throw exception, but constrain to the limits.
  if (probability < 0)
  {
  logic_error<RealType>(BOOST_CURRENT_FUNCTION, "probability %1% is < 0, so has been constrained to zero !", probability);
  return static_cast<RealType>(0.); // Constrain to zero if logic_error does not throw.
  }
  if(probability > 1)
  {
  logic_error<RealType>(BOOST_CURRENT_FUNCTION, "probability %1% is > 1, so has been constrained to unity!", probability);
  return static_cast<RealType>(1.); // Constrain to unity if logic_error does not throw.
  }
  return probability;
  } // binomial_imp

  template <class ArithmeticType, class RealType> // Binomial distribution (k, n, x)
  // Probability of number of events between 0 and k-1 inclusive, if expected probability of success events is success_fraction.
  inline typename tools::promote_arg3<ArithmeticType, ArithmeticType, RealType>::type
  // Return type is the wider of the two (perhaps promoted) floating-point types.
  binomial(ArithmeticType k, ArithmeticType n, RealType success_fraction)
  {
  typedef typename tools::promote_arg3<ArithmeticType, ArithmeticType, RealType>::type promote_type; // Arguments type.
  return detail::binomial_imp(static_cast<promote_type>(k), static_cast<promote_type>(n), static_cast<promote_type>(success_fraction));
  } // binomial

  template <class ArithmeticType, class RealType> // Binomial distribution complement (k, n, x)
  // Probability of number of events between 0 and k-1 inclusive, if expected mean is x.
  inline typename tools::promote_arg3<ArithmeticType, ArithmeticType, RealType>::type
  // Return type is the wider of the two (perhaps promoted) floating-point types.
  binomial_c(ArithmeticType k, ArithmeticType n, RealType success_fraction)
  {
  typedef typename tools::promote_arg3<ArithmeticType, RealType, ArithmeticType>::type promote_type; // Arguments type.
  return detail::binomial_c_imp(static_cast<promote_type>(k), static_cast<promote_type>(n), static_cast<promote_type>(success_fraction));
  } // binomial

  template <class ArithmeticType, class RealType>
  // success_fraction if number of events between 0 and k-1 inclusive, and probability p.
  inline typename tools::promote_arg3<ArithmeticType, ArithmeticType, RealType>::type
  // Return type is the wider of the two (perhaps promoted) floating point types.
  binomial_inv(ArithmeticType k,  ArithmeticType n, RealType probability) // Binomial distribution inverse (k, n, p)
  {
  typedef typename tools::promote_arg3<ArithmeticType, ArithmeticType, RealType>::type promote_type; // Arguments type.
  return detail::binomial_inv_imp(static_cast<promote_type>(k), static_cast<promote_type>(n),static_cast<promote_type>(probability));
  } // binomial_inv



  * Returns the sum of the terms 0 through k of the Binomial
  * probability density:
  *
  *   k
  *   --  ( n )   j      n-j
  *   >   (   )  p  (1-p)
  *   --  ( j )
  *  j=0
  *
  * The terms are not summed directly; instead the incomplete
  * beta integral is employed, according to the formula
  *
  * y = bdtr( k, n, p ) = incbet( n-k, k+1, 1-p ).
  *
  * The arguments must be positive, with p ranging from 0 to 1.
  *   message         condition      value returned
  * bdtr domain         k < 0            0.0
  *                     n < k
  *                     x < 0, x > 1


  * Complement
  * Returns the sum of the terms k+1 through n of the Binomial
  * probability density:
  *
  *   n
  *   --  ( n )   j      n-j
  *   >   (   )  p  (1-p)
  *   --  ( j )
  *  j=k+1
  *
  * The terms are not summed directly; instead the incomplete
  * beta integral is employed, according to the formula
  *
  * y = bdtrc( k, n, p ) = incbet( k+1, n-k, p ).
  *
  * The arguments must be positive, with p ranging from 0 to 1.
  *
  *   message         condition      value returned
  * bdtrc domain      x<0, x>1, n<k       0.0


  * p = bdtr( k, n, y );
  *
  * DESCRIPTION:
  *
  * Finds the event probability p such that the sum of the
  * terms 0 through k of the Binomial probability density
  * is equal to the given cumulative probability y.
  *
  * This is accomplished using the inverse beta integral
  * function and the relation
  *
  * 1 - p = incbi( n-k, k+1, y ).
  *
  *   message         condition      value returned
  * bdtri domain     k < 0, n <= k         0.0
  *                  x < 0, x > 1
  */