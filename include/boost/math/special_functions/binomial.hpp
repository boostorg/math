// boost\math\special_functions\binomial.hpp

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

// Also could get k from probability and mean x???  TODO

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


#ifndef BOOST_MATH_SPECIAL_BINOMIAL_HPP
#define BOOST_MATH_SPECIAL_BINOMIAL_HPP

#include <boost/math/special_functions/beta.hpp> // for incomplete beta.
#include <boost/math/tools/roots.hpp> // for domain_error & logic_error.
#include <boost/math/tools/promotion.hpp> // for argument promotion.

#include <boost/type_traits/is_floating_point.hpp>
#include <boost/type_traits/is_integral.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/mpl/if.hpp>

#include <limits> // 	using std::numeric_limits;

namespace boost
{
  namespace math
  {
    namespace detail
    { // Implementations called by actual functions
      // with arguments that, if necessary,
      // have been promoted from ArithmeticType to RealType.

      template <class RealType>
      RealType binomial_imp(RealType k, RealType n, RealType success_fraction)
      {
        using boost::math::ibeta; // Regularized incomplete beta function.
        using boost::math::tools::domain_error;
        using boost::math::tools::logic_error;
        using std::numeric_limits;

        // k argument may be integral, signed, or unsigned, or floating point.
        // If necessary, it has already been promoted from an integral type.
        if(k < 0)
        { // k must be >= 0!
          return domain_error<RealType>(BOOST_CURRENT_FUNCTION, "k argument is %1%, but must be >= 0 !", k);
        }
        if(n < 0)
        { // n must be >= 0!
          return domain_error<RealType>(BOOST_CURRENT_FUNCTION, "n argument is %1%, but must be >= 0 !", n);
        }

        if(k > n)
        { // n must be < k!
          return domain_error<RealType>(BOOST_CURRENT_FUNCTION, "n argument is %1%, but n must be >= k !", n);
        }
        // TODO Could use a 4 arg here for
        //if(n >= k)
        //{ // n must be < k!
        //  return domain_error<RealType>(BOOST_CURRENT_FUNCTION, "n argument is %1%, but must be < %2% !", n, k);
        //}
        if ((success_fraction < 0) || (success_fraction > 1))
        { // Check 0 <= probability probability <= 1.
          logic_error<RealType>(BOOST_CURRENT_FUNCTION, "success fraction is %1%, but must be >= 0 and <= 1 !", success_fraction);
          return static_cast<RealType>(0.); // Constrain to zero if logic_error does not throw.
        }

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
        return probability;
      } // binomial_imp

       template <class RealType>
      RealType binomial_c_imp(RealType k, RealType n, RealType success_fraction)
      {
        using boost::math::ibeta; // Regularized incomplete beta function.
        using boost::math::tools::domain_error;
        using boost::math::tools::logic_error;
        using std::numeric_limits;

        // k argument may be integral, signed, or unsigned, or floating point.
        // If necessary, it has already been promoted from an integral type.
        if(k < 0)
        { // k must be >= 0!
          return domain_error<RealType>(BOOST_CURRENT_FUNCTION, "k argument is %1%, but must be >= 0 !", k);
        }
        if(n < 0)
        { // n must be >= 0!
          return domain_error<RealType>(BOOST_CURRENT_FUNCTION, "n argument is %1%, but must be >= 0 !", n);
        }

        if(n <= k)
        { // n must be < k!
          return domain_error<RealType>(BOOST_CURRENT_FUNCTION, "n argument is %1%, but must be > k !", n);
        }
        // TODO Could use a 4 arg here for
        //if(n >= k)
        //{ // n must be < k!
        //  return domain_error<RealType>(BOOST_CURRENT_FUNCTION, "n argument is %1%, but must be < %2% !", n, k);
        //}
        if ((success_fraction < 0) || (success_fraction > 1))
        { // Check 0 <= probability probability <= 1.
          logic_error<RealType>(BOOST_CURRENT_FUNCTION, "success fraction is %1%, but must be >= 0 and <= 1 !", success_fraction);
          return static_cast<RealType>(0.); // Constrain to zero if logic_error does not throw.
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

     template <class RealType>
      RealType binomial_inv_imp(RealType k, RealType n, RealType probability)
      { // Inverse cumulative Distribution Function or Quintile (percentile / 100) of binomial Probability.
        // k argument may be integral, signed, or unsigned, or floating point.
        // If necessary, it has already been promoted from an integral type.

        using boost::math::gamma_Q_inv; // gamma_Q_inv
        using boost::math::tools::domain_error;
        using std::numeric_limits;

        // k argument may be integral, signed, or unsigned, or floating-point.
        if(k < 0)
        { // k must be >= 0!
          return domain_error<RealType>(BOOST_CURRENT_FUNCTION, "k argument is %1%, but must be >= 0 !", k);
        }
        if(n < 0)
        { // n must be >= 0!
          return domain_error<RealType>(BOOST_CURRENT_FUNCTION, "n argument is %1%, but must be >= 0 !", n);
        }
        if(n <= k)
        { // n must be <= k!
          return domain_error<RealType>(BOOST_CURRENT_FUNCTION, "n argument is %1%, but must be > k !", n);
        }

        if((probability < 0) || (probability > 1))
        { // probability must be >= 0 and <= 1!
          return domain_error<RealType>(BOOST_CURRENT_FUNCTION, "probability argument is %1%, but must be >= 0 and <= 1 !", probability);
        }
        // Special cases, regardless of k.  ?? TODO ???????????
        if (probability == 0)
          return -numeric_limits<RealType>::infinity();
        if (probability == 1)
          return +numeric_limits<RealType>::infinity();
        // Calculate quantile of binomial using the inverse incomplete beta function.
        return 1 - ibeta_inv(n - k, k + 1, probability);
      } // binomial_inv_imp
    } // namespace detail

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

  } // namespace math
} // namespace boost

#endif // BOOST_MATH_SPECIAL_BINOMIAL_HPP

/*

TODO remove this


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