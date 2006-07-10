// boost\math\special_functions\negative_binomial.hpp

// Copyright Paul A. Bristow 2006.

// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

// http://en.wikipedia.org/wiki/negative_binomial_distribution

// negative_binomial distribution is a discrete probability distribution.

// In a sequence of Bernoulli  trials or events
// (independent, yes or no, succeed or fail) with success probability p,
// negative_binomial is the probability that k or fewer failures
// preceed the nth trial's success.

#ifndef BOOST_MATH_SPECIAL_NEGATIVE_BINOMIAL_HPP
#define BOOST_MATH_SPECIAL_NEGATIVE_BINOMIAL_HPP

#include <boost/math/special_functions/beta.hpp> // for ibeta(a, b, x).
#include <boost/math/tools/roots.hpp> // for domain_error & logic_error.
#include <boost/math/tools/promotion.hpp> // for promotion.

#include <boost/type_traits/is_floating_point.hpp>
#include <boost/type_traits/is_integral.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/mpl/if.hpp>

#include <limits> // using std::numeric_limits;

namespace boost
{
	namespace math
	{
      namespace detail
			{ // Implementations called by actual functions 
        // with arguments that, if necessary,
        // have been promoted from ArithmeticType to RealType.
      template <class RealType>
      RealType negative_binomial_imp(RealType k, RealType n, RealType success_probability) 
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

        if ((success_probability < 0) || (success_probability > 1))
        { // Check 0 <= probability probability <= 1.  
          logic_error<RealType>(BOOST_CURRENT_FUNCTION, "success fraction is %1%, but must be >= 0 and <= 1 !", success_probability);
          return static_cast<RealType>(0.); // Constrain to zero if logic_error does not throw.
        }

        RealType probability = ibeta(static_cast<RealType>(k), static_cast<RealType>(n), 1 - success_probability);
        // Cephes nbdtr incbet(k, n, p);
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
      } // negative_binomial_imp

       template <class RealType>
      RealType negative_binomial_c_imp(RealType k, RealType n, RealType success_probability) 
      { 
        using boost::math::ibeta; // Regularized incomplete beta function.
        using boost::math::tools::domain_error;
        using boost::math::tools::logic_error;
        using std::numeric_limits;

        // k argument may be integral, signed, or unsigned, or floating point.
        // If necessary, it has already been promoted from an integral type.
        if (k < 0)
        { // k must be >= 0!
          return domain_error<RealType>(BOOST_CURRENT_FUNCTION, "k argument is %1%, but must be >= 0 !", k);
        }
        if (n < 0)
        { // n must be >= 0!
          return domain_error<RealType>(BOOST_CURRENT_FUNCTION, "n argument is %1%, but must be >= 0 !", n);
        }

        if (n <= k)
        { // n must be < k!
          return domain_error<RealType>(BOOST_CURRENT_FUNCTION, "n argument is %1%, but must be > k !", n);
        }
        // TODO Could use a 4 arg here for 
        //if(n >= k)
        //{ // n must be < k!
        //  return domain_error<RealType>(BOOST_CURRENT_FUNCTION, "n argument is %1%, but must be < %2% !", n, k);
        //}
        if ((success_probability < 0) || (success_probability > 1))
        { // Check 0 <= probability probability <= 1.  
          logic_error<RealType>(BOOST_CURRENT_FUNCTION, "success fraction is %1%, but must be >= 0 and <= 1 !", success_probability);
          return static_cast<RealType>(0.); // Constrain to zero if logic_error does not throw.
        }

        RealType probability = ibeta(static_cast<RealType>(k)+1, static_cast<RealType>(n), 1 - success_probability);
        // Cephes nbdtrc incbet(k+1, n, 1-p);
        // Numerical errors might cause probability to be slightly outside the range < 0 or > 1.
        // This might cause trouble downstream, so warn, possibly throw exception, but constrain to the limits.
        if (probability < 0)
        {
          logic_error<RealType>(BOOST_CURRENT_FUNCTION, "probability %1% is < 0, so has been constrained to zero !", probability);
          return static_cast<RealType>(0.); // Constrain to zero if logic_error does not throw.
        }
        if (probability > 1)
        {
          logic_error<RealType>(BOOST_CURRENT_FUNCTION, "probability %1% is > 1, so has been constrained to unity!", probability);
          return static_cast<RealType>(1.); // Constrain to unity if logic_error does not throw.
        }
        return probability;
      } // negative_binomial_imp

     template <class RealType>
      RealType negative_binomial_inv_imp(RealType k, RealType n, RealType probability)
      { // Inverse cumulative Distribution Function or Quintile (percentile / 100) of negative_binomial Probability.
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
        // Calculate quantile of negative_binomial using the inverse incomplete beta function.
        return ibeta_inv(n, k + 1, probability); // returns success_probability.
      } // negative_binomial_inv_imp
    } // namespace detail

    template <class ArithmeticType, class RealType> // negative_binomial distribution (k, n, x)
    // Probability of number of events between 0 and k-1 inclusive, if expected probability of success events is success_probability.
      inline typename tools::promote_arg3<ArithmeticType, ArithmeticType, RealType>::type
      // Return type is the wider of the two (perhaps promoted) floating-point types.
      negative_binomial(ArithmeticType k, ArithmeticType n, RealType success_probability)
    { 
      typedef typename tools::promote_arg3<ArithmeticType, ArithmeticType, RealType>::type promote_type; // Arguments type.
      return detail::negative_binomial_imp(static_cast<promote_type>(k), static_cast<promote_type>(n), static_cast<promote_type>(success_probability));
    } // negative_binomial

    template <class ArithmeticType, class RealType> // negative_binomial distribution complement (k, n, x)
    // Probability of number of events between 0 and k-1 inclusive, if expected mean is x.
      inline typename tools::promote_arg3<ArithmeticType, ArithmeticType, RealType>::type
      // Return type is the wider of the two (perhaps promoted) floating-point types.
      negative_binomial_c(ArithmeticType k, ArithmeticType n, RealType success_probability)
    { 
      typedef typename tools::promote_arg3<ArithmeticType, RealType, ArithmeticType>::type promote_type; // Arguments type.
      return detail::negative_binomial_c_imp(static_cast<promote_type>(k), static_cast<promote_type>(n), static_cast<promote_type>(success_probability));
    } // negative_binomial

    template <class ArithmeticType, class RealType>
    // success_probability if number of events between 0 and k-1 inclusive, and probability p.
    inline typename tools::promote_arg3<ArithmeticType, ArithmeticType, RealType>::type
    // Return type is the wider of the two (perhaps promoted) floating point types.
      negative_binomial_inv(ArithmeticType k,  ArithmeticType n, RealType probability)
      // negative_binomial distribution inverse (k, n, p)
    { 
      typedef typename tools::promote_arg3<ArithmeticType, ArithmeticType, RealType>::type promote_type; // Arguments type.
      return detail::negative_binomial_inv_imp(static_cast<promote_type>(k), static_cast<promote_type>(n),static_cast<promote_type>(probability));
    } // negative_binomial_inv


    // http://en.wikipedia.org/wiki/negative_binomial_distribution

    // Given a sample of N measured values k[i],
    // we wish to estimate the value of the parameter x (mean)
    // of the negative_binomial population from which the sample was drawn. 
    // To calculate the maximum likelihood value = 1/N sum i = 1 to N of k[i] 

    // Also could get k from probability and mean x???  TODO

    // Also may want a function for EXACTLY k.

    // And probability that there are EXACTLY k occurrences is
    // exp(-x) * pow(x, k) / factorial(k)
    // where x is expected occurrences (mean) during the given interval.
    // For example, if events occur, on average, every 4 min,
    // and we are interested in number of events occurring in 10 min,
    // then x = 10/4 = 2.5

	} // namespace math
} // namespace boost

#endif // BOOST_MATH_SPECIAL_NEGATIVE_BINOMIAL_HPP
