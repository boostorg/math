// boost\math\special_functions\poisson.hpp

// Copyright John Maddock & Paul A. Bristow 2006.

// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

// http://en.wikipedia.org/wiki/Poisson_distribution

// Poisson distribution is a discrete probability distribution.
// It expresses the probability of a number of events occurring in a fixed time
// if these events occur with a known average rate,
// and are independent of the time since the last event.

// The distribution was discovered by Siméon-Denis Poisson (1781–1840).

// The number of cars that pass through a certain point on a road during a given period of time.
// The number of spelling mistakes a secretary makes while typing a single page.
// The number of phone calls at a call center per minute.
// The number of times a web server is accessed per minute.
// The number of light bulbs that burn out in a certain amount of time.  
// The number of roadkill found per unit length of road

#ifndef BOOST_MATH_SPECIAL_POISSON_HPP
#define BOOST_MATH_SPECIAL_POISSON_HPP

#include <boost/math/special_functions/beta.hpp> // for ibeta(a, b, x).
#include <boost/math/tools/roots.hpp> // for domain_error & logic_error.
#include <boost/math/tools/promotion.hpp> // for promotion.

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
		    RealType poisson_imp(RealType k, RealType x)
		    { 
	        using boost::math::gamma_Q; // regularized gamma function.
	        using boost::math::tools::domain_error;
	        using boost::math::tools::logic_error;
          using std::numeric_limits;

			    // k argument may be integral, signed, or unsigned, or floating point.
          // If necessary, it has already been promoted from an integral type.
			    if(k <= 0 )
			    { // k must be >=1 !
				    return domain_error<RealType>(BOOST_CURRENT_FUNCTION, "k argument is %1%, but must be >= 1 !", k);
			    }
         RealType probability = gamma_Q(static_cast<RealType>(k), x);
			    // Check 0 <= probability probability <= 1.  
	  	    // Numerical errors might cause probability to be slightly outside the range < 0 or > 1.
	  	    // This might cause trouble downstream, so warn, possibly throw exception, but constrain to the limits.
			    if (probability < static_cast<RealType>(0.))
			    {
				    logic_error<RealType>(BOOST_CURRENT_FUNCTION, "probability %1% is < 0, so has been constrained to zero !", probability);
				    return static_cast<RealType>(0.); // Constrain to zero if logic_error does not throw.
			    }
			    if(probability > static_cast<RealType>(1.))
			    {
				    logic_error<RealType>(BOOST_CURRENT_FUNCTION, "probability %1% is > 1, so has been constrained to unity!", probability);
				    return static_cast<RealType>(1.); // Constrain to unity if logic_error does not throw.
			    }
			    return probability;
		    } // poisson_imp

 		  template <class RealType>
		  RealType poisson_inv_imp(RealType k, RealType probability)
		  { // Inverse cumulative Distribution Function or Quintile (percentile / 100) of Poisson Probability.
        // k argument may be integral, signed, or unsigned, or floating point.
        // If necessary, it has already been promoted from an integral type.

			  using boost::math::gamma_Q_inv; // gamma_Q_inv
			  using boost::math::tools::domain_error;
        using std::numeric_limits;

			  // k argument may be integral, signed, or unsigned, or floating-point.
			  if(k < 1)
			  { // Degrees of freedom must be > 0!
				  return domain_error<RealType>(BOOST_CURRENT_FUNCTION, "k argument is %1%, but must be >= 1 !", k);
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
			  // Calculate quantile of Poisson using the incomplete gamma function.
        return gamma_Q_inv(k, probability);
	    } // poisson_inv_imp
    } // namespace detail

    template <class ArithmeticType, class RealType> // Poisson distribution (k, x)
    inline // Probability of number of events between 0 and k-1 inclusive, if expected mean is x.
    typename tools::promote_arg2<RealType, ArithmeticType>::type
      // return type is the wider of the two ( perhaps promoted) floating-point types.
    poisson(ArithmeticType k, RealType x)
		{ 
         typedef typename tools::promote_arg2<RealType, ArithmeticType>::type promote_type; // Arguments type.
         return detail::poisson_imp(static_cast<promote_type>(k), static_cast<promote_type>(x));
		} // poisson_inv

    template <class ArithmeticType, class RealType>
    inline  // Mean if number of events between 0 and k-1 inclusive, and probability.
    typename tools::promote_arg2<RealType, ArithmeticType>::type
      // return type is the wider of the two ( perhaps promoted) floating point types.
    poisson_inv(ArithmeticType k, RealType probability) // Poisson distribution (k, p)
		{ 
         typedef typename tools::promote_arg2<RealType, ArithmeticType>::type promote_type; // Arguments type.
         return detail::poisson_inv_imp(static_cast<promote_type>(k), static_cast<promote_type>(probability));
		} // poisson_inv

    // http://en.wikipedia.org/wiki/Poisson_distribution

    // Given a sample of N measured values k[i],
    // we wish to estimate the value of the parameter x (mean)
    // of the Poisson population from which the sample was drawn. 
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

#endif // BOOST_MATH_SPECIAL_POISSON_HPP
