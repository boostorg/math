// boost\math\special_functions\students_t.hpp

// Copyright John Maddock & Paul A. Bristow 2006.

// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

// This statistical distribution was published by W. Gosset in 1908 [23]. His employer,
// Guinness Breweries, required him to publish under a pseudonym, so he chose “Student.”
// "Student" (W.S. Gosset) (1908) The probable error of a mean. Biometrika 6(1):1--25

// http://en.wikipedia.org/wiki/Student%27s_t-distribution

// The cumulative distribution function is given by an incomplete beta function.
// Abramowitz and Stegun, formula 26.5.27 (1966)

// A recent review of calculating t quantiles is at:
// http://www.mth.kcl.ac.uk/~shaww/web_page/papers/Tdistribution06.pdf

// A lookup table of quantiles of the t distribution
// for 1 to 25 in steps of 0.1 is provided in CSV form at:
// www.mth.kcl.ac.uk/~shaww/web_page/papers/Tsupp/tquantiles.csv

#ifndef BOOST_MATH_SPECIAL_STUDENTS_T_HPP
#define BOOST_MATH_SPECIAL_STUDENTS_T_HPP

#include <boost/math/special_functions/beta.hpp> // for ibeta(a, b, x).
#include <boost/math/tools/roots.hpp> // for domain_error & logic_error.
#include <boost/math/tools/promotion.hpp> // for promotion.

namespace boost
{
	namespace math
	{ // Forward declaration of students_t.
    template <class ArithmeticType, class RealType>
    typename tools::promote_arg2<RealType, ArithmeticType>::type
     // return type is the wider of the two (?promoted) floating point types.
    students_t(ArithmeticType degrees_of_freedom, RealType t);
  }
}

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
		    RealType students_t_imp(RealType degrees_of_freedom, RealType t)
		    { // Implementation of Probability of Student's t.
			    // Abramowitz and Stegun, formula 26.5.27 (1966)

	        using boost::math::ibeta; // ibeta(a, b, x)
	        using boost::math::tools::domain_error;
	        using boost::math::tools::logic_error;
          using std::numeric_limits;

			    // Degrees_of_freedom argument may be integral, signed, or unsigned, or floating point.
          // If necessary, it has already been promoted from an integral type.
			    if(degrees_of_freedom <= 0)
			    { // Degrees of freedom must be > 0!
				    return domain_error<RealType>(BOOST_CURRENT_FUNCTION, "degrees of freedom argument is %1%, but must be > 0 !", degrees_of_freedom);
			    }
			    RealType z = degrees_of_freedom / (degrees_of_freedom + t * t);
			    // Calculate probability of Student's t using the incomplete beta function.
			    // probability = ibeta(degrees_of_freedom/2, 1/2, degrees_of_freedom/ (degrees_of_freedom + t*t))
          RealType probability = ibeta(degrees_of_freedom / 2, static_cast<RealType>(0.5), z) / 2;
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
			    return (t > 0 ? 1	- probability : probability);
		    } // students_t_imp

 		  template <class RealType>
		  RealType students_t_inv_imp(RealType degrees_of_freedom, RealType probability)
		  { // Inverse cumulative Distribution Function or Quintile (percentile / 100) of Student's t Probability.
           using boost::math::ibeta_inv; // ibeta(a, b, x)
           using boost::math::tools::domain_error;
           using std::numeric_limits;
           using namespace std;  // for fabs

			  // Degrees of freedom argument may be integral, signed, or unsigned, or floating point.
			  if(degrees_of_freedom <= 0)
			  { // Degrees of freedom must be > 0!
				  return domain_error<RealType>(BOOST_CURRENT_FUNCTION, "degrees of freedom argument is %1%, but must be > 0 !", degrees_of_freedom);
			  }
			  if((probability < 0) || (probability > 1))
			  { // probability must be >= 0 and <= 1!
				  return domain_error<RealType>(BOOST_CURRENT_FUNCTION, "probability argument is %1%, but must be >= 0 and <= 1 !", probability);
			  }
			  // Special cases, regardless of degrees_of_freedom.
			  if (probability == 0)
				  return -numeric_limits<RealType>::infinity();
			  if (probability == 1)
				  return +numeric_limits<RealType>::infinity();
			  if (probability == static_cast<RealType>(0.5))
				  return 0; 
			  // Calculate quantile of Student's t using the incomplete beta function.
			  if ((probability > RealType(0.25)) && (probability < RealType(0.75)) )
			  { // probability is middling.
				  RealType z = 1 - 2 * probability;
				  z = ibeta_inv(RealType(0.5), degrees_of_freedom / 2 , fabs(z));
				  RealType t = sqrt(degrees_of_freedom * z / (1 - z));
				  return (probability < RealType(0.5)) ? -t : t;
			  }
			  else
			  { // probability is small or large.
				  int sign = -1;
				  if (probability >= 0.5)
				  { // large and > 0.75
					  probability = 1 - probability;
					  sign = +1;
				  }
				  RealType z = ibeta_inv(degrees_of_freedom / 2, static_cast<RealType>(0.5), 2 * probability);
				  if (((numeric_limits<RealType>::max)() * z) < degrees_of_freedom)
				  {
					  return sign * (numeric_limits<RealType>::max)();
				  }
				  else
				  {
					  return sign * sqrt(degrees_of_freedom/ z - degrees_of_freedom);
				  }
			  }
		  } // students_t_inv_imp
    } // namespace detail

    template <class ArithmeticType, class RealType>
      inline typename tools::promote_arg2<RealType, ArithmeticType>::type // return type is the wider of the two (?promoted) floating point types.
         students_t(ArithmeticType degrees_of_freedom, RealType t)
		{ 
         typedef typename tools::promote_arg2<RealType, ArithmeticType>::type promote_type; // Arguments type.
         return detail::students_t_imp(static_cast<promote_type>(degrees_of_freedom), static_cast<promote_type>(t));
		} // students_t

    template <class ArithmeticType, class RealType>
    inline typename tools::promote_arg2<RealType, ArithmeticType>::type // return type is the wider of the two (?promoted) floating-point types.
    students_t_inv(ArithmeticType degrees_of_freedom, RealType probability)
		{ 
         typedef typename tools::promote_arg2<RealType, ArithmeticType>::type promote_type; // Arguments type.
         return detail::students_t_inv_imp(static_cast<promote_type>(degrees_of_freedom), static_cast<promote_type>(probability));
		} // students_t_inv

	} // namespace math
} // namespace boost

#endif // BOOST_MATH_SPECIAL_STUDENTS_T_HPP
