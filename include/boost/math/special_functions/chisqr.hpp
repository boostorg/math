// chisqr.hpp    chisqr cumulative distribution function.

// Copyright John Maddock & Paul A. Bristow 2006.

// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

// The Chisqr Cumulative Distribution Function (CDF) is given
// by an incomplete gamma integral function.

// http://en.wikipedia.org/wiki/Chi-squared_distribution

// The chisqr-distribution is used for statistical testing:

// http://en.wikipedia.org/wiki/Pearson%27s_chi-square_test

// "Pearson's is one of a variety of chi-square tests –
// statistical procedures whose results are evaluated
// by reference to the chi-square distribution.
// It tests a null hypothesis that the relative frequencies 
// of occurrence of observed events follow a specified frequency distribution.
// The events must be mutually exclusive.
// One of the simplest examples is the hypothesis that 
// an ordinary six-sided die is "fair": all six outcomes occur equally often.
// Chi-square is calculated by finding the difference between
// each observed and theoretical frequency, squaring them,
// dividing each by the theoretical frequency, and taking the sum of the results:
//   chi^2 = sum_{i=1 to 6} {(O_i - E_i)^2 / E_i}
// where:
//    O_i = an observed frequency
//    E_i = an expected (theoretical) frequency, asserted by the null hypothesis.

// For very large values of both degrees_of_freedom1 and degrees_of_freedom2,
// greater than 10^5, a normal approximation is used.
// If only one of degrees_of_freedom1 or degrees_of_freedom2 is
// greater than 10^5 then a chisqr approximation is used,
// see Abramowitz and Stegun (1965).
// TODO?

#ifndef BOOST_MATH_SPECIAL_CHISQR_HPP
#define BOOST_MATH_SPECIAL_CHISQR_HPP

#include <boost/math/special_functions/gamma.hpp> // for gamma.
#include <boost/math/tools/roots.hpp>
#include <boost/math/tools/error_handling.hpp> // for domain_error, logic_error.
#include <boost/math/tools/promotion.hpp> // for promotion

#include <boost/type_traits/is_floating_point.hpp>
#include <boost/type_traits/is_integral.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/mpl/if.hpp>

namespace boost
{
	namespace math
	{
      namespace detail
			{ // Implementations called by actual functions 
        // with arguments that, if necessary,
        // have been promoted from ArithmeticType to RealType.

		   template <class RealType>
		   RealType chisqr_imp(RealType degrees_of_freedom, RealType chisqr)
		   { 
         // Implementation of probability of chisqr.
         // Returns the area under the left hand tail (from 0 to x)
         // of the Chi square probability density function
         // with v degrees of freedom.

			   using boost::math::gamma_Q; // gamma_Q(degrees_of_freedom/2, chisqr/2)
			   using boost::math::tools::domain_error;
			   using boost::math::tools::logic_error;

			   // Degrees of freedom argument may be integral, signed or unsigned, or floating-point, or User Defined real.
			   if(degrees_of_freedom <= 0)
			   { // Degrees of freedom must be > 0!
				   return domain_error<RealType>(BOOST_CURRENT_FUNCTION, "degrees of freedom argument is %1%, but must be > 0 !", degrees_of_freedom);
			   }
			   if(chisqr < 0)
			   { // chisqr must be > 0!
				   return domain_error<RealType>(BOOST_CURRENT_FUNCTION, "chisqr argument is %1%, but must be > 0 !", chisqr);
			   }

			   // Calculate probability of chisqr using the incomplete gamma integral function.
			   RealType probability = gamma_Q(degrees_of_freedom / 2, chisqr / 2);
			   // Expect 0 <= probability <= 1.
	  	   // Numerical errors might cause probability to be slightly outside the range < 0 or > 1.
	  	   // This might cause trouble downstream, so warn, possibly throw exception, but constrain to the limits.
			   if (probability < 0)
			   {
				   logic_error<RealType>(BOOST_CURRENT_FUNCTION, "probability %1% is < 0, so has been constrained to zero !", probability);
				   return 0; // Constrain to zero, if logic_error does not throw.
			   }
			   if(probability > 1)
			   {
				   logic_error<RealType>(BOOST_CURRENT_FUNCTION, "probability %1% is > 1, so has been constrained to unity!", probability);
				   return 1; // Constrain to unity, if logic_error does not throw.
			   }
			   return probability;
		   } // chisqr_imp

		   template <class RealType>
		   RealType chisqr_c_imp(RealType degrees_of_freedom, RealType chisqr)
		   { 
         // Implementation of probability of chisqr complemented.
         // Returns the area under the right hand tail (from x to infinity)
         // of the chisqr probability density function
         // with v degrees of freedom:

			   using boost::math::gamma_Q; // gamma_Q(degrees_of_freedom/2, chisqr/2)
			   using boost::math::tools::domain_error;
			   using boost::math::tools::logic_error;

			   // Degrees of freedom argument may be integral, signed or unsigned, or floating-point, or User Defined real.
			   if(degrees_of_freedom <= 0)
			   { // Degrees of freedom must be > 0!
				   return domain_error<RealType>(BOOST_CURRENT_FUNCTION, "degrees of freedom argument is %1%, but must be > 0 !", degrees_of_freedom);
			   }
			   if(chisqr < 0)
			   { // chisqr must be > 0!
				   return domain_error<RealType>(BOOST_CURRENT_FUNCTION, "chisqr argument is %1%, but must be > 0 !", chisqr);
			   }

			   // Calculate probability of chisqr using the incomplete beta function.
			   RealType probability = gamma_Q(degrees_of_freedom / 2, chisqr / 2);
			   // Expect 0 <= probability <= 1.
	  	   // Numerical errors might cause probability to be slightly outside the range < 0 or > 1.
	  	   // This might cause trouble downstream, so warn, possibly throw exception, but constrain to the limits.
			   if (probability < 0)
			   {
				   logic_error<RealType>(BOOST_CURRENT_FUNCTION, "probability %1% is < 0, so has been constrained to zero !", probability);
				   return 0; // Constrain to zero, if logic_error does not throw.
			   }
			   if(probability > 1)
			   {
				   logic_error<RealType>(BOOST_CURRENT_FUNCTION, "probability %1% is > 1, so has been constrained to unity!", probability);
				   return 1; // Constrain to unity, if logic_error does not throw.
			   }
			   return probability;
		   } // chisqr_imp

		   template <class RealType>
		   RealType chisqr_inv_imp(RealType degrees_of_freedom, RealType probability)
		   { 
         // Implementation of inverse of chisqr distribution.
         // Finds the chisqr argument x such that the integral
         // from x to infinity of the chisqr density 
         // is equal to the given cumulative probability y.

			   using boost::math::gamma_Q; // gamma_Q(degrees_of_freedom/2, chisqr/2)
			   using boost::math::tools::domain_error;
			   using boost::math::tools::logic_error;

			   // Degrees of freedom argument may be integral, signed or unsigned, or floating-point, or User Defined real.
			   if(degrees_of_freedom <= 0)
			   { // Degrees of freedom must be > 0!
				   return domain_error<RealType>(BOOST_CURRENT_FUNCTION, "degrees of freedom argument is %1%, but must be > 0 !", degrees_of_freedom);
			   }
			   if((probability < 0) || (probability > 1))
			   { // probability must be > 0!
				   return domain_error<RealType>(BOOST_CURRENT_FUNCTION, "probability argument is %1%, but must be >= 0 and =< 1 !", probability);
			   }

			   // Calculate chisqr from probability & degrees_of_freedom using the inverse gamma integral function.
			   RealType chisqr = gamma_Q_inv(degrees_of_freedom / 2, probability) * 2;
			   // Expect chisqr > 0.
	  	   // Numerical errors might cause probability to be slightly outside the range < 0 or > 1.
	  	   // This might cause trouble downstream, so warn, possibly throw exception, but constrain to the limits.
			   if (chisqr < 0)
			   {
				   logic_error<RealType>(BOOST_CURRENT_FUNCTION, "chisqr %1% is < 0, so has been constrained to zero !", probability);
				   return 0; // Constrain to zero, if logic_error does not throw.
			   }
			   return probability;
		   } // chisqr_inv_imp

		   template <class RealType>
		   RealType chisqr_df_inv_imp(RealType chisqr, RealType probability)
		   { 
         // Implementation of inverse (complemented?) chisqr distribution.
         // Finds the degress_fo_freedom argument x such that the 
         // integral from x to infinity of the chisqr density 
         // is equal to the given cumulative probability y.

			   using boost::math::gamma_Q_inv; // gamma_Q(degrees_of_freedom/2, chisqr/2)
			   using boost::math::tools::domain_error;
			   using boost::math::tools::logic_error;

			   // Degrees of freedom argument may be integral, signed or unsigned, or floating-point, or User Defined real.
			   if(chisqr <= 0)
			   { // chisqr must be > 0!
				   return domain_error<RealType>(BOOST_CURRENT_FUNCTION, "chisqr argument is %1%, but must be > 0 !", chisqr);
			   }
			   if ((probability < 0) || (probability > 1))
			   { // chisqr must be > 0!
				   return domain_error<RealType>(BOOST_CURRENT_FUNCTION, "probability argument is %1%, but must be >= 0 and =< 1 !", chisqr);
			   }

			   // Calculate degrees_of_freedom from probability & chisqr using the inverse gamma integral function??
			   RealType degrees_of_freedom = gamma_Q_inv(chisqr / 2, probability) * 2;// TODO <<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
			   // Expect degrees_of_freedom >= 0.
			   if (degrees_of_freedom < 0)
			   {
				   logic_error<RealType>(BOOST_CURRENT_FUNCTION, "degrees_of_freedom %1% is < 0, so has been constrained to zero !", probability);
				   return 0; // Constrain to zero, if logic_error does not throw.
			   }
			   return probability;
		   } // chisqr_inv_imp
     } // namespace detail

  	template <class ArithmeticType, class RealType> // probability chisqr(degrees_of_freedom, chisqr)
    inline typename tools::promote_arg2<RealType, ArithmeticType>::type
      // return type is the wider of the two (?promoted) floating point types.
         chisqr(ArithmeticType degrees_of_freedom, RealType chisqr)
		{ 
         typedef typename tools::promote_arg2<RealType, ArithmeticType>::type promote_type; // Arguments type.
         return detail::chisqr_imp(static_cast<promote_type>(degrees_of_freedom), static_cast<promote_type>(chisqr));
		} // chisqr

		template <class ArithmeticType, class RealType> // complement probability chisqr_c(degrees_of_freedom, chisqr)
      inline typename tools::promote_arg2<RealType, ArithmeticType>::type
         chisqr_c(ArithmeticType degrees_of_freedom, RealType chisqr)
		{ 
         typedef typename tools::promote_arg2<RealType, ArithmeticType>::type promote_type; // Arguments type.
         return detail::chisqr_c_imp(static_cast<promote_type>(degrees_of_freedom), static_cast<promote_type>(chisqr));
		} // chisqr_c

		template <class ArithmeticType, class RealType> // chisqr = chisqr_inv(degrees_of_freedom,  probability)
      inline typename tools::promote_arg2<RealType, ArithmeticType>::type
         chisqr_inv(ArithmeticType degrees_of_freedom, RealType probability)
		{ 
         typedef typename tools::promote_arg2<RealType, ArithmeticType>::type promote_type; // Arguments type.
         return detail::chisqr_inv_imp(static_cast<promote_type>(degrees_of_freedom), static_cast<promote_type>(probability));
		} // chisqr_inv

		template <class ArithmeticType, class RealType> // degrees_of_freedom = chisqr_inv_df(chisqr, probability)
      inline typename tools::promote_arg2<RealType, RealType>::type // <RealType, ArithmeticType> but both RealType.
         chisqr_df_inv(RealType chisqr, ArithmeticType probability)
		{ 
         typedef typename tools::promote_arg2<RealType, ArithmeticType>::type promote_type; // Arguments type.
         return detail::chisqr_df_inv_imp(static_cast<promote_type>(chisqr), static_cast<promote_type>(probability));
		} // chisqr_inv_df

	} // namespace math
} // namespace boost

#endif // BOOST_MATH_SPECIAL_CHISQR_HPP
