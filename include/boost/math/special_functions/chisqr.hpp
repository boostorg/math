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

#include <boost/type_traits/is_floating_point.hpp>
#include <boost/type_traits/is_integral.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/mpl/if.hpp>

namespace boost
{
	namespace math
	{
      namespace detail
			{
				// If either T or U is an integer type, 
				// pretend it was a double (for the purposes of further analysis).
				// Then pick the wider of the two floating-point types
				// as the actual signature to forward to.
				// For example:
				// foo(int, short) -> foo(double, double);
				// foo(int, float) -> foo(double, double);
				// foo(int, double) -> foo(double, double);
				// foo(double, float) -> foo(double, double);
				// foo(any-int-or-float-type, long double) -> foo(long double, long double);
				// but float foo(float, float) is unchanged.

         template <class T>
         struct promote_arg
         { // if T integral, then promote to double.
            typedef typename mpl::if_<is_integral<T>, double, T>::type type;
         };

         template <class T, class U>
         struct promote_arg2 
         { // Promote, if necessary, & pick the wider of the two floating-point types.
           // for both parameter types, if integral promote to double.
            typedef typename promote_arg<T>::type TP; // Perhaps promoted.
            typedef typename promote_arg<U>::type UP; // Perhaps promoted.

            typedef typename mpl::if_c<
               is_floating_point<TP>::value && is_floating_point<UP>::value,
               typename mpl::if_c<
                  is_same<long double, TP>::value || is_same<long double, UP>::value,
                  long double,
                  typename mpl::if_c<
                     is_same<double, TP>::value || is_same<double, UP>::value,
                     double,
                     float
                  >::type
               >::type,
               typename mpl::if_<is_convertible<TP, UP>, UP, TP>::type
            >::type type;
         }; // promote_arg2

		   template <class FPT>
		   FPT chisqr_imp(FPT degrees_of_freedom, FPT chisqr)
		   { 
        // Implementation of Probability of CHISQR chisqr.

			   using boost::math::gamma_Q; // gamma_Q(degrees_of_freedom/2, chisqr/2)
			   using boost::math::tools::domain_error;
			   using boost::math::tools::logic_error;

			   // Degrees of freedom argument may be integral, signed or unsigned, or floating-point, or User Defined real.
			   if(degrees_of_freedom <= 0)
			   { // Degrees of freedom must be > 0!
				   return domain_error<FPT>(BOOST_CURRENT_FUNCTION, "degrees of freedom argument is %1%, but must be > 0 !", degrees_of_freedom);
			   }
			   if(chisqr < 0)
			   { // chisqr must be > 0!
				   return domain_error<FPT>(BOOST_CURRENT_FUNCTION, "chisqr argument is %1%, but must be > 0 !", chisqr);
			   }

			   // Calculate probability of chisqr using the incomplete beta function.
			   FPT probability = gamma_Q(degrees_of_freedom / 2, chisqr / 2);
			   // Expect 0 <= probability <= 1.
	  	   // Numerical errors might cause probability to be slightly outside the range < 0 or > 1.
	  	   // This might cause trouble downstream, so warn, possibly throw exception, but constrain to the limits.
			   if (probability < 0)
			   {
				   logic_error<FPT>(BOOST_CURRENT_FUNCTION, "probability %1% is < 0, so has been constrained to zero !", probability);
				   return 0; // Constrain to zero, if logic_error does not throw.
			   }
			   if(probability > 1)
			   {
				   logic_error<FPT>(BOOST_CURRENT_FUNCTION, "probability %1% is > 1, so has been constrained to unity!", probability);
				   return 1; // Constrain to unity, if logic_error does not throw.
			   }
			   return probability;
		   } // chisqr_imp
      } // namespace detail

		template <class DFT, class FPT>
      inline typename detail::promote_arg2<FPT, DFT>::type // return type is the wider of the two (?promoted) floating point types.
         chisqr(DFT degrees_of_freedom, FPT chisqr)
		{ 
         typedef typename detail::promote_arg2<FPT, DFT>::type promote_type; // Arguments type.
         return detail::chisqr_imp(static_cast<promote_type>(degrees_of_freedom), static_cast<promote_type>(chisqr));
		} // chisqr
	} // namespace math
} // namespace boost

#endif // BOOST_MATH_SPECIAL_CHISQR_HPP
