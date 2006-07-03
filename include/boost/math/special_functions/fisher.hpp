// Copyright Paul A. Bristow 2006.

// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

// boost\math\special_functions\fisher.hpp

// Fisher-Snedecor distribution
// (named after Sir R.A. Fisher and George W. Snedecor).

// http://en.wikipedia.org/wiki/Fisher-Snedecor_distribution

// The Fisher Cumulative Distribution Function (CDF) is given
// by an incomplete beta function.

#ifndef BOOST_MATH_SPECIAL_FISHER_HPP
#define BOOST_MATH_SPECIAL_FISHER_HPP

#include <boost/math/special_functions/beta.hpp> // for ibeta(a, b, x)
#include <boost/math/tools/roots.hpp> // for domain_error.
#include <boost/math/tools/promotion.hpp> // for promotion.

// The F-distribution is relevant when we try to calculate the ratios of variances 
// (assuming normally distributed data) the ratio F = variance1/variance2.
// By convention, variances 1 and 2 are chosen so that F >= 1.

// For very large values of both degrees_of_freedom1 and degrees_of_freedom2,
// greater than 10^5, a normal approximation is used.
// If only one of degrees_of_freedom1 or degrees_of_freedom2
// is greater than 10^5 then a chi sq 2 approximation is used,
// see Abramowitz and Stegun (1965) Eq 26.6.12 and 26.6.13 on p 947.
// TODO?

namespace boost
{
	namespace math
	{
    namespace detail
    {  // Implementations called by actual functions 
       // with arguments that, if necessary,
       // have been promoted from ArithmeticType to RealType.

		  template <class ArithmeticType, class RealType>
		  RealType fisher_imp(ArithmeticType degrees_of_freedom1, ArithmeticType degrees_of_freedom2, RealType x)
		  { // Implementation of Probability of Fisher x.
        // Returns the area from zero to x under the F density function
        // (also known as Snedcor's density or the variance ratio density).
        // This is the density of x = (u1/df1)/(u2/df2),
        // where u1 and u2 are random variables having Chi square distributions
        // with df1 and df2 degrees of freedom, respectively.

        using boost::math::tools::domain_error;
        using boost::math::tools::logic_error;
	      using boost::math::ibeta; // ibeta(a, b, x)
			  // ArithmeticType degrees of freedom arguments may be positive integral, signed, or unsigned, or floating-point.
			  if(degrees_of_freedom1 <= ArithmeticType(0))
			  { // Degrees of freedom must not be negative!
				  domain_error<ArithmeticType>(BOOST_CURRENT_FUNCTION, "degrees of freedom 1 argument is %1%, but must be > 0 !");
			  }
			  if(degrees_of_freedom2 <= ArithmeticType(0))
			  { // Degrees of freedom must not be negative!
				  domain_error<ArithmeticType>(BOOST_CURRENT_FUNCTION, "degrees of freedom 2 argument is %1%, but must be > 0 !");
			  }

			  if(x < ArithmeticType(0))
			  { // x must not be negative!
				  domain_error<ArithmeticType>(BOOST_CURRENT_FUNCTION, "x argument is %1%, but must be >= 0!");
			  }

			  RealType z = (degrees_of_freedom1 * x)/(degrees_of_freedom1 * x + degrees_of_freedom2);

			  // Calculate probability of Fisher x using the incomplete beta function.
			  RealType probability = ibeta(degrees_of_freedom1 / 2, degrees_of_freedom2 / 2, z);
			  // Expect 0 <= probability result <= 1.
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
	
			  return (z > 0) ? 1	- probability : probability;
		  } // fisher_imp

		  template <class ArithmeticType, class RealType>
		  RealType fisher_c_imp(ArithmeticType degrees_of_freedom1, ArithmeticType degrees_of_freedom2, RealType x)
		  { // Implementation of Probability of complemented Fisher x.

        using boost::math::tools::domain_error;
        using boost::math::tools::logic_error;
	      using boost::math::ibeta; // ibeta(a, b, x)

			  // ArithmeticType degrees of freedom arguments may be positive integral, signed, or unsigned, or floating-point.
			  if(degrees_of_freedom1 <= 0)
			  { // Degrees of freedom must not be negative!
				  domain_error<ArithmeticType>(BOOST_CURRENT_FUNCTION, "degrees of freedom 1 argument is %1%, but must be > 0 !");
			  }
			  if(degrees_of_freedom2 <= 0)
			  { // Degrees of freedom must not be negative!
				  domain_error<ArithmeticType>(BOOST_CURRENT_FUNCTION, "degrees of freedom 2 argument is %1%, but must be > 0 !");
			  }

			  if(x < 0)
			  { // x must not be negative!
				  domain_error<ArithmeticType>(BOOST_CURRENT_FUNCTION, "x argument is %1%, but must be >= 0!");
			  }

			  RealType z = (degrees_of_freedom2 * x)/(degrees_of_freedom2 * x + degrees_of_freedom1);

			  // Calculate probability of Fisher x using the incomplete beta function.
			  RealType probability = ibeta(degrees_of_freedom2 / 2, degrees_of_freedom1 / 2, z);
			  // Expect 0 <= probability result <= 1.
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
	
			  return (z > 0) ? 1	- probability : probability;
		  } // fisher_c_imp

		  template <class ArithmeticType, class RealType>
		  RealType fisher_inv_imp(ArithmeticType degrees_of_freedom1, ArithmeticType degrees_of_freedom2, RealType probability)
		  { // Implementation of inverse Probability of Fisher x.
        // Finds the F density argument x & degrees_of_freedom 
        // such that the integral from x to infinity of the F density
        // is equal to the given probability p.
        // This is accomplished using the inverse beta integral function and the relations:
        // For the inverse of the uncomplemented F distribution:        
        //    z = incbi( df1/2, df2/2, p )
        //    x = df2 z / (df1 (1-z)).

        using boost::math::tools::domain_error;
        using boost::math::tools::logic_error;
	      using boost::math::ibeta; // ibeta(a, b, x)

			  // ArithmeticType degrees of freedom arguments may be positive integral, signed, or unsigned, or floating-point.
			  if(degrees_of_freedom1 <= 0)
			  { // Degrees of freedom must not be negative!
				  domain_error<ArithmeticType>(BOOST_CURRENT_FUNCTION, "degrees of freedom 1 argument is %1%, but must be > 0 !");
			  }
			  if(degrees_of_freedom2 <= 0)
			  { // Degrees of freedom must not be negative!
				  domain_error<ArithmeticType>(BOOST_CURRENT_FUNCTION, "degrees of freedom 2 argument is %1%, but must be > 0 !");
			  }
        if((probability < 0) || (probability > 1))
			  { // chisqr must be > 0!
				  return domain_error<RealType>(BOOST_CURRENT_FUNCTION, "probability argument is %1%, but must be >= 0 and =< 1 !", probability);
			  }

			  // Calculate Fisher from probability & degrees_of_freedom using the incomplete beta function.
			  RealType z = ibeta(degrees_of_freedom1 / 2, degrees_of_freedom2 / 2, probability);
			  return degrees_of_freedom2 * z / (degrees_of_freedom1 * ( 1 - z));
		  } // fisher_inv_imp

		  template <class ArithmeticType, class RealType>
		  RealType fisher_c_inv_imp(ArithmeticType degrees_of_freedom1, ArithmeticType degrees_of_freedom2, RealType probability)
		  { // Implementation of inverse Probability of complemented Fisher x.
        // Finds the F density argument x & degrees_of_freedom 
        // such that the integral from x to infinity of the F density
        // is equal to the given probability p.
        // This is accomplished using the inverse beta integral function and the relations
        // For the inverse of the complemented F distribution:        
        //    z = incbi( df2/2, df1/2, p )
        //    x = df2 (1-z) / (df1 z).

        using boost::math::tools::domain_error;
        using boost::math::tools::logic_error;
	      using boost::math::ibeta; // ibeta(a, b, x)

			  // ArithmeticType degrees of freedom arguments may be positive integral, signed, or unsigned, or floating-point.
			  if(degrees_of_freedom1 <= 0)
			  { // Degrees of freedom must not be negative!
				  domain_error<ArithmeticType>(BOOST_CURRENT_FUNCTION, "degrees of freedom 1 argument is %1%, but must be > 0 !");
			  }
			  if(degrees_of_freedom2 <= 0)
			  { // Degrees of freedom must not be negative!
				  domain_error<ArithmeticType>(BOOST_CURRENT_FUNCTION, "degrees of freedom 2 argument is %1%, but must be > 0 !");
			  }
        if((probability < 0) || (probability > 1))
			  { // chisqr must be > 0!
				  return domain_error<RealType>(BOOST_CURRENT_FUNCTION, "probability argument is %1%, but must be >= 0 and =< 1 !", probability);
			  }

			  // Calculate Fisher from probability & degrees_of_freedom using the incomplete beta function.
			  RealType z = ibeta(degrees_of_freedom2 / 2, degrees_of_freedom1 / 2, probability);
	
			  return degrees_of_freedom2 * (1 - z) / (degrees_of_freedom1 * z);
		  } // fisher_c_inv_imp

    } //namespace detail

		template <class ArithmeticType, class RealType>
      inline typename tools::promote_arg3<ArithmeticType, ArithmeticType, RealType>::type
      // return type is the wider of the three ( possibly promoted) floating point types.
      fisher(ArithmeticType degrees_of_freedom1, ArithmeticType degrees_of_freedom2, RealType fisher)
		{ 
         typedef typename tools::promote_arg3<ArithmeticType, ArithmeticType, RealType>::type promote_type; // Arguments type.
         return detail::fisher_imp(static_cast<promote_type>(degrees_of_freedom1), static_cast<promote_type>(degrees_of_freedom2),static_cast<promote_type>(fisher));
		} // fisher(df1, df2, x)

		template <class ArithmeticType, class RealType>
      inline typename tools::promote_arg3<ArithmeticType, ArithmeticType, RealType>::type
      fisher_c(ArithmeticType degrees_of_freedom1, ArithmeticType degrees_of_freedom2, RealType fisher)
		{ 
         typedef typename tools::promote_arg3<ArithmeticType, ArithmeticType, RealType>::type promote_type; // Arguments type.
         return detail::fisher_imp_c(static_cast<promote_type>(degrees_of_freedom1), static_cast<promote_type>(degrees_of_freedom2),static_cast<promote_type>(fisher));
		} // fisher_c(df1, df2, x)

		template <class ArithmeticType, class RealType>
      inline typename tools::promote_arg3<ArithmeticType, ArithmeticType, RealType>::type
      fisher_inv(ArithmeticType degrees_of_freedom1, ArithmeticType degrees_of_freedom2, RealType fisher)
		{ 
         typedef typename tools::promote_arg3<ArithmeticType, ArithmeticType, RealType>::type promote_type; // Arguments type.
         return detail::fisher_inv_imp(static_cast<promote_type>(degrees_of_freedom1), static_cast<promote_type>(degrees_of_freedom2),static_cast<promote_type>(fisher));
		} // fisher_inv(df1, df2, probability)

		template <class ArithmeticType, class RealType>
      inline typename tools::promote_arg3<ArithmeticType, ArithmeticType, RealType>::type
      fisher_c_inv(ArithmeticType degrees_of_freedom1, ArithmeticType degrees_of_freedom2, RealType fisher)
		{ 
         typedef typename tools::promote_arg3<ArithmeticType, ArithmeticType, RealType>::type promote_type; // Arguments type.
         return detail::fisher_c_inv_imp(static_cast<promote_type>(degrees_of_freedom1), static_cast<promote_type>(degrees_of_freedom2),static_cast<promote_type>(fisher));
		} // fisher_c_inv(df1, df2, probability)

	} // namespace math
} // namespace boost

#endif // BOOST_MATH_SPECIAL_FISHER_HPP
