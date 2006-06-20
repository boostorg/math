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

#include <boost/math/special_functions/beta.hpp> // for ibeta.
#include <boost/math/tools/roots.hpp> // for domain_error.

// The F-distribution is relevant when we try to calculate the ratios of variances 
// (assuming normally distributed data) the ratio F = variance1/variance2.
// By convention, variances 1 and 2 are chosen so that F >= 1.

// For very large values of both degrees_of_freedom1 and degrees_of_freedom2, greater than 10^5, a normal approximation is used.
// If only one of degrees_of_freedom1 or degrees_of_freedom2 is greater than 10^5 then a chi sq 2 approximation is used,
// see Abramowitz and Stegun (1965).
// TODO???

namespace boost
{
	namespace math
	{
		//	template <class DFT, class FPT> FPT students_t(DFT degrees_of_freedom1, DFT degrees_of_freedom2, FPT x); // Declaration.

		template <class DFT, class FPT>
		FPT fisher(DFT degrees_of_freedom1, DFT degrees_of_freedom2, FPT x)
		{ // Implementation of Probability of Fisher x.
			using boost::math::ibeta; // ibeta(a, b, x)
			using boost::math::tools::domain_error;

			// Degrees of freedom argument may be integral, signed, or unsigned, or floating-point.
			if(degrees_of_freedom1 <= DFT(0))
			{ // Degrees of freedom must not be negative!
				domain_error<DFT>(BOOST_CURRENT_FUNCTION, "<= 0 degrees of freedom 1 argument to the Fisher function!");
			}
			if(degrees_of_freedom2 <= DFT(0))
			{ // Degrees of freedom must not be negative!
				domain_error<DFT>(BOOST_CURRENT_FUNCTION, "<= 0 degrees of freedom 2 argument to the Fisher function!");
			}

			if(x < DFT(0))
			{ // x must not be negative!
				domain_error<DFT>(BOOST_CURRENT_FUNCTION, "< 0 x argument to the Fisher function!");
			}

			FPT dd1 = static_cast<FPT>(degrees_of_freedom1); // In case degrees_of_freedom1 was an integral type.
			FPT dd2 = static_cast<FPT>(degrees_of_freedom2); // In case degrees_of_freedom2 was an integral type.
			FPT z = (dd1 * x)/(dd1 * x + dd2);

			// Calculate probability of Fisher x using the incomplete beta function.
			FPT result = ibeta(dd1 * static_cast<FPT>(0.5), dd2 * static_cast<FPT>(0.5), z);
			// Expect 0 <= probability result <= 1.
			BOOST_ASSERT(result >= static_cast<FPT>(0.));
			BOOST_ASSERT(result <= static_cast<FPT>(1.));
	  	// Numerical errors might cause it to drift slightly < 0 or > 1.
	  	// This might cause trouble later.
	  	// Might also consider constraining to fit the 0 to 1 range of probability
	  	// by changing any value just outside to the limit 0 or 1?
			return (z > static_cast<FPT>(0)) ? static_cast<FPT>(1.)	- result : result;
		} // fisher
	} // namespace math
} // namespace boost

#endif // BOOST_MATH_SPECIAL_FISHER_HPP
