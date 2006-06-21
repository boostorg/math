// Copyright Paul A. Bristow 2006.

// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

// boost\math\special_functions\students_t.hpp

// This statistical distribution was published by W. Gosset in 1908 [23]. His employer,
// Guinness Breweries, required him to publish under a pseudonym, so he chose “Student.”

// http://en.wikipedia.org/wiki/Student%27s_t-distribution

// The cumulative distribution function is given by an incomplete beta function,

// A recent review of calculating t quantiles is at:
// http://www.mth.kcl.ac.uk/~shaww/web_page/papers/Tdistribution06.pdf

#ifndef BOOST_MATH_SPECIAL_STUDENTS_T_HPP
#define BOOST_MATH_SPECIAL_STUDENTS_T_HPP

#include <boost/math/special_functions/beta.hpp> // for ibeta.
#include <boost/math/tools/roots.hpp> // for domain_error.

namespace boost
{
	namespace math
	{
		template <class DFT, class FPT>
		FPT students_t(DFT degrees_of_freedom, FPT t); // Declaration.

		template <class DFT, class FPT>
		FPT students_t(DFT degrees_of_freedom, FPT t)
		{ // Probability of Student's t implementation.
			using boost::math::ibeta; // ibeta(a, b, x)
			using boost::math::tools::domain_error;

			// Degrees of freedom argument may be integral, signed, or unsigned, or floating point.
			if(degrees_of_freedom <= DFT(0))
			{ // Degrees of freedom must not be negative!
				domain_error<DFT>(BOOST_CURRENT_FUNCTION, "Negative degrees of freedom argument to the students_t function!");
			}
			FPT dd = static_cast<FPT>(degrees_of_freedom); // In case degrees_of_freedom was an integral type.
			FPT z = dd/(dd + t * t);
			// Calculate probability of Student's t using the incomplete beta function.
			FPT result = (static_cast<FPT>(0.5) * ibeta(dd * static_cast<FPT>(0.5), static_cast<FPT>(0.5), z));
			// Check 0 <= probability result <= 1.
			BOOST_ASSERT(result >= static_cast<FPT>(0));
			BOOST_ASSERT(result <= static_cast<FPT>(1.));
	  	// Numerical errors might cause it to drift slightly < 0 or > 1.
	  	// This might cause trouble later.
	  	// Might also consider constraining to fit the range 0 to 1
	  	// by changing any value just outside to the limit 0 or 1?
			return (t > static_cast<FPT>(0)) ? static_cast<FPT>(1.)	- result : result;
		} // students_t
	} // namespace math
} // namespace boost

#endif // BOOST_MATH_SPECIAL_STUDENTS_T_HPP
