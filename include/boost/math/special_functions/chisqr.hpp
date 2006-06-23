// chisqr.hpp    chisqr cumulative distribution function.

// Copyright Paul A. Bristow 2006.

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
// statistical procedures whose results are evaluated by referenceto the chi-square distribution.
// It tests a null hypothesis that the relative frequencies of occurrence of observed events
// follow a specified frequency distribution.
// The events must be mutually exclusive
// One of the simplest examples is the hypothesis that an ordinary six-sided die is "fair",
//	i.e., all six outcomes occur equally often.
//	Chi-square is calculated by finding the difference between
// each observed and theoretical frequency, squaring them, dividing each by the theoretical frequency,
// and taking the sum of the results:
//   chi^2 = sum_{i=1 to 6} {(O_i - E_i)^2 / E_i}
// where:
//    O_i = an observed frequency
//    E_i = an expected (theoretical) frequency, asserted by the null hypothesis.

// For very large values of both degrees_of_freedom1 and degrees_of_freedom2,
// greater than 10^5, a normal approximation is used.
// If only one of degrees_of_freedom1 or degrees_of_freedom2 is
// greater than 10^5 then a chisqr approximation is used,
// see Abramowitz and Stegun (1965).
// TODO???

#ifndef BOOST_MATH_SPECIAL_CHISQR_HPP
#define BOOST_MATH_SPECIAL_CHISQR_HPP

#include <boost/math/special_functions/gamma.hpp> // for gamma.
	using boost::math::gamma_Q; //

#include <boost/math/tools/roots.hpp>
#include <boost/math/tools/error_handling.hpp> // for domain_error.
	using boost::math::tools::domain_error;
	using boost::math::tools::overflow_error;
	using boost::math::tools::underflow_error;

namespace boost
{
	namespace math
	{
		template <class DFT, class FPT>
		FPT chisqr(DFT degrees_of_freedom, FPT chisqr); // Declaration.
		// degrees_of_freedom DFT may be integral or floating-point.

		template <class DFT, class FPT>
		FPT chisqr(DFT degrees_of_freedom, FPT chisqr)
		{ // Implementation of Probability of CHISQR chisqr.
			using boost::math::gamma_Q; // gamma_Q(degrees_of_freedom/2, chisqr/2)
			using boost::math::tools::domain_error;

			// Degrees of freedom argument may be integral, signed, or unsigned, or floating-point.
			if(degrees_of_freedom <= DFT(0.))
			{ // Degrees of freedom must be > 0!
				return domain_error<DFT>(BOOST_CURRENT_FUNCTION, "degrees of freedom argument is %1%, but must be > 0 !", degrees_of_freedom);
			}
			if(chisqr < FPT(0.))
			{ // chisqr must be > 0!
				return domain_error<FPT>(BOOST_CURRENT_FUNCTION, "chisqr argument is %1%, but must be > 0 !", chisqr);
			}

			// FPT degrees_of_freedom = static_cast<FPT>(degrees_of_freedom); // In case degrees_of_freedom was an integral type.

			// Calculate probability of chisqr using the incomplete beta function.
			FPT probability = gamma_Q(static_cast<FPT>(degrees_of_freedom) * static_cast<FPT>(0.5), chisqr * static_cast<FPT>(0.5));
			// Expect 0 <= probability <= 1.

	  	// Numerical errors might cause probability to be slightly outside the range < 0 or > 1.
	  	// This might cause trouble downstream, so warn, possibly throw exception, but constrain to the limits.
			if (probability < static_cast<FPT>(0.))
			{
				domain_error<FPT>(BOOST_CURRENT_FUNCTION, "probability %1% is < 0, so has been constrained to zero !", probability);
				return static_cast<FPT>(0.); // Constrain to zero.
			}
			if(probability > static_cast<FPT>(1.))
			{
				domain_error<FPT>(BOOST_CURRENT_FUNCTION, "probability %1% is > 1, so has been constrained to unity!", probability);
				return static_cast<FPT>(1.); // Constrain to unity.
			}

			return probability;
		//	return (chisqr > static_cast<FPT>(0)) ? static_cast<FPT>(1.)	- probability : probability;
		} // chisqr
	} // namespace math
} // namespace boost

#endif // BOOST_MATH_SPECIAL_CHISQR_HPP
