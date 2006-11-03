//  Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/math/concepts/std_real_concept.hpp>
#include <boost/math/concepts/distributions.hpp>

#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/students_t.hpp>
#include <boost/math/distributions/binomial.hpp>
#include <boost/math/distributions/cauchy.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/exponential.hpp>
#include <boost/math/distributions/extreme_value.hpp>
#include <boost/math/distributions/fisher_f.hpp>
#include <boost/math/distributions/weibull.hpp>
#include <boost/math/distributions/lognormal.hpp>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/cbrt.hpp>

//
// The purpose of this test is to verify that our code compiles
// cleanly with a type whose std lib functions are in namespace
// std and can *not* be found by ADL.  This verifies that we're
// not finding std lib functions that are in the global namespace
// for example calling ::pow(double) rather than std::pow(long double).
// This is a silent error that does the wrong thing at runtime, and
// of course we can't call std::pow() directly because we want
// the functions to be found by ADL when that's appropriate.
//
// Furthermore our code does different things internally depending
// on numeric_limits<>::digits, so there are some macros that can
// be defined that cause our concept-archetype to emulate various
// floating point types:
//
// EMULATE32: 32-bit float
// EMULATE64: 64-bit double
// EMULATE80: 80-bit long double
// EMULATE128: 128-bit long double
//
// In order to ensure total code coverage this file must be
// compiled with each of the above macros in turn, and then 
// without any of the above as well!
//

#define NULL_MACRO /**/
#ifdef EMULATE32
namespace std{
template<>
struct numeric_limits<boost::math::concepts::std_real_concept>
{
   static const bool is_specialized = true;
   static boost::math::concepts::std_real_concept min NULL_MACRO() throw();
   static boost::math::concepts::std_real_concept max NULL_MACRO() throw();
   static const int digits = 24;
   static const int digits10 = 6;
   static const bool is_signed = true;
   static const bool is_integer = false;
   static const bool is_exact = false;
   static const int radix = 2;
   static boost::math::concepts::std_real_concept epsilon() throw();
   static boost::math::concepts::std_real_concept round_error() throw();
   static const int min_exponent = -125;
   static const int min_exponent10 = -37;
   static const int max_exponent = 128;
   static const int max_exponent10 = 38;
   static const bool has_infinity = true;
   static const bool has_quiet_NaN = true;
   static const bool has_signaling_NaN = true;
   static const float_denorm_style has_denorm = denorm_absent;
   static const bool has_denorm_loss = false;
   static boost::math::concepts::std_real_concept infinity() throw();
   static boost::math::concepts::std_real_concept quiet_NaN() throw();
   static boost::math::concepts::std_real_concept signaling_NaN() throw();
   static boost::math::concepts::std_real_concept denorm_min() throw();
   static const bool is_iec559 = true;
   static const bool is_bounded = false;
   static const bool is_modulo = false;
   static const bool traps = false;
   static const bool tinyness_before = false;
   static const float_round_style round_style = round_toward_zero;
};
}
#endif
#ifdef EMULATE64
namespace std{
template<>
struct numeric_limits<boost::math::concepts::std_real_concept>
{
   static const bool is_specialized = true;
   static boost::math::concepts::std_real_concept min NULL_MACRO() throw();
   static boost::math::concepts::std_real_concept max NULL_MACRO() throw();
   static const int digits = 53;
   static const int digits10 = 15;
   static const bool is_signed = true;
   static const bool is_integer = false;
   static const bool is_exact = false;
   static const int radix = 2;
   static boost::math::concepts::std_real_concept epsilon() throw();
   static boost::math::concepts::std_real_concept round_error() throw();
   static const int min_exponent = -1021;
   static const int min_exponent10 = -307;
   static const int max_exponent = 1024;
   static const int max_exponent10 = 308;
   static const bool has_infinity = true;
   static const bool has_quiet_NaN = true;
   static const bool has_signaling_NaN = true;
   static const float_denorm_style has_denorm = denorm_absent;
   static const bool has_denorm_loss = false;
   static boost::math::concepts::std_real_concept infinity() throw();
   static boost::math::concepts::std_real_concept quiet_NaN() throw();
   static boost::math::concepts::std_real_concept signaling_NaN() throw();
   static boost::math::concepts::std_real_concept denorm_min() throw();
   static const bool is_iec559 = true;
   static const bool is_bounded = false;
   static const bool is_modulo = false;
   static const bool traps = false;
   static const bool tinyness_before = false;
   static const float_round_style round_style = round_toward_zero;
};
}
#endif
#ifdef EMULATE80
namespace std{
template<>
struct numeric_limits<boost::math::concepts::std_real_concept>
{
   static const bool is_specialized = true;
   static boost::math::concepts::std_real_concept min NULL_MACRO() throw();
   static boost::math::concepts::std_real_concept max NULL_MACRO() throw();
   static const int digits = 64;
   static const int digits10 = 18;
   static const bool is_signed = true;
   static const bool is_integer = false;
   static const bool is_exact = false;
   static const int radix = 2;
   static boost::math::concepts::std_real_concept epsilon() throw();
   static boost::math::concepts::std_real_concept round_error() throw();
   static const int min_exponent = -16381;
   static const int min_exponent10 = -4931;
   static const int max_exponent = 16384;
   static const int max_exponent10 = 4932;
   static const bool has_infinity = true;
   static const bool has_quiet_NaN = true;
   static const bool has_signaling_NaN = true;
   static const float_denorm_style has_denorm = denorm_absent;
   static const bool has_denorm_loss = false;
   static boost::math::concepts::std_real_concept infinity() throw();
   static boost::math::concepts::std_real_concept quiet_NaN() throw();
   static boost::math::concepts::std_real_concept signaling_NaN() throw();
   static boost::math::concepts::std_real_concept denorm_min() throw();
   static const bool is_iec559 = true;
   static const bool is_bounded = false;
   static const bool is_modulo = false;
   static const bool traps = false;
   static const bool tinyness_before = false;
   static const float_round_style round_style = round_toward_zero;
};
}
#endif
#ifdef EMULATE128
namespace std{
template<>
struct numeric_limits<boost::math::concepts::std_real_concept>
{
   static const bool is_specialized = true;
   static boost::math::concepts::std_real_concept min NULL_MACRO() throw();
   static boost::math::concepts::std_real_concept max NULL_MACRO() throw();
   static const int digits = 113;
   static const int digits10 = 33;
   static const bool is_signed = true;
   static const bool is_integer = false;
   static const bool is_exact = false;
   static const int radix = 2;
   static boost::math::concepts::std_real_concept epsilon() throw();
   static boost::math::concepts::std_real_concept round_error() throw();
   static const int min_exponent = -16381;
   static const int min_exponent10 = -4931;
   static const int max_exponent = 16384;
   static const int max_exponent10 = 4932;
   static const bool has_infinity = true;
   static const bool has_quiet_NaN = true;
   static const bool has_signaling_NaN = true;
   static const float_denorm_style has_denorm = denorm_absent;
   static const bool has_denorm_loss = false;
   static boost::math::concepts::std_real_concept infinity() throw();
   static boost::math::concepts::std_real_concept quiet_NaN() throw();
   static boost::math::concepts::std_real_concept signaling_NaN() throw();
   static boost::math::concepts::std_real_concept denorm_min() throw();
   static const bool is_iec559 = true;
   static const bool is_bounded = false;
   static const bool is_modulo = false;
   static const bool traps = false;
   static const bool tinyness_before = false;
   static const float_round_style round_style = round_toward_zero;
};
}
#endif

template <class RealType>
void instantiate(RealType)
{
   using namespace boost;
   using namespace boost::math;
   using namespace boost::math::concepts;

   function_requires<DistributionConcept<normal_distribution<RealType> > >();
   function_requires<DistributionConcept<binomial_distribution<RealType> > >();
   function_requires<DistributionConcept<cauchy_distribution<RealType> > >();
   function_requires<DistributionConcept<chi_squared_distribution<RealType> > >();
   function_requires<DistributionConcept<exponential_distribution<RealType> > >();
   function_requires<DistributionConcept<extreme_value_distribution<RealType> > >();
   function_requires<DistributionConcept<fisher_f_distribution<RealType> > >();
   function_requires<DistributionConcept<students_t_distribution<RealType> > >();
   function_requires<DistributionConcept<weibull_distribution<RealType> > >();
   function_requires<DistributionConcept<lognormal_distribution<RealType> > >();
}

void test_functions()
{
   int i;
   boost::math::concepts::std_real_concept v1(0.5), v2(0.5), v3(0.5);
   boost::math::tgamma(v1);
   boost::math::tgamma1pm1(v1);
   boost::math::lgamma(v1);
   boost::math::lgamma(v1, &i);
   boost::math::digamma(v1);
   boost::math::tgamma_ratio(v1, v2);
   boost::math::tgamma_delta_ratio(v1, v2);
   boost::math::factorial<boost::math::concepts::std_real_concept>(i);
   boost::math::unchecked_factorial<boost::math::concepts::std_real_concept>(i);
   i = boost::math::max_factorial<boost::math::concepts::std_real_concept>::value;
   boost::math::double_factorial<boost::math::concepts::std_real_concept>(i);
   boost::math::rising_factorial(v1, i);
   boost::math::falling_factorial(v1, i);
   boost::math::tgamma(v1, v2);
   boost::math::tgamma_lower(v1, v2);
   boost::math::gamma_P(v1, v2);
   boost::math::gamma_Q(v1, v2);
   boost::math::gamma_P_inv(v1, v2);
   boost::math::gamma_Q_inv(v1, v2);
   boost::math::gamma_P_inva(v1, v2);
   boost::math::gamma_Q_inva(v1, v2);
   boost::math::erf(v1);
   boost::math::erfc(v1);
   boost::math::erf_inv(v1);
   boost::math::erfc_inv(v1);
   boost::math::beta(v1, v2);
   boost::math::beta(v1, v2, v3);
   boost::math::betac(v1, v2, v3);
   boost::math::ibeta(v1, v2, v3);
   boost::math::ibetac(v1, v2, v3);
   boost::math::ibeta_inv(v1, v2, v3);
   boost::math::ibetac_inv(v1, v2, v3);
   boost::math::ibeta_inva(v1, v2, v3);
   boost::math::ibetac_inva(v1, v2, v3);
   boost::math::ibeta_invb(v1, v2, v3);
   boost::math::ibetac_invb(v1, v2, v3);
   boost::math::gamma_P_derivative(v2, v3);
   boost::math::ibeta_derivative(v1, v2, v3);
   boost::math::fpclassify(v1);
   boost::math::isfinite(v1);
   boost::math::isnormal(v1);
   boost::math::isnan(v1);
   boost::math::isinf(v1);
   boost::math::log1p(v1);
   boost::math::expm1(v1);
   boost::math::cbrt(v1);
   boost::math::sqrt1pm1(v1);
   boost::math::powm1(v1, v2);

}


int main()
{
   instantiate(boost::math::concepts::std_real_concept(0));
}

