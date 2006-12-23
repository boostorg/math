
#ifndef BOOST_LIBS_MATH_TEST_INSTANTIATE_HPP
#define BOOST_LIBS_MATH_TEST_INSTANTIATE_HPP

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
#include <boost/math/special_functions/legendre.hpp>
#include <boost/math/special_functions/laguerre.hpp>
#include <boost/math/special_functions/hermite.hpp>
#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/ellint_2.hpp>
#include <boost/math/special_functions/ellint_3.hpp>
#include <boost/math/concepts/distributions.hpp>

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

   int i;
   RealType v1(0.5), v2(0.5), v3(0.5);
   boost::math::tgamma(v1);
   boost::math::tgamma1pm1(v1);
   boost::math::lgamma(v1);
   boost::math::lgamma(v1, &i);
   boost::math::digamma(v1);
   boost::math::tgamma_ratio(v1, v2);
   boost::math::tgamma_delta_ratio(v1, v2);
   boost::math::factorial<RealType>(i);
   boost::math::unchecked_factorial<RealType>(i);
   i = boost::math::max_factorial<RealType>::value;
   boost::math::double_factorial<RealType>(i);
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
   boost::math::legendre_p(1, v1);
   boost::math::legendre_p(1, 0, v1);
   boost::math::legendre_q(1, v1);
   boost::math::legendre_next(2, v1, v2, v3);
   boost::math::legendre_next(2, 2, v1, v2, v3);
   boost::math::laguerre(1, v1);
   boost::math::laguerre(2, 1, v1);
   boost::math::laguerre_next(2, v1, v2, v3);
   boost::math::laguerre_next(2, 1, v1, v2, v3);
   boost::math::hermite(1, v1);
   boost::math::hermite_next(2, v1, v2, v3);
   boost::math::spherical_harmonic_r(2, 1, v1, v2);
   boost::math::spherical_harmonic_i(2, 1, v1, v2);
   boost::math::ellint_1(v1);
   boost::math::ellint_1(v1, v2);
   boost::math::ellint_2(v1);
   boost::math::ellint_2(v1, v2);
   boost::math::ellint_3(v1, v2);
   boost::math::ellint_3(v1, v2, v3);
}


#endif // BOOST_LIBS_MATH_TEST_INSTANTIATE_HPP
