//  Copyright John Maddock 2006.
//  Copyright Paul A. Bristow 2007.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_LIBS_MATH_TEST_INSTANTIATE_HPP
#define BOOST_LIBS_MATH_TEST_INSTANTIATE_HPP

#ifndef BOOST_MATH_ASSERT_UNDEFINED_POLICY
#  define BOOST_MATH_ASSERT_UNDEFINED_POLICY false
#endif

#include <boost/math/distributions.hpp>

#include <boost/math/special_functions.hpp>
#include <boost/math/concepts/distributions.hpp>

typedef boost::math::policies::policy<> test_policy;

namespace test{

BOOST_MATH_DECLARE_SPECIAL_FUNCTIONS(test_policy)

}

namespace dist_test{

BOOST_MATH_DECLARE_DISTRIBUTIONS(double, test_policy)

}

template <class RealType>
void instantiate(RealType)
{
   using namespace boost;
   using namespace boost::math;
   using namespace boost::math::concepts;

   function_requires<DistributionConcept<bernoulli_distribution<RealType> > >();
   function_requires<DistributionConcept<beta_distribution<RealType> > >();
   function_requires<DistributionConcept<binomial_distribution<RealType> > >();
   function_requires<DistributionConcept<cauchy_distribution<RealType> > >();
   function_requires<DistributionConcept<chi_squared_distribution<RealType> > >();
   function_requires<DistributionConcept<exponential_distribution<RealType> > >();
   function_requires<DistributionConcept<extreme_value_distribution<RealType> > >();
   function_requires<DistributionConcept<fisher_f_distribution<RealType> > >();
   function_requires<DistributionConcept<gamma_distribution<RealType> > >();
   function_requires<DistributionConcept<lognormal_distribution<RealType> > >();
   function_requires<DistributionConcept<negative_binomial_distribution<RealType> > >();
   function_requires<DistributionConcept<normal_distribution<RealType> > >();
   function_requires<DistributionConcept<rayleigh_distribution<RealType> > >();
   function_requires<DistributionConcept<pareto_distribution<RealType> > >();
   function_requires<DistributionConcept<poisson_distribution<RealType> > >();
   function_requires<DistributionConcept<students_t_distribution<RealType> > >();
   function_requires<DistributionConcept<triangular_distribution<RealType> > >();
   function_requires<DistributionConcept<uniform_distribution<RealType> > >();
   function_requires<DistributionConcept<weibull_distribution<RealType> > >();

   function_requires<DistributionConcept<bernoulli_distribution<RealType, test_policy> > >();
   function_requires<DistributionConcept<beta_distribution<RealType, test_policy> > >();
   function_requires<DistributionConcept<binomial_distribution<RealType, test_policy> > >();
   function_requires<DistributionConcept<cauchy_distribution<RealType, test_policy> > >();
   function_requires<DistributionConcept<chi_squared_distribution<RealType, test_policy> > >();
   function_requires<DistributionConcept<exponential_distribution<RealType, test_policy> > >();
   function_requires<DistributionConcept<extreme_value_distribution<RealType, test_policy> > >();
   function_requires<DistributionConcept<fisher_f_distribution<RealType, test_policy> > >();
   function_requires<DistributionConcept<gamma_distribution<RealType, test_policy> > >();
   function_requires<DistributionConcept<lognormal_distribution<RealType, test_policy> > >();
   function_requires<DistributionConcept<negative_binomial_distribution<RealType, test_policy> > >();
   function_requires<DistributionConcept<normal_distribution<RealType, test_policy> > >();
   function_requires<DistributionConcept<rayleigh_distribution<RealType, test_policy> > >();
   function_requires<DistributionConcept<pareto_distribution<RealType, test_policy> > >();
   function_requires<DistributionConcept<poisson_distribution<RealType, test_policy> > >();
   function_requires<DistributionConcept<students_t_distribution<RealType, test_policy> > >();
   function_requires<DistributionConcept<triangular_distribution<RealType, test_policy> > >();
   function_requires<DistributionConcept<uniform_distribution<RealType, test_policy> > >();
   function_requires<DistributionConcept<weibull_distribution<RealType, test_policy> > >();

   function_requires<DistributionConcept<dist_test::bernoulli > >();
   function_requires<DistributionConcept<dist_test::beta > >();
   function_requires<DistributionConcept<dist_test::binomial > >();
   function_requires<DistributionConcept<dist_test::cauchy > >();
   function_requires<DistributionConcept<dist_test::chi_squared > >();
   function_requires<DistributionConcept<dist_test::exponential > >();
   function_requires<DistributionConcept<dist_test::extreme_value > >();
   function_requires<DistributionConcept<dist_test::fisher_f > >();
   function_requires<DistributionConcept<dist_test::gamma > >();
   function_requires<DistributionConcept<dist_test::lognormal > >();
   function_requires<DistributionConcept<dist_test::negative_binomial > >();
   function_requires<DistributionConcept<dist_test::normal > >();
   function_requires<DistributionConcept<dist_test::rayleigh > >();
   function_requires<DistributionConcept<dist_test::pareto > >();
   function_requires<DistributionConcept<dist_test::poisson > >();
   function_requires<DistributionConcept<dist_test::students_t > >();
   function_requires<DistributionConcept<dist_test::triangular > >();
   function_requires<DistributionConcept<dist_test::uniform > >();
   function_requires<DistributionConcept<dist_test::weibull > >();

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
   boost::math::gamma_p(v1, v2);
   boost::math::gamma_q(v1, v2);
   boost::math::gamma_p_inv(v1, v2);
   boost::math::gamma_q_inv(v1, v2);
   boost::math::gamma_p_inva(v1, v2);
   boost::math::gamma_q_inva(v1, v2);
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
   boost::math::gamma_p_derivative(v2, v3);
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
   boost::math::ellint_rc(v1, v2);
   boost::math::ellint_rd(v1, v2, v3);
   boost::math::ellint_rf(v1, v2, v3);
   boost::math::ellint_rj(v1, v2, v3, v1);
   boost::math::hypot(v1, v2);
   boost::math::sinc_pi(v1);
   boost::math::sinhc_pi(v1);
   boost::math::asinh(v1);
   boost::math::acosh(v1);
   boost::math::atanh(v1);
   boost::math::sin_pi(v1);
   boost::math::cos_pi(v1);
   boost::math::cyl_neumann(v1, v2);
   boost::math::cyl_neumann(i, v2);
   boost::math::cyl_bessel_j(v1, v2);
   boost::math::cyl_bessel_j(i, v2);
   boost::math::cyl_bessel_i(v1, v2);
   boost::math::cyl_bessel_i(i, v2);
   boost::math::cyl_bessel_k(v1, v2);
   boost::math::cyl_bessel_k(i, v2);
   boost::math::sph_bessel(i, v2);
   boost::math::sph_bessel(i, 1);
   boost::math::sph_neumann(i, v2);
   boost::math::sph_neumann(i, i);
   //
   // All over again, with a policy this time:
   //
   test_policy pol;
   boost::math::tgamma(v1, pol);
   boost::math::tgamma1pm1(v1, pol);
   boost::math::lgamma(v1, pol);
   boost::math::lgamma(v1, &i, pol);
   boost::math::digamma(v1, pol);
   boost::math::tgamma_ratio(v1, v2, pol);
   boost::math::tgamma_delta_ratio(v1, v2, pol);
   boost::math::factorial<RealType>(i, pol);
   boost::math::unchecked_factorial<RealType>(i);
   i = boost::math::max_factorial<RealType>::value;
   boost::math::double_factorial<RealType>(i, pol);
   boost::math::rising_factorial(v1, i, pol);
   boost::math::falling_factorial(v1, i, pol);
   boost::math::tgamma(v1, v2, pol);
   boost::math::tgamma_lower(v1, v2, pol);
   boost::math::gamma_p(v1, v2, pol);
   boost::math::gamma_q(v1, v2, pol);
   boost::math::gamma_p_inv(v1, v2, pol);
   boost::math::gamma_q_inv(v1, v2, pol);
   boost::math::gamma_p_inva(v1, v2, pol);
   boost::math::gamma_q_inva(v1, v2, pol);
   boost::math::erf(v1, pol);
   boost::math::erfc(v1, pol);
   boost::math::erf_inv(v1, pol);
   boost::math::erfc_inv(v1, pol);
   boost::math::beta(v1, v2, pol);
   boost::math::beta(v1, v2, v3, pol);
   boost::math::betac(v1, v2, v3, pol);
   boost::math::ibeta(v1, v2, v3, pol);
   boost::math::ibetac(v1, v2, v3, pol);
   boost::math::ibeta_inv(v1, v2, v3, pol);
   boost::math::ibetac_inv(v1, v2, v3, pol);
   boost::math::ibeta_inva(v1, v2, v3, pol);
   boost::math::ibetac_inva(v1, v2, v3, pol);
   boost::math::ibeta_invb(v1, v2, v3, pol);
   boost::math::ibetac_invb(v1, v2, v3, pol);
   boost::math::gamma_p_derivative(v2, v3, pol);
   boost::math::ibeta_derivative(v1, v2, v3, pol);
   boost::math::fpclassify(v1);
   boost::math::isfinite(v1);
   boost::math::isnormal(v1);
   boost::math::isnan(v1);
   boost::math::isinf(v1);
   boost::math::log1p(v1, pol);
   boost::math::expm1(v1, pol);
   boost::math::cbrt(v1, pol);
   boost::math::sqrt1pm1(v1, pol);
   boost::math::powm1(v1, v2, pol);
   boost::math::legendre_p(1, v1, pol);
   boost::math::legendre_p(1, 0, v1, pol);
   boost::math::legendre_q(1, v1, pol);
   boost::math::legendre_next(2, v1, v2, v3);
   boost::math::legendre_next(2, 2, v1, v2, v3);
   boost::math::laguerre(1, v1, pol);
   boost::math::laguerre(2, 1, v1, pol);
   boost::math::laguerre_next(2, v1, v2, v3);
   boost::math::laguerre_next(2, 1, v1, v2, v3);
   boost::math::hermite(1, v1, pol);
   boost::math::hermite_next(2, v1, v2, v3);
   boost::math::spherical_harmonic_r(2, 1, v1, v2, pol);
   boost::math::spherical_harmonic_i(2, 1, v1, v2, pol);
   boost::math::ellint_1(v1, pol);
   boost::math::ellint_1(v1, v2, pol);
   boost::math::ellint_2(v1, pol);
   boost::math::ellint_2(v1, v2, pol);
   boost::math::ellint_3(v1, v2, pol);
   boost::math::ellint_3(v1, v2, v3, pol);
   boost::math::ellint_rc(v1, v2, pol);
   boost::math::ellint_rd(v1, v2, v3, pol);
   boost::math::ellint_rf(v1, v2, v3, pol);
   boost::math::ellint_rj(v1, v2, v3, v1, pol);
   boost::math::hypot(v1, v2, pol);
   boost::math::sinc_pi(v1, pol);
   boost::math::sinhc_pi(v1, pol);
   boost::math::asinh(v1, pol);
   boost::math::acosh(v1, pol);
   boost::math::atanh(v1, pol);
   boost::math::sin_pi(v1, pol);
   boost::math::cos_pi(v1, pol);
   boost::math::cyl_neumann(v1, v2, pol);
   boost::math::cyl_neumann(i, v2, pol);
   boost::math::cyl_bessel_j(v1, v2, pol);
   boost::math::cyl_bessel_j(i, v2, pol);
   boost::math::cyl_bessel_i(v1, v2, pol);
   boost::math::cyl_bessel_i(i, v2, pol);
   boost::math::cyl_bessel_k(v1, v2, pol);
   boost::math::cyl_bessel_k(i, v2, pol);
   boost::math::sph_bessel(i, v2, pol);
   boost::math::sph_bessel(i, 1, pol);
   boost::math::sph_neumann(i, v2, pol);
   boost::math::sph_neumann(i, i, pol);
   //
   // All over again with the versions in test::
   //
   test::tgamma(v1);
   test::tgamma1pm1(v1);
   test::lgamma(v1);
   test::lgamma(v1, &i);
   test::digamma(v1);
   test::tgamma_ratio(v1, v2);
   test::tgamma_delta_ratio(v1, v2);
   test::factorial<RealType>(i);
   test::unchecked_factorial<RealType>(i);
   i = test::max_factorial<RealType>::value;
   test::double_factorial<RealType>(i);
   test::rising_factorial(v1, i);
   test::falling_factorial(v1, i);
   test::tgamma(v1, v2);
   test::tgamma_lower(v1, v2);
   test::gamma_p(v1, v2);
   test::gamma_q(v1, v2);
   test::gamma_p_inv(v1, v2);
   test::gamma_q_inv(v1, v2);
   test::gamma_p_inva(v1, v2);
   test::gamma_q_inva(v1, v2);
   test::erf(v1);
   test::erfc(v1);
   test::erf_inv(v1);
   test::erfc_inv(v1);
   test::beta(v1, v2);
   test::beta(v1, v2, v3);
   test::betac(v1, v2, v3);
   test::ibeta(v1, v2, v3);
   test::ibetac(v1, v2, v3);
   test::ibeta_inv(v1, v2, v3);
   test::ibetac_inv(v1, v2, v3);
   test::ibeta_inva(v1, v2, v3);
   test::ibetac_inva(v1, v2, v3);
   test::ibeta_invb(v1, v2, v3);
   test::ibetac_invb(v1, v2, v3);
   test::gamma_p_derivative(v2, v3);
   test::ibeta_derivative(v1, v2, v3);
   test::fpclassify(v1);
   test::isfinite(v1);
   test::isnormal(v1);
   test::isnan(v1);
   test::isinf(v1);
   test::log1p(v1);
   test::expm1(v1);
   test::cbrt(v1);
   test::sqrt1pm1(v1);
   test::powm1(v1, v2);
   test::legendre_p(1, v1);
   test::legendre_p(1, 0, v1);
   test::legendre_q(1, v1);
   test::legendre_next(2, v1, v2, v3);
   test::legendre_next(2, 2, v1, v2, v3);
   test::laguerre(1, v1);
   test::laguerre(2, 1, v1);
   test::laguerre_next(2, v1, v2, v3);
   test::laguerre_next(2, 1, v1, v2, v3);
   test::hermite(1, v1);
   test::hermite_next(2, v1, v2, v3);
   test::spherical_harmonic_r(2, 1, v1, v2);
   test::spherical_harmonic_i(2, 1, v1, v2);
   test::ellint_1(v1);
   test::ellint_1(v1, v2);
   test::ellint_2(v1);
   test::ellint_2(v1, v2);
   test::ellint_3(v1, v2);
   test::ellint_3(v1, v2, v3);
   test::ellint_rc(v1, v2);
   test::ellint_rd(v1, v2, v3);
   test::ellint_rf(v1, v2, v3);
   test::ellint_rj(v1, v2, v3, v1);
   test::hypot(v1, v2);
   test::sinc_pi(v1);
   test::sinhc_pi(v1);
   test::asinh(v1);
   test::acosh(v1);
   test::atanh(v1);
   test::sin_pi(v1);
   test::cos_pi(v1);
   test::cyl_neumann(v1, v2);
   test::cyl_neumann(i, v2);
   test::cyl_bessel_j(v1, v2);
   test::cyl_bessel_j(i, v2);
   test::cyl_bessel_i(v1, v2);
   test::cyl_bessel_i(i, v2);
   test::cyl_bessel_k(v1, v2);
   test::cyl_bessel_k(i, v2);
   test::sph_bessel(i, v2);
   test::sph_bessel(i, 1);
   test::sph_neumann(i, v2);
   test::sph_neumann(i, i);
}


#endif // BOOST_LIBS_MATH_TEST_INSTANTIATE_HPP
