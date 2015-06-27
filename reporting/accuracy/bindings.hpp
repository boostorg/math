//  Copyright John Maddock 2015.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_BINDINGS
#define BOOST_MATH_BINDINGS

#define ERROR_REPORTING_MODE

#if TEST_LIBSTDCXX

#include <tr1/cmath>
#include <stdexcept>

#define TEST_LIBRARY_NAME "<tr1/cmath>"

#define CBRT_FUNCTION_TO_TEST std::tr1::cbrt
#define ERF_FUNCTION_TO_TEST std::tr1::erf
#define ERFC_FUNCTION_TO_TEST std::tr1::erfc

#define LGAMMA_FUNCTION_TO_TEST std::tr1::lgamma
#define TGAMMA_FUNCTION_TO_TEST std::tr1::tgamma

#define BESSEL_I_FUNCTION_TO_TEST std::tr1::cyl_bessel_i
#define BESSEL_IN_FUNCTION_TO_TEST std::tr1::cyl_bessel_i
#define BESSEL_J_FUNCTION_TO_TEST std::tr1::cyl_bessel_j
#define BESSEL_JN_FUNCTION_TO_TEST std::tr1::cyl_bessel_j
#define BESSEL_JS_FUNCTION_TO_TEST std::tr1::sph_bessel
#define BESSEL_K_FUNCTION_TO_TEST std::tr1::cyl_bessel_k
#define BESSEL_KN_FUNCTION_TO_TEST std::tr1::cyl_bessel_k
#define BESSEL_Y_FUNCTION_TO_TEST std::tr1::cyl_neumannl
#define BESSEL_YN_FUNCTION_TO_TEST std::tr1::cyl_neumannl
#define BESSEL_YS_FUNCTION_TO_TEST std::tr1::sph_neumannl

#define BETA_FUNCTION_TO_TEST std::tr1::betal

#define ELLINT_1_FUNCTION_TO_TEST std::tr1::ellint_1l
#define ELLINT_1C_FUNCTION_TO_TEST std::tr1::comp_ellint_1l
#define ELLINT_2_FUNCTION_TO_TEST std::tr1::ellint_2l
#define ELLINT_2C_FUNCTION_TO_TEST std::tr1::comp_ellint_2l
#define ELLINT_3_FUNCTION_TO_TEST std::tr1::ellint_3l
#define ELLINT_3C_FUNCTION_TO_TEST std::tr1::comp_ellint_3l

#define EI_FUNCTION_TO_TEST std::tr1::expintl

#define LAGUERRE_FUNCTION_TO_TEST std::tr1::laguerrel
#define ASSOC_LAGUERRE_FUNCTION_TO_TEST std::tr1::assoc_laguerrel

inline long double legendre_p_binder(int i, long double d)
{
   if(i < 0)
      throw std::domain_error("order parameters less than 0 not supported in TR1");
   return std::tr1::legendre(i, d);
}
inline long double assoc_legendre_p_binder(int i, int j, long double d)
{
   if((i < 0) || (j < 0))
      throw std::domain_error("order parameters less than 0 not supported in TR1");
   return std::tr1::assoc_legendre(i, j, d);
}

#define LEGENDRE_P_FUNCTION_TO_TEST legendre_p_binder
#define LEGENDRE_PA_FUNCTION_TO_TEST assoc_legendre_p_binder
#define ZETA_FUNCTION_TO_TEST std::tr1::riemann_zeta

#define TYPE_TO_TEST long double

#elif defined(TEST_C99)

#include <math.h>

#define TEST_LIBRARY_NAME "<math.h>"

#ifdef _MSC_VER

#define CBRT_FUNCTION_TO_TEST ::cbrt
#define ERF_FUNCTION_TO_TEST ::erf
#define ERFC_FUNCTION_TO_TEST ::erfc

#define LGAMMA_FUNCTION_TO_TEST ::lgamma
#define TGAMMA_FUNCTION_TO_TEST ::tgamma
#define BESSEL_JN_FUNCTION_TO_TEST ::jn
#define BESSEL_YN_FUNCTION_TO_TEST ::yn

#define TYPE_TO_TEST double

#else

#define CBRT_FUNCTION_TO_TEST ::cbrt
#define ERF_FUNCTION_TO_TEST ::erfl
#define ERFC_FUNCTION_TO_TEST ::erfcl

#define LGAMMA_FUNCTION_TO_TEST ::lgammal
#define TGAMMA_FUNCTION_TO_TEST ::tgammal
//#define BESSEL_JN_FUNCTION_TO_TEST ::jnl
//#define BESSEL_JN_FUNCTION_TO_TEST ::ynl

#define TYPE_TO_TEST long double
#endif

#elif defined(TEST_GSL)

//#define CBRT_FUNCTION_TO_TEST boost::cbrt
#define ERF_FUNCTION_TO_TEST gsl_sf_erf
#define ERFC_FUNCTION_TO_TEST gsl_sf_erfc
//#define ERF_INV_FUNCTION_TO_TEST boost::math::erf_inv
//#define ERFC_INV_FUNCTION_TO_TEST boost::math::erfc_inv

#define LGAMMA_FUNCTION_TO_TEST gsl_sf_lngamma
#define TGAMMA_FUNCTION_TO_TEST gsl_sf_gamma
//#define TGAMMA1PM1_FUNCTION_TO_TEST boost::math::tgamma1pm1

#define TYPE_TO_TEST double

#else

#define TEST_LIBRARY_NAME "boost"

#define CBRT_FUNCTION_TO_TEST boost::cbrt
#define ERF_FUNCTION_TO_TEST boost::math::erf
#define ERFC_FUNCTION_TO_TEST boost::math::erfc
#define ERF_INV_FUNCTION_TO_TEST boost::math::erf_inv
#define ERFC_INV_FUNCTION_TO_TEST boost::math::erfc_inv

#define LGAMMA_FUNCTION_TO_TEST boost::math::lgamma
#define TGAMMA_FUNCTION_TO_TEST boost::math::tgamma
#define TGAMMA1PM1_FUNCTION_TO_TEST boost::math::tgamma1pm1

#define BESSEL_I_FUNCTION_TO_TEST boost::math::cyl_bessel_i
#define BESSEL_IN_FUNCTION_TO_TEST boost::math::cyl_bessel_i
#define BESSEL_IP_FUNCTION_TO_TEST boost::math::cyl_bessel_i_prime
#define BESSEL_IPN_FUNCTION_TO_TEST boost::math::cyl_bessel_i_prime
#define BESSEL_J_FUNCTION_TO_TEST boost::math::cyl_bessel_j
#define BESSEL_JN_FUNCTION_TO_TEST boost::math::cyl_bessel_j
#define BESSEL_JS_FUNCTION_TO_TEST boost::math::sph_bessel
#define BESSEL_JP_FUNCTION_TO_TEST boost::math::cyl_bessel_j_prime
#define BESSEL_JPN_FUNCTION_TO_TEST boost::math::cyl_bessel_j_prime
#define BESSEL_JPS_FUNCTION_TO_TEST boost::math::sph_bessel_prime
#define BESSEL_K_FUNCTION_TO_TEST boost::math::cyl_bessel_k
#define BESSEL_KN_FUNCTION_TO_TEST boost::math::cyl_bessel_k
#define BESSEL_KP_FUNCTION_TO_TEST boost::math::cyl_bessel_k_prime
#define BESSEL_KPN_FUNCTION_TO_TEST boost::math::cyl_bessel_k_prime
#define BESSEL_Y_FUNCTION_TO_TEST boost::math::cyl_neumann
#define BESSEL_YN_FUNCTION_TO_TEST boost::math::cyl_neumann
#define BESSEL_YS_FUNCTION_TO_TEST boost::math::sph_neumann
#define BESSEL_YP_FUNCTION_TO_TEST boost::math::cyl_neumann_prime
#define BESSEL_YNP_FUNCTION_TO_TEST boost::math::cyl_neumann_prime
#define BESSEL_YSP_FUNCTION_TO_TEST boost::math::sph_neumann_prime

#define BETA_FUNCTION_TO_TEST boost::math::beta
#define BINOMIAL_FUNCTION_TO_TEST boost::math::binomial_coefficient<T>

#define ELLINT_RC_FUNCTION_TO_TEST boost::math::ellint_rc
#define ELLINT_RD_FUNCTION_TO_TEST boost::math::ellint_rd
#define ELLINT_RF_FUNCTION_TO_TEST boost::math::ellint_rf
#define ELLINT_RG_FUNCTION_TO_TEST boost::math::ellint_rg
#define ELLINT_RJ_FUNCTION_TO_TEST boost::math::ellint_rj

#define DIGAMMA_FUNCTION_TO_TEST boost::math::digamma

#define ELLINT_1_FUNCTION_TO_TEST boost::math::ellint_1
#define ELLINT_1C_FUNCTION_TO_TEST boost::math::ellint_1
#define ELLINT_2_FUNCTION_TO_TEST boost::math::ellint_2
#define ELLINT_2C_FUNCTION_TO_TEST boost::math::ellint_2
#define ELLINT_3_FUNCTION_TO_TEST boost::math::ellint_3
#define ELLINT_3C_FUNCTION_TO_TEST boost::math::ellint_3
#define ELLINT_D2_FUNCTION_TO_TEST boost::math::ellint_d
#define ELLINT_D1_FUNCTION_TO_TEST boost::math::ellint_d

#define EI_FUNCTION_TO_TEST boost::math::expint
#define EN_FUNCTION_TO_TEST boost::math::expint

#define HERMITE_FUNCTION_TO_TEST boost::math::hermite
#define HEUMAN_LAMBDA_FUNCTION_TO_TEST boost::math::heuman_lambda

#define BETA_INC_FUNCTION_TO_TEST boost::math::beta
#define BETAC_INC_FUNCTION_TO_TEST boost::math::betac
#define IBETA_FUNCTION_TO_TEST boost::math::ibeta
#define IBETAC_FUNCTION_TO_TEST boost::math::ibetac
#define IBETA_INV_FUNCTION_TO_TEST boost::math::ibeta_inv
#define IBETAC_INV_FUNCTION_TO_TEST boost::math::ibetac_inv
#define IBETA_INVA_FUNCTION_TO_TEST boost::math::ibeta_inva
#define IBETAC_INVA_FUNCTION_TO_TEST boost::math::ibetac_inva
#define IBETA_INVB_FUNCTION_TO_TEST boost::math::ibeta_invb
#define IBETAC_INVB_FUNCTION_TO_TEST boost::math::ibetac_invb

#define IGAMMA_FUNCTION_TO_TEST boost::math::tgamma
#define IGAMMAL_FUNCTION_TO_TEST boost::math::tgamma_lower
#define GAMMAP_FUNCTION_TO_TEST boost::math::gamma_p
#define GAMMAQ_FUNCTION_TO_TEST boost::math::gamma_q
#define GAMMAP_INV_FUNCTION_TO_TEST boost::math::gamma_p_inv
#define GAMMAQ_INV_FUNCTION_TO_TEST boost::math::gamma_q_inv
#define GAMMAP_INVA_FUNCTION_TO_TEST boost::math::gamma_p_inva
#define GAMMAQ_INVA_FUNCTION_TO_TEST boost::math::gamma_q_inva

#define SN_FUNCTION_TO_TEST boost::math::jacobi_sn
#define CN_FUNCTION_TO_TEST boost::math::jacobi_cn
#define DN_FUNCTION_TO_TEST boost::math::jacobi_dn
#define JACOBI_ZETA_FUNCTION_TO_TEST boost::math::jacobi_zeta

#define LAGUERRE_FUNCTION_TO_TEST boost::math::laguerre
#define ASSOC_LAGUERRE_FUNCTION_TO_TEST boost::math::laguerre

#define LEGENDRE_P_FUNCTION_TO_TEST boost::math::legendre_p
#define LEGENDRE_Q_FUNCTION_TO_TEST boost::math::legendre_q
#define LEGENDRE_PA_FUNCTION_TO_TEST boost::math::legendre_p

#define POLYGAMMA_FUNCTION_TO_TEST boost::math::polygamma
#define TGAMMA_RATIO_FUNCTION_TO_TEST boost::math::tgamma_ratio
#define TGAMMA_DELTA_RATIO_FUNCTION_TO_TEST boost::math::tgamma_delta_ratio
#define SIN_PI_RATIO_FUNCTION_TO_TEST boost::math::sin_pi
#define COS_PI_RATIO_FUNCTION_TO_TEST boost::math::cos_pi
#define TRIGAMMA_RATIO_FUNCTION_TO_TEST boost::math::trigamma
#define ZETA_FUNCTION_TO_TEST boost::math::zeta

#endif

#if defined(TYPE_TO_TEST) && !defined(NAME_OF_TYPE_TO_TEST)
#define NAME_OF_TYPE_TO_TEST BOOST_STRINGIZE(TYPE_TO_TEST)
#endif

//
// This include has to come at the end after all the setup is done:
//
#include "handle_test_result.hpp"


#endif

