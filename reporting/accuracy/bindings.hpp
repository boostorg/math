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

#define LOG1P_FUNCTION_TO_TEST std::tr1::log1p
#define EXPM1_FUNCTION_TO_TEST std::tr1::log1p

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

#define LOG1P_FUNCTION_TO_TEST ::log1p
#define EXPM1_FUNCTION_TO_TEST ::expm1

#define CBRT_FUNCTION_TO_TEST ::cbrt
#define ERF_FUNCTION_TO_TEST ::erf
#define ERFC_FUNCTION_TO_TEST ::erfc

#define LGAMMA_FUNCTION_TO_TEST ::lgamma
#define TGAMMA_FUNCTION_TO_TEST ::tgamma
#define BESSEL_JN_FUNCTION_TO_TEST ::jn
#define BESSEL_YN_FUNCTION_TO_TEST ::yn

#define TYPE_TO_TEST double

#else

#define LOG1P_FUNCTION_TO_TEST ::log1pl
#define EXPM1_FUNCTION_TO_TEST ::expm1l

#define CBRT_FUNCTION_TO_TEST ::cbrtl
#define ERF_FUNCTION_TO_TEST ::erfl
#define ERFC_FUNCTION_TO_TEST ::erfcl

#define LGAMMA_FUNCTION_TO_TEST ::lgammal
#define TGAMMA_FUNCTION_TO_TEST ::tgammal
//#define BESSEL_JN_FUNCTION_TO_TEST ::jnl
//#define BESSEL_JN_FUNCTION_TO_TEST ::ynl

#define TYPE_TO_TEST long double
#endif

#elif defined(TEST_GSL)

#include <stdexcept>

#include <gsl/gsl_sf.h>
#include <gsl/gsl_errno.h>

#define TEST_LIBRARY_NAME "GSL"

void gsl_handler(const char * reason, const char * file, int line, int gsl_errno)
{
   if(gsl_errno == GSL_ERANGE) return; // handle zero or infinity in our test code.
   throw std::domain_error(reason);
}

struct gsl_error_handler_setter
{
   gsl_error_handler_t * old_handler;
   gsl_error_handler_setter()
   {
      old_handler = gsl_set_error_handler(gsl_handler);
   }
   ~gsl_error_handler_setter()
   {
      gsl_set_error_handler(old_handler);
   }
};

static const gsl_error_handler_setter handler;

inline double gsl_bessel_ys(unsigned i, double d)
{
   return gsl_sf_bessel_yl(i, d);
}

inline double gsl_bessel_js(unsigned i, double d)
{
   return gsl_sf_bessel_jl(i, d);
}

//#define CBRT_FUNCTION_TO_TEST boost::cbrt
#define ERF_FUNCTION_TO_TEST gsl_sf_erf
#define ERFC_FUNCTION_TO_TEST gsl_sf_erfc
//#define ERF_INV_FUNCTION_TO_TEST boost::math::erf_inv
//#define ERFC_INV_FUNCTION_TO_TEST boost::math::erfc_inv

#define LGAMMA_FUNCTION_TO_TEST gsl_sf_lngamma
#define TGAMMA_FUNCTION_TO_TEST gsl_sf_gamma
//#define TGAMMA1PM1_FUNCTION_TO_TEST boost::math::tgamma1pm1

#define BESSEL_I_FUNCTION_TO_TEST gsl_sf_bessel_Inu
#define BESSEL_IN_FUNCTION_TO_TEST gsl_sf_bessel_In
//#define BESSEL_IP_FUNCTION_TO_TEST boost::math::cyl_bessel_i_prime
//#define BESSEL_IPN_FUNCTION_TO_TEST boost::math::cyl_bessel_i_prime
#define BESSEL_J_FUNCTION_TO_TEST  gsl_sf_bessel_Jnu
#define BESSEL_JN_FUNCTION_TO_TEST  gsl_sf_bessel_Jn
#define BESSEL_JS_FUNCTION_TO_TEST  gsl_bessel_js
//#define BESSEL_JP_FUNCTION_TO_TEST boost::math::cyl_bessel_j_prime
//#define BESSEL_JPN_FUNCTION_TO_TEST boost::math::cyl_bessel_j_prime
//#define BESSEL_JPS_FUNCTION_TO_TEST boost::math::sph_bessel_prime
#define BESSEL_K_FUNCTION_TO_TEST  gsl_sf_bessel_Knu
#define BESSEL_KN_FUNCTION_TO_TEST gsl_sf_bessel_Kn
//#define BESSEL_KP_FUNCTION_TO_TEST boost::math::cyl_bessel_k_prime
//#define BESSEL_KPN_FUNCTION_TO_TEST boost::math::cyl_bessel_k_prime
#define BESSEL_Y_FUNCTION_TO_TEST gsl_sf_bessel_Ynu
#define BESSEL_YN_FUNCTION_TO_TEST gsl_sf_bessel_Yn
#define BESSEL_YS_FUNCTION_TO_TEST gsl_bessel_ys
//#define BESSEL_YP_FUNCTION_TO_TEST boost::math::cyl_neumann_prime
//#define BESSEL_YNP_FUNCTION_TO_TEST boost::math::cyl_neumann_prime
//#define BESSEL_YSP_FUNCTION_TO_TEST boost::math::sph_neumann_prime

#define BETA_FUNCTION_TO_TEST gsl_sf_beta
//#define BINOMIAL_FUNCTION_TO_TEST boost::math::binomial_coefficient<T>

inline double RC(double a, double b)
{
   return gsl_sf_ellint_RC(a, b, GSL_PREC_DOUBLE);
}
inline double RD(double a, double b, double c)
{
   return gsl_sf_ellint_RD(a, b, c, GSL_PREC_DOUBLE);
}
inline double RF(double a, double b, double c)
{
   return gsl_sf_ellint_RF(a, b, c, GSL_PREC_DOUBLE);
}
inline double RJ(double a, double b, double c, double d)
{
   return gsl_sf_ellint_RJ(a, b, c, d, GSL_PREC_DOUBLE);
}


#define ELLINT_RC_FUNCTION_TO_TEST RC
#define ELLINT_RD_FUNCTION_TO_TEST RD
#define ELLINT_RF_FUNCTION_TO_TEST RF
//#define ELLINT_RG_FUNCTION_TO_TEST boost::math::ellint_rg
#define ELLINT_RJ_FUNCTION_TO_TEST RJ

#define DIGAMMA_FUNCTION_TO_TEST  gsl_sf_psi

inline double ellintK(double a) { return gsl_sf_ellint_Kcomp(a, GSL_PREC_DOUBLE); }
inline double ellintE(double a) { return gsl_sf_ellint_Ecomp(a, GSL_PREC_DOUBLE); }
inline double ellintP(double a, double b) { return gsl_sf_ellint_Pcomp(a, -b, GSL_PREC_DOUBLE); }

inline double ellintF(double a, double b) { return gsl_sf_ellint_F(b, a, GSL_PREC_DOUBLE); }
inline double ellintE2(double a, double b) { return gsl_sf_ellint_E(b, a, GSL_PREC_DOUBLE); }
inline double ellintP3(double a, double b, double c) { return gsl_sf_ellint_P(c, a, -b, GSL_PREC_DOUBLE); }
inline double ellintD2(double a, double b) { return gsl_sf_ellint_D(b, a, 0.0, GSL_PREC_DOUBLE); }

#define ELLINT_1_FUNCTION_TO_TEST ellintF
#define ELLINT_1C_FUNCTION_TO_TEST ellintK
#define ELLINT_2_FUNCTION_TO_TEST ellintE2
#define ELLINT_2C_FUNCTION_TO_TEST ellintE
#define ELLINT_3_FUNCTION_TO_TEST ellintP3
#define ELLINT_3C_FUNCTION_TO_TEST ellintP
#define ELLINT_D2_FUNCTION_TO_TEST ellintD2
//#define ELLINT_D1_FUNCTION_TO_TEST boost::math::ellint_d

#define EI_FUNCTION_TO_TEST gsl_sf_expint_Ei
#define EN_FUNCTION_TO_TEST  gsl_sf_expint_En

//#define HERMITE_FUNCTION_TO_TEST boost::math::hermite
//#define HEUMAN_LAMBDA_FUNCTION_TO_TEST boost::math::heuman_lambda

//#define BETA_INC_FUNCTION_TO_TEST boost::math::beta
//#define BETAC_INC_FUNCTION_TO_TEST boost::math::betac
#define IBETA_FUNCTION_TO_TEST gsl_sf_beta_inc
//#define IBETAC_FUNCTION_TO_TEST boost::math::ibetac
//#define IBETA_INV_FUNCTION_TO_TEST boost::math::ibeta_inv
//#define IBETAC_INV_FUNCTION_TO_TEST boost::math::ibetac_inv
//#define IBETA_INVA_FUNCTION_TO_TEST boost::math::ibeta_inva
//#define IBETAC_INVA_FUNCTION_TO_TEST boost::math::ibetac_inva
//#define IBETA_INVB_FUNCTION_TO_TEST boost::math::ibeta_invb
//#define IBETAC_INVB_FUNCTION_TO_TEST boost::math::ibetac_invb

#define IGAMMA_FUNCTION_TO_TEST gsl_sf_gamma_inc
//#define IGAMMAL_FUNCTION_TO_TEST boost::math::tgamma_lower
#define GAMMAP_FUNCTION_TO_TEST gsl_sf_gamma_inc_P
#define GAMMAQ_FUNCTION_TO_TEST gsl_sf_gamma_inc_Q
//#define GAMMAP_INV_FUNCTION_TO_TEST boost::math::gamma_p_inv
//#define GAMMAQ_INV_FUNCTION_TO_TEST boost::math::gamma_q_inv
//#define GAMMAP_INVA_FUNCTION_TO_TEST boost::math::gamma_p_inva
//#define GAMMAQ_INVA_FUNCTION_TO_TEST boost::math::gamma_q_inva

inline double sn(double k, double u)
{
   double s, c, d;
   gsl_sf_elljac_e(u, k * k, &s, &c, &d);
   return s;
}
inline double cn(double k, double u)
{
   double s, c, d;
   gsl_sf_elljac_e(u, k * k, &s, &c, &d);
   return c;
}
inline double dn(double k, double u)
{
   double s, c, d;
   gsl_sf_elljac_e(u, k * k, &s, &c, &d);
   return d;
}

#define SN_FUNCTION_TO_TEST sn
#define CN_FUNCTION_TO_TEST cn
#define DN_FUNCTION_TO_TEST dn
//#define JACOBI_ZETA_FUNCTION_TO_TEST boost::math::jacobi_zeta

inline double laguerre(unsigned n, unsigned m, double x){ return gsl_sf_laguerre_n(n, m, x); }
inline double laguerre_0(unsigned n, double x){ return gsl_sf_laguerre_n(n, 0, x); }

#define LAGUERRE_FUNCTION_TO_TEST laguerre_0
#define ASSOC_LAGUERRE_FUNCTION_TO_TEST laguerre

inline double legendre_q(unsigned n, double x) { return gsl_sf_legendre_Ql(n, x); }

#define LEGENDRE_P_FUNCTION_TO_TEST gsl_sf_legendre_Pl
#define LEGENDRE_Q_FUNCTION_TO_TEST  legendre_q
#define LEGENDRE_PA_FUNCTION_TO_TEST  gsl_sf_legendre_Plm

#define POLYGAMMA_FUNCTION_TO_TEST  gsl_sf_psi_n
//#define TGAMMA_RATIO_FUNCTION_TO_TEST boost::math::tgamma_ratio
//#define TGAMMA_DELTA_RATIO_FUNCTION_TO_TEST boost::math::tgamma_delta_ratio
//#define SIN_PI_RATIO_FUNCTION_TO_TEST boost::math::sin_pi
//#define COS_PI_RATIO_FUNCTION_TO_TEST boost::math::cos_pi
#define TRIGAMMA_RATIO_FUNCTION_TO_TEST gsl_sf_psi_1
#define ZETA_FUNCTION_TO_TEST gsl_sf_zeta

#define TYPE_TO_TEST double

#else

#include <boost/math/distributions/non_central_beta.hpp>
#include <boost/math/distributions/non_central_chi_squared.hpp>
#include <boost/math/distributions/non_central_t.hpp>

#define TEST_LIBRARY_NAME "boost"

#define LOG1P_FUNCTION_TO_TEST boost::math::log1p
#define EXPM1_FUNCTION_TO_TEST boost::math::expm1

#define CBRT_FUNCTION_TO_TEST boost::math::cbrt
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

#define SQRT1PM1_FUNCTION_TO_TEST boost::math::sqrt1pm1
#define POWM1_FUNCTION_TO_TEST boost::math::powm1
#define OWENS_T_FUNCTION_TO_TEST boost::math::owens_t
#define SPHERICAL_HARMONIC_R_FUNCTION_TO_TEST boost::math::spherical_harmonic_r
#define SPHERICAL_HARMONIC_I_FUNCTION_TO_TEST boost::math::spherical_harmonic_i

template <class T> T do_nc_beta_cdf(T a, T b, T nc, T x){ return cdf(boost::math::non_central_beta_distribution<T>(a, b, nc), x); }
template <class T> T do_nc_beta_ccdf(T a, T b, T nc, T x){ return cdf(complement(boost::math::non_central_beta_distribution<T>(a, b, nc), x)); }
template <class T> T do_nc_chi_squared_cdf(T df, T nc, T x){ return cdf(boost::math::non_central_chi_squared_distribution<T>(df, nc), x); }
template <class T> T do_nc_chi_squared_ccdf(T df, T nc, T x){ return cdf(complement(boost::math::non_central_chi_squared_distribution<T>(df, nc), x)); }
template <class T> T do_nc_t_cdf(T df, T nc, T x){ return cdf(boost::math::non_central_t_distribution<T>(df, nc), x); }
template <class T> T do_nc_t_ccdf(T df, T nc, T x){ return cdf(complement(boost::math::non_central_t_distribution<T>(df, nc), x)); }

#define NC_BETA_CDF_FUNCTION_TO_TEST do_nc_beta_cdf
#define NC_BETA_CCDF_FUNCTION_TO_TEST do_nc_beta_ccdf
#define NC_CHI_SQUARED_CDF_FUNCTION_TO_TEST do_nc_chi_squared_cdf
#define NC_CHI_SQUARED_CCDF_FUNCTION_TO_TEST do_nc_chi_squared_ccdf
#define NC_T_CDF_FUNCTION_TO_TEST do_nc_t_cdf
#define NC_T_CCDF_FUNCTION_TO_TEST do_nc_t_ccdf


#endif

#if defined(TYPE_TO_TEST) && !defined(NAME_OF_TYPE_TO_TEST)
#define NAME_OF_TYPE_TO_TEST BOOST_STRINGIZE(TYPE_TO_TEST)
#endif

//
// This include has to come at the end after all the setup is done:
//
#include "handle_test_result.hpp"


#endif

