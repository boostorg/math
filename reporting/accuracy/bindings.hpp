//  Copyright John Maddock 2015.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_BINDINGS
#define BOOST_MATH_BINDINGS

#define ERROR_REPORTING_MODE

#if TEST_LIBSTDCXX

#include <tr1/cmath>

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
#define BESSEL_Y_FUNCTION_TO_TEST std::tr1::cyl_neumann
#define BESSEL_YN_FUNCTION_TO_TEST std::tr1::cyl_neumann
#define BESSEL_YS_FUNCTION_TO_TEST std::tr1::sph_neumann

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

#endif

#if defined(TYPE_TO_TEST) && !defined(NAME_OF_TYPE_TO_TEST)
#define NAME_OF_TYPE_TO_TEST BOOST_STRINGIZE(TYPE_TO_TEST)
#endif

//
// This include has to come at the end after all the setup is done:
//
#include "handle_test_result.hpp"


#endif

