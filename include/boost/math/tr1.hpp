// Copyright John Maddock 2008.
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_TR1_HPP
#define BOOST_MATH_TR1_HPP

#ifdef _MSC_VER
#pragma once
#endif

#ifdef __cplusplus

#include <boost/config.hpp>

namespace boost{ namespace math{ namespace tr1{ extern "C"{

#endif // __cplusplus

#ifdef BOOST_HAS_DECLSPEC // defined in config system
// we need to import/export our code only if the user has specifically
// asked for it by defining either BOOST_ALL_DYN_LINK if they want all boost
// libraries to be dynamically linked, or BOOST_MATH_TR1_DYN_LINK
// if they want just this one to be dynamically liked:
#if defined(BOOST_ALL_DYN_LINK) || defined(BOOST_MATH_TR1_DYN_LINK)
// export if this is our own source, otherwise import:
#ifdef BOOST_MATH_TR1_SOURCE
# define BOOST_MATH_TR1_DECL __declspec(dllexport)
#else
# define BOOST_MATH_TR1_DECL __declspec(dllimport)
#endif  // BOOST_MATH_TR1_SOURCE
#endif  // DYN_LINK
#endif  // BOOST_HAS_DECLSPEC
//
// if BOOST_MATH_TR1_DECL isn't defined yet define it now:
#ifndef BOOST_MATH_TR1_DECL
#define BOOST_MATH_TR1_DECL
#endif

//
// Now set up the libraries to link against:
//
#if !defined(BOOST_MATH_TR1_NO_LIB) && !defined(BOOST_MATH_TR1_SOURCE) \
   && !defined(BOOST_ALL_NO_LIB) && defined(__cplusplus)
#  define BOOST_LIB_NAME boost_math_c99
#  if defined(BOOST_MATH_TR1_DYN_LINK) || defined(BOOST_ALL_DYN_LINK)
#     define BOOST_DYN_LINK
#  endif
#  include <boost/config/auto_link.hpp>
#endif
#if !defined(BOOST_MATH_TR1_NO_LIB) && !defined(BOOST_MATH_TR1_SOURCE) \
   && !defined(BOOST_ALL_NO_LIB) && defined(__cplusplus)
#  define BOOST_LIB_NAME boost_math_c99f
#  if defined(BOOST_MATH_TR1_DYN_LINK) || defined(BOOST_ALL_DYN_LINK)
#     define BOOST_DYN_LINK
#  endif
#  include <boost/config/auto_link.hpp>
#endif
#if !defined(BOOST_MATH_TR1_NO_LIB) && !defined(BOOST_MATH_TR1_SOURCE) \
   && !defined(BOOST_ALL_NO_LIB) && defined(__cplusplus) \
   && !defined(BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS)
#  define BOOST_LIB_NAME boost_math_c99l
#  if defined(BOOST_MATH_TR1_DYN_LINK) || defined(BOOST_ALL_DYN_LINK)
#     define BOOST_DYN_LINK
#  endif
#  include <boost/config/auto_link.hpp>
#endif
#if !defined(BOOST_MATH_TR1_NO_LIB) && !defined(BOOST_MATH_TR1_SOURCE) \
   && !defined(BOOST_ALL_NO_LIB) && defined(__cplusplus)
#  define BOOST_LIB_NAME boost_math_tr1
#  if defined(BOOST_MATH_TR1_DYN_LINK) || defined(BOOST_ALL_DYN_LINK)
#     define BOOST_DYN_LINK
#  endif
#  include <boost/config/auto_link.hpp>
#endif
#if !defined(BOOST_MATH_TR1_NO_LIB) && !defined(BOOST_MATH_TR1_SOURCE) \
   && !defined(BOOST_ALL_NO_LIB) && defined(__cplusplus)
#  define BOOST_LIB_NAME boost_math_tr1f
#  if defined(BOOST_MATH_TR1_DYN_LINK) || defined(BOOST_ALL_DYN_LINK)
#     define BOOST_DYN_LINK
#  endif
#  include <boost/config/auto_link.hpp>
#endif
#if !defined(BOOST_MATH_TR1_NO_LIB) && !defined(BOOST_MATH_TR1_SOURCE) \
   && !defined(BOOST_ALL_NO_LIB) && defined(__cplusplus) \
   && !defined(BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS)
#  define BOOST_LIB_NAME boost_math_tr1l
#  if defined(BOOST_MATH_TR1_DYN_LINK) || defined(BOOST_ALL_DYN_LINK)
#     define BOOST_DYN_LINK
#  endif
#  include <boost/config/auto_link.hpp>
#endif

#ifndef FLT_EVAL_METHOD
typedef float float_t;
typedef double double_t;
#elif FLT_EVAL_METHOD == 0
typedef float float_t;
typedef double double_t;
#elif FLT_EVAL_METHOD == 1
typedef double float_t;
typedef double double_t;
#else
typedef long double float_t;
typedef long double double_t;
#endif

// C99 Functions:
double BOOST_MATH_TR1_DECL acosh(double x);
float BOOST_MATH_TR1_DECL acoshf(float x);
long double BOOST_MATH_TR1_DECL acoshl(long double x);

double BOOST_MATH_TR1_DECL asinh(double x);
float BOOST_MATH_TR1_DECL asinhf(float x);
long double BOOST_MATH_TR1_DECL asinhl(long double x);

double BOOST_MATH_TR1_DECL atanh(double x);
float BOOST_MATH_TR1_DECL atanhf(float x);
long double BOOST_MATH_TR1_DECL atanhl(long double x);

double BOOST_MATH_TR1_DECL cbrt(double x);
float BOOST_MATH_TR1_DECL cbrtf(float x);
long double BOOST_MATH_TR1_DECL cbrtl(long double x);

double BOOST_MATH_TR1_DECL copysign(double x, double y);
float BOOST_MATH_TR1_DECL copysignf(float x, float y);
long double BOOST_MATH_TR1_DECL copysignl(long double x, long double y);

double BOOST_MATH_TR1_DECL erf(double x);
float BOOST_MATH_TR1_DECL erff(float x);
long double BOOST_MATH_TR1_DECL erfl(long double x);

double BOOST_MATH_TR1_DECL erfc(double x);
float BOOST_MATH_TR1_DECL erfcf(float x);
long double BOOST_MATH_TR1_DECL erfcl(long double x);
#if 0
double BOOST_MATH_TR1_DECL exp2(double x);
float BOOST_MATH_TR1_DECL exp2f(float x);
long double BOOST_MATH_TR1_DECL exp2l(long double x);
#endif
double BOOST_MATH_TR1_DECL boost_expm1(double x);
float BOOST_MATH_TR1_DECL boost_expm1f(float x);
long double BOOST_MATH_TR1_DECL boost_expm1l(long double x);
#if 0
double BOOST_MATH_TR1_DECL fdim(double x, double y);
float BOOST_MATH_TR1_DECL fdimf(float x, float y);
long double BOOST_MATH_TR1_DECL fdiml(long double x, long double y);
double BOOST_MATH_TR1_DECL fma(double x, double y, double z);
float BOOST_MATH_TR1_DECL fmaf(float x, float y, float z);
long double BOOST_MATH_TR1_DECL fmal(long double x, long double y, long double z);
#endif
double BOOST_MATH_TR1_DECL fmax(double x, double y);
float BOOST_MATH_TR1_DECL fmaxf(float x, float y);
long double BOOST_MATH_TR1_DECL fmaxl(long double x, long double y);

double BOOST_MATH_TR1_DECL fmin(double x, double y);
float BOOST_MATH_TR1_DECL fminf(float x, float y);
long double BOOST_MATH_TR1_DECL fminl(long double x, long double y);

double BOOST_MATH_TR1_DECL hypot(double x, double y);
float BOOST_MATH_TR1_DECL hypotf(float x, float y);
long double BOOST_MATH_TR1_DECL hypotl(long double x, long double y);
#if 0
int BOOST_MATH_TR1_DECL ilogb(double x);
int BOOST_MATH_TR1_DECL ilogbf(float x);
int BOOST_MATH_TR1_DECL ilogbl(long double x);
#endif
double BOOST_MATH_TR1_DECL lgamma(double x);
float BOOST_MATH_TR1_DECL lgammaf(float x);
long double BOOST_MATH_TR1_DECL lgammal(long double x);
#if 0
long long BOOST_MATH_TR1_DECL llrint(double x);
long long BOOST_MATH_TR1_DECL llrintf(float x);
long long BOOST_MATH_TR1_DECL llrintl(long double x);
#endif
long long BOOST_MATH_TR1_DECL llround(double x);
long long BOOST_MATH_TR1_DECL llroundf(float x);
long long BOOST_MATH_TR1_DECL llroundl(long double x);

double BOOST_MATH_TR1_DECL boost_log1p(double x);
float BOOST_MATH_TR1_DECL boost_log1pf(float x);
long double BOOST_MATH_TR1_DECL boost_log1pl(long double x);
#if 0
double BOOST_MATH_TR1_DECL log2(double x);
float BOOST_MATH_TR1_DECL log2f(float x);
long double BOOST_MATH_TR1_DECL log2l(long double x);

double BOOST_MATH_TR1_DECL logb(double x);
float BOOST_MATH_TR1_DECL logbf(float x);
long double BOOST_MATH_TR1_DECL logbl(long double x);
long BOOST_MATH_TR1_DECL lrint(double x);
long BOOST_MATH_TR1_DECL lrintf(float x);
long BOOST_MATH_TR1_DECL lrintl(long double x);
#endif
long BOOST_MATH_TR1_DECL lround(double x);
long BOOST_MATH_TR1_DECL lroundf(float x);
long BOOST_MATH_TR1_DECL lroundl(long double x);
#if 0
double BOOST_MATH_TR1_DECL nan(const char *str);
float BOOST_MATH_TR1_DECL nanf(const char *str);
long double BOOST_MATH_TR1_DECL nanl(const char *str);
double BOOST_MATH_TR1_DECL nearbyint(double x);
float BOOST_MATH_TR1_DECL nearbyintf(float x);
long double BOOST_MATH_TR1_DECL nearbyintl(long double x);
#endif
double BOOST_MATH_TR1_DECL boost_nextafter(double x, double y);
float BOOST_MATH_TR1_DECL boost_nextafterf(float x, float y);
long double BOOST_MATH_TR1_DECL boost_nextafterl(long double x, long double y);

double BOOST_MATH_TR1_DECL nexttoward(double x, long double y);
float BOOST_MATH_TR1_DECL nexttowardf(float x, long double y);
long double BOOST_MATH_TR1_DECL nexttowardl(long double x, long double y);
#if 0
double BOOST_MATH_TR1_DECL remainder(double x, double y);
float BOOST_MATH_TR1_DECL remainderf(float x, float y);
long double BOOST_MATH_TR1_DECL remainderl(long double x, long double y);
double BOOST_MATH_TR1_DECL remquo(double x, double y, int *pquo);
float BOOST_MATH_TR1_DECL remquof(float x, float y, int *pquo);
long double BOOST_MATH_TR1_DECL remquol(long double x, long double y, int *pquo);
double BOOST_MATH_TR1_DECL rint(double x);
float BOOST_MATH_TR1_DECL rintf(float x);
long double BOOST_MATH_TR1_DECL rintl(long double x);
#endif
double BOOST_MATH_TR1_DECL round(double x);
float BOOST_MATH_TR1_DECL roundf(float x);
long double BOOST_MATH_TR1_DECL roundl(long double x);
#if 0
double BOOST_MATH_TR1_DECL scalbln(double x, long ex);
float BOOST_MATH_TR1_DECL scalblnf(float x, long ex);
long double BOOST_MATH_TR1_DECL scalblnl(long double x, long ex);
double BOOST_MATH_TR1_DECL scalbn(double x, int ex);
float BOOST_MATH_TR1_DECL scalbnf(float x, int ex);
long double BOOST_MATH_TR1_DECL scalbnl(long double x, int ex);
#endif
double BOOST_MATH_TR1_DECL tgamma(double x);
float BOOST_MATH_TR1_DECL tgammaf(float x);
long double BOOST_MATH_TR1_DECL tgammal(long double x);

double BOOST_MATH_TR1_DECL trunc(double x);
float BOOST_MATH_TR1_DECL truncf(float x);
long double BOOST_MATH_TR1_DECL truncl(long double x);

// [5.2.1.1] associated Laguerre polynomials:
double BOOST_MATH_TR1_DECL assoc_laguerre(unsigned n, unsigned m, double x);
float BOOST_MATH_TR1_DECL assoc_laguerref(unsigned n, unsigned m, float x);
long double BOOST_MATH_TR1_DECL assoc_laguerrel(unsigned n, unsigned m, long double x);

// [5.2.1.2] associated Legendre functions:
double BOOST_MATH_TR1_DECL assoc_legendre(unsigned l, unsigned m, double x);
float BOOST_MATH_TR1_DECL assoc_legendref(unsigned l, unsigned m, float x);
long double BOOST_MATH_TR1_DECL assoc_legendrel(unsigned l, unsigned m, long double x);

// [5.2.1.3] beta function:
double BOOST_MATH_TR1_DECL beta(double x, double y);
float BOOST_MATH_TR1_DECL betaf(float x, float y);
long double BOOST_MATH_TR1_DECL betal(long double x, long double y);

// [5.2.1.4] (complete) elliptic integral of the first kind:
double BOOST_MATH_TR1_DECL comp_ellint_1(double k);
float BOOST_MATH_TR1_DECL comp_ellint_1f(float k);
long double BOOST_MATH_TR1_DECL comp_ellint_1l(long double k);

// [5.2.1.5] (complete) elliptic integral of the second kind:
double BOOST_MATH_TR1_DECL comp_ellint_2(double k);
float BOOST_MATH_TR1_DECL comp_ellint_2f(float k);
long double BOOST_MATH_TR1_DECL comp_ellint_2l(long double k);

// [5.2.1.6] (complete) elliptic integral of the third kind:
double BOOST_MATH_TR1_DECL comp_ellint_3(double k, double nu);
float BOOST_MATH_TR1_DECL comp_ellint_3f(float k, float nu);
long double BOOST_MATH_TR1_DECL comp_ellint_3l(long double k, long double nu);
#if 0
// [5.2.1.7] confluent hypergeometric functions:
double BOOST_MATH_TR1_DECL conf_hyperg(double a, double c, double x);
float BOOST_MATH_TR1_DECL conf_hypergf(float a, float c, float x);
long double BOOST_MATH_TR1_DECL conf_hypergl(long double a, long double c, long double x);
#endif
// [5.2.1.8] regular modified cylindrical Bessel functions:
double BOOST_MATH_TR1_DECL cyl_bessel_i(double nu, double x);
float BOOST_MATH_TR1_DECL cyl_bessel_if(float nu, float x);
long double BOOST_MATH_TR1_DECL cyl_bessel_il(long double nu, long double x);

// [5.2.1.9] cylindrical Bessel functions (of the first kind):
double BOOST_MATH_TR1_DECL cyl_bessel_j(double nu, double x);
float BOOST_MATH_TR1_DECL cyl_bessel_jf(float nu, float x);
long double BOOST_MATH_TR1_DECL cyl_bessel_jl(long double nu, long double x);

// [5.2.1.10] irregular modified cylindrical Bessel functions:
double BOOST_MATH_TR1_DECL cyl_bessel_k(double nu, double x);
float BOOST_MATH_TR1_DECL cyl_bessel_kf(float nu, float x);
long double BOOST_MATH_TR1_DECL cyl_bessel_kl(long double nu, long double x);

// [5.2.1.11] cylindrical Neumann functions;
// cylindrical Bessel functions (of the second kind):
double BOOST_MATH_TR1_DECL cyl_neumann(double nu, double x);
float BOOST_MATH_TR1_DECL cyl_neumannf(float nu, float x);
long double BOOST_MATH_TR1_DECL cyl_neumannl(long double nu, long double x);

// [5.2.1.12] (incomplete) elliptic integral of the first kind:
double BOOST_MATH_TR1_DECL ellint_1(double k, double phi);
float BOOST_MATH_TR1_DECL ellint_1f(float k, float phi);
long double BOOST_MATH_TR1_DECL ellint_1l(long double k, long double phi);

// [5.2.1.13] (incomplete) elliptic integral of the second kind:
double BOOST_MATH_TR1_DECL ellint_2(double k, double phi);
float BOOST_MATH_TR1_DECL ellint_2f(float k, float phi);
long double BOOST_MATH_TR1_DECL ellint_2l(long double k, long double phi);

// [5.2.1.14] (incomplete) elliptic integral of the third kind:
double BOOST_MATH_TR1_DECL ellint_3(double k, double nu, double phi);
float BOOST_MATH_TR1_DECL ellint_3f(float k, float nu, float phi);
long double BOOST_MATH_TR1_DECL ellint_3l(long double k, long double nu, long double phi);

// [5.2.1.15] exponential integral:
double BOOST_MATH_TR1_DECL expint(double x);
float BOOST_MATH_TR1_DECL expintf(float x);
long double BOOST_MATH_TR1_DECL expintl(long double x);

// [5.2.1.16] Hermite polynomials:
double BOOST_MATH_TR1_DECL hermite(unsigned n, double x);
float BOOST_MATH_TR1_DECL hermitef(unsigned n, float x);
long double BOOST_MATH_TR1_DECL hermitel(unsigned n, long double x);

#if 0
// [5.2.1.17] hypergeometric functions:
double BOOST_MATH_TR1_DECL hyperg(double a, double b, double c, double x);
float BOOST_MATH_TR1_DECL hypergf(float a, float b, float c, float x);
long double BOOST_MATH_TR1_DECL hypergl(long double a, long double b, long double c,
long double x);
#endif

// [5.2.1.18] Laguerre polynomials:
double BOOST_MATH_TR1_DECL laguerre(unsigned n, double x);
float BOOST_MATH_TR1_DECL laguerref(unsigned n, float x);
long double BOOST_MATH_TR1_DECL laguerrel(unsigned n, long double x);

// [5.2.1.19] Legendre polynomials:
double BOOST_MATH_TR1_DECL legendre(unsigned l, double x);
float BOOST_MATH_TR1_DECL legendref(unsigned l, float x);
long double BOOST_MATH_TR1_DECL legendrel(unsigned l, long double x);

// [5.2.1.20] Riemann zeta function:
double BOOST_MATH_TR1_DECL riemann_zeta(double);
float BOOST_MATH_TR1_DECL riemann_zetaf(float);
long double BOOST_MATH_TR1_DECL riemann_zetal(long double);

// [5.2.1.21] spherical Bessel functions (of the first kind):
double BOOST_MATH_TR1_DECL sph_bessel(unsigned n, double x);
float BOOST_MATH_TR1_DECL sph_besself(unsigned n, float x);
long double BOOST_MATH_TR1_DECL sph_bessell(unsigned n, long double x);

// [5.2.1.22] spherical associated Legendre functions:
double BOOST_MATH_TR1_DECL sph_legendre(unsigned l, unsigned m, double theta);
float BOOST_MATH_TR1_DECL sph_legendref(unsigned l, unsigned m, float theta);
long double BOOST_MATH_TR1_DECL sph_legendrel(unsigned l, unsigned m, long double theta);

// [5.2.1.23] spherical Neumann functions;
// spherical Bessel functions (of the second kind):
double BOOST_MATH_TR1_DECL sph_neumann(unsigned n, double x);
float BOOST_MATH_TR1_DECL sph_neumannf(unsigned n, float x);
long double BOOST_MATH_TR1_DECL sph_neumannl(unsigned n, long double x);

#ifdef __cplusplus

}}}}  // namespaces

#include <boost/math/tools/promotion.hpp>

namespace boost{ namespace math{ namespace tr1{
//
// Declare overload of the functions which forward to the
// C interfaces:
//
// C99 Functions:
inline float acosh(float x)
{ return boost::math::tr1::acoshf(x); }
inline long double acosh(long double x)
{ return boost::math::tr1::acoshl(x); }
template <class T>
inline typename tools::promote_args<T>::type acosh(T x)
{ return boost::math::tr1::acosh(static_cast<typename tools::promote_args<T>::type>(x)); }

inline float asinh(float x){ return boost::math::tr1::asinhf(x); }
inline long double asinh(long double x){ return boost::math::tr1::asinhl(x); }
template <class T>
inline typename tools::promote_args<T>::type asinh(T x)
{ return boost::math::tr1::asinh(static_cast<typename tools::promote_args<T>::type>(x)); }

inline float atanh(float x){ return boost::math::tr1::atanhf(x); }
inline long double atanh(long double x){ return boost::math::tr1::atanhl(x); }
template <class T>
inline typename tools::promote_args<T>::type atanh(T x)
{ return boost::math::tr1::atanh(static_cast<typename tools::promote_args<T>::type>(x)); }

inline float cbrt(float x)
{ return boost::math::tr1::cbrtf(x); }
inline long double cbrt(long double x)
{ return boost::math::tr1::cbrtl(x); }
template <class T>
inline typename tools::promote_args<T>::type cbrt(T x)
{ return boost::math::tr1::cbrt(static_cast<typename tools::promote_args<T>::type>(x)); }

inline float copysign(float x, float y)
{ return boost::math::tr1::copysignf(x, y); }
inline long double copysign(long double x, long double y)
{ return boost::math::tr1::copysignl(x, y); }
template <class T1, class T2>
inline typename tools::promote_args<T1, T2>::type copysign(T1 x, T2 y)
{ return boost::math::tr1::copysign(static_cast<typename tools::promote_args<T1, T2>::type>(x), static_cast<typename tools::promote_args<T1, T2>::type>(y)); }

inline float erf(float x)
{ return boost::math::tr1::erff(x); }
inline long double erf(long double x)
{ return boost::math::tr1::erfl(x); }
template <class T>
inline typename tools::promote_args<T>::type erf(T x)
{ return boost::math::tr1::erf(static_cast<typename tools::promote_args<T>::type>(x)); }

inline float erfc(float x)
{ return boost::math::tr1::erfcf(x); }
inline long double erfc(long double x)
{ return boost::math::tr1::erfcl(x); }
template <class T>
inline typename tools::promote_args<T>::type erfc(T x)
{ return boost::math::tr1::erfc(static_cast<typename tools::promote_args<T>::type>(x)); }
#if 0
double exp2(double x);
float exp2f(float x);
long double exp2l(long double x);
#endif
inline float expm1f(float x)
{ return boost::math::tr1::boost_expm1f(x); }
inline double expm1(double x)
{ return boost::math::tr1::boost_expm1(x); }
inline long double expm1l(long double x)
{ return boost::math::tr1::boost_expm1l(x); }
inline float expm1(float x)
{ return boost::math::tr1::expm1f(x); }
inline long double expm1(long double x)
{ return boost::math::tr1::expm1l(x); }
template <class T>
inline typename tools::promote_args<T>::type expm1(T x)
{ return boost::math::tr1::expm1(static_cast<typename tools::promote_args<T>::type>(x)); }
#if 0
double fdim(double x, double y);
float fdimf(float x, float y);
long double fdiml(long double x, long double y);
double fma(double x, double y, double z);
float fmaf(float x, float y, float z);
long double fmal(long double x, long double y, long double z);
#endif
inline float fmax(float x, float y)
{ return boost::math::tr1::fmaxf(x, y); }
inline long double fmax(long double x, long double y)
{ return boost::math::tr1::fmaxl(x, y); }
template <class T1, class T2>
inline typename tools::promote_args<T1, T2>::type fmax(T1 x, T2 y)
{ return boost::math::tr1::fmax(static_cast<typename tools::promote_args<T1, T2>::type>(x), static_cast<typename tools::promote_args<T1, T2>::type>(y)); }

inline float fmin(float x, float y)
{ return boost::math::tr1::fminf(x, y); }
inline long double fmin(long double x, long double y)
{ return boost::math::tr1::fminl(x, y); }
template <class T1, class T2>
inline typename tools::promote_args<T1, T2>::type fmin(T1 x, T2 y)
{ return boost::math::tr1::fmin(static_cast<typename tools::promote_args<T1, T2>::type>(x), static_cast<typename tools::promote_args<T1, T2>::type>(y)); }

inline float hypot(float x, float y)
{ return boost::math::tr1::hypotf(x, y); }
inline long double hypot(long double x, long double y)
{ return boost::math::tr1::hypotl(x, y); }
template <class T1, class T2>
inline typename tools::promote_args<T1, T2>::type hypot(T1 x, T2 y)
{ return boost::math::tr1::hypot(static_cast<typename tools::promote_args<T1, T2>::type>(x), static_cast<typename tools::promote_args<T1, T2>::type>(y)); }
#if 0
int ilogb(double x);
int ilogbf(float x);
int ilogbl(long double x);
#endif
inline float lgamma(float x)
{ return boost::math::tr1::lgammaf(x); }
inline long double lgamma(long double x)
{ return boost::math::tr1::lgammal(x); }
template <class T>
inline typename tools::promote_args<T>::type lgamma(T x)
{ return boost::math::tr1::lgamma(static_cast<typename tools::promote_args<T>::type>(x)); }
#if 0
long long llrint(double x);
long long llrintf(float x);
long long llrintl(long double x);
#endif
inline long long llround(float x)
{ return boost::math::tr1::llroundf(x); }
inline long long llround(long double x)
{ return boost::math::tr1::llroundl(x); }
template <class T>
inline long long llround(T x)
{ return llround(static_cast<double>(x)); }

inline float log1pf(float x)
{ return boost::math::tr1::boost_log1pf(x); }
inline double log1p(double x)
{ return boost::math::tr1::boost_log1p(x); }
inline long double log1pl(long double x)
{ return boost::math::tr1::boost_log1pl(x); }
inline float log1p(float x)
{ return boost::math::tr1::log1pf(x); }
inline long double log1p(long double x)
{ return boost::math::tr1::log1pl(x); }
template <class T>
inline typename tools::promote_args<T>::type log1p(T x)
{ return boost::math::tr1::log1p(static_cast<typename tools::promote_args<T>::type>(x)); }
#if 0
double log2(double x);
float log2f(float x);
long double log2l(long double x);

double logb(double x);
float logbf(float x);
long double logbl(long double x);
long lrint(double x);
long lrintf(float x);
long lrintl(long double x);
#endif
inline long lround(float x)
{ return boost::math::tr1::lroundf(x); }
inline long lround(long double x)
{ return boost::math::tr1::lroundl(x); }
template <class T>
long lround(T x)
{ return boost::math::tr1::lround(static_cast<double>(x)); }
#if 0
double nan(const char *str);
float nanf(const char *str);
long double nanl(const char *str);
double nearbyint(double x);
float nearbyintf(float x);
long double nearbyintl(long double x);
#endif
inline float nextafterf(float x, float y)
{ return boost::math::tr1::boost_nextafterf(x, y); }
inline double nextafter(double x, double y)
{ return boost::math::tr1::boost_nextafter(x, y); }
inline long double nextafterl(long double x, long double y)
{ return boost::math::tr1::boost_nextafterl(x, y); }
inline float nextafter(float x, float y)
{ return boost::math::tr1::nextafterf(x, y); }
inline long double nextafter(long double x, long double y)
{ return boost::math::tr1::nextafterl(x, y); }
template <class T1, class T2>
inline typename tools::promote_args<T1, T2>::type nextafter(T1 x, T2 y)
{ return boost::math::tr1::nextafter(static_cast<typename tools::promote_args<T1, T2>::type>(x), static_cast<typename tools::promote_args<T1, T2>::type>(y)); }

inline float nexttoward(float x, long double y)
{ return boost::math::tr1::nexttowardf(x, y); }
inline long double nexttoward(long double x, long double y)
{ return boost::math::tr1::nexttowardl(x, y); }
template <class T1, class T2>
inline typename tools::promote_args<T1, T2>::type nexttoward(T1 x, T2 y)
{ return boost::math::tr1::nexttoward(static_cast<typename tools::promote_args<T1, T2>::type>(x), static_cast<long double>(y)); }
#if 0
double remainder(double x, double y);
float remainderf(float x, float y);
long double remainderl(long double x, long double y);
double remquo(double x, double y, int *pquo);
float remquof(float x, float y, int *pquo);
long double remquol(long double x, long double y, int *pquo);
double rint(double x);
float rintf(float x);
long double rintl(long double x);
#endif
inline float round(float x)
{ return boost::math::tr1::roundf(x); }
inline long double round(long double x)
{ return boost::math::tr1::roundl(x); }
template <class T>
inline typename tools::promote_args<T>::type round(T x)
{ return boost::math::tr1::round(static_cast<typename tools::promote_args<T>::type>(x)); }
#if 0
double scalbln(double x, long ex);
float scalblnf(float x, long ex);
long double scalblnl(long double x, long ex);
double scalbn(double x, int ex);
float scalbnf(float x, int ex);
long double scalbnl(long double x, int ex);
#endif
inline float tgamma(float x)
{ return boost::math::tr1::tgammaf(x); }
inline long double tgamma(long double x)
{ return boost::math::tr1::tgammal(x); }
template <class T>
inline typename tools::promote_args<T>::type tgamma(T x)
{ return boost::math::tr1::tgamma(static_cast<typename tools::promote_args<T>::type>(x)); }

inline float trunc(float x)
{ return boost::math::tr1::truncf(x); }
inline long double trunc(long double x)
{ return boost::math::tr1::truncl(x); }
template <class T>
inline typename tools::promote_args<T>::type trunc(T x)
{ return boost::math::tr1::trunc(static_cast<typename tools::promote_args<T>::type>(x)); }

# define NO_MACRO_EXPAND /**/
// C99 macros defined as C++ templates
template<class T> bool signbit NO_MACRO_EXPAND(T x)
{ BOOST_STATIC_ASSERT(sizeof(T) == 0); return false; } // must not be instantiated
template<> bool BOOST_MATH_TR1_DECL signbit<float> NO_MACRO_EXPAND(float x);
template<> bool BOOST_MATH_TR1_DECL signbit<double> NO_MACRO_EXPAND(double x);
template<> bool BOOST_MATH_TR1_DECL signbit<long double> NO_MACRO_EXPAND(long double x);

template<class T> int fpclassify NO_MACRO_EXPAND(T x)
{ BOOST_STATIC_ASSERT(sizeof(T) == 0); return false; } // must not be instantiated
template<> int BOOST_MATH_TR1_DECL fpclassify<float> NO_MACRO_EXPAND(float x);
template<> int BOOST_MATH_TR1_DECL fpclassify<double> NO_MACRO_EXPAND(double x);
template<> int BOOST_MATH_TR1_DECL fpclassify<long double> NO_MACRO_EXPAND(long double x);

template<class T> bool isfinite NO_MACRO_EXPAND(T x)
{ BOOST_STATIC_ASSERT(sizeof(T) == 0); return false; } // must not be instantiated
template<> bool BOOST_MATH_TR1_DECL isfinite<float> NO_MACRO_EXPAND(float x);
template<> bool BOOST_MATH_TR1_DECL isfinite<double> NO_MACRO_EXPAND(double x);
template<> bool BOOST_MATH_TR1_DECL isfinite<long double> NO_MACRO_EXPAND(long double x);

template<class T> bool isinf NO_MACRO_EXPAND(T x)
{ BOOST_STATIC_ASSERT(sizeof(T) == 0); return false; } // must not be instantiated
template<> bool BOOST_MATH_TR1_DECL isinf<float> NO_MACRO_EXPAND(float x);
template<> bool BOOST_MATH_TR1_DECL isinf<double> NO_MACRO_EXPAND(double x);
template<> bool BOOST_MATH_TR1_DECL isinf<long double> NO_MACRO_EXPAND(long double x);

template<class T> bool isnan NO_MACRO_EXPAND(T x)
{ BOOST_STATIC_ASSERT(sizeof(T) == 0); return false; } // must not be instantiated
template<> bool BOOST_MATH_TR1_DECL isnan<float> NO_MACRO_EXPAND(float x);
template<> bool BOOST_MATH_TR1_DECL isnan<double> NO_MACRO_EXPAND(double x);
template<> bool BOOST_MATH_TR1_DECL isnan<long double> NO_MACRO_EXPAND(long double x);

template<class T> bool isnormal NO_MACRO_EXPAND(T x)
{ BOOST_STATIC_ASSERT(sizeof(T) == 0); return false; } // must not be instantiated
template<> bool BOOST_MATH_TR1_DECL isnormal<float> NO_MACRO_EXPAND(float x);
template<> bool BOOST_MATH_TR1_DECL isnormal<double> NO_MACRO_EXPAND(double x);
template<> bool BOOST_MATH_TR1_DECL isnormal<long double> NO_MACRO_EXPAND(long double x);

#undef NO_MACRO_EXPAND   
   
// [5.2.1.1] associated Laguerre polynomials:
inline float assoc_laguerre(unsigned n, unsigned m, float x)
{ return boost::math::tr1::assoc_laguerref(n, m, x); }
inline long double assoc_laguerre(unsigned n, unsigned m, long double x)
{ return boost::math::tr1::assoc_laguerrel(n, m, x); }
template <class T> 
inline typename tools::promote_args<T>::type assoc_laguerre(unsigned n, unsigned m, T x)
{ return boost::math::tr1::assoc_laguerre(n, m, static_cast<typename tools::promote_args<T>::type>(x)); }

// [5.2.1.2] associated Legendre functions:
inline float assoc_legendre(unsigned l, unsigned m, float x)
{ return boost::math::tr1::assoc_legendref(l, m, x); }
inline long double assoc_legendre(unsigned l, unsigned m, long double x)
{ return boost::math::tr1::assoc_legendrel(l, m, x); }
template <class T>
inline typename tools::promote_args<T>::type assoc_legendre(unsigned l, unsigned m, T x)
{ return boost::math::tr1::assoc_legendre(l, m, static_cast<typename tools::promote_args<T>::type>(x)); }

// [5.2.1.3] beta function:
inline float beta(float x, float y)
{ return boost::math::tr1::betaf(x, y); }
inline long double beta(long double x, long double y)
{ return boost::math::tr1::betal(x, y); }
template <class T1, class T2>
inline typename tools::promote_args<T1, T2>::type beta(T2 x, T1 y)
{ return boost::math::tr1::beta(static_cast<typename tools::promote_args<T1, T2>::type>(x), static_cast<typename tools::promote_args<T1, T2>::type>(y)); }

// [5.2.1.4] (complete) elliptic integral of the first kind:
inline float comp_ellint_1(float k)
{ return boost::math::tr1::comp_ellint_1f(k); }
inline long double comp_ellint_1(long double k)
{ return boost::math::tr1::comp_ellint_1l(k); }
template <class T>
inline typename tools::promote_args<T>::type comp_ellint_1(T k)
{ return boost::math::tr1::comp_ellint_1(static_cast<typename tools::promote_args<T>::type>(k)); }

// [5.2.1.5] (complete) elliptic integral of the second kind:
inline float comp_ellint_2(float k)
{ return boost::math::tr1::comp_ellint_2f(k); }
inline long double comp_ellint_2(long double k)
{ return boost::math::tr1::comp_ellint_2l(k); }
template <class T>
inline typename tools::promote_args<T>::type comp_ellint_2(T k)
{ return boost::math::tr1::comp_ellint_2(static_cast<typename tools::promote_args<T>::type>(k)); }

// [5.2.1.6] (complete) elliptic integral of the third kind:
inline float comp_ellint_3(float k, float nu)
{ return boost::math::tr1::comp_ellint_3f(k, nu); }
inline long double comp_ellint_3(long double k, long double nu)
{ return boost::math::tr1::comp_ellint_3l(k, nu); }
template <class T1, class T2>
inline typename tools::promote_args<T1, T2>::type comp_ellint_3(T1 k, T2 nu)
{ return boost::math::tr1::comp_ellint_3(static_cast<typename tools::promote_args<T1, T2>::type>(k), static_cast<typename tools::promote_args<T1, T2>::type>(nu)); }

#if 0
// [5.2.1.7] confluent hypergeometric functions:
double conf_hyperg(double a, double c, double x);
float conf_hypergf(float a, float c, float x);
long double conf_hypergl(long double a, long double c, long double x);
#endif

// [5.2.1.8] regular modified cylindrical Bessel functions:
inline float cyl_bessel_i(float nu, float x)
{ return boost::math::tr1::cyl_bessel_if(nu, x); }
inline long double cyl_bessel_i(long double nu, long double x)
{ return boost::math::tr1::cyl_bessel_il(nu, x); }
template <class T1, class T2>
inline typename tools::promote_args<T1, T2>::type cyl_bessel_i(T1 nu, T2 x)
{ return boost::math::tr1::cyl_bessel_i(static_cast<typename tools::promote_args<T1, T2>::type>(nu), static_cast<typename tools::promote_args<T1, T2>::type>(x)); }

// [5.2.1.9] cylindrical Bessel functions (of the first kind):
inline float cyl_bessel_j(float nu, float x)
{ return boost::math::tr1::cyl_bessel_jf(nu, x); }
inline long double cyl_bessel_j(long double nu, long double x)
{ return boost::math::tr1::cyl_bessel_jl(nu, x); }
template <class T1, class T2>
inline typename tools::promote_args<T1, T2>::type cyl_bessel_j(T1 nu, T2 x)
{ return boost::math::tr1::cyl_bessel_j(static_cast<typename tools::promote_args<T1, T2>::type>(nu), static_cast<typename tools::promote_args<T1, T2>::type>(x)); }

// [5.2.1.10] irregular modified cylindrical Bessel functions:
inline float cyl_bessel_k(float nu, float x)
{ return boost::math::tr1::cyl_bessel_kf(nu, x); }
inline long double cyl_bessel_k(long double nu, long double x)
{ return boost::math::tr1::cyl_bessel_kl(nu, x); }
template <class T1, class T2>
inline typename tools::promote_args<T1, T2>::type cyl_bessel_k(T1 nu, T2 x)
{ return boost::math::tr1::cyl_bessel_k(static_cast<typename tools::promote_args<T1, T2>::type>(nu), static_cast<typename tools::promote_args<T1, T2>::type>(x)); }

// [5.2.1.11] cylindrical Neumann functions;
// cylindrical Bessel functions (of the second kind):
inline float cyl_neumann(float nu, float x)
{ return boost::math::tr1::cyl_neumannf(nu, x); }
inline long double cyl_neumann(long double nu, long double x)
{ return boost::math::tr1::cyl_neumannl(nu, x); }
template <class T1, class T2>
inline typename tools::promote_args<T1, T2>::type cyl_neumann(T1 nu, T2 x)
{ return boost::math::tr1::cyl_neumann(static_cast<typename tools::promote_args<T1, T2>::type>(nu), static_cast<typename tools::promote_args<T1, T2>::type>(x)); }

// [5.2.1.12] (incomplete) elliptic integral of the first kind:
inline float ellint_1(float k, float phi)
{ return boost::math::tr1::ellint_1f(k, phi); }
inline long double ellint_1(long double k, long double phi)
{ return boost::math::tr1::ellint_1l(k, phi); }
template <class T1, class T2>
inline typename tools::promote_args<T1, T2>::type ellint_1(T1 k, T2 phi)
{ return boost::math::tr1::ellint_1(static_cast<typename tools::promote_args<T1, T2>::type>(k), static_cast<typename tools::promote_args<T1, T2>::type>(phi)); }

// [5.2.1.13] (incomplete) elliptic integral of the second kind:
inline float ellint_2(float k, float phi)
{ return boost::math::tr1::ellint_2f(k, phi); }
inline long double ellint_2(long double k, long double phi)
{ return boost::math::tr1::ellint_2l(k, phi); }
template <class T1, class T2>
inline typename tools::promote_args<T1, T2>::type ellint_2(T1 k, T2 phi)
{ return boost::math::tr1::ellint_2(static_cast<typename tools::promote_args<T1, T2>::type>(k), static_cast<typename tools::promote_args<T1, T2>::type>(phi)); }

// [5.2.1.14] (incomplete) elliptic integral of the third kind:
inline float ellint_3(float k, float nu, float phi)
{ return boost::math::tr1::ellint_3f(k, nu, phi); }
inline long double ellint_3(long double k, long double nu, long double phi)
{ return boost::math::tr1::ellint_3l(k, nu, phi); }
template <class T1, class T2, class T3>
inline typename tools::promote_args<T1, T2, T3>::type ellint_3(T1 k, T2 nu, T3 phi)
{ return boost::math::tr1::ellint_3(static_cast<typename tools::promote_args<T1, T2, T3>::type>(k), static_cast<typename tools::promote_args<T1, T2, T3>::type>(nu), static_cast<typename tools::promote_args<T1, T2, T3>::type>(phi)); }

// [5.2.1.15] exponential integral:
inline float expint(float x)
{ return boost::math::tr1::expintf(x); }
inline long double expint(long double x)
{ return boost::math::tr1::expintl(x); }
template <class T>
inline typename tools::promote_args<T>::type expint(T x)
{ return boost::math::tr1::expint(static_cast<typename tools::promote_args<T>::type>(x)); }

// [5.2.1.16] Hermite polynomials:
inline float hermite(unsigned n, float x)
{ return boost::math::tr1::hermitef(n, x); }
inline long double hermite(unsigned n, long double x)
{ return boost::math::tr1::hermitel(n, x); }
template <class T>
inline typename tools::promote_args<T>::type hermite(unsigned n, T x)
{ return boost::math::tr1::hermite(n, static_cast<typename tools::promote_args<T>::type>(x)); }

#if 0
// [5.2.1.17] hypergeometric functions:
double hyperg(double a, double b, double c, double x);
float hypergf(float a, float b, float c, float x);
long double hypergl(long double a, long double b, long double c,
long double x);
#endif

// [5.2.1.18] Laguerre polynomials:
inline float laguerre(unsigned n, float x)
{ return boost::math::tr1::laguerref(n, x); }
inline long double laguerre(unsigned n, long double x)
{ return boost::math::tr1::laguerrel(n, x); }
template <class T>
inline typename tools::promote_args<T>::type laguerre(unsigned n, T x)
{ return boost::math::tr1::laguerre(n, static_cast<typename tools::promote_args<T>::type>(x)); }

// [5.2.1.19] Legendre polynomials:
inline float legendre(unsigned l, float x)
{ return boost::math::tr1::legendref(l, x); }
inline long double legendre(unsigned l, long double x)
{ return boost::math::tr1::legendrel(l, x); }
template <class T>
inline typename tools::promote_args<T>::type legendre(unsigned l, T x)
{ return boost::math::tr1::legendre(l, static_cast<typename tools::promote_args<T>::type>(x)); }

// [5.2.1.20] Riemann zeta function:
inline float riemann_zeta(float z)
{ return boost::math::tr1::riemann_zetaf(z); }
inline long double riemann_zeta(long double z)
{ return boost::math::tr1::riemann_zetal(z); }
template <class T>
inline typename tools::promote_args<T>::type riemann_zeta(T z)
{ return boost::math::tr1::riemann_zeta(static_cast<typename tools::promote_args<T>::type>(z)); }

// [5.2.1.21] spherical Bessel functions (of the first kind):
inline float sph_bessel(unsigned n, float x)
{ return boost::math::tr1::sph_besself(n, x); }
inline long double sph_bessel(unsigned n, long double x)
{ return boost::math::tr1::sph_bessell(n, x); }
template <class T>
inline typename tools::promote_args<T>::type sph_bessel(unsigned n, T x)
{ return boost::math::tr1::sph_bessel(n, static_cast<typename tools::promote_args<T>::type>(x)); }

// [5.2.1.22] spherical associated Legendre functions:
inline float sph_legendre(unsigned l, unsigned m, float theta)
{ return boost::math::tr1::sph_legendref(l, m, theta); }
inline long double sph_legendre(unsigned l, unsigned m, long double theta)
{ return boost::math::tr1::sph_legendrel(l, m, theta); }
template <class T>
inline typename tools::promote_args<T>::type sph_legendre(unsigned l, unsigned m, T theta)
{ return boost::math::tr1::sph_legendre(l, m, static_cast<typename tools::promote_args<T>::type>(theta)); }

// [5.2.1.23] spherical Neumann functions;
// spherical Bessel functions (of the second kind):
inline float sph_neumann(unsigned n, float x)
{ return boost::math::tr1::sph_neumannf(n, x); }
inline long double sph_neumann(unsigned n, long double x)
{ return boost::math::tr1::sph_neumannl(n, x); }
template <class T>
inline typename tools::promote_args<T>::type sph_neumann(unsigned n, T x)
{ return boost::math::tr1::sph_neumann(n, static_cast<typename tools::promote_args<T>::type>(x)); }

}}} // namespaces

#endif // __cplusplus

#endif // BOOST_MATH_TR1_HPP

