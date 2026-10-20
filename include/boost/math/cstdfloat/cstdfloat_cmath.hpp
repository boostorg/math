///////////////////////////////////////////////////////////////////////////////
// Copyright Christopher Kormanyos 2014.
// Copyright John Maddock 2014.
// Copyright Paul Bristow 2014.
// Distributed under the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)
//

// Implement quadruple-precision <cmath> support.

#ifndef BOOST_MATH_CSTDFLOAT_CMATH_2014_02_15_HPP_
#define BOOST_MATH_CSTDFLOAT_CMATH_2014_02_15_HPP_

#include <boost/math/cstdfloat/cstdfloat_types.hpp>
#include <boost/math/cstdfloat/cstdfloat_limits.hpp>

#if (defined(__GNUC__) && defined(BOOST_HAS_FLOAT128))
//
// This is the only way we can avoid
// warning: non-standard suffix on floating constant [-Wpedantic]
// when building with -Wall -pedantic.  Neither __extension__
// nor #pragma diagnostic ignored work :(
//
#pragma GCC system_header
#endif

#if defined(BOOST_CSTDFLOAT_HAS_INTERNAL_FLOAT128_T) && defined(BOOST_MATH_USE_FLOAT128) && !defined(BOOST_CSTDFLOAT_NO_LIBQUADMATH_SUPPORT)

#include <cstdint>
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <type_traits>
#include <memory>
#include <boost/math/tools/assert.hpp>
#include <boost/math/tools/nothrow.hpp>
#include <boost/math/tools/throw_exception.hpp>

#if defined(_WIN32) && defined(__GNUC__)
  // Several versions of Mingw and probably cygwin too have broken
  // libquadmath implementations that segfault as soon as you call
  // expq or any function that depends on it.
#define BOOST_CSTDFLOAT_BROKEN_FLOAT128_MATH_FUNCTIONS
#endif

// Here is a helper function used for raising the value of a given
// floating-point type to the power of n, where n has integral type.
BOOST_MATH_NAMESPACE_BEGIN
      namespace cstdfloat {
         namespace detail {

            template<class float_type, class integer_type>
            inline float_type pown(const float_type& x, const integer_type p)
            {
               const bool isneg = (x < 0);
               const bool isnan = (x != x);
               const bool isinf = ((!isneg) ? bool(+x > (std::numeric_limits<float_type>::max)())
                  : bool(-x > (std::numeric_limits<float_type>::max)()));

               if (isnan) { return x; }

               if (isinf) { return std::numeric_limits<float_type>::quiet_NaN(); }

               const bool       x_is_neg = (x < 0);
               const float_type abs_x = (x_is_neg ? -x : x);

               if (p < static_cast<integer_type>(0))
               {
                  if (abs_x < (std::numeric_limits<float_type>::min)())
                  {
                     return (x_is_neg ? -std::numeric_limits<float_type>::infinity()
                        : +std::numeric_limits<float_type>::infinity());
                  }
                  else
                  {
                     return float_type(1) / pown(x, static_cast<integer_type>(-p));
                  }
               }

               if (p == static_cast<integer_type>(0))
               {
                  return float_type(1);
               }
               else
               {
                  if (p == static_cast<integer_type>(1)) { return x; }

                  if (abs_x > (std::numeric_limits<float_type>::max)())
                  {
                     return (x_is_neg ? -std::numeric_limits<float_type>::infinity()
                        : +std::numeric_limits<float_type>::infinity());
                  }

                  if      (p == static_cast<integer_type>(2)) { return  (x * x); }
                  else if (p == static_cast<integer_type>(3)) { return ((x * x) * x); }
                  else if (p == static_cast<integer_type>(4)) { const float_type x2 = (x * x); return (x2 * x2); }
                  else
                  {
                     // The variable xn stores the binary powers of x.
                     float_type result(((p % integer_type(2)) != integer_type(0)) ? x : float_type(1));
                     float_type xn(x);

                     integer_type p2 = p;

                     while (integer_type(p2 /= 2) != integer_type(0))
                     {
                        // Square xn for each binary power.
                        xn *= xn;

                        const bool has_binary_power = (integer_type(p2 % integer_type(2)) != integer_type(0));

                        if (has_binary_power)
                        {
                           // Multiply the result with each binary power contained in the exponent.
                           result *= xn;
                        }
                     }

                     return result;
                  }
               }
            }

         }
      }
   BOOST_MATH_NAMESPACE_END // boost::math::cstdfloat::detail

// We will now define preprocessor symbols representing quadruple-precision <cmath> functions.
#if defined(__INTEL_COMPILER)
#define BOOST_CSTDFLOAT_FLOAT128_LDEXP  __ldexpq
#define BOOST_CSTDFLOAT_FLOAT128_FREXP  __frexpq
#define BOOST_CSTDFLOAT_FLOAT128_FABS   __fabsq
#define BOOST_CSTDFLOAT_FLOAT128_FLOOR  __floorq
#define BOOST_CSTDFLOAT_FLOAT128_CEIL   __ceilq
#if !defined(BOOST_CSTDFLOAT_FLOAT128_SQRT)
#define BOOST_CSTDFLOAT_FLOAT128_SQRT   __sqrtq
#endif
#define BOOST_CSTDFLOAT_FLOAT128_TRUNC  __truncq
#define BOOST_CSTDFLOAT_FLOAT128_EXP    __expq
#define BOOST_CSTDFLOAT_FLOAT128_EXPM1  __expm1q
#define BOOST_CSTDFLOAT_FLOAT128_POW    __powq
#define BOOST_CSTDFLOAT_FLOAT128_LOG    __logq
#define BOOST_CSTDFLOAT_FLOAT128_LOG10  __log10q
#define BOOST_CSTDFLOAT_FLOAT128_SIN    __sinq
#define BOOST_CSTDFLOAT_FLOAT128_COS    __cosq
#define BOOST_CSTDFLOAT_FLOAT128_TAN    __tanq
#define BOOST_CSTDFLOAT_FLOAT128_ASIN   __asinq
#define BOOST_CSTDFLOAT_FLOAT128_ACOS   __acosq
#define BOOST_CSTDFLOAT_FLOAT128_ATAN   __atanq
#define BOOST_CSTDFLOAT_FLOAT128_SINH   __sinhq
#define BOOST_CSTDFLOAT_FLOAT128_COSH   __coshq
#define BOOST_CSTDFLOAT_FLOAT128_TANH   __tanhq
#define BOOST_CSTDFLOAT_FLOAT128_ASINH  __asinhq
#define BOOST_CSTDFLOAT_FLOAT128_ACOSH  __acoshq
#define BOOST_CSTDFLOAT_FLOAT128_ATANH  __atanhq
#define BOOST_CSTDFLOAT_FLOAT128_FMOD   __fmodq
#define BOOST_CSTDFLOAT_FLOAT128_ATAN2  __atan2q
#define BOOST_CSTDFLOAT_FLOAT128_LGAMMA __lgammaq
#define BOOST_CSTDFLOAT_FLOAT128_TGAMMA __tgammaq
//   begin more functions
#define BOOST_CSTDFLOAT_FLOAT128_REMAINDER   __remainderq
#define BOOST_CSTDFLOAT_FLOAT128_REMQUO      __remquoq
#define BOOST_CSTDFLOAT_FLOAT128_FMA         __fmaq
#define BOOST_CSTDFLOAT_FLOAT128_FMAX        __fmaxq
#define BOOST_CSTDFLOAT_FLOAT128_FMIN        __fminq
#define BOOST_CSTDFLOAT_FLOAT128_FDIM        __fdimq
#define BOOST_CSTDFLOAT_FLOAT128_NAN         __nanq
//#define BOOST_CSTDFLOAT_FLOAT128_EXP2      __exp2q
#define BOOST_CSTDFLOAT_FLOAT128_LOG2        __log2q
#define BOOST_CSTDFLOAT_FLOAT128_LOG1P       __log1pq
#define BOOST_CSTDFLOAT_FLOAT128_CBRT        __cbrtq
#define BOOST_CSTDFLOAT_FLOAT128_HYPOT       __hypotq
#define BOOST_CSTDFLOAT_FLOAT128_ERF         __erfq
#define BOOST_CSTDFLOAT_FLOAT128_ERFC        __erfcq
#define BOOST_CSTDFLOAT_FLOAT128_LLROUND     __llroundq
#define BOOST_CSTDFLOAT_FLOAT128_LROUND      __lroundq
#define BOOST_CSTDFLOAT_FLOAT128_ROUND       __roundq
#define BOOST_CSTDFLOAT_FLOAT128_NEARBYINT   __nearbyintq
#define BOOST_CSTDFLOAT_FLOAT128_LLRINT      __llrintq
#define BOOST_CSTDFLOAT_FLOAT128_LRINT       __lrintq
#define BOOST_CSTDFLOAT_FLOAT128_RINT        __rintq
#define BOOST_CSTDFLOAT_FLOAT128_MODF        __modfq
#define BOOST_CSTDFLOAT_FLOAT128_SCALBLN     __scalblnq
#define BOOST_CSTDFLOAT_FLOAT128_SCALBN      __scalbnq
#define BOOST_CSTDFLOAT_FLOAT128_ILOGB       __ilogbq
#define BOOST_CSTDFLOAT_FLOAT128_LOGB        __logbq
#define BOOST_CSTDFLOAT_FLOAT128_NEXTAFTER   __nextafterq
//#define BOOST_CSTDFLOAT_FLOAT128_NEXTTOWARD  __nexttowardq
#define BOOST_CSTDFLOAT_FLOAT128_COPYSIGN     __copysignq
#define BOOST_CSTDFLOAT_FLOAT128_SIGNBIT      __signbitq
//#define BOOST_CSTDFLOAT_FLOAT128_FPCLASSIFY __fpclassifyq
//#define BOOST_CSTDFLOAT_FLOAT128_ISFINITE   __isfiniteq
#define BOOST_CSTDFLOAT_FLOAT128_ISINF        __isinfq
#define BOOST_CSTDFLOAT_FLOAT128_ISNAN        __isnanq
//#define BOOST_CSTDFLOAT_FLOAT128_ISNORMAL   __isnormalq
//#define BOOST_CSTDFLOAT_FLOAT128_ISGREATER  __isgreaterq
//#define BOOST_CSTDFLOAT_FLOAT128_ISGREATEREQUAL __isgreaterequalq
//#define BOOST_CSTDFLOAT_FLOAT128_ISLESS         __islessq
//#define BOOST_CSTDFLOAT_FLOAT128_ISLESSEQUAL    __islessequalq
//#define BOOST_CSTDFLOAT_FLOAT128_ISLESSGREATER  __islessgreaterq
//#define BOOST_CSTDFLOAT_FLOAT128_ISUNORDERED    __isunorderedq
//   end more functions
#elif defined(__GNUC__)
#define BOOST_CSTDFLOAT_FLOAT128_LDEXP  ldexpq
#define BOOST_CSTDFLOAT_FLOAT128_FREXP  frexpq
#define BOOST_CSTDFLOAT_FLOAT128_FABS   fabsq
#define BOOST_CSTDFLOAT_FLOAT128_FLOOR  floorq
#define BOOST_CSTDFLOAT_FLOAT128_CEIL   ceilq
#if !defined(BOOST_CSTDFLOAT_FLOAT128_SQRT)
#define BOOST_CSTDFLOAT_FLOAT128_SQRT   sqrtq
#endif
#define BOOST_CSTDFLOAT_FLOAT128_TRUNC  truncq
#define BOOST_CSTDFLOAT_FLOAT128_POW    powq
#define BOOST_CSTDFLOAT_FLOAT128_LOG    logq
#define BOOST_CSTDFLOAT_FLOAT128_LOG10  log10q
#define BOOST_CSTDFLOAT_FLOAT128_SIN    sinq
#define BOOST_CSTDFLOAT_FLOAT128_COS    cosq
#define BOOST_CSTDFLOAT_FLOAT128_TAN    tanq
#define BOOST_CSTDFLOAT_FLOAT128_ASIN   asinq
#define BOOST_CSTDFLOAT_FLOAT128_ACOS   acosq
#define BOOST_CSTDFLOAT_FLOAT128_ATAN   atanq
#define BOOST_CSTDFLOAT_FLOAT128_FMOD   fmodq
#define BOOST_CSTDFLOAT_FLOAT128_ATAN2  atan2q
#define BOOST_CSTDFLOAT_FLOAT128_LGAMMA lgammaq
#if !defined(BOOST_CSTDFLOAT_BROKEN_FLOAT128_MATH_FUNCTIONS)
#define BOOST_CSTDFLOAT_FLOAT128_EXP    expq
#define BOOST_CSTDFLOAT_FLOAT128_EXPM1  expm1q
#define BOOST_CSTDFLOAT_FLOAT128_SINH   sinhq
#define BOOST_CSTDFLOAT_FLOAT128_COSH   coshq
#define BOOST_CSTDFLOAT_FLOAT128_TANH   tanhq
#define BOOST_CSTDFLOAT_FLOAT128_ASINH  asinhq
#define BOOST_CSTDFLOAT_FLOAT128_ACOSH  acoshq
#define BOOST_CSTDFLOAT_FLOAT128_ATANH  atanhq
#define BOOST_CSTDFLOAT_FLOAT128_TGAMMA tgammaq
#else // BOOST_CSTDFLOAT_BROKEN_FLOAT128_MATH_FUNCTIONS
#define BOOST_CSTDFLOAT_FLOAT128_EXP    expq_patch
#define BOOST_CSTDFLOAT_FLOAT128_SINH   sinhq_patch
#define BOOST_CSTDFLOAT_FLOAT128_COSH   coshq_patch
#define BOOST_CSTDFLOAT_FLOAT128_TANH   tanhq_patch
#define BOOST_CSTDFLOAT_FLOAT128_ASINH  asinhq_patch
#define BOOST_CSTDFLOAT_FLOAT128_ACOSH  acoshq_patch
#define BOOST_CSTDFLOAT_FLOAT128_ATANH  atanhq_patch
#define BOOST_CSTDFLOAT_FLOAT128_TGAMMA tgammaq_patch
#endif // BOOST_CSTDFLOAT_BROKEN_FLOAT128_MATH_FUNCTIONS
//   begin more functions
#define BOOST_CSTDFLOAT_FLOAT128_REMAINDER   remainderq
#define BOOST_CSTDFLOAT_FLOAT128_REMQUO      remquoq
#define BOOST_CSTDFLOAT_FLOAT128_FMA         fmaq
#define BOOST_CSTDFLOAT_FLOAT128_FMAX        fmaxq
#define BOOST_CSTDFLOAT_FLOAT128_FMIN        fminq
#define BOOST_CSTDFLOAT_FLOAT128_FDIM        fdimq
#define BOOST_CSTDFLOAT_FLOAT128_NAN         nanq
//#define BOOST_CSTDFLOAT_FLOAT128_EXP2      exp2q
#define BOOST_CSTDFLOAT_FLOAT128_LOG2        log2q
#define BOOST_CSTDFLOAT_FLOAT128_LOG1P       log1pq
#define BOOST_CSTDFLOAT_FLOAT128_CBRT        cbrtq
#define BOOST_CSTDFLOAT_FLOAT128_HYPOT       hypotq
#define BOOST_CSTDFLOAT_FLOAT128_ERF         erfq
#define BOOST_CSTDFLOAT_FLOAT128_ERFC        erfcq
#define BOOST_CSTDFLOAT_FLOAT128_LLROUND     llroundq
#define BOOST_CSTDFLOAT_FLOAT128_LROUND      lroundq
#define BOOST_CSTDFLOAT_FLOAT128_ROUND       roundq
#define BOOST_CSTDFLOAT_FLOAT128_NEARBYINT   nearbyintq
#define BOOST_CSTDFLOAT_FLOAT128_LLRINT      llrintq
#define BOOST_CSTDFLOAT_FLOAT128_LRINT       lrintq
#define BOOST_CSTDFLOAT_FLOAT128_RINT        rintq
#define BOOST_CSTDFLOAT_FLOAT128_MODF        modfq
#define BOOST_CSTDFLOAT_FLOAT128_SCALBLN     scalblnq
#define BOOST_CSTDFLOAT_FLOAT128_SCALBN      scalbnq
#define BOOST_CSTDFLOAT_FLOAT128_ILOGB       ilogbq
#define BOOST_CSTDFLOAT_FLOAT128_LOGB        logbq
#define BOOST_CSTDFLOAT_FLOAT128_NEXTAFTER   nextafterq
//#define BOOST_CSTDFLOAT_FLOAT128_NEXTTOWARD nexttowardq
#define BOOST_CSTDFLOAT_FLOAT128_COPYSIGN    copysignq
#define BOOST_CSTDFLOAT_FLOAT128_SIGNBIT     signbitq
//#define BOOST_CSTDFLOAT_FLOAT128_FPCLASSIFY fpclassifyq
//#define BOOST_CSTDFLOAT_FLOAT128_ISFINITE   isfiniteq
#define BOOST_CSTDFLOAT_FLOAT128_ISINF        isinfq
#define BOOST_CSTDFLOAT_FLOAT128_ISNAN        isnanq
//#define BOOST_CSTDFLOAT_FLOAT128_ISNORMAL   isnormalq
//#define BOOST_CSTDFLOAT_FLOAT128_ISGREATER  isgreaterq
//#define BOOST_CSTDFLOAT_FLOAT128_ISGREATEREQUAL isgreaterequalq
//#define BOOST_CSTDFLOAT_FLOAT128_ISLESS         islessq
//#define BOOST_CSTDFLOAT_FLOAT128_ISLESSEQUAL    islessequalq
//#define BOOST_CSTDFLOAT_FLOAT128_ISLESSGREATER  islessgreaterq
//#define BOOST_CSTDFLOAT_FLOAT128_ISUNORDERED    isunorderedq
//   end more functions
#endif

// Implement quadruple-precision <cmath> functions in the namespace
// boost::math::cstdfloat::detail. Subsequently inject these into the
// std namespace via *using* directive.

// Begin with some forward function declarations. Also implement patches
// for compilers that have broken float128 exponential functions.

extern "C" int quadmath_snprintf(char*, std::size_t, const char*, ...) BOOST_MATH_NOTHROW;

extern "C" BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_LDEXP(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t, int) BOOST_MATH_NOTHROW;
extern "C" BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_FREXP(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t, int*) BOOST_MATH_NOTHROW;
extern "C" BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_FABS(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t) BOOST_MATH_NOTHROW;
extern "C" BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_FLOOR(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t) BOOST_MATH_NOTHROW;
extern "C" BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_CEIL(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t) BOOST_MATH_NOTHROW;
extern "C" BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_SQRT(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t) BOOST_MATH_NOTHROW;
extern "C" BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_TRUNC(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t) BOOST_MATH_NOTHROW;
extern "C" BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_POW(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t) BOOST_MATH_NOTHROW;
extern "C" BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_LOG(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t) BOOST_MATH_NOTHROW;
extern "C" BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_LOG10(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t) BOOST_MATH_NOTHROW;
extern "C" BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_SIN(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t) BOOST_MATH_NOTHROW;
extern "C" BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_COS(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t) BOOST_MATH_NOTHROW;
extern "C" BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_TAN(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t) BOOST_MATH_NOTHROW;
extern "C" BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_ASIN(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t) BOOST_MATH_NOTHROW;
extern "C" BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_ACOS(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t) BOOST_MATH_NOTHROW;
extern "C" BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_ATAN(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t) BOOST_MATH_NOTHROW;
extern "C" BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_FMOD(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t) BOOST_MATH_NOTHROW;
extern "C" BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_ATAN2(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t) BOOST_MATH_NOTHROW;
extern "C" BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_LGAMMA(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t) BOOST_MATH_NOTHROW;

//   begin more functions
extern "C" BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t  BOOST_CSTDFLOAT_FLOAT128_REMAINDER(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t) BOOST_MATH_NOTHROW;
extern "C" BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t  BOOST_CSTDFLOAT_FLOAT128_REMQUO(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t, int*) BOOST_MATH_NOTHROW;
extern "C" BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t  BOOST_CSTDFLOAT_FLOAT128_FMA(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t) BOOST_MATH_NOTHROW;
extern "C" BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t  BOOST_CSTDFLOAT_FLOAT128_FMAX(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t) BOOST_MATH_NOTHROW;
extern "C" BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t  BOOST_CSTDFLOAT_FLOAT128_FMIN(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t) BOOST_MATH_NOTHROW;
extern "C" BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t  BOOST_CSTDFLOAT_FLOAT128_FDIM(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t) BOOST_MATH_NOTHROW;
extern "C" BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t  BOOST_CSTDFLOAT_FLOAT128_NAN(const char*) BOOST_MATH_NOTHROW;
//extern "C" boost::math::cstdfloat::detail::float_internal128_t  BOOST_CSTDFLOAT_FLOAT128_EXP2         (boost::math::cstdfloat::detail::float_internal128_t) BOOST_MATH_NOTHROW;
extern "C" BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t  BOOST_CSTDFLOAT_FLOAT128_LOG2(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t) BOOST_MATH_NOTHROW;
extern "C" BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t  BOOST_CSTDFLOAT_FLOAT128_LOG1P(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t) BOOST_MATH_NOTHROW;
extern "C" BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t  BOOST_CSTDFLOAT_FLOAT128_CBRT(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t) BOOST_MATH_NOTHROW;
extern "C" BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t  BOOST_CSTDFLOAT_FLOAT128_HYPOT(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t) BOOST_MATH_NOTHROW;
extern "C" BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t  BOOST_CSTDFLOAT_FLOAT128_ERF(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t) BOOST_MATH_NOTHROW;
extern "C" BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t  BOOST_CSTDFLOAT_FLOAT128_ERFC(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t) BOOST_MATH_NOTHROW;
extern "C" long long int                                        BOOST_CSTDFLOAT_FLOAT128_LLROUND(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t) BOOST_MATH_NOTHROW;
extern "C" long int                                             BOOST_CSTDFLOAT_FLOAT128_LROUND(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t) BOOST_MATH_NOTHROW;
extern "C" BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t  BOOST_CSTDFLOAT_FLOAT128_ROUND(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t) BOOST_MATH_NOTHROW;
extern "C" BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t  BOOST_CSTDFLOAT_FLOAT128_NEARBYINT(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t) BOOST_MATH_NOTHROW;
extern "C" long long int                                        BOOST_CSTDFLOAT_FLOAT128_LLRINT(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t) BOOST_MATH_NOTHROW;
extern "C" long int                                             BOOST_CSTDFLOAT_FLOAT128_LRINT(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t) BOOST_MATH_NOTHROW;
extern "C" BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t  BOOST_CSTDFLOAT_FLOAT128_RINT(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t) BOOST_MATH_NOTHROW;
extern "C" BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t  BOOST_CSTDFLOAT_FLOAT128_MODF(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t*) BOOST_MATH_NOTHROW;
extern "C" BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t  BOOST_CSTDFLOAT_FLOAT128_SCALBLN(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t, long int) BOOST_MATH_NOTHROW;
extern "C" BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t  BOOST_CSTDFLOAT_FLOAT128_SCALBN(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t, int) BOOST_MATH_NOTHROW;
extern "C" int                                                  BOOST_CSTDFLOAT_FLOAT128_ILOGB(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t) BOOST_MATH_NOTHROW;
extern "C" BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t  BOOST_CSTDFLOAT_FLOAT128_LOGB(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t) BOOST_MATH_NOTHROW;
extern "C" BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t  BOOST_CSTDFLOAT_FLOAT128_NEXTAFTER(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t) BOOST_MATH_NOTHROW;
//extern "C" boost::math::cstdfloat::detail::float_internal128_t  BOOST_CSTDFLOAT_FLOAT128_NEXTTOWARD   (boost::math::cstdfloat::detail::float_internal128_t, boost::math::cstdfloat::detail::float_internal128_t) BOOST_MATH_NOTHROW;
extern "C" BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t  BOOST_CSTDFLOAT_FLOAT128_COPYSIGN(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t) BOOST_MATH_NOTHROW;
extern "C" int                                                  BOOST_CSTDFLOAT_FLOAT128_SIGNBIT(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t) BOOST_MATH_NOTHROW;
//extern "C" int                                                BOOST_CSTDFLOAT_FLOAT128_FPCLASSIFY   (boost::math::cstdfloat::detail::float_internal128_t) BOOST_MATH_NOTHROW;
//extern "C" int                                                BOOST_CSTDFLOAT_FLOAT128_ISFINITE      (boost::math::cstdfloat::detail::float_internal128_t) BOOST_MATH_NOTHROW;
extern "C" int                                                  BOOST_CSTDFLOAT_FLOAT128_ISINF(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t) BOOST_MATH_NOTHROW;
extern "C" int                                                  BOOST_CSTDFLOAT_FLOAT128_ISNAN(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t) BOOST_MATH_NOTHROW;
//extern "C" boost::math::cstdfloat::detail::float_internal128_t  BOOST_CSTDFLOAT_FLOAT128_ISNORMAL   (boost::math::cstdfloat::detail::float_internal128_t) BOOST_MATH_NOTHROW;
//extern "C" int                                                BOOST_CSTDFLOAT_FLOAT128_ISGREATER   (boost::math::cstdfloat::detail::float_internal128_t, boost::math::cstdfloat::detail::float_internal128_t) BOOST_MATH_NOTHROW;
//extern "C" int                                                BOOST_CSTDFLOAT_FLOAT128_ISGREATEREQUAL(boost::math::cstdfloat::detail::float_internal128_t, boost::math::cstdfloat::detail::float_internal128_t) BOOST_MATH_NOTHROW;
//extern "C" int                                                BOOST_CSTDFLOAT_FLOAT128_ISLESS      (boost::math::cstdfloat::detail::float_internal128_t, boost::math::cstdfloat::detail::float_internal128_t) BOOST_MATH_NOTHROW;
//extern "C" int                                                BOOST_CSTDFLOAT_FLOAT128_ISLESSEQUAL   (boost::math::cstdfloat::detail::float_internal128_t, boost::math::cstdfloat::detail::float_internal128_t) BOOST_MATH_NOTHROW;
//extern "C" int                                                BOOST_CSTDFLOAT_FLOAT128_ISLESSGREATER(boost::math::cstdfloat::detail::float_internal128_t, boost::math::cstdfloat::detail::float_internal128_t) BOOST_MATH_NOTHROW;
//extern "C" int                                                BOOST_CSTDFLOAT_FLOAT128_ISUNORDERED   (boost::math::cstdfloat::detail::float_internal128_t, boost::math::cstdfloat::detail::float_internal128_t) BOOST_MATH_NOTHROW;
 //   end more functions

#if !defined(BOOST_CSTDFLOAT_BROKEN_FLOAT128_MATH_FUNCTIONS)

extern "C" BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_EXP(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x) BOOST_MATH_NOTHROW;
extern "C" BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_EXPM1(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x) BOOST_MATH_NOTHROW;
extern "C" BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_SINH(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x) BOOST_MATH_NOTHROW;
extern "C" BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_COSH(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x) BOOST_MATH_NOTHROW;
extern "C" BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_TANH(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x) BOOST_MATH_NOTHROW;
extern "C" BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_ASINH(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x) BOOST_MATH_NOTHROW;
extern "C" BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_ACOSH(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x) BOOST_MATH_NOTHROW;
extern "C" BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_ATANH(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x) BOOST_MATH_NOTHROW;
extern "C" BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_TGAMMA(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x) BOOST_MATH_NOTHROW;
 
#else // BOOST_CSTDFLOAT_BROKEN_FLOAT128_MATH_FUNCTIONS

// Forward declaration of the patched exponent function, exp(x).
inline BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_EXP(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x);

inline BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_EXPM1(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x)
{
   // Compute exp(x) - 1 for x small.

   // Use an order-12 Pade approximation of the exponential function.
   // PadeApproximant[Exp[x] - 1, {x, 0, 12, 12}].

   typedef BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t float_type;

   float_type sum;

   if (x > BOOST_FLOAT128_C(0.693147180559945309417232121458176568075500134360255))
   {
      sum = ::BOOST_CSTDFLOAT_FLOAT128_EXP(x) - float_type(1);
   }
   else
   {
      const float_type x2 = (x * x);

      const float_type top = (((((  float_type(BOOST_FLOAT128_C(2.4087176110456818621091195109360728010934088788572E-13))  * x2
                                  + float_type(BOOST_FLOAT128_C(9.2735628025258751691201101171038802842096241836000E-10))) * x2
                                  + float_type(BOOST_FLOAT128_C(9.0806726962333369656024118266681195742980640005812E-07))) * x2
                                  + float_type(BOOST_FLOAT128_C(3.1055900621118012422360248447204968944099378881988E-04))) * x2
                                  + float_type(BOOST_FLOAT128_C(3.6231884057971014492753623188405797101449275362319E-02))) * x2
                                  + float_type(BOOST_FLOAT128_C(1.00000000000000000000000000000000000000000000000000000)))
                                  ;

      const float_type bot = ((((((((((((  float_type(BOOST_FLOAT128_C(+7.7202487533515444298369215094104897470942592271063E-16))  * x
                                         + float_type(BOOST_FLOAT128_C(-1.2043588055228409310545597554680364005467044394286E-13))) * x
                                         + float_type(BOOST_FLOAT128_C(+9.2735628025258751691201101171038802842096241836000E-12))) * x
                                         + float_type(BOOST_FLOAT128_C(-4.6367814012629375845600550585519401421048120918000E-10))) * x
                                         + float_type(BOOST_FLOAT128_C(+1.6692413044546575304416198210786984511577323530480E-08))) * x
                                         + float_type(BOOST_FLOAT128_C(-4.5403363481166684828012059133340597871490320002906E-07))) * x
                                         + float_type(BOOST_FLOAT128_C(+9.5347063310450038138825324180015255530129672006102E-06))) * x
                                         + float_type(BOOST_FLOAT128_C(-1.5527950310559006211180124223602484472049689440994E-04))) * x
                                         + float_type(BOOST_FLOAT128_C(+1.9409937888198757763975155279503105590062111801242E-03))) * x
                                         + float_type(BOOST_FLOAT128_C(-1.8115942028985507246376811594202898550724637681159E-02))) * x
                                         + float_type(BOOST_FLOAT128_C(+1.1956521739130434782608695652173913043478260869565E-01))) * x
                                         + float_type(BOOST_FLOAT128_C(-0.50000000000000000000000000000000000000000000000000000))) * x
                                         + float_type(BOOST_FLOAT128_C(+1.00000000000000000000000000000000000000000000000000000)))
                                         ;

      sum = (x * top) / bot;
   }

   return sum;
}
inline BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_EXP(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x)
{
   // Patch the expq() function for a subset of broken GCC compilers
   // like GCC 4.7, 4.8 on MinGW.

   typedef BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t float_type;

   // Scale the argument x to the range (-ln2 < x < ln2).
   constexpr float_type one_over_ln2 = float_type(BOOST_FLOAT128_C(1.44269504088896340735992468100189213742664595415299));
   const float_type x_over_ln2 = x * one_over_ln2;

   int n;

   if (x != x)
   {
      // The argument is NaN.
      return std::numeric_limits<float_type>::quiet_NaN();
   }
   else if (::BOOST_CSTDFLOAT_FLOAT128_FABS(x) > BOOST_FLOAT128_C(+0.693147180559945309417232121458176568075500134360255))
   {
      // The absolute value of the argument exceeds ln2.
      n = static_cast<int>(::BOOST_CSTDFLOAT_FLOAT128_FLOOR(x_over_ln2));
   }
   else if (::BOOST_CSTDFLOAT_FLOAT128_FABS(x) < BOOST_FLOAT128_C(+0.693147180559945309417232121458176568075500134360255))
   {
      // The absolute value of the argument is less than ln2.
      n = 0;
   }
   else
   {
      // The absolute value of the argument is exactly equal to ln2 (in the sense of floating-point equality).
      return float_type(2);
   }

   // Check if the argument is very near an integer.
   const float_type floor_of_x = ::BOOST_CSTDFLOAT_FLOAT128_FLOOR(x);

   if (::BOOST_CSTDFLOAT_FLOAT128_FABS(x - floor_of_x) < float_type(BOOST_CSTDFLOAT_FLOAT128_EPS))
   {
      // Return e^n for arguments very near an integer.
      return BOOST_MATH_NAMESPACE::cstdfloat::detail::pown(BOOST_FLOAT128_C(2.71828182845904523536028747135266249775724709369996), static_cast<std::int_fast32_t>(floor_of_x));
   }

   // Compute the scaled argument alpha.
   const float_type alpha = x - (n * BOOST_FLOAT128_C(0.693147180559945309417232121458176568075500134360255));

   // Compute the polynomial approximation of expm1(alpha) and add to it
   // in order to obtain the scaled result.
   const float_type scaled_result = ::BOOST_CSTDFLOAT_FLOAT128_EXPM1(alpha) + float_type(1);

   // Rescale the result and return it.
   return scaled_result * BOOST_MATH_NAMESPACE::cstdfloat::detail::pown(float_type(2), n);
}
inline BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_SINH(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x)
{
   // Patch the sinhq() function for a subset of broken GCC compilers
   // like GCC 4.7, 4.8 on MinGW.
   typedef BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t float_type;

   // Here, we use the following:
   // Set: ex  = exp(x)
   // Set: em1 = expm1(x)
   // Then
   // sinh(x) = (ex - 1/ex) / 2         ; for |x| >= 1
   // sinh(x) = (2em1 + em1^2) / (2ex)  ; for |x| < 1

   const float_type ex = ::BOOST_CSTDFLOAT_FLOAT128_EXP(x);

   if (::BOOST_CSTDFLOAT_FLOAT128_FABS(x) < float_type(+1))
   {
      const float_type em1 = ::BOOST_CSTDFLOAT_FLOAT128_EXPM1(x);

      return ((em1 * 2) + (em1 * em1)) / (ex * 2);
   }
   else
   {
      return (ex - (float_type(1) / ex)) / 2;
   }
}
inline BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_COSH(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x)
{
   // Patch the coshq() function for a subset of broken GCC compilers
   // like GCC 4.7, 4.8 on MinGW.
   typedef BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t float_type;
   const float_type ex = ::BOOST_CSTDFLOAT_FLOAT128_EXP(x);
   return (ex + (float_type(1) / ex)) / 2;
}
inline BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_TANH(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x)
{
   // Patch the tanhq() function for a subset of broken GCC compilers
   // like GCC 4.7, 4.8 on MinGW.
   typedef BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t float_type;
   const float_type ex_plus = ::BOOST_CSTDFLOAT_FLOAT128_EXP(x);
   const float_type ex_minus = (float_type(1) / ex_plus);
   return (ex_plus - ex_minus) / (ex_plus + ex_minus);
}
inline BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_ASINH(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x) BOOST_MATH_NOTHROW
{
   // Patch the asinh() function since quadmath does not have it.
   typedef BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t float_type;
   return ::BOOST_CSTDFLOAT_FLOAT128_LOG(x + ::BOOST_CSTDFLOAT_FLOAT128_SQRT((x * x) + float_type(1)));
}
inline BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_ACOSH(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x) BOOST_MATH_NOTHROW
{
   // Patch the acosh() function since quadmath does not have it.
   typedef BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t float_type;
   const float_type zp(x + float_type(1));
   const float_type zm(x - float_type(1));

   return ::BOOST_CSTDFLOAT_FLOAT128_LOG(x + (zp * ::BOOST_CSTDFLOAT_FLOAT128_SQRT(zm / zp)));
}
inline BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_ATANH(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x) BOOST_MATH_NOTHROW
{
   // Patch the atanh() function since quadmath does not have it.
   typedef BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t float_type;
   return (::BOOST_CSTDFLOAT_FLOAT128_LOG(float_type(1) + x)
      - ::BOOST_CSTDFLOAT_FLOAT128_LOG(float_type(1) - x)) / 2;
}
inline BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_TGAMMA(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x) BOOST_MATH_NOTHROW
{
   // Patch the tgammaq() function for a subset of broken GCC compilers
   // like GCC 4.7, 4.8 on MinGW.
   typedef BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t float_type;

   if (x > float_type(0))
   {
      return ::BOOST_CSTDFLOAT_FLOAT128_EXP(::BOOST_CSTDFLOAT_FLOAT128_LGAMMA(x));
   }
   else if (x < float_type(0))
   {
      // For x < 0, compute tgamma(-x) and use the reflection formula.
      const float_type positive_x = -x;
      float_type gamma_value = ::BOOST_CSTDFLOAT_FLOAT128_TGAMMA(positive_x);
      const float_type floor_of_positive_x = ::BOOST_CSTDFLOAT_FLOAT128_FLOOR(positive_x);

      // Take the reflection checks (slightly adapted) from <boost/math/gamma.hpp>.
      const bool floor_of_z_is_equal_to_z = (positive_x == ::BOOST_CSTDFLOAT_FLOAT128_FLOOR(positive_x));

      constexpr float_type my_pi = BOOST_FLOAT128_C(3.14159265358979323846264338327950288419716939937511);

      if (floor_of_z_is_equal_to_z)
      {
         const bool is_odd = ((std::int32_t(floor_of_positive_x) % std::int32_t(2)) != std::int32_t(0));

         return (is_odd ? -std::numeric_limits<float_type>::infinity()
            : +std::numeric_limits<float_type>::infinity());
      }

      const float_type sinpx_value = x * ::BOOST_CSTDFLOAT_FLOAT128_SIN(my_pi * x);

      gamma_value *= sinpx_value;

      const bool result_is_too_large_to_represent = ((::BOOST_CSTDFLOAT_FLOAT128_FABS(gamma_value) < float_type(1))
         && (((std::numeric_limits<float_type>::max)() * ::BOOST_CSTDFLOAT_FLOAT128_FABS(gamma_value)) < my_pi));

      if (result_is_too_large_to_represent)
      {
         const bool is_odd = ((std::int32_t(floor_of_positive_x) % std::int32_t(2)) != std::int32_t(0));

         return (is_odd ? -std::numeric_limits<float_type>::infinity()
            : +std::numeric_limits<float_type>::infinity());
      }

      gamma_value = -my_pi / gamma_value;

      if ((gamma_value > float_type(0)) || (gamma_value < float_type(0)))
      {
         return gamma_value;
      }
      else
      {
         // The value of gamma is too small to represent. Return 0.0 here.
         return float_type(0);
      }
   }
   else
   {
      // Gamma of zero is complex infinity. Return NaN here.
      return std::numeric_limits<float_type>::quiet_NaN();
   }
}
#endif // BOOST_CSTDFLOAT_BROKEN_FLOAT128_MATH_FUNCTIONS

// Define the quadruple-precision <cmath> functions in the namespace boost::math::cstdfloat::detail.

BOOST_MATH_NAMESPACE_BEGIN
      namespace cstdfloat {
         namespace detail {
            inline   BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t ldexp(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x, int n) { return ::BOOST_CSTDFLOAT_FLOAT128_LDEXP(x, n); }
            inline   BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t frexp(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x, int* pn) { return ::BOOST_CSTDFLOAT_FLOAT128_FREXP(x, pn); }
            inline   BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t fabs(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x) { return ::BOOST_CSTDFLOAT_FLOAT128_FABS(x); }
            inline   BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t abs(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x) { return ::BOOST_CSTDFLOAT_FLOAT128_FABS(x); }
            inline   BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t floor(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x) { return ::BOOST_CSTDFLOAT_FLOAT128_FLOOR(x); }
            inline   BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t ceil(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x) { return ::BOOST_CSTDFLOAT_FLOAT128_CEIL(x); }
            inline   BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t sqrt(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x) { return ::BOOST_CSTDFLOAT_FLOAT128_SQRT(x); }
            inline   BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t trunc(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x) { return ::BOOST_CSTDFLOAT_FLOAT128_TRUNC(x); }
            inline   BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t exp(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x) { return ::BOOST_CSTDFLOAT_FLOAT128_EXP(x); }
            inline   BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t expm1(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x) { return ::BOOST_CSTDFLOAT_FLOAT128_EXPM1(x); }
            inline   BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t pow(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t a) { return ::BOOST_CSTDFLOAT_FLOAT128_POW(x, a); }
            inline   BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t pow(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x, int a) { return ::BOOST_CSTDFLOAT_FLOAT128_POW(x, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t(a)); }
            inline   BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t log(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x) { return ::BOOST_CSTDFLOAT_FLOAT128_LOG(x); }
            inline   BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t log10(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x) { return ::BOOST_CSTDFLOAT_FLOAT128_LOG10(x); }
            inline   BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t sin(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x) { return ::BOOST_CSTDFLOAT_FLOAT128_SIN(x); }
            inline   BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t cos(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x) { return ::BOOST_CSTDFLOAT_FLOAT128_COS(x); }
            inline   BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t tan(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x) { return ::BOOST_CSTDFLOAT_FLOAT128_TAN(x); }
            inline   BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t asin(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x) { return ::BOOST_CSTDFLOAT_FLOAT128_ASIN(x); }
            inline   BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t acos(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x) { return ::BOOST_CSTDFLOAT_FLOAT128_ACOS(x); }
            inline   BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t atan(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x) { return ::BOOST_CSTDFLOAT_FLOAT128_ATAN(x); }
            inline   BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t sinh(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x) { return ::BOOST_CSTDFLOAT_FLOAT128_SINH(x); }
            inline   BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t cosh(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x) { return ::BOOST_CSTDFLOAT_FLOAT128_COSH(x); }
            inline   BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t tanh(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x) { return ::BOOST_CSTDFLOAT_FLOAT128_TANH(x); }
            inline   BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t asinh(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x) { return ::BOOST_CSTDFLOAT_FLOAT128_ASINH(x); }
            inline   BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t acosh(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x) { return ::BOOST_CSTDFLOAT_FLOAT128_ACOSH(x); }
            inline   BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t atanh(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x) { return ::BOOST_CSTDFLOAT_FLOAT128_ATANH(x); }
            inline   BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t fmod(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t a, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t b) { return ::BOOST_CSTDFLOAT_FLOAT128_FMOD(a, b); }
            inline   BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t atan2(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t y, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x) { return ::BOOST_CSTDFLOAT_FLOAT128_ATAN2(y, x); }
            inline   BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t lgamma(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x) { return ::BOOST_CSTDFLOAT_FLOAT128_LGAMMA(x); }
            inline   BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t tgamma(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x) { return ::BOOST_CSTDFLOAT_FLOAT128_TGAMMA(x); }
            //   begin more functions
            inline BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t  remainder(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t y) { return ::BOOST_CSTDFLOAT_FLOAT128_REMAINDER(x, y); }
            inline BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t  remquo(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t y, int* z) { return ::BOOST_CSTDFLOAT_FLOAT128_REMQUO(x, y, z); }
            inline BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t  fma(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t y, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t z) { return BOOST_CSTDFLOAT_FLOAT128_FMA(x, y, z); }

            inline BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t  fmax(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t y) { return ::BOOST_CSTDFLOAT_FLOAT128_FMAX(x, y); }
            template <class T>
            inline typename std::enable_if<
               std::is_convertible<T, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t>::value
               && !std::is_same<T, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t>::value, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t>::type
               fmax(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x, T y) { return ::BOOST_CSTDFLOAT_FLOAT128_FMAX(x, y); }
            template <class T>
            inline typename std::enable_if<
               std::is_convertible<T, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t>::value
               && !std::is_same<T, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t>::value, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t>::type
               fmax(T x, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t y) { return ::BOOST_CSTDFLOAT_FLOAT128_FMAX(x, y); }
            inline BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t  fmin(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t y) { return ::BOOST_CSTDFLOAT_FLOAT128_FMIN(x, y); }
            template <class T>
            inline typename std::enable_if<
               std::is_convertible<T, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t>::value
               && !std::is_same<T, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t>::value, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t>::type
               fmin(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x, T y) { return ::BOOST_CSTDFLOAT_FLOAT128_FMIN(x, y); }
            template <class T>
            inline typename std::enable_if<
               std::is_convertible<T, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t>::value
               && !std::is_same<T, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t>::value, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t>::type
               fmin(T x, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t y) { return ::BOOST_CSTDFLOAT_FLOAT128_FMIN(x, y); }

            inline BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t  fdim(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t y) { return ::BOOST_CSTDFLOAT_FLOAT128_FDIM(x, y); }
            inline BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t  nanq(const char* x) { return ::BOOST_CSTDFLOAT_FLOAT128_NAN(x); }
            inline BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t  exp2(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x)
            {
               return ::BOOST_CSTDFLOAT_FLOAT128_POW(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t(2), x);
            }
            inline BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t  log2(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x) { return ::BOOST_CSTDFLOAT_FLOAT128_LOG2(x); }
            inline BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t  log1p(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x) { return ::BOOST_CSTDFLOAT_FLOAT128_LOG1P(x); }
            inline BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t  cbrt(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x) { return ::BOOST_CSTDFLOAT_FLOAT128_CBRT(x); }
            inline BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t  hypot(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t y, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t z) { return ::BOOST_CSTDFLOAT_FLOAT128_SQRT(x*x + y * y + z * z); }
            inline BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t  hypot(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t y) { return ::BOOST_CSTDFLOAT_FLOAT128_HYPOT(x, y); }
            template <class T>
            inline typename std::enable_if<
               std::is_convertible<T, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t>::value
               && !std::is_same<T, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t>::value, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t>::type
               hypot(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x, T y) { return ::BOOST_CSTDFLOAT_FLOAT128_HYPOT(x, y); }
            template <class T>
            inline typename std::enable_if<
               std::is_convertible<T, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t>::value
               && !std::is_same<T, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t>::value, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t>::type
               hypot(T x, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t y) { return ::BOOST_CSTDFLOAT_FLOAT128_HYPOT(x, y); }


            inline BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t  erf(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x) { return ::BOOST_CSTDFLOAT_FLOAT128_ERF(x); }
            inline BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t  erfc(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x) { return ::BOOST_CSTDFLOAT_FLOAT128_ERFC(x); }
            inline long long int                                        llround(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x) { return ::BOOST_CSTDFLOAT_FLOAT128_LLROUND(x); }
            inline long int                                             lround(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x) { return ::BOOST_CSTDFLOAT_FLOAT128_LROUND(x); }
            inline BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t  round(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x) { return ::BOOST_CSTDFLOAT_FLOAT128_ROUND(x); }
            inline BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t  nearbyint(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x) { return ::BOOST_CSTDFLOAT_FLOAT128_NEARBYINT(x); }
            inline long long int                                        llrint(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x) { return ::BOOST_CSTDFLOAT_FLOAT128_LLRINT(x); }
            inline long int                                             lrint(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x) { return ::BOOST_CSTDFLOAT_FLOAT128_LRINT(x); }
            inline BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t  rint(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x) { return ::BOOST_CSTDFLOAT_FLOAT128_RINT(x); }
            inline BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t  modf(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t* y) { return ::BOOST_CSTDFLOAT_FLOAT128_MODF(x, y); }
            inline BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t  scalbln(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x, long int y) { return ::BOOST_CSTDFLOAT_FLOAT128_SCALBLN(x, y); }
            inline BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t  scalbn(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x, int y) { return ::BOOST_CSTDFLOAT_FLOAT128_SCALBN(x, y); }
            inline int                                                  ilogb(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x) { return ::BOOST_CSTDFLOAT_FLOAT128_ILOGB(x); }
            inline BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t  logb(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x) { return ::BOOST_CSTDFLOAT_FLOAT128_LOGB(x); }
            inline BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t  nextafter(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t y) { return ::BOOST_CSTDFLOAT_FLOAT128_NEXTAFTER(x, y); }
            inline BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t  nexttoward(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t y) { return -(::BOOST_CSTDFLOAT_FLOAT128_NEXTAFTER(-x, -y)); }
            inline BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t  copysign   BOOST_MATH_PREVENT_MACRO_SUBSTITUTION(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t y) { return ::BOOST_CSTDFLOAT_FLOAT128_COPYSIGN(x, y); }
            inline bool                                                 signbit   BOOST_MATH_PREVENT_MACRO_SUBSTITUTION(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x) { return ::BOOST_CSTDFLOAT_FLOAT128_SIGNBIT(x); }
            inline int                                                  fpclassify BOOST_MATH_PREVENT_MACRO_SUBSTITUTION(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x)
            {
               if (::BOOST_CSTDFLOAT_FLOAT128_ISNAN(x))
                  return FP_NAN;
               else if (::BOOST_CSTDFLOAT_FLOAT128_ISINF(x))
                  return FP_INFINITE;
               else if (x == BOOST_FLOAT128_C(0.0))
                  return FP_ZERO;

               if (::BOOST_CSTDFLOAT_FLOAT128_FABS(x) < BOOST_CSTDFLOAT_FLOAT128_MIN)
                  return FP_SUBNORMAL;
               else
                  return FP_NORMAL;
            }
            inline bool                                      isfinite   BOOST_MATH_PREVENT_MACRO_SUBSTITUTION(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x)
            {
               return !::BOOST_CSTDFLOAT_FLOAT128_ISNAN(x) && !::BOOST_CSTDFLOAT_FLOAT128_ISINF(x);
            }
            inline bool                                      isinf      BOOST_MATH_PREVENT_MACRO_SUBSTITUTION(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x) { return ::BOOST_CSTDFLOAT_FLOAT128_ISINF(x); }
            inline bool                                      isnan      BOOST_MATH_PREVENT_MACRO_SUBSTITUTION(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x) { return ::BOOST_CSTDFLOAT_FLOAT128_ISNAN(x); }
            inline bool                                      isnormal   BOOST_MATH_PREVENT_MACRO_SUBSTITUTION(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x) { return BOOST_MATH_NAMESPACE::cstdfloat::detail::fpclassify BOOST_MATH_PREVENT_MACRO_SUBSTITUTION(x) == FP_NORMAL; }
            inline bool                                      isgreater      BOOST_MATH_PREVENT_MACRO_SUBSTITUTION(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t y)
            {
               if (isnan BOOST_MATH_PREVENT_MACRO_SUBSTITUTION(x) || isnan BOOST_MATH_PREVENT_MACRO_SUBSTITUTION(y))
                  return false;
               return x > y;
            }
            template <class T>
            inline typename std::enable_if<
               std::is_convertible<T, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t>::value
               && !std::is_same<T, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t>::value, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t>::type
               isgreater BOOST_MATH_PREVENT_MACRO_SUBSTITUTION(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x, T y) { return isgreater BOOST_MATH_PREVENT_MACRO_SUBSTITUTION(x, (BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t)y); }
            template <class T>
            inline typename std::enable_if<
               std::is_convertible<T, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t>::value
               && !std::is_same<T, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t>::value, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t>::type
               isgreater BOOST_MATH_PREVENT_MACRO_SUBSTITUTION(T x, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t y) { return isgreater BOOST_MATH_PREVENT_MACRO_SUBSTITUTION((BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t)x, y); }

            inline bool                                      isgreaterequal BOOST_MATH_PREVENT_MACRO_SUBSTITUTION(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t y)
            {
               if (isnan BOOST_MATH_PREVENT_MACRO_SUBSTITUTION(x) || isnan BOOST_MATH_PREVENT_MACRO_SUBSTITUTION(y))
                  return false;
               return x >= y;
            }
            template <class T>
            inline typename std::enable_if<
               std::is_convertible<T, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t>::value
               && !std::is_same<T, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t>::value, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t>::type
               isgreaterequal BOOST_MATH_PREVENT_MACRO_SUBSTITUTION(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x, T y) { return isgreaterequal BOOST_MATH_PREVENT_MACRO_SUBSTITUTION(x, (BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t)y); }
            template <class T>
            inline typename std::enable_if<
               std::is_convertible<T, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t>::value
               && !std::is_same<T, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t>::value, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t>::type
               isgreaterequal BOOST_MATH_PREVENT_MACRO_SUBSTITUTION(T x, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t y) { return isgreaterequal BOOST_MATH_PREVENT_MACRO_SUBSTITUTION((BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t)x, y); }

            inline bool                                      isless      BOOST_MATH_PREVENT_MACRO_SUBSTITUTION(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t y)
            {
               if (isnan BOOST_MATH_PREVENT_MACRO_SUBSTITUTION(x) || isnan BOOST_MATH_PREVENT_MACRO_SUBSTITUTION(y))
                  return false;
               return x < y;
            }
            template <class T>
            inline typename std::enable_if<
               std::is_convertible<T, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t>::value
               && !std::is_same<T, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t>::value, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t>::type
               isless BOOST_MATH_PREVENT_MACRO_SUBSTITUTION(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x, T y) { return isless BOOST_MATH_PREVENT_MACRO_SUBSTITUTION(x, (BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t)y); }
            template <class T>
            inline typename std::enable_if<
               std::is_convertible<T, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t>::value
               && !std::is_same<T, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t>::value, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t>::type
               isless BOOST_MATH_PREVENT_MACRO_SUBSTITUTION(T x, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t y) { return isless BOOST_MATH_PREVENT_MACRO_SUBSTITUTION((BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t)x, y); }


            inline bool                                      islessequal   BOOST_MATH_PREVENT_MACRO_SUBSTITUTION(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t y)
            {
               if (isnan BOOST_MATH_PREVENT_MACRO_SUBSTITUTION(x) || isnan BOOST_MATH_PREVENT_MACRO_SUBSTITUTION(y))
                  return false;
               return x <= y;
            }
            template <class T>
            inline typename std::enable_if<
               std::is_convertible<T, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t>::value
               && !std::is_same<T, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t>::value, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t>::type
               islessequal BOOST_MATH_PREVENT_MACRO_SUBSTITUTION(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x, T y) { return islessequal BOOST_MATH_PREVENT_MACRO_SUBSTITUTION(x, (BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t)y); }
            template <class T>
            inline typename std::enable_if<
               std::is_convertible<T, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t>::value
               && !std::is_same<T, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t>::value, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t>::type
               islessequal BOOST_MATH_PREVENT_MACRO_SUBSTITUTION(T x, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t y) { return islessequal BOOST_MATH_PREVENT_MACRO_SUBSTITUTION((BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t)x, y); }


            inline bool                                      islessgreater   BOOST_MATH_PREVENT_MACRO_SUBSTITUTION(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t y)
            {
               if (isnan BOOST_MATH_PREVENT_MACRO_SUBSTITUTION(x) || isnan BOOST_MATH_PREVENT_MACRO_SUBSTITUTION(y))
                  return false;
               return (x < y) || (x > y);
            }
            template <class T>
            inline typename std::enable_if<
               std::is_convertible<T, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t>::value
               && !std::is_same<T, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t>::value, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t>::type
               islessgreater BOOST_MATH_PREVENT_MACRO_SUBSTITUTION(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x, T y) { return islessgreater BOOST_MATH_PREVENT_MACRO_SUBSTITUTION(x, (BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t)y); }
            template <class T>
            inline typename std::enable_if<
               std::is_convertible<T, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t>::value
               && !std::is_same<T, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t>::value, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t>::type
               islessgreater BOOST_MATH_PREVENT_MACRO_SUBSTITUTION(T x, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t y) { return islessgreater BOOST_MATH_PREVENT_MACRO_SUBSTITUTION((BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t)x, y); }


            inline bool                                      isunordered   BOOST_MATH_PREVENT_MACRO_SUBSTITUTION(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t y) { return ::BOOST_CSTDFLOAT_FLOAT128_ISNAN(x) || ::BOOST_CSTDFLOAT_FLOAT128_ISNAN(y); }
            template <class T>
            inline typename std::enable_if<
               std::is_convertible<T, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t>::value
               && !std::is_same<T, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t>::value, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t>::type
               isunordered BOOST_MATH_PREVENT_MACRO_SUBSTITUTION(BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t x, T y) { return isunordered BOOST_MATH_PREVENT_MACRO_SUBSTITUTION(x, (BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t)y); }
            template <class T>
            inline typename std::enable_if<
               std::is_convertible<T, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t>::value
               && !std::is_same<T, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t>::value, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t>::type
               isunordered BOOST_MATH_PREVENT_MACRO_SUBSTITUTION(T x, BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t y) { return isunordered BOOST_MATH_PREVENT_MACRO_SUBSTITUTION((BOOST_MATH_NAMESPACE::cstdfloat::detail::float_internal128_t)x, y); }


            //   end more functions
         }
      }
   BOOST_MATH_NAMESPACE_END // boost::math::cstdfloat::detail

// We will now inject the quadruple-precision <cmath> functions
// into the std namespace. This is done via *using* directive.
namespace std
{
   using BOOST_MATH_NAMESPACE::cstdfloat::detail::ldexp;
   using BOOST_MATH_NAMESPACE::cstdfloat::detail::frexp;
   using BOOST_MATH_NAMESPACE::cstdfloat::detail::fabs;

#if !(defined(_GLIBCXX_USE_FLOAT128) && defined(__GNUC__) && (__GNUC__ >= 7))
#if (defined(__clang__) && !(!defined(__STRICT_ANSI__) && defined(_GLIBCXX_USE_FLOAT128))) || (__GNUC__ <= 6 && !defined(__clang__)) 
   // workaround for clang using libstdc++ and old GCC
   using BOOST_MATH_NAMESPACE::cstdfloat::detail::abs;
#endif
#endif

   using BOOST_MATH_NAMESPACE::cstdfloat::detail::floor;
   using BOOST_MATH_NAMESPACE::cstdfloat::detail::ceil;
   using BOOST_MATH_NAMESPACE::cstdfloat::detail::sqrt;
   using BOOST_MATH_NAMESPACE::cstdfloat::detail::trunc;
   using BOOST_MATH_NAMESPACE::cstdfloat::detail::exp;
   using BOOST_MATH_NAMESPACE::cstdfloat::detail::expm1;
   using BOOST_MATH_NAMESPACE::cstdfloat::detail::pow;
   using BOOST_MATH_NAMESPACE::cstdfloat::detail::log;
   using BOOST_MATH_NAMESPACE::cstdfloat::detail::log10;
   using BOOST_MATH_NAMESPACE::cstdfloat::detail::sin;
   using BOOST_MATH_NAMESPACE::cstdfloat::detail::cos;
   using BOOST_MATH_NAMESPACE::cstdfloat::detail::tan;
   using BOOST_MATH_NAMESPACE::cstdfloat::detail::asin;
   using BOOST_MATH_NAMESPACE::cstdfloat::detail::acos;
   using BOOST_MATH_NAMESPACE::cstdfloat::detail::atan;
   using BOOST_MATH_NAMESPACE::cstdfloat::detail::sinh;
   using BOOST_MATH_NAMESPACE::cstdfloat::detail::cosh;
   using BOOST_MATH_NAMESPACE::cstdfloat::detail::tanh;
   using BOOST_MATH_NAMESPACE::cstdfloat::detail::asinh;
   using BOOST_MATH_NAMESPACE::cstdfloat::detail::acosh;
   using BOOST_MATH_NAMESPACE::cstdfloat::detail::atanh;
   using BOOST_MATH_NAMESPACE::cstdfloat::detail::fmod;
   using BOOST_MATH_NAMESPACE::cstdfloat::detail::atan2;
   using BOOST_MATH_NAMESPACE::cstdfloat::detail::lgamma;
   using BOOST_MATH_NAMESPACE::cstdfloat::detail::tgamma;

   //   begin more functions
   using BOOST_MATH_NAMESPACE::cstdfloat::detail::remainder;
   using BOOST_MATH_NAMESPACE::cstdfloat::detail::remquo;
   using BOOST_MATH_NAMESPACE::cstdfloat::detail::fma;
   using BOOST_MATH_NAMESPACE::cstdfloat::detail::fmax;
   using BOOST_MATH_NAMESPACE::cstdfloat::detail::fmin;
   using BOOST_MATH_NAMESPACE::cstdfloat::detail::fdim;
   using BOOST_MATH_NAMESPACE::cstdfloat::detail::nanq;
   using BOOST_MATH_NAMESPACE::cstdfloat::detail::exp2;
   using BOOST_MATH_NAMESPACE::cstdfloat::detail::log2;
   using BOOST_MATH_NAMESPACE::cstdfloat::detail::log1p;
   using BOOST_MATH_NAMESPACE::cstdfloat::detail::cbrt;
   using BOOST_MATH_NAMESPACE::cstdfloat::detail::hypot;
   using BOOST_MATH_NAMESPACE::cstdfloat::detail::erf;
   using BOOST_MATH_NAMESPACE::cstdfloat::detail::erfc;
   using BOOST_MATH_NAMESPACE::cstdfloat::detail::llround;
   using BOOST_MATH_NAMESPACE::cstdfloat::detail::lround;
   using BOOST_MATH_NAMESPACE::cstdfloat::detail::round;
   using BOOST_MATH_NAMESPACE::cstdfloat::detail::nearbyint;
   using BOOST_MATH_NAMESPACE::cstdfloat::detail::llrint;
   using BOOST_MATH_NAMESPACE::cstdfloat::detail::lrint;
   using BOOST_MATH_NAMESPACE::cstdfloat::detail::rint;
   using BOOST_MATH_NAMESPACE::cstdfloat::detail::modf;
   using BOOST_MATH_NAMESPACE::cstdfloat::detail::scalbln;
   using BOOST_MATH_NAMESPACE::cstdfloat::detail::scalbn;
   using BOOST_MATH_NAMESPACE::cstdfloat::detail::ilogb;
   using BOOST_MATH_NAMESPACE::cstdfloat::detail::logb;
   using BOOST_MATH_NAMESPACE::cstdfloat::detail::nextafter;
   using BOOST_MATH_NAMESPACE::cstdfloat::detail::nexttoward;
   using BOOST_MATH_NAMESPACE::cstdfloat::detail::copysign;
   using BOOST_MATH_NAMESPACE::cstdfloat::detail::signbit;
   using BOOST_MATH_NAMESPACE::cstdfloat::detail::fpclassify;
   using BOOST_MATH_NAMESPACE::cstdfloat::detail::isfinite;
   using BOOST_MATH_NAMESPACE::cstdfloat::detail::isinf;
   using BOOST_MATH_NAMESPACE::cstdfloat::detail::isnan;
   using BOOST_MATH_NAMESPACE::cstdfloat::detail::isnormal;
   using BOOST_MATH_NAMESPACE::cstdfloat::detail::isgreater;
   using BOOST_MATH_NAMESPACE::cstdfloat::detail::isgreaterequal;
   using BOOST_MATH_NAMESPACE::cstdfloat::detail::isless;
   using BOOST_MATH_NAMESPACE::cstdfloat::detail::islessequal;
   using BOOST_MATH_NAMESPACE::cstdfloat::detail::islessgreater;
   using BOOST_MATH_NAMESPACE::cstdfloat::detail::isunordered;
   //   end more functions

} // namespace std

// We will now remove the preprocessor symbols representing quadruple-precision <cmath>
// functions from the preprocessor.

#undef BOOST_CSTDFLOAT_FLOAT128_LDEXP
#undef BOOST_CSTDFLOAT_FLOAT128_FREXP
#undef BOOST_CSTDFLOAT_FLOAT128_FABS
#undef BOOST_CSTDFLOAT_FLOAT128_FLOOR
#undef BOOST_CSTDFLOAT_FLOAT128_CEIL
#undef BOOST_CSTDFLOAT_FLOAT128_SQRT
#undef BOOST_CSTDFLOAT_FLOAT128_TRUNC
#undef BOOST_CSTDFLOAT_FLOAT128_EXP
#undef BOOST_CSTDFLOAT_FLOAT128_EXPM1
#undef BOOST_CSTDFLOAT_FLOAT128_POW
#undef BOOST_CSTDFLOAT_FLOAT128_LOG
#undef BOOST_CSTDFLOAT_FLOAT128_LOG10
#undef BOOST_CSTDFLOAT_FLOAT128_SIN
#undef BOOST_CSTDFLOAT_FLOAT128_COS
#undef BOOST_CSTDFLOAT_FLOAT128_TAN
#undef BOOST_CSTDFLOAT_FLOAT128_ASIN
#undef BOOST_CSTDFLOAT_FLOAT128_ACOS
#undef BOOST_CSTDFLOAT_FLOAT128_ATAN
#undef BOOST_CSTDFLOAT_FLOAT128_SINH
#undef BOOST_CSTDFLOAT_FLOAT128_COSH
#undef BOOST_CSTDFLOAT_FLOAT128_TANH
#undef BOOST_CSTDFLOAT_FLOAT128_ASINH
#undef BOOST_CSTDFLOAT_FLOAT128_ACOSH
#undef BOOST_CSTDFLOAT_FLOAT128_ATANH
#undef BOOST_CSTDFLOAT_FLOAT128_FMOD
#undef BOOST_CSTDFLOAT_FLOAT128_ATAN2
#undef BOOST_CSTDFLOAT_FLOAT128_LGAMMA
#undef BOOST_CSTDFLOAT_FLOAT128_TGAMMA

//   begin more functions
#undef BOOST_CSTDFLOAT_FLOAT128_REMAINDER
#undef BOOST_CSTDFLOAT_FLOAT128_REMQUO
#undef BOOST_CSTDFLOAT_FLOAT128_FMA
#undef BOOST_CSTDFLOAT_FLOAT128_FMAX
#undef BOOST_CSTDFLOAT_FLOAT128_FMIN
#undef BOOST_CSTDFLOAT_FLOAT128_FDIM
#undef BOOST_CSTDFLOAT_FLOAT128_NAN
#undef BOOST_CSTDFLOAT_FLOAT128_EXP2
#undef BOOST_CSTDFLOAT_FLOAT128_LOG2
#undef BOOST_CSTDFLOAT_FLOAT128_LOG1P
#undef BOOST_CSTDFLOAT_FLOAT128_CBRT
#undef BOOST_CSTDFLOAT_FLOAT128_HYPOT
#undef BOOST_CSTDFLOAT_FLOAT128_ERF
#undef BOOST_CSTDFLOAT_FLOAT128_ERFC
#undef BOOST_CSTDFLOAT_FLOAT128_LLROUND
#undef BOOST_CSTDFLOAT_FLOAT128_LROUND
#undef BOOST_CSTDFLOAT_FLOAT128_ROUND
#undef BOOST_CSTDFLOAT_FLOAT128_NEARBYINT
#undef BOOST_CSTDFLOAT_FLOAT128_LLRINT
#undef BOOST_CSTDFLOAT_FLOAT128_LRINT
#undef BOOST_CSTDFLOAT_FLOAT128_RINT
#undef BOOST_CSTDFLOAT_FLOAT128_MODF
#undef BOOST_CSTDFLOAT_FLOAT128_SCALBLN
#undef BOOST_CSTDFLOAT_FLOAT128_SCALBN
#undef BOOST_CSTDFLOAT_FLOAT128_ILOGB
#undef BOOST_CSTDFLOAT_FLOAT128_LOGB
#undef BOOST_CSTDFLOAT_FLOAT128_NEXTAFTER
#undef BOOST_CSTDFLOAT_FLOAT128_NEXTTOWARD
#undef BOOST_CSTDFLOAT_FLOAT128_COPYSIGN
#undef BOOST_CSTDFLOAT_FLOAT128_SIGNBIT
#undef BOOST_CSTDFLOAT_FLOAT128_FPCLASSIFY
#undef BOOST_CSTDFLOAT_FLOAT128_ISFINITE
#undef BOOST_CSTDFLOAT_FLOAT128_ISINF
#undef BOOST_CSTDFLOAT_FLOAT128_ISNAN
#undef BOOST_CSTDFLOAT_FLOAT128_ISNORMAL
#undef BOOST_CSTDFLOAT_FLOAT128_ISGREATER
#undef BOOST_CSTDFLOAT_FLOAT128_ISGREATEREQUAL
#undef BOOST_CSTDFLOAT_FLOAT128_ISLESS
#undef BOOST_CSTDFLOAT_FLOAT128_ISLESSEQUAL
#undef BOOST_CSTDFLOAT_FLOAT128_ISLESSGREATER
#undef BOOST_CSTDFLOAT_FLOAT128_ISUNORDERED
//   end more functions

#endif // Not BOOST_CSTDFLOAT_NO_LIBQUADMATH_SUPPORT (i.e., the user would like to have libquadmath support)

#endif // BOOST_MATH_CSTDFLOAT_CMATH_2014_02_15_HPP_

