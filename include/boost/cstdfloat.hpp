///////////////////////////////////////////////////////////////////////////////
// Copyright Christopher Kormanyos 2014.
// Copyright John Maddock 2014.
// Copyright Paul Bristow 2014.
// Distributed under the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)
//

// <boost/cstdfloat.hpp> implements floating-point typedefs having
// specified widths, as described in N3626 (proposed for C++14).
// See: http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2013/n3626.pdf

#ifndef _BOOST_CSTDFLOAT_2014_01_09_HPP_
  #define _BOOST_CSTDFLOAT_2014_01_09_HPP_

  #include <float.h>
  #include <cstddef>
  #include <limits>
  #include <stdexcept>
  #include <boost/static_assert.hpp>
  #include <boost/throw_exception.hpp>
  #include <boost/math/tools/config.hpp>

  // This is the beginning of the preamble.

  // In this preamble, the preprocessor is used to query certain
  // preprocessor definitions from <float.h>. Based on the results
  // of these queries, an attempt is made to automatically detect
  // the presence of built-in floating-point types having specified
  // widths. These are *thought* to be conformant with IEEE-754,
  // whereby an unequivocal test based on numeric_limits follows below.

  // In addition, macros that are used for initializing floating-point
  // literal values and some basic min/max values are defined.

  // First, we will pre-load certain preprocessor definitions
  // with a dummy value.

  #define BOOST_CSTDFLOAT_MAXIMUM_AVAILABLE_WIDTH  0

  #define BOOST_CSTDFLOAT_HAS_FLOAT16_NATIVE_TYPE  0
  #define BOOST_CSTDFLOAT_HAS_FLOAT32_NATIVE_TYPE  0
  #define BOOST_CSTDFLOAT_HAS_FLOAT64_NATIVE_TYPE  0
  #define BOOST_CSTDFLOAT_HAS_FLOAT80_NATIVE_TYPE  0
  #define BOOST_CSTDFLOAT_HAS_FLOAT128_NATIVE_TYPE 0

  // Ensure that the compiler has a radix-2 floating-point representation.
  #if (!defined(FLT_RADIX) || ((defined(FLT_RADIX) && (FLT_RADIX != 2))))
    #error The compiler does not support radix-2 floating-point types required for <boost/cstdfloat.hpp>.
  #endif

  // Check if built-in float is equivalent to float16_t, float32_t, float64_t, float80_t, or float128_t.
  #if(defined(FLT_MANT_DIG) && defined(FLT_MAX_EXP))
    #if  ((FLT_MANT_DIG == 11) && (FLT_MAX_EXP == 16) && (BOOST_CSTDFLOAT_HAS_FLOAT16_NATIVE_TYPE == 0))
      #define BOOST_CSTDFLOAT_FLOAT16_NATIVE_TYPE float
      #undef  BOOST_CSTDFLOAT_MAXIMUM_AVAILABLE_WIDTH
      #define BOOST_CSTDFLOAT_MAXIMUM_AVAILABLE_WIDTH 16
      #undef  BOOST_CSTDFLOAT_HAS_FLOAT16_NATIVE_TYPE
      #define BOOST_CSTDFLOAT_HAS_FLOAT16_NATIVE_TYPE  1
      #define BOOST_FLOAT16_C(x)  (x ## F)
      #define BOOST_FLOAT_16_MIN  FLT_MIN
      #define BOOST_FLOAT_16_MAX  FLT_MAX
    #elif((FLT_MANT_DIG == 24) && (FLT_MAX_EXP == 128) && (BOOST_CSTDFLOAT_HAS_FLOAT32_NATIVE_TYPE == 0))
      #define BOOST_CSTDFLOAT_FLOAT32_NATIVE_TYPE float
      #undef  BOOST_CSTDFLOAT_MAXIMUM_AVAILABLE_WIDTH
      #define BOOST_CSTDFLOAT_MAXIMUM_AVAILABLE_WIDTH 32
      #undef  BOOST_CSTDFLOAT_HAS_FLOAT32_NATIVE_TYPE
      #define BOOST_CSTDFLOAT_HAS_FLOAT32_NATIVE_TYPE  1
      #define BOOST_FLOAT32_C(x)  (x ## F)
      #define BOOST_FLOAT_32_MIN  FLT_MIN
      #define BOOST_FLOAT_32_MAX  FLT_MAX
    #elif((FLT_MANT_DIG == 53) && (FLT_MAX_EXP == 1024) && (BOOST_CSTDFLOAT_HAS_FLOAT64_NATIVE_TYPE == 0))
      #define BOOST_CSTDFLOAT_FLOAT64_NATIVE_TYPE float
      #undef  BOOST_CSTDFLOAT_MAXIMUM_AVAILABLE_WIDTH
      #define BOOST_CSTDFLOAT_MAXIMUM_AVAILABLE_WIDTH 64
      #undef  BOOST_CSTDFLOAT_HAS_FLOAT64_NATIVE_TYPE
      #define BOOST_CSTDFLOAT_HAS_FLOAT64_NATIVE_TYPE  1
      #define BOOST_FLOAT64_C(x)  (x ## F)
      #define BOOST_FLOAT_64_MIN  FLT_MIN
      #define BOOST_FLOAT_64_MAX  FLT_MAX
    #elif((FLT_MANT_DIG == 64) && (FLT_MAX_EXP == 16384) && (BOOST_CSTDFLOAT_HAS_FLOAT80_NATIVE_TYPE == 0))
      #define BOOST_CSTDFLOAT_FLOAT80_NATIVE_TYPE float
      #undef  BOOST_CSTDFLOAT_MAXIMUM_AVAILABLE_WIDTH
      #define BOOST_CSTDFLOAT_MAXIMUM_AVAILABLE_WIDTH 80
      #undef  BOOST_CSTDFLOAT_HAS_FLOAT80_NATIVE_TYPE
      #define BOOST_CSTDFLOAT_HAS_FLOAT80_NATIVE_TYPE  1
      #define BOOST_FLOAT80_C(x)  (x ## F)
      #define BOOST_FLOAT_80_MIN  FLT_MIN
      #define BOOST_FLOAT_80_MAX  FLT_MAX
    #elif((FLT_MANT_DIG == 113) && (FLT_MAX_EXP == 16384) && (BOOST_CSTDFLOAT_HAS_FLOAT128_NATIVE_TYPE == 0))
      #define BOOST_CSTDFLOAT_FLOAT128_NATIVE_TYPE float
      #undef  BOOST_CSTDFLOAT_MAXIMUM_AVAILABLE_WIDTH
      #define BOOST_CSTDFLOAT_MAXIMUM_AVAILABLE_WIDTH 128
      #undef  BOOST_CSTDFLOAT_HAS_FLOAT128_NATIVE_TYPE
      #define BOOST_CSTDFLOAT_HAS_FLOAT128_NATIVE_TYPE  1
      #define BOOST_FLOAT128_C(x)  (x ## F)
      #define BOOST_FLOAT_128_MIN  FLT_MIN
      #define BOOST_FLOAT_128_MAX  FLT_MAX
    #endif
  #endif

  // Check if built-in double is equivalent to float16_t, float32_t, float64_t, float80_t, or float128_t.
  #if(defined(DBL_MANT_DIG) && defined(DBL_MAX_EXP))
    #if  ((DBL_MANT_DIG == 11) && (DBL_MAX_EXP == 16) && (BOOST_CSTDFLOAT_HAS_FLOAT16_NATIVE_TYPE == 0))
      #define BOOST_CSTDFLOAT_FLOAT16_NATIVE_TYPE double
      #undef  BOOST_CSTDFLOAT_MAXIMUM_AVAILABLE_WIDTH
      #define BOOST_CSTDFLOAT_MAXIMUM_AVAILABLE_WIDTH 16
      #undef  BOOST_CSTDFLOAT_HAS_FLOAT16_NATIVE_TYPE
      #define BOOST_CSTDFLOAT_HAS_FLOAT16_NATIVE_TYPE  1
      #define BOOST_FLOAT16_C(x)  (x)
      #define BOOST_FLOAT_16_MIN  DBL_MIN
      #define BOOST_FLOAT_16_MAX  DBL_MAX
    #elif((DBL_MANT_DIG == 24) && (DBL_MAX_EXP == 128) && (BOOST_CSTDFLOAT_HAS_FLOAT32_NATIVE_TYPE == 0))
      #define BOOST_CSTDFLOAT_FLOAT32_NATIVE_TYPE double
      #undef  BOOST_CSTDFLOAT_MAXIMUM_AVAILABLE_WIDTH
      #define BOOST_CSTDFLOAT_MAXIMUM_AVAILABLE_WIDTH 32
      #undef  BOOST_CSTDFLOAT_HAS_FLOAT32_NATIVE_TYPE
      #define BOOST_CSTDFLOAT_HAS_FLOAT32_NATIVE_TYPE  1
      #define BOOST_FLOAT32_C(x)  (x)
      #define BOOST_FLOAT_32_MIN  DBL_MIN
      #define BOOST_FLOAT_32_MAX  DBL_MAX
    #elif((DBL_MANT_DIG == 53) && (DBL_MAX_EXP == 1024) && (BOOST_CSTDFLOAT_HAS_FLOAT64_NATIVE_TYPE == 0))
      #define BOOST_CSTDFLOAT_FLOAT64_NATIVE_TYPE double
      #undef  BOOST_CSTDFLOAT_MAXIMUM_AVAILABLE_WIDTH
      #define BOOST_CSTDFLOAT_MAXIMUM_AVAILABLE_WIDTH 64
      #undef  BOOST_CSTDFLOAT_HAS_FLOAT64_NATIVE_TYPE
      #define BOOST_CSTDFLOAT_HAS_FLOAT64_NATIVE_TYPE  1
      #define BOOST_FLOAT64_C(x)  (x)
      #define BOOST_FLOAT_64_MIN  DBL_MIN
      #define BOOST_FLOAT_64_MAX  DBL_MAX
    #elif((DBL_MANT_DIG == 64) && (DBL_MAX_EXP == 16384) && (BOOST_CSTDFLOAT_HAS_FLOAT80_NATIVE_TYPE == 0))
      #define BOOST_CSTDFLOAT_FLOAT80_NATIVE_TYPE double
      #undef  BOOST_CSTDFLOAT_MAXIMUM_AVAILABLE_WIDTH
      #define BOOST_CSTDFLOAT_MAXIMUM_AVAILABLE_WIDTH 80
      #undef  BOOST_CSTDFLOAT_HAS_FLOAT80_NATIVE_TYPE
      #define BOOST_CSTDFLOAT_HAS_FLOAT80_NATIVE_TYPE  1
      #define BOOST_FLOAT80_C(x)  (x)
      #define BOOST_FLOAT_80_MIN  DBL_MIN
      #define BOOST_FLOAT_80_MAX  DBL_MAX
    #elif((DBL_MANT_DIG == 113) && (DBL_MAX_EXP == 16384) && (BOOST_CSTDFLOAT_HAS_FLOAT128_NATIVE_TYPE == 0))
      #define BOOST_CSTDFLOAT_FLOAT128_NATIVE_TYPE double
      #undef  BOOST_CSTDFLOAT_MAXIMUM_AVAILABLE_WIDTH
      #define BOOST_CSTDFLOAT_MAXIMUM_AVAILABLE_WIDTH 128
      #undef  BOOST_CSTDFLOAT_HAS_FLOAT128_NATIVE_TYPE
      #define BOOST_CSTDFLOAT_HAS_FLOAT128_NATIVE_TYPE  1
      #define BOOST_FLOAT128_C(x)  (x)
      #define BOOST_FLOAT_128_MIN  DBL_MIN
      #define BOOST_FLOAT_128_MAX  DBL_MAX
    #endif
  #endif

  // Check if built-in long double is equivalent to float16_t, float32_t, float64_t, float80_t, or float128_t.
  #if(defined(LDBL_MANT_DIG) && defined(LDBL_MAX_EXP))
    #if  ((LDBL_MANT_DIG == 11) && (LDBL_MAX_EXP == 16) && (BOOST_CSTDFLOAT_HAS_FLOAT16_NATIVE_TYPE == 0))
      #define BOOST_CSTDFLOAT_FLOAT16_NATIVE_TYPE long double
      #undef  BOOST_CSTDFLOAT_MAXIMUM_AVAILABLE_WIDTH
      #define BOOST_CSTDFLOAT_MAXIMUM_AVAILABLE_WIDTH 16
      #undef  BOOST_CSTDFLOAT_HAS_FLOAT16_NATIVE_TYPE
      #define BOOST_CSTDFLOAT_HAS_FLOAT16_NATIVE_TYPE  1
      #define BOOST_FLOAT16_C(x)  (x ## L)
      #define BOOST_FLOAT_16_MIN  LDBL_MIN
      #define BOOST_FLOAT_16_MAX  LDBL_MAX
    #elif((LDBL_MANT_DIG == 24) && (LDBL_MAX_EXP == 128) && (BOOST_CSTDFLOAT_HAS_FLOAT32_NATIVE_TYPE == 0))
      #define BOOST_CSTDFLOAT_FLOAT32_NATIVE_TYPE long double
      #undef  BOOST_CSTDFLOAT_MAXIMUM_AVAILABLE_WIDTH
      #define BOOST_CSTDFLOAT_MAXIMUM_AVAILABLE_WIDTH 32
      #undef  BOOST_CSTDFLOAT_HAS_FLOAT32_NATIVE_TYPE
      #define BOOST_CSTDFLOAT_HAS_FLOAT32_NATIVE_TYPE  1
      #define BOOST_FLOAT32_C(x)  (x ## L)
      #define BOOST_FLOAT_32_MIN  LDBL_MIN
      #define BOOST_FLOAT_32_MAX  LDBL_MAX
    #elif((LDBL_MANT_DIG == 53) && (LDBL_MAX_EXP == 1024) && (BOOST_CSTDFLOAT_HAS_FLOAT64_NATIVE_TYPE == 0))
      #define BOOST_CSTDFLOAT_FLOAT64_NATIVE_TYPE long double
      #undef  BOOST_CSTDFLOAT_MAXIMUM_AVAILABLE_WIDTH
      #define BOOST_CSTDFLOAT_MAXIMUM_AVAILABLE_WIDTH 64
      #undef  BOOST_CSTDFLOAT_HAS_FLOAT64_NATIVE_TYPE
      #define BOOST_CSTDFLOAT_HAS_FLOAT64_NATIVE_TYPE  1
      #define BOOST_FLOAT64_C(x)  (x ## L)
      #define BOOST_FLOAT_64_MIN  LDBL_MIN
      #define BOOST_FLOAT_64_MAX  LDBL_MAX
    #elif((LDBL_MANT_DIG == 64) && (LDBL_MAX_EXP == 16384) && (BOOST_CSTDFLOAT_HAS_FLOAT80_NATIVE_TYPE == 0))
      #define BOOST_CSTDFLOAT_FLOAT80_NATIVE_TYPE long double
      #undef  BOOST_CSTDFLOAT_MAXIMUM_AVAILABLE_WIDTH
      #define BOOST_CSTDFLOAT_MAXIMUM_AVAILABLE_WIDTH 80
      #undef  BOOST_CSTDFLOAT_HAS_FLOAT80_NATIVE_TYPE
      #define BOOST_CSTDFLOAT_HAS_FLOAT80_NATIVE_TYPE  1
      #define BOOST_FLOAT80_C(x)  (x ## L)
      #define BOOST_FLOAT_80_MIN  LDBL_MIN
      #define BOOST_FLOAT_80_MAX  LDBL_MAX
    #elif((LDBL_MANT_DIG == 113) && (LDBL_MAX_EXP == 16384) && (BOOST_CSTDFLOAT_HAS_FLOAT128_NATIVE_TYPE == 0))
      #define BOOST_CSTDFLOAT_FLOAT128_NATIVE_TYPE long double
      #undef  BOOST_CSTDFLOAT_MAXIMUM_AVAILABLE_WIDTH
      #define BOOST_CSTDFLOAT_MAXIMUM_AVAILABLE_WIDTH 128
      #undef  BOOST_CSTDFLOAT_HAS_FLOAT128_NATIVE_TYPE
      #define BOOST_CSTDFLOAT_HAS_FLOAT128_NATIVE_TYPE  1
      #define BOOST_FLOAT128_C(x)  (x ## L)
      #define BOOST_FLOAT_128_MIN  LDBL_MIN
      #define BOOST_FLOAT_128_MAX  LDBL_MAX
    #endif
  #endif

  // Check if float_internal128_t from GCC's libquadmath or if (potentially)
  // ICC's /Qlong-double flag is supported.
  // TODO: Should we allow BOOST_MATH_USE_FLOAT128 for ICC?
  // Here, we use the BOOST_MATH_USE_FLOAT128 pre-processor
  // definition from <boost/math/tools/config.hpp>.
  #if (BOOST_CSTDFLOAT_HAS_FLOAT128_NATIVE_TYPE == 0) && defined(BOOST_MATH_USE_FLOAT128)

    namespace boost { namespace cstdfloat { namespace detail {
    #if defined(BOOST_INTEL)
      typedef _Quad      float_internal128_t;

      #define BOOST_CSTDFLOAT_FLOAT128_LDEXP  __ldexpq
      #define BOOST_CSTDFLOAT_FLOAT128_FREXP  __frexpq
      #define BOOST_CSTDFLOAT_FLOAT128_FABS   __fabsq
      #define BOOST_CSTDFLOAT_FLOAT128_FLOOR  __floorq
      #define BOOST_CSTDFLOAT_FLOAT128_CEIL   __ceilq
      #define BOOST_CSTDFLOAT_FLOAT128_SQRT   __sqrtq
      #define BOOST_CSTDFLOAT_FLOAT128_TRUND  __truncq
      #define BOOST_CSTDFLOAT_FLOAT128_EXP    __expq
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
      #define BOOST_CSTDFLOAT_FLOAT128_FMOD   __fmodq
      #define BOOST_CSTDFLOAT_FLOAT128_ATAN2  __atan2q
      #define BOOST_CSTDFLOAT_FLOAT128_LGAMMA __lgammaq
      #define BOOST_CSTDFLOAT_FLOAT128_TGAMMA __tgammaq

    #elif defined(__GNUC__)
      typedef __float128 float_internal128_t;

      #define BOOST_CSTDFLOAT_FLOAT128_LDEXP  ldexpq
      #define BOOST_CSTDFLOAT_FLOAT128_FREXP  frexpq
      #define BOOST_CSTDFLOAT_FLOAT128_FABS   fabsq
      #define BOOST_CSTDFLOAT_FLOAT128_FLOOR  floorq
      #define BOOST_CSTDFLOAT_FLOAT128_CEIL   ceilq
      #define BOOST_CSTDFLOAT_FLOAT128_SQRT   sqrtq
      #define BOOST_CSTDFLOAT_FLOAT128_TRUND  truncq
      #define BOOST_CSTDFLOAT_FLOAT128_EXP    expq
      #define BOOST_CSTDFLOAT_FLOAT128_POW    powq
      #define BOOST_CSTDFLOAT_FLOAT128_LOG    logq
      #define BOOST_CSTDFLOAT_FLOAT128_LOG10  log10q
      #define BOOST_CSTDFLOAT_FLOAT128_SIN    sinq
      #define BOOST_CSTDFLOAT_FLOAT128_COS    cosq
      #define BOOST_CSTDFLOAT_FLOAT128_TAN    tanq
      #define BOOST_CSTDFLOAT_FLOAT128_ASIN   asinq
      #define BOOST_CSTDFLOAT_FLOAT128_ACOS   acosq
      #define BOOST_CSTDFLOAT_FLOAT128_ATAN   atanq
      #define BOOST_CSTDFLOAT_FLOAT128_SINH   sinhq
      #define BOOST_CSTDFLOAT_FLOAT128_COSH   coshq
      #define BOOST_CSTDFLOAT_FLOAT128_TANH   tanhq
      #define BOOST_CSTDFLOAT_FLOAT128_FMOD   fmodq
      #define BOOST_CSTDFLOAT_FLOAT128_ATAN2  atan2q
      #define BOOST_CSTDFLOAT_FLOAT128_LGAMMA lgammaq
      #define BOOST_CSTDFLOAT_FLOAT128_TGAMMA tgammaq

    #else
      #error "Sorry compiler is neither GCC, nor Intel, don't know how to configure this header."
    #endif
    } } } // boost::cstdfloat::detail

    #define BOOST_CSTDFLOAT_FLOAT128_NATIVE_TYPE boost::cstdfloat::detail::float_internal128_t
    #undef  BOOST_CSTDFLOAT_MAXIMUM_AVAILABLE_WIDTH
    #define BOOST_CSTDFLOAT_MAXIMUM_AVAILABLE_WIDTH 128
    #undef  BOOST_CSTDFLOAT_HAS_FLOAT128_NATIVE_TYPE
    #define BOOST_CSTDFLOAT_HAS_FLOAT128_NATIVE_TYPE  1
    #define BOOST_CSTDFLOAT_FLOAT128_MIN  3.36210314311209350626267781732175260e-4932Q
    #define BOOST_CSTDFLOAT_FLOAT128_MAX  1.18973149535723176508575932662800702e+4932Q
    #define BOOST_CSTDFLOAT_FLOAT128_EPS  1.92592994438723585305597794258492732e-0034Q
    #define BOOST_FLOAT128_C(x)  (x ## Q)

    #if !defined(BOOST_CSTDFLOAT_NO_GCC_FLOAT128_LIMITS)

    // For float_internal128_t, implement a specialization of std::numeric_limits<>.

    // Forward declaration of quadruple-precisdion square root function.
    extern "C" boost::cstdfloat::detail::float_internal128_t sqrtq(boost::cstdfloat::detail::float_internal128_t);

    namespace std
    {
      template<>
      class numeric_limits<boost::cstdfloat::detail::float_internal128_t>
      {
      public:
        BOOST_STATIC_CONSTEXPR bool                                           is_specialized           = true;
        static                 boost::cstdfloat::detail::float_internal128_t  (min) () BOOST_NOEXCEPT  { return BOOST_CSTDFLOAT_FLOAT128_MIN; }
        static                 boost::cstdfloat::detail::float_internal128_t  (max) () BOOST_NOEXCEPT  { return BOOST_CSTDFLOAT_FLOAT128_MAX; }
        static                 boost::cstdfloat::detail::float_internal128_t  lowest() BOOST_NOEXCEPT  { return -(max)(); }
        BOOST_STATIC_CONSTEXPR int                                            digits                   = 113;
        BOOST_STATIC_CONSTEXPR int                                            digits10                 = 34;
        BOOST_STATIC_CONSTEXPR int                                            max_digits10             = 36;
        BOOST_STATIC_CONSTEXPR bool                                           is_signed                = true;
        BOOST_STATIC_CONSTEXPR bool                                           is_integer               = false;
        BOOST_STATIC_CONSTEXPR bool                                           is_exact                 = false;
        BOOST_STATIC_CONSTEXPR int                                            radix                    = 2;
        static                 boost::cstdfloat::detail::float_internal128_t  epsilon    ()            { return BOOST_CSTDFLOAT_FLOAT128_EPS; }
        static                 boost::cstdfloat::detail::float_internal128_t  round_error()            { return BOOST_FLOAT128_C(0.5); }
        BOOST_STATIC_CONSTEXPR int                                            min_exponent             = -16381;
        BOOST_STATIC_CONSTEXPR int                                            min_exponent10           = static_cast<int>((min_exponent * 301L) / 1000L);
        BOOST_STATIC_CONSTEXPR int                                            max_exponent             = +16384;
        BOOST_STATIC_CONSTEXPR int                                            max_exponent10           = static_cast<int>((max_exponent * 301L) / 1000L);
        BOOST_STATIC_CONSTEXPR bool                                           has_infinity             = true;
        BOOST_STATIC_CONSTEXPR bool                                           has_quiet_NaN            = true;
        BOOST_STATIC_CONSTEXPR bool                                           has_signaling_NaN        = false;
        BOOST_STATIC_CONSTEXPR float_denorm_style                             has_denorm               = denorm_absent;
        BOOST_STATIC_CONSTEXPR bool                                           has_denorm_loss          = false;
        static                 boost::cstdfloat::detail::float_internal128_t  infinity     ()          { return BOOST_FLOAT128_C(1.0) / BOOST_FLOAT128_C(0.0); }
        static                 boost::cstdfloat::detail::float_internal128_t  quiet_NaN    ()          { return ::sqrtq(BOOST_FLOAT128_C(-1.0)); }
        static                 boost::cstdfloat::detail::float_internal128_t  signaling_NaN()          { return BOOST_FLOAT128_C(0.0); }
        static                 boost::cstdfloat::detail::float_internal128_t  denorm_min   ()          { return BOOST_FLOAT128_C(0.0); }
        BOOST_STATIC_CONSTEXPR bool                                           is_iec559                = true;
        BOOST_STATIC_CONSTEXPR bool                                           is_bounded               = false;
        BOOST_STATIC_CONSTEXPR bool                                           is_modulo                = false;
        BOOST_STATIC_CONSTEXPR bool                                           traps                    = false;
        BOOST_STATIC_CONSTEXPR bool                                           tinyness_before          = false;
        BOOST_STATIC_CONSTEXPR float_round_style                              round_style              = round_to_nearest;
      };
    }
    #endif

    #if !defined(BOOST_CSTDFLOAT_NO_GCC_FLOAT128_CMATH)

      // For float_internal128_t, implement <math.h> functions in the namespace
      // boost::cstdfloat::detail. Subsequently *use* these in the global namespace.

      // Begin with some forward function declarations.

      // Forward declaration of quad string print function.
      extern "C" int quadmath_snprintf(char *str, size_t size, const char *format, ...);
      extern "C" boost::cstdfloat::detail::float_internal128_t strtoflt128(const char*, char **);

      // Forward declarations of quad elementary functions.
      extern "C" boost::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_LDEXP (boost::cstdfloat::detail::float_internal128_t, int);
      extern "C" boost::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_FREXP (boost::cstdfloat::detail::float_internal128_t, int*);
      extern "C" boost::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_FABS  (boost::cstdfloat::detail::float_internal128_t);
      extern "C" boost::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_FLOOR (boost::cstdfloat::detail::float_internal128_t);
      extern "C" boost::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_CEIL  (boost::cstdfloat::detail::float_internal128_t);
      extern "C" boost::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_SQRT  (boost::cstdfloat::detail::float_internal128_t);
      extern "C" boost::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_TRUND (boost::cstdfloat::detail::float_internal128_t);
      extern "C" boost::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_EXP   (boost::cstdfloat::detail::float_internal128_t);
      extern "C" boost::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_POW   (boost::cstdfloat::detail::float_internal128_t, boost::cstdfloat::detail::float_internal128_t);
      extern "C" boost::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_LOG   (boost::cstdfloat::detail::float_internal128_t);
      extern "C" boost::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_LOG10 (boost::cstdfloat::detail::float_internal128_t);
      extern "C" boost::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_SIN   (boost::cstdfloat::detail::float_internal128_t);
      extern "C" boost::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_COS   (boost::cstdfloat::detail::float_internal128_t);
      extern "C" boost::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_TAN   (boost::cstdfloat::detail::float_internal128_t);
      extern "C" boost::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_ASIN  (boost::cstdfloat::detail::float_internal128_t);
      extern "C" boost::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_ACOS  (boost::cstdfloat::detail::float_internal128_t);
      extern "C" boost::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_ATAN  (boost::cstdfloat::detail::float_internal128_t);
      extern "C" boost::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_SINH  (boost::cstdfloat::detail::float_internal128_t);
      extern "C" boost::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_COSH  (boost::cstdfloat::detail::float_internal128_t);
      extern "C" boost::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_TANH  (boost::cstdfloat::detail::float_internal128_t);
      extern "C" boost::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_FMOD  (boost::cstdfloat::detail::float_internal128_t, boost::cstdfloat::detail::float_internal128_t);
      extern "C" boost::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_ATAN2 (boost::cstdfloat::detail::float_internal128_t, boost::cstdfloat::detail::float_internal128_t);
      extern "C" boost::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_LGAMMA(boost::cstdfloat::detail::float_internal128_t);
      extern "C" boost::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_TGAMMA(boost::cstdfloat::detail::float_internal128_t);

      // Put the float_internal128_t <math.h> functions in the namespace
      // boost::cstdfloat::detail.

      namespace boost { namespace cstdfloat { namespace detail {
      inline   boost::cstdfloat::detail::float_internal128_t ldexp (boost::cstdfloat::detail::float_internal128_t x, int n)                                           { return ::BOOST_CSTDFLOAT_FLOAT128_LDEXP (x, n); }
      inline   boost::cstdfloat::detail::float_internal128_t frexp (boost::cstdfloat::detail::float_internal128_t x, int* pn)                                         { return ::BOOST_CSTDFLOAT_FLOAT128_FREXP (x, pn); }
      inline   boost::cstdfloat::detail::float_internal128_t fabs  (boost::cstdfloat::detail::float_internal128_t x)                                                  { return ::BOOST_CSTDFLOAT_FLOAT128_FABS  (x); }
      inline   boost::cstdfloat::detail::float_internal128_t floor (boost::cstdfloat::detail::float_internal128_t x)                                                  { return ::BOOST_CSTDFLOAT_FLOAT128_FLOOR (x); }
      inline   boost::cstdfloat::detail::float_internal128_t ceil  (boost::cstdfloat::detail::float_internal128_t x)                                                  { return ::BOOST_CSTDFLOAT_FLOAT128_CEIL  (x); }
      inline   boost::cstdfloat::detail::float_internal128_t sqrt  (boost::cstdfloat::detail::float_internal128_t x)                                                  { return ::BOOST_CSTDFLOAT_FLOAT128_SQRT  (x); }
      inline   boost::cstdfloat::detail::float_internal128_t trunc (boost::cstdfloat::detail::float_internal128_t x)                                                  { return ::BOOST_CSTDFLOAT_FLOAT128_TRUND (x); }
      inline   boost::cstdfloat::detail::float_internal128_t exp   (boost::cstdfloat::detail::float_internal128_t x)                                                  { return ::BOOST_CSTDFLOAT_FLOAT128_EXP   (x); }
      inline   boost::cstdfloat::detail::float_internal128_t pow   (boost::cstdfloat::detail::float_internal128_t x, boost::cstdfloat::detail::float_internal128_t a) { return ::BOOST_CSTDFLOAT_FLOAT128_POW   (x, a); }
      inline   boost::cstdfloat::detail::float_internal128_t log   (boost::cstdfloat::detail::float_internal128_t x)                                                  { return ::BOOST_CSTDFLOAT_FLOAT128_LOG   (x); }
      inline   boost::cstdfloat::detail::float_internal128_t log10 (boost::cstdfloat::detail::float_internal128_t x)                                                  { return ::BOOST_CSTDFLOAT_FLOAT128_LOG10 (x); }
      inline   boost::cstdfloat::detail::float_internal128_t sin   (boost::cstdfloat::detail::float_internal128_t x)                                                  { return ::BOOST_CSTDFLOAT_FLOAT128_SIN   (x); }
      inline   boost::cstdfloat::detail::float_internal128_t cos   (boost::cstdfloat::detail::float_internal128_t x)                                                  { return ::BOOST_CSTDFLOAT_FLOAT128_COS   (x); }
      inline   boost::cstdfloat::detail::float_internal128_t tan   (boost::cstdfloat::detail::float_internal128_t x)                                                  { return ::BOOST_CSTDFLOAT_FLOAT128_TAN   (x); }
      inline   boost::cstdfloat::detail::float_internal128_t asin  (boost::cstdfloat::detail::float_internal128_t x)                                                  { return ::BOOST_CSTDFLOAT_FLOAT128_ASIN  (x); }
      inline   boost::cstdfloat::detail::float_internal128_t acos  (boost::cstdfloat::detail::float_internal128_t x)                                                  { return ::BOOST_CSTDFLOAT_FLOAT128_ACOS  (x); }
      inline   boost::cstdfloat::detail::float_internal128_t atan  (boost::cstdfloat::detail::float_internal128_t x)                                                  { return ::BOOST_CSTDFLOAT_FLOAT128_ATAN  (x); }
      inline   boost::cstdfloat::detail::float_internal128_t sinh  (boost::cstdfloat::detail::float_internal128_t x)                                                  { return ::BOOST_CSTDFLOAT_FLOAT128_SINH  (x); }
      inline   boost::cstdfloat::detail::float_internal128_t cosh  (boost::cstdfloat::detail::float_internal128_t x)                                                  { return ::BOOST_CSTDFLOAT_FLOAT128_COSH  (x); }
      inline   boost::cstdfloat::detail::float_internal128_t tanh  (boost::cstdfloat::detail::float_internal128_t x)                                                  { return ::BOOST_CSTDFLOAT_FLOAT128_TANH  (x); }
      inline   boost::cstdfloat::detail::float_internal128_t fmod  (boost::cstdfloat::detail::float_internal128_t a, boost::cstdfloat::detail::float_internal128_t b) { return ::BOOST_CSTDFLOAT_FLOAT128_FMOD  (a, b); }
      inline   boost::cstdfloat::detail::float_internal128_t atan2 (boost::cstdfloat::detail::float_internal128_t y, boost::cstdfloat::detail::float_internal128_t x) { return ::BOOST_CSTDFLOAT_FLOAT128_ATAN2 (y, x); }
      inline   boost::cstdfloat::detail::float_internal128_t lgamma(boost::cstdfloat::detail::float_internal128_t x)                                                  { return ::BOOST_CSTDFLOAT_FLOAT128_LGAMMA(x); }
      inline   boost::cstdfloat::detail::float_internal128_t tgamma(boost::cstdfloat::detail::float_internal128_t x)                                                  { return ::BOOST_CSTDFLOAT_FLOAT128_TGAMMA(x); }
      } } }  // boost::cstdfloat::detail

      #undef BOOST_CSTDFLOAT_FLOAT128_LDEXP
      #undef BOOST_CSTDFLOAT_FLOAT128_FREXP
      #undef BOOST_CSTDFLOAT_FLOAT128_FABS
      #undef BOOST_CSTDFLOAT_FLOAT128_FLOOR
      #undef BOOST_CSTDFLOAT_FLOAT128_CEIL
      #undef BOOST_CSTDFLOAT_FLOAT128_SQRT
      #undef BOOST_CSTDFLOAT_FLOAT128_TRUND
      #undef BOOST_CSTDFLOAT_FLOAT128_EXP
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
      #undef BOOST_CSTDFLOAT_FLOAT128_FMOD
      #undef BOOST_CSTDFLOAT_FLOAT128_ATAN2
      #undef BOOST_CSTDFLOAT_FLOAT128_LGAMMA
      #undef BOOST_CSTDFLOAT_FLOAT128_TGAMMA

      using boost::cstdfloat::detail::ldexp;
      using boost::cstdfloat::detail::frexp;
      using boost::cstdfloat::detail::fabs;
      using boost::cstdfloat::detail::floor;
      using boost::cstdfloat::detail::ceil;
      using boost::cstdfloat::detail::sqrt;
      using boost::cstdfloat::detail::trunc;
      using boost::cstdfloat::detail::exp;
      using boost::cstdfloat::detail::pow;
      using boost::cstdfloat::detail::log;
      using boost::cstdfloat::detail::log10;
      using boost::cstdfloat::detail::sin;
      using boost::cstdfloat::detail::cos;
      using boost::cstdfloat::detail::tan;
      using boost::cstdfloat::detail::asin;
      using boost::cstdfloat::detail::acos;
      using boost::cstdfloat::detail::atan;
      using boost::cstdfloat::detail::sinh;
      using boost::cstdfloat::detail::cosh;
      using boost::cstdfloat::detail::tanh;
      using boost::cstdfloat::detail::fmod;
      using boost::cstdfloat::detail::atan2;
      using boost::cstdfloat::detail::lgamma;
      using boost::cstdfloat::detail::tgamma;

      // For float_internal128_t, implement <cmath> functions in the std namespace.
      namespace std
      {
        inline boost::cstdfloat::detail::float_internal128_t ldexp (boost::cstdfloat::detail::float_internal128_t x, int n)                                           { return boost::cstdfloat::detail::ldexp (x, n); }
        inline boost::cstdfloat::detail::float_internal128_t frexp (boost::cstdfloat::detail::float_internal128_t x, int* pn)                                         { return boost::cstdfloat::detail::frexp (x, pn); }
        inline boost::cstdfloat::detail::float_internal128_t fabs  (boost::cstdfloat::detail::float_internal128_t x)                                                  { return boost::cstdfloat::detail::fabs  (x); }
        inline boost::cstdfloat::detail::float_internal128_t floor (boost::cstdfloat::detail::float_internal128_t x)                                                  { return boost::cstdfloat::detail::floor (x); }
        inline boost::cstdfloat::detail::float_internal128_t ceil  (boost::cstdfloat::detail::float_internal128_t x)                                                  { return boost::cstdfloat::detail::ceil  (x); }
        inline boost::cstdfloat::detail::float_internal128_t sqrt  (boost::cstdfloat::detail::float_internal128_t x)                                                  { return boost::cstdfloat::detail::sqrt  (x); }
        inline boost::cstdfloat::detail::float_internal128_t trunc (boost::cstdfloat::detail::float_internal128_t x)                                                  { return boost::cstdfloat::detail::trunc (x); }
        inline boost::cstdfloat::detail::float_internal128_t exp   (boost::cstdfloat::detail::float_internal128_t x)                                                  { return boost::cstdfloat::detail::exp   (x); }
        inline boost::cstdfloat::detail::float_internal128_t pow   (boost::cstdfloat::detail::float_internal128_t x, boost::cstdfloat::detail::float_internal128_t a) { return boost::cstdfloat::detail::pow   (x, a); }
        inline boost::cstdfloat::detail::float_internal128_t log   (boost::cstdfloat::detail::float_internal128_t x)                                                  { return boost::cstdfloat::detail::log   (x); }
        inline boost::cstdfloat::detail::float_internal128_t log10 (boost::cstdfloat::detail::float_internal128_t x)                                                  { return boost::cstdfloat::detail::log10 (x); }
        inline boost::cstdfloat::detail::float_internal128_t sin   (boost::cstdfloat::detail::float_internal128_t x)                                                  { return boost::cstdfloat::detail::sin   (x); }
        inline boost::cstdfloat::detail::float_internal128_t cos   (boost::cstdfloat::detail::float_internal128_t x)                                                  { return boost::cstdfloat::detail::cos   (x); }
        inline boost::cstdfloat::detail::float_internal128_t tan   (boost::cstdfloat::detail::float_internal128_t x)                                                  { return boost::cstdfloat::detail::tan   (x); }
        inline boost::cstdfloat::detail::float_internal128_t asin  (boost::cstdfloat::detail::float_internal128_t x)                                                  { return boost::cstdfloat::detail::asin  (x); }
        inline boost::cstdfloat::detail::float_internal128_t acos  (boost::cstdfloat::detail::float_internal128_t x)                                                  { return boost::cstdfloat::detail::acos  (x); }
        inline boost::cstdfloat::detail::float_internal128_t atan  (boost::cstdfloat::detail::float_internal128_t x)                                                  { return boost::cstdfloat::detail::atan  (x); }
        inline boost::cstdfloat::detail::float_internal128_t sinh  (boost::cstdfloat::detail::float_internal128_t x)                                                  { return boost::cstdfloat::detail::sinh  (x); }
        inline boost::cstdfloat::detail::float_internal128_t cosh  (boost::cstdfloat::detail::float_internal128_t x)                                                  { return boost::cstdfloat::detail::cosh  (x); }
        inline boost::cstdfloat::detail::float_internal128_t tanh  (boost::cstdfloat::detail::float_internal128_t x)                                                  { return boost::cstdfloat::detail::tanh  (x); }
        inline boost::cstdfloat::detail::float_internal128_t fmod  (boost::cstdfloat::detail::float_internal128_t a, boost::cstdfloat::detail::float_internal128_t b) { return boost::cstdfloat::detail::fmod  (a, b); }
        inline boost::cstdfloat::detail::float_internal128_t atan2 (boost::cstdfloat::detail::float_internal128_t y, boost::cstdfloat::detail::float_internal128_t x) { return boost::cstdfloat::detail::atan2 (y, x); }
        inline boost::cstdfloat::detail::float_internal128_t lgamma(boost::cstdfloat::detail::float_internal128_t x )                                                 { return boost::cstdfloat::detail::lgamma(x); }
        inline boost::cstdfloat::detail::float_internal128_t tgamma(boost::cstdfloat::detail::float_internal128_t x )                                                 { return boost::cstdfloat::detail::tgamma(x); }
      }
    #endif

    #if !defined(BOOST_CSTDFLOAT_NO_GCC_FLOAT128_IOSTREAM)

      // For float_internal128_t, implement I/O stream operations.

      #include <istream>
      #include <ostream>
      #include <string>

      inline std::ostream& operator<<(std::ostream& os, const boost::cstdfloat::detail::float_internal128_t& x)
      {
        #if defined(__GNUC__)

        char my_buffer[128U];

        const int my_prec   = static_cast<int>(os.precision());
        const int my_digits = ((my_prec == 0) ? 36 : my_prec);

        const std::ios_base::fmtflags my_flags  = os.flags();

        char my_format_string[8U];

        std::size_t my_format_string_index = 0U;

        my_format_string[my_format_string_index] = '%';
        ++my_format_string_index;

        if(my_flags & std::ios_base::showpos)    { my_format_string[my_format_string_index] = '+'; ++my_format_string_index; }
        if(my_flags & std::ios_base::showpoint)  { my_format_string[my_format_string_index] = '#'; ++my_format_string_index; }

        my_format_string[my_format_string_index + 0U] = '.';
        my_format_string[my_format_string_index + 1U] = '*';
        my_format_string[my_format_string_index + 2U] = 'Q';

        my_format_string_index += 3U;

        char the_notation_char;

        if     (my_flags & std::ios_base::scientific) { the_notation_char = 'e'; }
        else if(my_flags & std::ios_base::fixed)      { the_notation_char = 'f'; }
        else                                          { the_notation_char = 'g'; }

        my_format_string[my_format_string_index + 0U] = the_notation_char;
        my_format_string[my_format_string_index + 1U] = 0;

        const int v = ::quadmath_snprintf(my_buffer,
                                          static_cast<int>(sizeof(my_buffer)),
                                          my_format_string,
                                          my_digits,
                                          x);

        if(v < 0) { BOOST_THROW_EXCEPTION(std::runtime_error("Formatting of boost::float128_t failed.")); }

        if(v >= static_cast<int>(sizeof(my_buffer) - 1U))
        {
          // Evidently there is a really long floating-point string here.
          // So we finally have to resort to dynamic memory allocation.
          char* my_buffer2 = static_cast<char*>(0U);

          try
          {
            my_buffer2 = new char[v + 3];
          }
          catch(const std::bad_alloc&)
          {
            BOOST_THROW_EXCEPTION(std::runtime_error("Formatting of boost::float128_t failed while allocating memory."));
          }

          const int v2 = ::quadmath_snprintf(my_buffer2,
                                             v + 3,
                                             my_format_string,
                                             my_digits,
                                             x);

          if(v2 >= v + 3)
          {
            BOOST_THROW_EXCEPTION(std::runtime_error("Formatting of boost::float128_t failed."));
          }

          os << my_buffer2;

          delete [] my_buffer2;

          return os;
        }
        else
        {
          return (os << my_buffer);
        }

        #elif defined(BOOST_INTEL)

        BOOST_THROW_EXCEPTION(std::runtime_error("<boost/cstdfloat.hpp> does not yet support output stream operation for ICC quadruple-precision type"));
        return os;

        #endif
      }

      inline std::istream& operator>>(std::istream& is, boost::cstdfloat::detail::float_internal128_t& x)
      {
        #if defined(__GNUC__)

        std::string str;

        is >> str;

        char* p_end;

        x = strtoflt128(str.c_str(), &p_end);

        if(static_cast<std::ptrdiff_t>(p_end - str.c_str()) != static_cast<std::ptrdiff_t>(str.length()))
        {
          BOOST_THROW_EXCEPTION(std::runtime_error("Unable to interpret input string as a boost::float128_t"));
        }

        return is;

        #elif defined(BOOST_INTEL)

        BOOST_THROW_EXCEPTION(std::runtime_error("<boost/cstdfloat.hpp> does not yet support input stream operation for ICC quadruple-precision type"));
        return is;

        #endif
      }

    #endif
  #endif

  // This is the end of the preamble. Now we use the results
  // of the queries that have been obtained in the preamble.

  // Make sure that the compiler has any floating-point type whatsoever.
  #if (   (BOOST_CSTDFLOAT_HAS_FLOAT16_NATIVE_TYPE  == 0)  \
       && (BOOST_CSTDFLOAT_HAS_FLOAT32_NATIVE_TYPE  == 0)  \
       && (BOOST_CSTDFLOAT_HAS_FLOAT64_NATIVE_TYPE  == 0)  \
       && (BOOST_CSTDFLOAT_HAS_FLOAT80_NATIVE_TYPE  == 0)  \
       && (BOOST_CSTDFLOAT_HAS_FLOAT128_NATIVE_TYPE == 0))
    #error The compiler does not support any of the floating-point types required for <boost/cstdfloat.hpp>.
  #endif

  // The following section contains the various min/max macros
  // for the *leastN and *fastN types.

  #if(BOOST_CSTDFLOAT_HAS_FLOAT16_NATIVE_TYPE == 1)
    #define BOOST_FLOAT_FAST16_MIN   BOOST_FLOAT_16_MIN
    #define BOOST_FLOAT_LEAST16_MIN  BOOST_FLOAT_16_MIN
    #define BOOST_FLOAT_FAST16_MAX   BOOST_FLOAT_16_MAX
    #define BOOST_FLOAT_LEAST16_MAX  BOOST_FLOAT_16_MAX
  #endif

  #if(BOOST_CSTDFLOAT_HAS_FLOAT32_NATIVE_TYPE == 1)
    #define BOOST_FLOAT_FAST32_MIN   BOOST_FLOAT_32_MIN
    #define BOOST_FLOAT_LEAST32_MIN  BOOST_FLOAT_32_MIN
    #define BOOST_FLOAT_FAST32_MAX   BOOST_FLOAT_32_MAX
    #define BOOST_FLOAT_LEAST32_MAX  BOOST_FLOAT_32_MAX
  #endif

  #if(BOOST_CSTDFLOAT_HAS_FLOAT64_NATIVE_TYPE == 1)
    #define BOOST_FLOAT_FAST64_MIN   BOOST_FLOAT_64_MIN
    #define BOOST_FLOAT_LEAST64_MIN  BOOST_FLOAT_64_MIN
    #define BOOST_FLOAT_FAST64_MAX   BOOST_FLOAT_64_MAX
    #define BOOST_FLOAT_LEAST64_MAX  BOOST_FLOAT_64_MAX
  #endif

  #if(BOOST_CSTDFLOAT_HAS_FLOAT80_NATIVE_TYPE == 1)
    #define BOOST_FLOAT_FAST80_MIN   BOOST_FLOAT_80_MIN
    #define BOOST_FLOAT_LEAST80_MIN  BOOST_FLOAT_80_MIN
    #define BOOST_FLOAT_FAST80_MAX   BOOST_FLOAT_80_MAX
    #define BOOST_FLOAT_LEAST80_MAX  BOOST_FLOAT_80_MAX
  #endif

  #define BOOST_NO_FLOAT128_T

  #if(BOOST_CSTDFLOAT_HAS_FLOAT128_NATIVE_TYPE == 1)
    #undef BOOST_NO_FLOAT128_T
    #define BOOST_FLOAT_FAST128_MIN   BOOST_FLOAT_128_MIN
    #define BOOST_FLOAT_LEAST128_MIN  BOOST_FLOAT_128_MIN
    #define BOOST_FLOAT_FAST128_MAX   BOOST_FLOAT_128_MAX
    #define BOOST_FLOAT_LEAST128_MAX  BOOST_FLOAT_128_MAX
  #endif

  // The following section contains the various min/max macros
  // for the *floatmax types.

  #if  (BOOST_CSTDFLOAT_MAXIMUM_AVAILABLE_WIDTH == 16)
    #define BOOST_FLOATMAX_C(x) BOOST_FLOAT16_C(x)
    #define BOOST_FLOATMAX_MIN  BOOST_FLOAT_16_MIN
    #define BOOST_FLOATMAX_MAX  BOOST_FLOAT_16_MAX
  #elif(BOOST_CSTDFLOAT_MAXIMUM_AVAILABLE_WIDTH == 32)
    #define BOOST_FLOATMAX_C(x) BOOST_FLOAT32_C(x)
    #define BOOST_FLOATMAX_MIN  BOOST_FLOAT_32_MIN
    #define BOOST_FLOATMAX_MAX  BOOST_FLOAT_32_MAX
  #elif(BOOST_CSTDFLOAT_MAXIMUM_AVAILABLE_WIDTH == 64)
    #define BOOST_FLOATMAX_C(x) BOOST_FLOAT64_C(x)
    #define BOOST_FLOATMAX_MIN  BOOST_FLOAT_64_MIN
    #define BOOST_FLOATMAX_MAX  BOOST_FLOAT_64_MAX
  #elif(BOOST_CSTDFLOAT_MAXIMUM_AVAILABLE_WIDTH == 80)
    #define BOOST_FLOATMAX_C(x) BOOST_FLOAT80_C(x)
    #define BOOST_FLOATMAX_MIN  BOOST_FLOAT_80_MIN
    #define BOOST_FLOATMAX_MAX  BOOST_FLOAT_80_MAX
  #elif(BOOST_CSTDFLOAT_MAXIMUM_AVAILABLE_WIDTH == 128)
    #define BOOST_FLOATMAX_C(x) BOOST_FLOAT128_C(x)
    #define BOOST_FLOATMAX_MIN  BOOST_FLOAT_128_MIN
    #define BOOST_FLOATMAX_MAX  BOOST_FLOAT_128_MAX
  #else
    #error The maximum available floating-point width for <boost/cstdfloat.hpp> is undefined.
  #endif

  // Here, we define the floating-point typedefs having specified widths.
  // The types are defined in the namespace boost.

  // For simplicity, the least and fast types are type defined identically
  // as the corresponding fixed-width type. This behavior can, however,
  // be modified in order to be optimized for a given compiler implementation.

  // In addition, a clear assessment of IEEE-754 comformance is carried out
  // using compile-time assertion.

  namespace boost
  {
    #if(BOOST_CSTDFLOAT_HAS_FLOAT16_NATIVE_TYPE == 1)
      typedef BOOST_CSTDFLOAT_FLOAT16_NATIVE_TYPE float16_t;
      typedef boost::float16_t float_fast16_t;
      typedef boost::float16_t float_least16_t;

      BOOST_STATIC_ASSERT(std::numeric_limits<boost::float16_t>::is_iec559    == true);
      BOOST_STATIC_ASSERT(std::numeric_limits<boost::float16_t>::radix        ==    2);
      BOOST_STATIC_ASSERT(std::numeric_limits<boost::float16_t>::digits       ==   11);
      BOOST_STATIC_ASSERT(std::numeric_limits<boost::float16_t>::max_exponent ==   16);
    #endif

    #if(BOOST_CSTDFLOAT_HAS_FLOAT32_NATIVE_TYPE == 1)
      typedef BOOST_CSTDFLOAT_FLOAT32_NATIVE_TYPE float32_t;
      typedef boost::float32_t float_fast32_t;
      typedef boost::float32_t float_least32_t;

      BOOST_STATIC_ASSERT(std::numeric_limits<boost::float32_t>::is_iec559    == true);
      BOOST_STATIC_ASSERT(std::numeric_limits<boost::float32_t>::radix        ==    2);
      BOOST_STATIC_ASSERT(std::numeric_limits<boost::float32_t>::digits       ==   24);
      BOOST_STATIC_ASSERT(std::numeric_limits<boost::float32_t>::max_exponent ==  128);
    #endif

    #if(BOOST_CSTDFLOAT_HAS_FLOAT64_NATIVE_TYPE == 1)
      typedef BOOST_CSTDFLOAT_FLOAT64_NATIVE_TYPE float64_t;
      typedef boost::float64_t float_fast64_t;
      typedef boost::float64_t float_least64_t;

      BOOST_STATIC_ASSERT(std::numeric_limits<boost::float64_t>::is_iec559    == true);
      BOOST_STATIC_ASSERT(std::numeric_limits<boost::float64_t>::radix        ==    2);
      BOOST_STATIC_ASSERT(std::numeric_limits<boost::float64_t>::digits       ==   53);
      BOOST_STATIC_ASSERT(std::numeric_limits<boost::float64_t>::max_exponent == 1024);
    #endif

    #if(BOOST_CSTDFLOAT_HAS_FLOAT80_NATIVE_TYPE == 1)
      typedef BOOST_CSTDFLOAT_FLOAT80_NATIVE_TYPE float80_t;
      typedef boost::float80_t float_fast80_t;
      typedef boost::float80_t float_least80_t;

      BOOST_STATIC_ASSERT(std::numeric_limits<boost::float80_t>::is_iec559    ==  true);
      BOOST_STATIC_ASSERT(std::numeric_limits<boost::float80_t>::radix        ==     2);
      BOOST_STATIC_ASSERT(std::numeric_limits<boost::float80_t>::digits       ==    64);
      BOOST_STATIC_ASSERT(std::numeric_limits<boost::float80_t>::max_exponent == 16384);
    #endif

    #if(BOOST_CSTDFLOAT_HAS_FLOAT128_NATIVE_TYPE == 1)
      typedef BOOST_CSTDFLOAT_FLOAT128_NATIVE_TYPE float128_t;
      typedef boost::float128_t float_fast128_t;
      typedef boost::float128_t float_least128_t;

      BOOST_STATIC_ASSERT(std::numeric_limits<boost::float128_t>::is_iec559    ==  true);
      BOOST_STATIC_ASSERT(std::numeric_limits<boost::float128_t>::radix        ==     2);
      BOOST_STATIC_ASSERT(std::numeric_limits<boost::float128_t>::digits       ==   113);
      BOOST_STATIC_ASSERT(std::numeric_limits<boost::float128_t>::max_exponent == 16384);
    #endif

    #if  (BOOST_CSTDFLOAT_MAXIMUM_AVAILABLE_WIDTH ==  16)
      typedef boost::float16_t  floatmax_t;
    #elif(BOOST_CSTDFLOAT_MAXIMUM_AVAILABLE_WIDTH ==  32)
      typedef boost::float32_t  floatmax_t;
    #elif(BOOST_CSTDFLOAT_MAXIMUM_AVAILABLE_WIDTH ==  64)
      typedef boost::float64_t  floatmax_t;
    #elif(BOOST_CSTDFLOAT_MAXIMUM_AVAILABLE_WIDTH ==  80)
      typedef boost::float80_t  floatmax_t;
    #elif(BOOST_CSTDFLOAT_MAXIMUM_AVAILABLE_WIDTH == 128)
      typedef boost::float128_t floatmax_t;
    #else
      #error The maximum available floating-point width for <boost/cstdfloat.hpp> is undefined.
    #endif
  }
  // namespace boost

#endif // _BOOST_CSTDFLOAT_2014_01_09_HPP_
