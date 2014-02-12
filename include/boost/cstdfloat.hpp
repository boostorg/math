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

  // Check if float_internal128_t is supported (i.e., __float128 from
  // GCC's quadmath or _Quad from ICC's /Qlong-double flag).
  // Here, we use the BOOST_MATH_USE_FLOAT128 pre-processor definition
  // from <boost/math/tools/config.hpp> to query for libquadmath.

  // What now follows are some rather long sections, each one of which
  // optionally implements part of the C++ standard library for float128_t.
  // These parts of the C++ standard library include <limits>, <cmath>,
  // I/O stream support, and <complex> for boost::float128_t.

  #if (BOOST_CSTDFLOAT_HAS_FLOAT128_NATIVE_TYPE == 0) && defined(BOOST_MATH_USE_FLOAT128) && !defined(BOOST_CSTDFLOAT_NO_LIBQUADMATH_SUPPORT)

    namespace boost { namespace cstdfloat { namespace detail {
    #if defined(BOOST_INTEL)
      typedef _Quad      float_internal128_t;

      #define BOOST_CSTDFLOAT_FLOAT128_LDEXP  __ldexpq
      #define BOOST_CSTDFLOAT_FLOAT128_FREXP  __frexpq
      #define BOOST_CSTDFLOAT_FLOAT128_FABS   __fabsq
      #define BOOST_CSTDFLOAT_FLOAT128_FLOOR  __floorq
      #define BOOST_CSTDFLOAT_FLOAT128_CEIL   __ceilq
      #define BOOST_CSTDFLOAT_FLOAT128_SQRT   __sqrtq
      #define BOOST_CSTDFLOAT_FLOAT128_TRUNC  __truncq
      #define BOOST_CSTDFLOAT_FLOAT128_POW    __powq
      #define BOOST_CSTDFLOAT_FLOAT128_EXP    __expq_patch
      #define BOOST_CSTDFLOAT_FLOAT128_LOG    __logq
      #define BOOST_CSTDFLOAT_FLOAT128_LOG10  __log10q
      #define BOOST_CSTDFLOAT_FLOAT128_SIN    __sinq
      #define BOOST_CSTDFLOAT_FLOAT128_COS    __cosq
      #define BOOST_CSTDFLOAT_FLOAT128_TAN    __tanq
      #define BOOST_CSTDFLOAT_FLOAT128_ASIN   __asinq
      #define BOOST_CSTDFLOAT_FLOAT128_ACOS   __acosq
      #define BOOST_CSTDFLOAT_FLOAT128_ATAN   __atanq
      #define BOOST_CSTDFLOAT_FLOAT128_SINH   __sinhq_patch
      #define BOOST_CSTDFLOAT_FLOAT128_COSH   __coshq_patch
      #define BOOST_CSTDFLOAT_FLOAT128_TANH   __tanhq_patch
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
      #define BOOST_CSTDFLOAT_FLOAT128_TRUNC  truncq
      #define BOOST_CSTDFLOAT_FLOAT128_POW    powq
      #define BOOST_CSTDFLOAT_FLOAT128_EXP    expq_patch
      #define BOOST_CSTDFLOAT_FLOAT128_LOG    logq
      #define BOOST_CSTDFLOAT_FLOAT128_LOG10  log10q
      #define BOOST_CSTDFLOAT_FLOAT128_SIN    sinq
      #define BOOST_CSTDFLOAT_FLOAT128_COS    cosq
      #define BOOST_CSTDFLOAT_FLOAT128_TAN    tanq
      #define BOOST_CSTDFLOAT_FLOAT128_ASIN   asinq
      #define BOOST_CSTDFLOAT_FLOAT128_ACOS   acosq
      #define BOOST_CSTDFLOAT_FLOAT128_ATAN   atanq
      #define BOOST_CSTDFLOAT_FLOAT128_SINH   sinhq_patch
      #define BOOST_CSTDFLOAT_FLOAT128_COSH   coshq_patch
      #define BOOST_CSTDFLOAT_FLOAT128_TANH   tanhq_patch
      #define BOOST_CSTDFLOAT_FLOAT128_FMOD   fmodq
      #define BOOST_CSTDFLOAT_FLOAT128_ATAN2  atan2q
      #define BOOST_CSTDFLOAT_FLOAT128_LGAMMA lgammaq
      #define BOOST_CSTDFLOAT_FLOAT128_TGAMMA tgammaq

    #else
      #error "Sorry, the compiler is neither GCC, nor Intel, I don't know how to configure <boost/cstdfloat.hpp>."
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

    #if !defined(BOOST_CSTDFLOAT_NO_LIBQUADMATH_LIMITS)

    // Implement a specialization of std::numeric_limits<> for quadruple-precision.

    // Forward declaration of quadruple-precision square root function.
    extern "C" boost::cstdfloat::detail::float_internal128_t sqrtq(boost::cstdfloat::detail::float_internal128_t) throw();

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
    #endif // Not BOOST_CSTDFLOAT_NO_LIBQUADMATH_LIMITS (i.e., has libquadmath <limits> support)

    #if !defined(BOOST_CSTDFLOAT_NO_LIBQUADMATH_CMATH)

      #include <boost/cstdint.hpp>

      // Implement quadruple-precision <math.h> functions in the namespace
      // boost::cstdfloat::detail. Subsequently *use* these in the global
      // namespace and in the std namespace.

      // Begin with some forward function declarations.

      // Forward declarations of quadruple-precision string functions.
      extern "C" int quadmath_snprintf(char *str, size_t size, const char *format, ...) throw();
      extern "C" boost::cstdfloat::detail::float_internal128_t strtoflt128(const char*, char **) throw();

      // Forward declarations of quadruple-precision elementary functions.
      extern "C" boost::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_LDEXP (boost::cstdfloat::detail::float_internal128_t, int)  throw();
      extern "C" boost::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_FREXP (boost::cstdfloat::detail::float_internal128_t, int*) throw();
      extern "C" boost::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_FABS  (boost::cstdfloat::detail::float_internal128_t) throw();
      extern "C" boost::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_FLOOR (boost::cstdfloat::detail::float_internal128_t) throw();
      extern "C" boost::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_CEIL  (boost::cstdfloat::detail::float_internal128_t) throw();
      extern "C" boost::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_SQRT  (boost::cstdfloat::detail::float_internal128_t) throw();
      extern "C" boost::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_TRUNC (boost::cstdfloat::detail::float_internal128_t) throw();
      extern "C" boost::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_POW   (boost::cstdfloat::detail::float_internal128_t, boost::cstdfloat::detail::float_internal128_t) throw();
      inline     boost::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_EXP   (boost::cstdfloat::detail::float_internal128_t x)
      {
        // Patch the expq() function for a subset of broken GCC compilers
        // like GCC 4.7, 4.8 on MinGW.
        typedef boost::cstdfloat::detail::float_internal128_t float_type;

        // Use an order-36 polynomial approximation of the exponential function
        // in the range of (-ln2 < x < ln2). Scale the argument to this range
        // and multiply the result by 2^n accordingly.

        // Generate the polynomial coefficients with Mathematica(R).

        // Table[{x, Exp[x] - 1}, {x, -Log[2], Log[2], 1/180}]
        // N[%, 120]
        // Fit[%, {x, x^2, x^3, x^4, x^5, x^6, x^7, x^8, x^9, x^10, x^11, x^12,
        //         x^13, x^14, x^15, x^16, x^17, x^18, x^19, x^20, x^21, x^22,
        //         x^23, x^24, x^25, x^26, x^27, x^28, x^29, x^30, x^31, x^32,
        //         x^33, x^34, x^35, x^36}, x]

        // Scale the argument x to the range (-ln2 < x < ln2).
        BOOST_CONSTEXPR_OR_CONST float_type one_over_ln2 = float_type(BOOST_FLOAT128_C(1.44269504088896340735992468100189213742664595415299));
        const float_type x_over_ln2   = x * one_over_ln2;

        boost::int_fast32_t n;

        if     (x < -1) { n = static_cast<boost::int_fast32_t>(BOOST_CSTDFLOAT_FLOAT128_CEIL (x_over_ln2)); }
        else if(x > +1) { n = static_cast<boost::int_fast32_t>(BOOST_CSTDFLOAT_FLOAT128_FLOOR(x_over_ln2)); }
        else            { n = static_cast<boost::int_fast32_t>(0); }

        // Here, alpha is the scaled argument.
        const float_type alpha = x - (n * float_type(BOOST_FLOAT128_C(0.693147180559945309417232121458176568075500134360255)));

        // Compute the polynomial approximation of the scaled-argument exponential function.
        const float_type sum =
          ((((((((((((((((((((((((((((((((((((  float_type(BOOST_FLOAT128_C(2.69291698127774166063293705964720493864630783729857438187365E-42))  * alpha
                                              + float_type(BOOST_FLOAT128_C(9.70937085471487654794114679403710456028986572118859594614033E-41))) * alpha
                                              + float_type(BOOST_FLOAT128_C(3.38715585158055097155585505318085512156885389014410753080500E-39))) * alpha
                                              + float_type(BOOST_FLOAT128_C(1.15162718532861050809222658798662695267019717760563645440433E-37))) * alpha
                                              + float_type(BOOST_FLOAT128_C(3.80039074689434663295873584133017767349635602413675471702393E-36))) * alpha
                                              + float_type(BOOST_FLOAT128_C(1.21612504934087520075905434734158045947460467096773246215239E-34))) * alpha
                                              + float_type(BOOST_FLOAT128_C(3.76998762883139753126119821241037824830069851253295480396224E-33))) * alpha
                                              + float_type(BOOST_FLOAT128_C(1.13099628863830344684998293828608215735777107850991029729440E-31))) * alpha
                                              + float_type(BOOST_FLOAT128_C(3.27988923706982293204067897468714277771890104022419696770352E-30))) * alpha
                                              + float_type(BOOST_FLOAT128_C(9.18368986379558482800593745627556950089950023355628325088207E-29))) * alpha
                                              + float_type(BOOST_FLOAT128_C(2.47959626322479746949155352659617642905315302382639380521497E-27))) * alpha
                                              + float_type(BOOST_FLOAT128_C(6.44695028438447337900255966737803112935639344283098705091949E-26))) * alpha
                                              + float_type(BOOST_FLOAT128_C(1.61173757109611834904452725462599961406036904573072897122957E-24))) * alpha
                                              + float_type(BOOST_FLOAT128_C(3.86817017063068403772269360016918092488847584660382953555804E-23))) * alpha
                                              + float_type(BOOST_FLOAT128_C(8.89679139245057328674891109315654704307721758924206107351744E-22))) * alpha
                                              + float_type(BOOST_FLOAT128_C(1.95729410633912612308475595397946731738088422488032228717097E-20))) * alpha
                                              + float_type(BOOST_FLOAT128_C(4.11031762331216485847799061511674191805055663711439605760231E-19))) * alpha
                                              + float_type(BOOST_FLOAT128_C(8.22063524662432971695598123977873600603370758794431071426640E-18))) * alpha
                                              + float_type(BOOST_FLOAT128_C(1.56192069685862264622163643500633782667263448653185159383285E-16))) * alpha
                                              + float_type(BOOST_FLOAT128_C(2.81145725434552076319894558300988749849555291507956994126835E-15))) * alpha
                                              + float_type(BOOST_FLOAT128_C(4.77947733238738529743820749111754320727153728139716409114011E-14))) * alpha
                                              + float_type(BOOST_FLOAT128_C(7.64716373181981647590113198578807092707697416852226691068627E-13))) * alpha
                                              + float_type(BOOST_FLOAT128_C(1.14707455977297247138516979786821056670509688396295740818677E-11))) * alpha
                                              + float_type(BOOST_FLOAT128_C(1.60590438368216145993923771701549479323291461578567184216302E-10))) * alpha
                                              + float_type(BOOST_FLOAT128_C(2.08767569878680989792100903212014323125428376052986408239620E-09))) * alpha
                                              + float_type(BOOST_FLOAT128_C(2.50521083854417187750521083854417187750523408006206780016659E-08))) * alpha
                                              + float_type(BOOST_FLOAT128_C(2.75573192239858906525573192239858906525573195144226062684604E-07))) * alpha
                                              + float_type(BOOST_FLOAT128_C(2.75573192239858906525573192239858906525573191310049321957902E-06))) * alpha
                                              + float_type(BOOST_FLOAT128_C(0.00002480158730158730158730158730158730158730158730149317774)))     * alpha
                                              + float_type(BOOST_FLOAT128_C(0.00019841269841269841269841269841269841269841269841293575920)))     * alpha
                                              + float_type(BOOST_FLOAT128_C(0.00138888888888888888888888888888888888888888888888889071045)))     * alpha
                                              + float_type(BOOST_FLOAT128_C(0.00833333333333333333333333333333333333333333333333332986595)))     * alpha
                                              + float_type(BOOST_FLOAT128_C(0.04166666666666666666666666666666666666666666666666666664876)))     * alpha
                                              + float_type(BOOST_FLOAT128_C(0.16666666666666666666666666666666666666666666666666666669048)))     * alpha
                                              + float_type(BOOST_FLOAT128_C(0.50000000000000000000000000000000000000000000000000000000006)))     * alpha
                                              + float_type(BOOST_FLOAT128_C(0.99999999999999999999999999999999999999999999999999999999995)))     * alpha
                                              + float_type(BOOST_FLOAT128_C(1.0)));

        // Rescale the result and return it.
        return sum * BOOST_CSTDFLOAT_FLOAT128_POW(float_type(2), float_type(n));
      }
      extern "C" boost::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_LOG   (boost::cstdfloat::detail::float_internal128_t) throw();
      extern "C" boost::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_LOG10 (boost::cstdfloat::detail::float_internal128_t) throw();
      extern "C" boost::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_SIN   (boost::cstdfloat::detail::float_internal128_t) throw();
      extern "C" boost::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_COS   (boost::cstdfloat::detail::float_internal128_t) throw();
      extern "C" boost::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_TAN   (boost::cstdfloat::detail::float_internal128_t) throw();
      extern "C" boost::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_ASIN  (boost::cstdfloat::detail::float_internal128_t) throw();
      extern "C" boost::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_ACOS  (boost::cstdfloat::detail::float_internal128_t) throw();
      extern "C" boost::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_ATAN  (boost::cstdfloat::detail::float_internal128_t) throw();
      inline     boost::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_SINH  (boost::cstdfloat::detail::float_internal128_t x)
      {
        // Patch the sinhq() function for a subset of broken GCC compilers
        // like GCC 4.7, 4.8 on MinGW.
        const boost::cstdfloat::detail::float_internal128_t ex = BOOST_CSTDFLOAT_FLOAT128_EXP(x);
        return (ex - (1 / ex)) / 2;
      }
      inline     boost::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_COSH  (boost::cstdfloat::detail::float_internal128_t x)
      {
        // Patch the coshq() function for a subset of broken GCC compilers
        // like GCC 4.7, 4.8 on MinGW.
        const boost::cstdfloat::detail::float_internal128_t ex = BOOST_CSTDFLOAT_FLOAT128_EXP(x);
        return (ex + (1 / ex)) / 2;
      }
      inline     boost::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_TANH  (boost::cstdfloat::detail::float_internal128_t x)
      {
        // Patch the tanhq() function for a subset of broken GCC compilers
        // like GCC 4.7, 4.8 on MinGW.
        return BOOST_CSTDFLOAT_FLOAT128_SINH(x) / BOOST_CSTDFLOAT_FLOAT128_COSH(x);
      }
      extern "C" boost::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_FMOD  (boost::cstdfloat::detail::float_internal128_t, boost::cstdfloat::detail::float_internal128_t) throw();
      extern "C" boost::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_ATAN2 (boost::cstdfloat::detail::float_internal128_t, boost::cstdfloat::detail::float_internal128_t) throw();
      extern "C" boost::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_LGAMMA(boost::cstdfloat::detail::float_internal128_t) throw();
      extern "C" boost::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_TGAMMA(boost::cstdfloat::detail::float_internal128_t) throw();

      // Define the quadruple-precision <math.h> functions in the namespace
      // boost::cstdfloat::detail.

      namespace boost { namespace cstdfloat { namespace detail {
      inline   boost::cstdfloat::detail::float_internal128_t ldexp (boost::cstdfloat::detail::float_internal128_t x, int n)                                           { return ::BOOST_CSTDFLOAT_FLOAT128_LDEXP (x, n); }
      inline   boost::cstdfloat::detail::float_internal128_t frexp (boost::cstdfloat::detail::float_internal128_t x, int* pn)                                         { return ::BOOST_CSTDFLOAT_FLOAT128_FREXP (x, pn); }
      inline   boost::cstdfloat::detail::float_internal128_t fabs  (boost::cstdfloat::detail::float_internal128_t x)                                                  { return ::BOOST_CSTDFLOAT_FLOAT128_FABS  (x); }
      inline   boost::cstdfloat::detail::float_internal128_t floor (boost::cstdfloat::detail::float_internal128_t x)                                                  { return ::BOOST_CSTDFLOAT_FLOAT128_FLOOR (x); }
      inline   boost::cstdfloat::detail::float_internal128_t ceil  (boost::cstdfloat::detail::float_internal128_t x)                                                  { return ::BOOST_CSTDFLOAT_FLOAT128_CEIL  (x); }
      inline   boost::cstdfloat::detail::float_internal128_t sqrt  (boost::cstdfloat::detail::float_internal128_t x)                                                  { return ::BOOST_CSTDFLOAT_FLOAT128_SQRT  (x); }
      inline   boost::cstdfloat::detail::float_internal128_t trunc (boost::cstdfloat::detail::float_internal128_t x)                                                  { return ::BOOST_CSTDFLOAT_FLOAT128_TRUNC (x); }
      inline   boost::cstdfloat::detail::float_internal128_t pow   (boost::cstdfloat::detail::float_internal128_t x, boost::cstdfloat::detail::float_internal128_t a) { return ::BOOST_CSTDFLOAT_FLOAT128_POW   (x, a); }
      inline   boost::cstdfloat::detail::float_internal128_t exp   (boost::cstdfloat::detail::float_internal128_t x)                                                  { return ::BOOST_CSTDFLOAT_FLOAT128_EXP   (x); }
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

      // Now *use* the quadruple-precision <math.h> functions in the global namespace.
      using boost::cstdfloat::detail::ldexp;
      using boost::cstdfloat::detail::frexp;
      using boost::cstdfloat::detail::fabs;
      using boost::cstdfloat::detail::floor;
      using boost::cstdfloat::detail::ceil;
      using boost::cstdfloat::detail::sqrt;
      using boost::cstdfloat::detail::trunc;
      using boost::cstdfloat::detail::pow;
      using boost::cstdfloat::detail::exp;
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

      // Now *use* the quadruple-precision <math.h> functions in the std namespace.
      namespace std
      {
        using boost::cstdfloat::detail::ldexp;
        using boost::cstdfloat::detail::frexp;
        using boost::cstdfloat::detail::fabs;
        using boost::cstdfloat::detail::floor;
        using boost::cstdfloat::detail::ceil;
        using boost::cstdfloat::detail::sqrt;
        using boost::cstdfloat::detail::trunc;
        using boost::cstdfloat::detail::pow;
        using boost::cstdfloat::detail::exp;
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
      }
    #endif // Not BOOST_CSTDFLOAT_NO_LIBQUADMATH_CMATH (i.e., has libquadmath <cmath> support)

    #if !defined(BOOST_CSTDFLOAT_NO_LIBQUADMATH_IOSTREAM)

      // Implement quadruple-precision I/O stream operations.

      #include <istream>
      #include <ostream>
      #include <sstream>
      #include <string>

      #if defined(__GNUC__)
//      #if 0

      namespace std
      {
        template<typename char_type, class traits_type>
        inline std::basic_ostream<char_type, traits_type>& operator<<(std::basic_ostream<char_type, traits_type>& os, const boost::cstdfloat::detail::float_internal128_t& x)
        {
          std::basic_ostringstream<char_type, traits_type> ostr;
          ostr.flags(os.flags());
          ostr.imbue(os.getloc());
          ostr.precision(os.precision());

          char my_buffer[64U];

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
            // Evidently there is a really long floating-point string here,
            // such as a small decimal representation. So we have to use
            // dynamic memory allocation.

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

            static_cast<void>(ostr << my_buffer2);

            delete [] my_buffer2;
          }
          else
          {
            static_cast<void>(ostr << my_buffer);
          }

          return (os << ostr.str());
        }

        template<typename char_type, class traits_type>
        inline std::basic_istream<char_type, traits_type>& operator>>(std::basic_istream<char_type, traits_type>& is, boost::cstdfloat::detail::float_internal128_t& x)
        {
          std::string str;

          is >> str;

          char* p_end;

          x = strtoflt128(str.c_str(), &p_end);

          if(static_cast<std::ptrdiff_t>(p_end - str.c_str()) != static_cast<std::ptrdiff_t>(str.length()))
          {
            BOOST_THROW_EXCEPTION(std::runtime_error("Unable to interpret input string as a boost::float128_t"));

            is.setstate(ios_base::failbit);
          }

          return is;
        }
      }

      #elif defined(BOOST_INTEL)
//      #elif defined(__GNUC__)

      // The following string-extraction routines are based on the methodology
      // used in Boost.Multiprecision by John Maddock and Christopher Kormanyos.
      // This methodology has been slightly modified here for boost::float128_t.

      #include <cstring>
      #include <cctype>
      #include <boost/lexical_cast.hpp>

      namespace boost { namespace cstdfloat { namespace detail {

      template<class string_type>
      void format_float_string(string_type& str,
                               int my_exp,
                               int digits,
                               const std::ios_base::fmtflags f,
                               const bool iszero)
      {
        typedef typename string_type::size_type size_type;

        const bool scientific = ((f & std::ios_base::scientific) == std::ios_base::scientific);
        const bool fixed      = ((f & std::ios_base::fixed)      == std::ios_base::fixed);
        const bool showpoint  = ((f & std::ios_base::showpoint)  == std::ios_base::showpoint);
        const bool showpos    = ((f & std::ios_base::showpos)    == std::ios_base::showpos);

        const bool neg = str.size() && (str[0] == '-');

        if(neg)
        {
          str.erase(0, 1);
        }

        if(digits == 0)
        {
          digits = static_cast<int>((std::max)(str.size(), size_type(16)));
        }

        if(iszero || str.empty() || (str.find_first_not_of('0') == string_type::npos))
        {
          // We will be printing zero, even though the value might not
          // actually be zero (it just may have been rounded to zero).
          str = "0";

          if(scientific || fixed)
          {
            str.append(1, '.');
            str.append(size_type(digits), '0');

            if(scientific)
            {
              str.append("e+00");
            }
          }
          else
          {
            if(showpoint)
            {
              str.append(1, '.');
              if(digits > 1)
              {
                str.append(size_type(digits - 1), '0');
              }
            }
          }
          if(neg)
          {
            str.insert(0U, 1U, '-');
          }
          else if(showpos)
          {
            str.insert(0U, 1U, '+');
          }

          return;
        }

        if(!fixed && !scientific && !showpoint)
        {
          // Suppress trailing zeros:
          typename string_type::iterator pos = str.end();

          while(pos != str.begin() && *--pos == '0') { ; }

          if(pos != str.end())
          {
            ++pos;
          }

          str.erase(pos, str.end());

          if(str.empty())
          {
            str = '0';
          }
        }
        else if(!fixed || (my_exp >= 0))
        {
          // Pad out the end with zero's if we need to.

          int chars = static_cast<int>(str.size());
          chars = digits - chars;

          if(scientific)
          {
            ++chars;
          }

          if(chars > 0)
          {
            str.append(static_cast<size_type>(chars), '0');
          }
        }

        if(fixed || (!scientific && (my_exp >= -4) && (my_exp < digits)))
        {
          if((1 + my_exp) > static_cast<int>(str.size()))
          {
            // Just pad out the end with zeros:
            str.append(static_cast<size_type>((1 + my_exp) - static_cast<int>(str.size())), '0');

            if(showpoint || fixed)
            {
              str.append(".");
            }
          }
          else if(my_exp + 1 < static_cast<int>(str.size()))
          {
            if(my_exp < 0)
            {
              str.insert(0U, static_cast<size_type>(-1 - my_exp), '0');
              str.insert(0U, "0.");
            }
            else
            {
              // Insert the decimal point:
              str.insert(static_cast<size_type>(my_exp + 1), 1, '.');
            }
          }
          else if(showpoint || fixed) // we have exactly the digits we require to left of the point
          {
            str += ".";
          }

          if(fixed)
          {
            // We may need to add trailing zeros.
            int l = static_cast<int>(str.find('.') + 1U);
            l = digits - (static_cast<int>(str.size()) - l);

            if(l > 0)
            {
              str.append(size_type(l), '0');
            }
          }
        }
        else
        {
          // Scientific format:
          if(showpoint || (str.size() > 1))
          {
            str.insert(1U, 1U, '.');
          }

          str.append(1U, 'e');
          string_type e = boost::lexical_cast<string_type>(std::abs(my_exp));

          if(e.size() < 2U)
          {
            e.insert(0U, 2U - e.size(), '0');
          }

          if(my_exp < 0)
          {
            e.insert(0U, 1U, '-');
          }
          else
          {
            e.insert(0U, 1U, '+');
          }

          str.append(e);
        }

        if(neg)
        {
          str.insert(0U, 1U, '-');
        }
        else if(showpos)
        {
          str.insert(0U, 1U, '+');
        }
      }

      template<class float_type, class type_a> inline void eval_convert_to(type_a* pa,    const float_type& cb)                        { *pa  = static_cast<type_a>(cb); }
      template<class float_type, class type_a> inline void eval_add       (float_type& b, const type_a& a)                             { b   += a; }
      template<class float_type, class type_a> inline void eval_subtract  (float_type& b, const type_a& a)                             { b   -= a; }
      template<class float_type, class type_a> inline void eval_multiply  (float_type& b, const type_a& a)                             { b   *= a; }
      template<class float_type>               inline void eval_multiply  (float_type& b, const float_type& cb, const float_type& cb2) { b    = (cb * cb2); }
      template<class float_type, class type_a> inline void eval_divide    (float_type& b, const type_a& a)                             { b   /= a; }
      template<class float_type>               inline void eval_log10     (float_type& b, const float_type& cb)                        { b    = std::log10(cb); }
      template<class float_type>               inline void eval_floor     (float_type& b, const float_type& cb)                        { b    = std::floor(cb); }

      template<class float_type, class type_n> inline float_type pown(const float_type& cb, const type_n p)
      {
        if     (p <  static_cast<type_n>(0)) { return 1 / pown(cb, static_cast<type_n>(-p)); }
        else if(p == static_cast<type_n>(0)) { return float_type(1); }
        else if(p == static_cast<type_n>(1)) { return  cb; }
        else if(p == static_cast<type_n>(2)) { return  cb * cb; }
        else if(p == static_cast<type_n>(3)) { return (cb * cb) * cb; }
        else
        {
          float_type value = cb;

          type_n n;

          for(n = static_cast<type_n>(1); n <= static_cast<type_n>(p / 2); n *= 2)
          {
            value *= value;
          }

          const type_n p_minus_n = static_cast<type_n>(p - n);

          // Call the function recursively for computing the remaining power of n.
          return ((p_minus_n == static_cast<type_n>(0)) ? value : (value * pown(cb, p_minus_n)));
        }
      }

      inline void round_string_up_at(std::string& s, int pos, int& expon)
      {
        // This subroutine rounds up a string representation of a
        // number at the given position pos.

        if(pos < 0)
        {
          s.insert(0U, 1U, '1');
          s.erase(s.size() - 1U);
          ++expon;
        }
        else if(s[pos] == '9')
        {
          s[pos] = '0';
          round_string_up_at(s, pos - 1, expon);
        }
        else
        {
          if((pos == 0) && (s[pos] == '0') && (s.size() == 1))
          {
            ++expon;
          }

          ++s[pos];
        }
      }

      template<class float_type>
      std::string convert_to_string(float_type& x,
                                    std::streamsize digits,
                                    const std::ios_base::fmtflags f)
      {
        const bool iszero = (std::fabs(x) < (std::numeric_limits<float_type>::min)());
        const bool isneg  = (x < 0);
        const bool isnan  = (x != x);
        const bool isinf  = (std::fabs(x) > (std::numeric_limits<float_type>::max)());

        int expon = 0;

        if(digits <= 0) { digits = std::numeric_limits<float_type>::max_digits10; }

        const int org_digits = static_cast<int>(digits);

        std::string result;

        if(iszero)
        {
          result = "0";
        }
        else if(isinf)
        {
          if(x < 0)
          {
            return "-inf";
          }
          else
          {
            return ((f & std::ios_base::showpos) == std::ios_base::showpos) ? "+inf" : "inf";
          }
        }
        else if(isnan)
        {
          return "nan";
        }
        else
        {
          // Start by figuring out the base-10 exponent.
          if(isneg) { x = -x; }

          float_type t;
          float_type ten = 10;

          eval_log10(t, x);
          eval_floor(t, t);
          eval_convert_to(&expon, t);

          if(-expon > std::numeric_limits<float_type>::max_exponent10 - 3)
          {
            int e = -expon / 2;

            const float_type t2 = pown(ten, e);

            eval_multiply(t, t2, x);
            eval_multiply(t, t2);

            if((expon & 1) != 0)
            {
              eval_multiply(t, ten);
            }
          }
          else
          {
            t = pown(ten, -expon);
            eval_multiply(t, x);
          }

          // Make sure that the value lies between [1, 10), and adjust if not.
          if(t < 1)
          {
            eval_multiply(t, 10);

            --expon;
          }
          else if(t >= 10)
          {
            eval_divide(t, 10);

            ++expon;
          }

          float_type digit;
          int        cdigit;

          // Adjust the number of digits required based on formatting options.
          if(((f & std::ios_base::fixed) == std::ios_base::fixed) && (expon != -1))
          {
            digits += (expon + 1);
          }

          if((f & std::ios_base::scientific) == std::ios_base::scientific)
          {
            ++digits;
          }

          // Extract the base-10 digits one at a time.
          for(int i = 0; i < digits; ++i)
          {
            eval_floor(digit, t);
            eval_convert_to(&cdigit, digit);

            result += static_cast<char>('0' + cdigit);

            eval_subtract(t, digit);
            eval_multiply(t, ten);
          }

          // Possibly round the result.
          if(digits >= 0)
          {
            eval_floor(digit, t);
            eval_convert_to(&cdigit, digit);
            eval_subtract(t, digit);

            if((cdigit == 5) && (t == 0))
            {
              // Use simple bankers rounding.

              if((static_cast<int>(*result.rbegin() - '0') & 1) != 0)
              {
                round_string_up_at(result, static_cast<int>(result.size() - 1U), expon);
              }
            }
            else if(cdigit >= 5)
            {
              round_string_up_at(result, static_cast<int>(result.size() - 1), expon);
            }
          }
        }

        while((result.size() > static_cast<std::string::size_type>(digits)) && result.size())
        {
          // We may get here as a result of rounding.

          if(result.size() > 1U)
          {
            result.erase(result.size() - 1U);
          }
          else
          {
            if(expon > 0)
            {
              --expon; // so we put less padding in the result.
            }
            else
            {
              ++expon;
            }

            ++digits;
          }
        }

        if(isneg)
        {
          result.insert(0U, 1U, '-');
        }

        format_float_string(result, expon, org_digits, f, iszero);

        return result;
      }

      template <class float_type>
      bool convert_from_string(float_type& value, const char* p)
      {
        value = 0;

        if((p == static_cast<const char*>(0U)) || (*p == static_cast<char>(0)))
        {
          return;
        }

        bool is_neg       = false;
        bool is_neg_expon = false;

        const int ten = 10;

        int expon       = 0;
        int digits_seen = 0;

        const int max_digits = (std::numeric_limits<float_type>::is_specialized ? (std::numeric_limits<float_type>::max_digits10 + 1) : 128);

        if(*p == static_cast<char>('+'))
        {
          ++p;
        }
        else if(*p == static_cast<char>('-'))
        {
          is_neg = true;
          ++p;
        }

        const bool isnan = ((std::strcmp(p, "nan") == 0) || (std::strcmp(p, "NaN") == 0) || (std::strcmp(p, "NAN") == 0));

        if(isnan)
        {
          eval_divide(value, 0);

          if(is_neg)
          {
            value = -value;
          }

          return true;
        }

        const bool isinf = ((std::strcmp(p, "inf") == 0) || (std::strcmp(p, "Inf") == 0) || (std::strcmp(p, "INF") == 0));

        if(isinf)
        {
          value = 1;
          eval_divide(value, 0);

          if(is_neg)
          {
            value = -value;
          }

          return true;
        }

        // Grab all the leading digits before the decimal point.
        while(std::isdigit(*p))
        {
          eval_multiply(value, ten);
          eval_add(value, static_cast<int>(*p - '0'));
          ++p;
          ++digits_seen;
        }

        if(*p == static_cast<char>('.'))
        {
          // Grab everything after the point, stop when we've seen
          // enough digits, even if there are actually more available.

          ++p;

          while(std::isdigit(*p))
          {
            eval_multiply(value, ten);
            eval_add(value, static_cast<int>(*p - '0'));
            ++p;
            --expon;

            if(++digits_seen > max_digits)
            {
              break;
            }
          }

          while(std::isdigit(*p))
          {
            ++p;
          }
        }

        // Parse the exponent.
        if((*p == static_cast<char>('e')) || (*p == static_cast<char>('E')))
        {
          ++p;

          if(*p == static_cast<char>('+'))
          {
            ++p;
          }
          else if(*p == static_cast<char>('-'))
          {
            is_neg_expon = true;
            ++p;
          }

          int e2 = 0;

          while(std::isdigit(*p))
          {
            e2 *= 10;
            e2 += (*p - '0');
            ++p;
          }

          if(is_neg_expon)
          {
            e2 = -e2;
          }

          expon += e2;
        }

        if(expon)
        {
          // Scale by 10^expon. Note that 10^expon can be outside the range
          // of our number type, even though the result is within range.
          // If that looks likely, then split the calculation in two parts.
          float_type t;
          t = ten;

          if(expon > (std::numeric_limits<float_type>::min_exponent10 + 2))
          {
            t = pown(t, expon);
            eval_multiply(value, t);
          }
          else
          {
            t = pown(t, (expon + digits_seen + 1));
            eval_multiply(value, t);
            t = ten;
            t = pown(t, (-digits_seen - 1));
            eval_multiply(value, t);
          }
        }

        if(is_neg)
        {
          value = -value;
        }

        return (*p == static_cast<char>(0));
      }
      } } } // boost::cstdfloat::detail

      namespace std
      {
        template<typename char_type, class traits_type>
        inline std::basic_ostream<char_type, traits_type>& operator<<(std::basic_ostream<char_type, traits_type>& os, const boost::cstdfloat::detail::float_internal128_t& x)
        {
          boost::cstdfloat::detail::float_internal128_t non_const_x = x;

          const std::string str = boost::cstdfloat::detail::convert_to_string(non_const_x,
                                                                              os.precision(),
                                                                              os.flags());

          std::basic_ostringstream<char_type, traits_type> ostr;
          ostr.flags(os.flags());
          ostr.imbue(os.getloc());
          ostr.precision(os.precision());

          static_cast<void>(ostr << str);

          return (os << ostr.str());
        }

        template<typename char_type, class traits_type>
        inline std::basic_istream<char_type, traits_type>& operator>>(std::basic_istream<char_type, traits_type>& is, boost::cstdfloat::detail::float_internal128_t& x)
        {
          std::string str;

          is >> str;

          const bool conversion_is_ok = boost::cstdfloat::detail::convert_from_string(x, str.c_str());

          if(false == conversion_is_ok)
          {
            BOOST_THROW_EXCEPTION(std::runtime_error("Unable to interpret input string as a boost::float128_t"));

            is.putback(str);

            is.setstate(ios_base::failbit);
          }

          return is;
        }
      }

      #endif // Use __GNUC__ or BOOST_INTEL libquadmath

    #endif // Not BOOST_CSTDFLOAT_NO_LIBQUADMATH_IOSTREAM (i.e., has libquadmath I/O stream support)

    #if !defined(BOOST_CSTDFLOAT_NO_LIBQUADMATH_COMPLEX)

    #include <complex>

    // Implement a specialization of std::complex<> for quadruple-precision.
    namespace std
    {
      // Forward template function declarations.
      #if defined(BOOST_NO_CXX11_CONSTEXPR)
      template<> boost::cstdfloat::detail::float_internal128_t& real<boost::cstdfloat::detail::float_internal128_t>(complex<boost::cstdfloat::detail::float_internal128_t>&);
      template<> boost::cstdfloat::detail::float_internal128_t& imag<boost::cstdfloat::detail::float_internal128_t>(complex<boost::cstdfloat::detail::float_internal128_t>&);
      #else
      template<> BOOST_CONSTEXPR boost::cstdfloat::detail::float_internal128_t real<boost::cstdfloat::detail::float_internal128_t>(const complex<boost::cstdfloat::detail::float_internal128_t>&);
      template<> BOOST_CONSTEXPR boost::cstdfloat::detail::float_internal128_t imag<boost::cstdfloat::detail::float_internal128_t>(const complex<boost::cstdfloat::detail::float_internal128_t>&);
      #endif
      template<> boost::cstdfloat::detail::float_internal128_t abs <boost::cstdfloat::detail::float_internal128_t>(const complex<boost::cstdfloat::detail::float_internal128_t>&);
      template<> boost::cstdfloat::detail::float_internal128_t arg <boost::cstdfloat::detail::float_internal128_t>(const complex<boost::cstdfloat::detail::float_internal128_t>&);
      template<> boost::cstdfloat::detail::float_internal128_t norm<boost::cstdfloat::detail::float_internal128_t>(const complex<boost::cstdfloat::detail::float_internal128_t>&);

      template<> complex<boost::cstdfloat::detail::float_internal128_t> sqrt (const complex<boost::cstdfloat::detail::float_internal128_t>&);
      template<> complex<boost::cstdfloat::detail::float_internal128_t> sin  (const complex<boost::cstdfloat::detail::float_internal128_t>&);
      template<> complex<boost::cstdfloat::detail::float_internal128_t> cos  (const complex<boost::cstdfloat::detail::float_internal128_t>&);
      template<> complex<boost::cstdfloat::detail::float_internal128_t> tan  (const complex<boost::cstdfloat::detail::float_internal128_t>&);
      #if !defined(BOOST_NO_CXX11_CONSTEXPR)
      template<> complex<boost::cstdfloat::detail::float_internal128_t> asin (const complex<boost::cstdfloat::detail::float_internal128_t>&);
      template<> complex<boost::cstdfloat::detail::float_internal128_t> acos (const complex<boost::cstdfloat::detail::float_internal128_t>&);
      template<> complex<boost::cstdfloat::detail::float_internal128_t> atan (const complex<boost::cstdfloat::detail::float_internal128_t>&);
      #endif
      template<> complex<boost::cstdfloat::detail::float_internal128_t> exp  (const complex<boost::cstdfloat::detail::float_internal128_t>&);
      template<> complex<boost::cstdfloat::detail::float_internal128_t> log  (const complex<boost::cstdfloat::detail::float_internal128_t>&);
      template<> complex<boost::cstdfloat::detail::float_internal128_t> log10(const complex<boost::cstdfloat::detail::float_internal128_t>&);
      template<> complex<boost::cstdfloat::detail::float_internal128_t> pow  (const complex<boost::cstdfloat::detail::float_internal128_t>&, const boost::cstdfloat::detail::float_internal128_t&);
      template<> complex<boost::cstdfloat::detail::float_internal128_t> pow  (const complex<boost::cstdfloat::detail::float_internal128_t>&, const complex<boost::cstdfloat::detail::float_internal128_t>&);
      template<> complex<boost::cstdfloat::detail::float_internal128_t> pow  (const boost::cstdfloat::detail::float_internal128_t&, const complex<boost::cstdfloat::detail::float_internal128_t>&);
      template<> complex<boost::cstdfloat::detail::float_internal128_t> sinh (const complex<boost::cstdfloat::detail::float_internal128_t>&);
      template<> complex<boost::cstdfloat::detail::float_internal128_t> cosh (const complex<boost::cstdfloat::detail::float_internal128_t>&);
      template<> complex<boost::cstdfloat::detail::float_internal128_t> tanh (const complex<boost::cstdfloat::detail::float_internal128_t>&);
      #if !defined(BOOST_NO_CXX11_CONSTEXPR)
      template<> complex<boost::cstdfloat::detail::float_internal128_t> asinh(const complex<boost::cstdfloat::detail::float_internal128_t>&);
      template<> complex<boost::cstdfloat::detail::float_internal128_t> acosh(const complex<boost::cstdfloat::detail::float_internal128_t>&);
      template<> complex<boost::cstdfloat::detail::float_internal128_t> atanh(const complex<boost::cstdfloat::detail::float_internal128_t>&);
      #endif

      template<>
      class complex<boost::cstdfloat::detail::float_internal128_t>
      {
      public:
        typedef boost::cstdfloat::detail::float_internal128_t value_type;

        #if defined(BOOST_NO_CXX11_CONSTEXPR)

        complex(const complex<value_type>& z) : re(z.re),
                                                im(z.im) { }

        complex(value_type Re = BOOST_FLOAT128_C(0.0),
                value_type Im = BOOST_FLOAT128_C(0.0)) : re(Re),
                                                         im(Im) { }

        explicit complex(const complex<float>&);
        explicit complex(const complex<double>&);
        explicit complex(const complex<long double>&);

        const value_type& real() const { return re; }
        const value_type& imag() const { return im; }
              value_type& real()       { return re; }
              value_type& imag()       { return im; }

        #else

        BOOST_CONSTEXPR complex(const complex<value_type>& z) : re(z.re),
                                                                im(z.im) { }

        BOOST_CONSTEXPR complex(value_type Re = BOOST_FLOAT128_C(0.0),
                                value_type Im = BOOST_FLOAT128_C(0.0)) : re(Re),
                                                                         im(Im) { }

        explicit BOOST_CONSTEXPR complex(const complex<float>&);
        explicit BOOST_CONSTEXPR complex(const complex<double>&);
        explicit BOOST_CONSTEXPR complex(const complex<long double>&);

        BOOST_CONSTEXPR value_type real() { return re; }
        BOOST_CONSTEXPR value_type imag() { return im; }

        #endif

        void real(value_type v) { re = v; }
        void imag(value_type v) { im = v; }

        complex<value_type>& operator+=(value_type v)
        {
          re += v;
          return *this;
        }

        complex<value_type>& operator-=(value_type v)
        {
          re -= v;
          return *this;
        }

        complex<value_type>& operator*=(value_type v)
        {
          re *= v;
          im *= v;
          return *this;
        }

        complex<value_type>& operator/=(value_type v)
        {
          re /= v;
          im /= v;
          return *this;
        }

        template<typename X> complex<value_type>& operator+=(const complex<X>& x)
        {
          re += value_type(x.re);
          im += value_type(x.im);
          return *this;
        }

        template<typename X> complex<value_type>& operator-=(const complex<X>& x)
        {
          re -= value_type(x.re);
          im -= value_type(x.im);
          return *this;
        }

        template<typename X> complex<value_type>& operator*=(const complex<X>& x)
        {
          const value_type re_x(x.re);
          const value_type im_x(x.im);
          const value_type tmp_re((re * re_x) - (im * im_x));
          const value_type tmp_im((re * im_x) + (im * re_x));
          re = tmp_re;
          im = tmp_im;
          return *this;
        }

        template<typename X> complex<value_type>& operator/=(const complex<X>& x)
        {
          const value_type re_x(x.re);
          const value_type im_x(x.im);
          const value_type one_over_denom = 1 / std::sqrt((re_x * re_x) + (im_x * im_x));
          const value_type tmp_re = ((re * re_x) + (im * im_x)) * one_over_denom;
          const value_type tmp_im = ((im * re_x) - (re * im_x)) * one_over_denom;
          re = tmp_re;
          im = tmp_im;
          return *this;
        }

        template<typename X>
        complex<value_type>& operator=(const complex<X>& z)
        {
          re = z.real();
          im = z.imag();
          return *this;
        }

        complex<value_type>& operator=(const value_type& v) { re = v; im = value_type(0); return *this; }

      private:
        value_type re;
        value_type im;
      };

      #if defined(BOOST_NO_CXX11_CONSTEXPR)
      complex<boost::cstdfloat::detail::float_internal128_t>::complex(const complex<float>&        f) : re(boost::cstdfloat::detail::float_internal128_t( f.real())), im(boost::cstdfloat::detail::float_internal128_t( f.imag())) { }
      complex<boost::cstdfloat::detail::float_internal128_t>::complex(const complex<double>&       d) : re(boost::cstdfloat::detail::float_internal128_t( d.real())), im(boost::cstdfloat::detail::float_internal128_t( d.imag())) { }
      complex<boost::cstdfloat::detail::float_internal128_t>::complex(const complex<long double>& ld) : re(boost::cstdfloat::detail::float_internal128_t(ld.real())), im(boost::cstdfloat::detail::float_internal128_t(ld.imag())) { }
      #else
      BOOST_CONSTEXPR complex<boost::cstdfloat::detail::float_internal128_t>::complex(const complex<float>&        f) : re(boost::cstdfloat::detail::float_internal128_t( f.real())), im(boost::cstdfloat::detail::float_internal128_t( f.imag())) { }
      BOOST_CONSTEXPR complex<boost::cstdfloat::detail::float_internal128_t>::complex(const complex<double>&       d) : re(boost::cstdfloat::detail::float_internal128_t( d.real())), im(boost::cstdfloat::detail::float_internal128_t( d.imag())) { }
      BOOST_CONSTEXPR complex<boost::cstdfloat::detail::float_internal128_t>::complex(const complex<long double>& ld) : re(boost::cstdfloat::detail::float_internal128_t(ld.real())), im(boost::cstdfloat::detail::float_internal128_t(ld.imag())) { }
      #endif
    } // namespace std

    namespace boost { namespace cstdfloat { namespace detail {
    template<class T> std::complex<T> iz_helper__x(const std::complex<T>& x)
    {
      const T tmp_r = x.real();
      return std::complex<T>(-x.imag(), tmp_r);
    }
    } } } // boost::cstdfloat::detail

    namespace std
    {
      // 26.4.7, specific values.
      #if defined(BOOST_NO_CXX11_CONSTEXPR)
      template<> boost::cstdfloat::detail::float_internal128_t& real<boost::cstdfloat::detail::float_internal128_t>(complex<boost::cstdfloat::detail::float_internal128_t>& x) { return x.real(); }
      template<> boost::cstdfloat::detail::float_internal128_t& imag<boost::cstdfloat::detail::float_internal128_t>(complex<boost::cstdfloat::detail::float_internal128_t>& x) { return x.imag(); }
      #else
      template<> BOOST_CONSTEXPR boost::cstdfloat::detail::float_internal128_t real<boost::cstdfloat::detail::float_internal128_t>(const complex<boost::cstdfloat::detail::float_internal128_t>& x) { return x.real(); }
      template<> BOOST_CONSTEXPR boost::cstdfloat::detail::float_internal128_t imag<boost::cstdfloat::detail::float_internal128_t>(const complex<boost::cstdfloat::detail::float_internal128_t>& x) { return x.imag(); }
      #endif
      template<> boost::cstdfloat::detail::float_internal128_t abs <boost::cstdfloat::detail::float_internal128_t>(const complex<boost::cstdfloat::detail::float_internal128_t>& x) { return std::sqrt((real(x) * real(x)) + (imag(x) * imag(x))); }
      template<> boost::cstdfloat::detail::float_internal128_t arg <boost::cstdfloat::detail::float_internal128_t>(const complex<boost::cstdfloat::detail::float_internal128_t>& x) { return std::atan2(x.imag(), x.real()); }
      template<> boost::cstdfloat::detail::float_internal128_t norm<boost::cstdfloat::detail::float_internal128_t>(const complex<boost::cstdfloat::detail::float_internal128_t>& x) { return (real(x) * real(x)) + (imag(x) * imag(x)); }

      template<> complex<boost::cstdfloat::detail::float_internal128_t> conj (const complex<boost::cstdfloat::detail::float_internal128_t>& x) { return complex<boost::cstdfloat::detail::float_internal128_t>(x.real(), -x.imag()); }
      #if !defined(BOOST_NO_CXX11_CONSTEXPR)
      template<> complex<boost::cstdfloat::detail::float_internal128_t> proj (const complex<boost::cstdfloat::detail::float_internal128_t>& x)
      {
        const boost::cstdfloat::detail::float_internal128_t two_over_denom = BOOST_FLOAT128_C(2.0) / (std::norm(x) + 1);

        return complex<boost::cstdfloat::detail::float_internal128_t>(x.real() * two_over_denom,
                                                                      x.imag() * two_over_denom);
      }
      #endif
      template<> complex<boost::cstdfloat::detail::float_internal128_t> polar(const boost::cstdfloat::detail::float_internal128_t& rho,
                                                                              const boost::cstdfloat::detail::float_internal128_t& theta)      { return complex<boost::cstdfloat::detail::float_internal128_t>(rho * std::cos(theta), rho * std::sin(theta)); }

      // Global add, sub, mul, div.
      template<> complex<boost::cstdfloat::detail::float_internal128_t> operator+(const complex<boost::cstdfloat::detail::float_internal128_t>& u, const complex<boost::cstdfloat::detail::float_internal128_t>& v) { return complex<boost::cstdfloat::detail::float_internal128_t>(u.real() + v.real(), u.imag() + v.imag()); }
      template<> complex<boost::cstdfloat::detail::float_internal128_t> operator-(const complex<boost::cstdfloat::detail::float_internal128_t>& u, const complex<boost::cstdfloat::detail::float_internal128_t>& v) { return complex<boost::cstdfloat::detail::float_internal128_t>(u.real() - v.real(), u.imag() - v.imag()); }

      template<> complex<boost::cstdfloat::detail::float_internal128_t> operator*(const complex<boost::cstdfloat::detail::float_internal128_t>& u, const complex<boost::cstdfloat::detail::float_internal128_t>& v)
      {
        const typename complex<boost::cstdfloat::detail::float_internal128_t>::value_type ur = u.real();
        const typename complex<boost::cstdfloat::detail::float_internal128_t>::value_type ui = u.imag();
        const typename complex<boost::cstdfloat::detail::float_internal128_t>::value_type vr = v.real();
        const typename complex<boost::cstdfloat::detail::float_internal128_t>::value_type vi = v.imag();

        return complex<boost::cstdfloat::detail::float_internal128_t>((ur * vr) - (ui * vi), (ur * vi) + (ui * vr));
      }

      template<> inline complex<boost::cstdfloat::detail::float_internal128_t> operator/(const complex<boost::cstdfloat::detail::float_internal128_t>& u, const complex<boost::cstdfloat::detail::float_internal128_t>& v)
      {
        const boost::cstdfloat::detail::float_internal128_t one_over_denom = 1 / std::norm(v);
        const boost::cstdfloat::detail::float_internal128_t tmp_re = ((u.real() * v.real()) + (u.imag() * v.imag())) * one_over_denom;
        const boost::cstdfloat::detail::float_internal128_t tmp_im = ((u.imag() * v.real()) - (u.real() * v.imag())) * one_over_denom;

        return complex<boost::cstdfloat::detail::float_internal128_t>(tmp_re, tmp_im);
      }

      template<> complex<boost::cstdfloat::detail::float_internal128_t> operator+(const complex<boost::cstdfloat::detail::float_internal128_t>& u, const boost::cstdfloat::detail::float_internal128_t& v) { return complex<boost::cstdfloat::detail::float_internal128_t>(u.real() + v, u.imag()); }
      template<> complex<boost::cstdfloat::detail::float_internal128_t> operator-(const complex<boost::cstdfloat::detail::float_internal128_t>& u, const boost::cstdfloat::detail::float_internal128_t& v) { return complex<boost::cstdfloat::detail::float_internal128_t>(u.real() - v, u.imag()); }
      template<> complex<boost::cstdfloat::detail::float_internal128_t> operator*(const complex<boost::cstdfloat::detail::float_internal128_t>& u, const boost::cstdfloat::detail::float_internal128_t& v) { return complex<boost::cstdfloat::detail::float_internal128_t>(u.real() * v, u.imag() * v); }
      template<> complex<boost::cstdfloat::detail::float_internal128_t> operator/(const complex<boost::cstdfloat::detail::float_internal128_t>& u, const boost::cstdfloat::detail::float_internal128_t& v) { return complex<boost::cstdfloat::detail::float_internal128_t>(u.real() / v, u.imag() / v); }

      template<> complex<boost::cstdfloat::detail::float_internal128_t> operator+(const boost::cstdfloat::detail::float_internal128_t& u, const complex<boost::cstdfloat::detail::float_internal128_t>& v) { return complex<boost::cstdfloat::detail::float_internal128_t>(u + v.real(), v.imag()); }
      template<> complex<boost::cstdfloat::detail::float_internal128_t> operator-(const boost::cstdfloat::detail::float_internal128_t& u, const complex<boost::cstdfloat::detail::float_internal128_t>& v) { return complex<boost::cstdfloat::detail::float_internal128_t>(u - v.real(), -v.imag()); }
      template<> complex<boost::cstdfloat::detail::float_internal128_t> operator*(const boost::cstdfloat::detail::float_internal128_t& u, const complex<boost::cstdfloat::detail::float_internal128_t>& v) { return complex<boost::cstdfloat::detail::float_internal128_t>(u * v.real(), u * v.imag()); }
      template<> complex<boost::cstdfloat::detail::float_internal128_t> operator/(const boost::cstdfloat::detail::float_internal128_t& u, const complex<boost::cstdfloat::detail::float_internal128_t>& v) { const boost::cstdfloat::detail::float_internal128_t v_norm = norm(v); return complex<boost::cstdfloat::detail::float_internal128_t>((u * v.real()) / v_norm, (-u * v.imag()) / v_norm); }

      // Unary plus / minus.
      template<> complex<boost::cstdfloat::detail::float_internal128_t> operator+(const complex<boost::cstdfloat::detail::float_internal128_t>& u) { return u; }
      template<> complex<boost::cstdfloat::detail::float_internal128_t> operator-(const complex<boost::cstdfloat::detail::float_internal128_t>& u) { return complex<boost::cstdfloat::detail::float_internal128_t>(-u.real(), -u.imag()); }

      // Equality and inequality.
      template<> BOOST_CONSTEXPR bool operator==(const complex<boost::cstdfloat::detail::float_internal128_t>& u, const complex<boost::cstdfloat::detail::float_internal128_t>& v) { return ((u.real() == v.real()) && (u.imag() == v.imag())); }
      template<> BOOST_CONSTEXPR bool operator!=(const complex<boost::cstdfloat::detail::float_internal128_t>& u, const complex<boost::cstdfloat::detail::float_internal128_t>& v) { return ((u.real() != v.real()) || (u.imag() != v.imag())); }

      template<> BOOST_CONSTEXPR bool operator==(const complex<boost::cstdfloat::detail::float_internal128_t>& u, const boost::cstdfloat::detail::float_internal128_t& v) { return ((u.real() == v) && (u.imag() == BOOST_FLOAT128_C(0.0))); }
      template<> BOOST_CONSTEXPR bool operator!=(const complex<boost::cstdfloat::detail::float_internal128_t>& u, const boost::cstdfloat::detail::float_internal128_t& v) { return ((u.real() != v) || (u.imag() != BOOST_FLOAT128_C(0.0))); }

      template<> BOOST_CONSTEXPR bool operator==(const boost::cstdfloat::detail::float_internal128_t& u, const complex<boost::cstdfloat::detail::float_internal128_t>& v) { return ((u == v.real()) && (v.imag() == BOOST_FLOAT128_C(0.0))); }
      template<> BOOST_CONSTEXPR bool operator!=(const boost::cstdfloat::detail::float_internal128_t& u, const complex<boost::cstdfloat::detail::float_internal128_t>& v) { return ((u != v.real()) || (v.imag() != BOOST_FLOAT128_C(0.0))); }

      // 26.4.8, transcendentals
      template<> complex<boost::cstdfloat::detail::float_internal128_t> sqrt(const complex<boost::cstdfloat::detail::float_internal128_t>& x)
      {
        // sqrt(*this) = (s, I / 2s) for R >= 0
        // (|I| / 2s, +-s) for R < 0
        // where s = sqrt{ [ |R| + sqrt(R^2 + I^2) ] / 2 },
        // and the +- sign is the same as the sign of I.

        const boost::cstdfloat::detail::float_internal128_t zr = x.real();
        const boost::cstdfloat::detail::float_internal128_t s  = std::sqrt((std::fabs(zr) + std::abs(x)) / 2);
        const boost::cstdfloat::detail::float_internal128_t zi = x.imag();

        if(zr >= 0)
        {
          return complex<boost::cstdfloat::detail::float_internal128_t>(s, (zi / s) / 2);
        }
        else
        {
          const bool imag_is_pos = (zi >= 0);

          return complex<boost::cstdfloat::detail::float_internal128_t>((std::fabs(zi) / s) / 2, (imag_is_pos ? s : -s));
        }
      }

      template<> complex<boost::cstdfloat::detail::float_internal128_t> sin(const complex<boost::cstdfloat::detail::float_internal128_t>& x)
      {
        const boost::cstdfloat::detail::float_internal128_t sin_x  = std::sin (x.real());
        const boost::cstdfloat::detail::float_internal128_t cos_x  = std::cos (x.real());
        const boost::cstdfloat::detail::float_internal128_t sinh_y = std::sinh(x.imag());
        const boost::cstdfloat::detail::float_internal128_t cosh_y = std::cosh(x.imag());

        return complex<boost::cstdfloat::detail::float_internal128_t>(sin_x * cosh_y, cos_x * sinh_y);
      }

      template<> complex<boost::cstdfloat::detail::float_internal128_t> cos(const complex<boost::cstdfloat::detail::float_internal128_t>& x)
      {
        const boost::cstdfloat::detail::float_internal128_t sin_x  = std::sin (x.real());
        const boost::cstdfloat::detail::float_internal128_t cos_x  = std::cos (x.real());
        const boost::cstdfloat::detail::float_internal128_t sinh_y = std::sinh(x.imag());
        const boost::cstdfloat::detail::float_internal128_t cosh_y = std::cosh(x.imag());

        return complex<boost::cstdfloat::detail::float_internal128_t>(cos_x * cosh_y, -(sin_x * sinh_y));
      }

      template<> complex<boost::cstdfloat::detail::float_internal128_t> tan(const complex<boost::cstdfloat::detail::float_internal128_t>& x)
      {
        return std::sin(x) / std::cos(x);
      }

      #if !defined(BOOST_NO_CXX11_CONSTEXPR)
      template<> complex<boost::cstdfloat::detail::float_internal128_t> asin(const complex<boost::cstdfloat::detail::float_internal128_t>& x)
      {
        return -boost::cstdfloat::detail::iz_helper__x(std::log(boost::cstdfloat::detail::iz_helper__x(x) + std::sqrt(BOOST_FLOAT128_C(1.0) - (x * x))));
      }

      template<> complex<boost::cstdfloat::detail::float_internal128_t> acos(const complex<boost::cstdfloat::detail::float_internal128_t>& x)
      {
        return BOOST_FLOAT128_C(1.57079632679489661923132169163975144209858469968755) - std::asin(x);
      }

      template<> complex<boost::cstdfloat::detail::float_internal128_t> atan(const complex<boost::cstdfloat::detail::float_internal128_t>& x)
      {
        const complex<boost::cstdfloat::detail::float_internal128_t> izz = boost::cstdfloat::detail::iz_helper__x(x);

        return boost::cstdfloat::detail::iz_helper__x(std::log(BOOST_FLOAT128_C(1.0) - izz) - std::log(BOOST_FLOAT128_C(1.0) + izz)) / BOOST_FLOAT128_C(2.0);
      }
      #endif

      template<> complex<boost::cstdfloat::detail::float_internal128_t> exp(const complex<boost::cstdfloat::detail::float_internal128_t>& x)
      {
        return std::polar(std::exp(x.real()), x.imag());
      }

      template<> complex<boost::cstdfloat::detail::float_internal128_t> log(const complex<boost::cstdfloat::detail::float_internal128_t>& x)
      {
        return complex<boost::cstdfloat::detail::float_internal128_t>(std::log(std::norm(x)) / 2, std::atan2(x.imag(), x.real()));
      }

      template<> complex<boost::cstdfloat::detail::float_internal128_t> log10(const complex<boost::cstdfloat::detail::float_internal128_t>& x)
      {
        return std::log(x) / BOOST_FLOAT128_C(2.30258509299404568401799145468436420760110148862877);
      }

      template<> complex<boost::cstdfloat::detail::float_internal128_t> pow(const complex<boost::cstdfloat::detail::float_internal128_t>& x,
                                                                            const boost::cstdfloat::detail::float_internal128_t& a)
      {
        std::exp(a * std::log(x));
      }

      template<> complex<boost::cstdfloat::detail::float_internal128_t> pow(const complex<boost::cstdfloat::detail::float_internal128_t>& x,
                                                                            const complex<boost::cstdfloat::detail::float_internal128_t>& a)
      {
        std::exp(a * std::log(x));
      }

      template<> complex<boost::cstdfloat::detail::float_internal128_t> pow(const boost::cstdfloat::detail::float_internal128_t& x,
                                                                            const complex<boost::cstdfloat::detail::float_internal128_t>& a)
      {
        std::exp(a * std::log(x));
      }

      template<> complex<boost::cstdfloat::detail::float_internal128_t> sinh(const complex<boost::cstdfloat::detail::float_internal128_t>& x)
      {
        const boost::cstdfloat::detail::float_internal128_t sin_y  = std::sin (x.imag());
        const boost::cstdfloat::detail::float_internal128_t cos_y  = std::cos (x.imag());
        const boost::cstdfloat::detail::float_internal128_t sinh_x = std::sinh(x.real());
        const boost::cstdfloat::detail::float_internal128_t cosh_x = std::cosh(x.real());

        return complex<boost::cstdfloat::detail::float_internal128_t>(cos_y * sinh_x, cosh_x * sin_y);
      }

      template<> complex<boost::cstdfloat::detail::float_internal128_t> cosh(const complex<boost::cstdfloat::detail::float_internal128_t>& x)
      {
        const boost::cstdfloat::detail::float_internal128_t sin_y  = std::sin (x.imag());
        const boost::cstdfloat::detail::float_internal128_t cos_y  = std::cos (x.imag());
        const boost::cstdfloat::detail::float_internal128_t sinh_x = std::sinh(x.real());
        const boost::cstdfloat::detail::float_internal128_t cosh_x = std::cosh(x.real());

        return complex<boost::cstdfloat::detail::float_internal128_t>(cos_y * cosh_x, sin_y * sinh_x);
      }

      template<> complex<boost::cstdfloat::detail::float_internal128_t> tanh(const complex<boost::cstdfloat::detail::float_internal128_t>& x)
      {
        return std::sinh(x) / std::cosh(x);
      }

      #if !defined(BOOST_NO_CXX11_CONSTEXPR)
      template<> complex<boost::cstdfloat::detail::float_internal128_t> asinh(const complex<boost::cstdfloat::detail::float_internal128_t>& x)
      {
        return std::log(x + std::sqrt((x * x) + BOOST_FLOAT128_C(1.0)));
      }

      template<> complex<boost::cstdfloat::detail::float_internal128_t> acosh(const complex<boost::cstdfloat::detail::float_internal128_t>& x)
      {
        const complex<boost::cstdfloat::detail::float_internal128_t> zp(x.real() + 1, x.imag());
        const complex<boost::cstdfloat::detail::float_internal128_t> zm(x.real() - 1, x.imag());

        return std::log(x + (zp * std::sqrt(zm / zp)));
      }

      template<> complex<boost::cstdfloat::detail::float_internal128_t> atanh(const complex<boost::cstdfloat::detail::float_internal128_t>& x)
      {
        return (std::log(BOOST_FLOAT128_C(1.0) + x) - std::log(BOOST_FLOAT128_C(1.0) - x)) / BOOST_FLOAT128_C(2.0);
      }
      #endif

      template<class char_type, class traits_type>
      inline std::basic_ostream<char_type, traits_type>& operator<< (std::basic_ostream<char_type, traits_type>& os, const std::complex<boost::cstdfloat::detail::float_internal128_t>&);

      template<class char_type, class traits_type>
      inline std::basic_ostream<char_type, traits_type>& operator<< (std::basic_ostream<char_type, traits_type>& os, const std::complex<boost::cstdfloat::detail::float_internal128_t>& x)
      {
        std::basic_ostringstream<char_type, traits_type> ostr;
        ostr.flags(os.flags());
        ostr.imbue(os.getloc());
        ostr.precision(os.precision());

        ostr << '(' << x.real() << ',' << x.imag() << ')';

        return (os << ostr.str());
      }
    } // namespace std

    #endif // Not BOOST_CSTDFLOAT_NO_LIBQUADMATH_COMPLEX (i.e., has libquadmath complex support)

  #endif // Not BOOST_CSTDFLOAT_NO_LIBQUADMATH_SUPPORT (i.e., has libquadmath support)

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
    #undef  BOOST_NO_FLOAT128_T

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
    #undef  BOOST_FLOAT_16_MIN
    #undef  BOOST_FLOAT_16_MAX
  #elif(BOOST_CSTDFLOAT_MAXIMUM_AVAILABLE_WIDTH == 32)
    #define BOOST_FLOATMAX_C(x) BOOST_FLOAT32_C(x)
    #define BOOST_FLOATMAX_MIN  BOOST_FLOAT_32_MIN
    #define BOOST_FLOATMAX_MAX  BOOST_FLOAT_32_MAX
    #undef  BOOST_FLOAT_32_MIN
    #undef  BOOST_FLOAT_32_MAX
  #elif(BOOST_CSTDFLOAT_MAXIMUM_AVAILABLE_WIDTH == 64)
    #define BOOST_FLOATMAX_C(x) BOOST_FLOAT64_C(x)
    #define BOOST_FLOATMAX_MIN  BOOST_FLOAT_64_MIN
    #define BOOST_FLOATMAX_MAX  BOOST_FLOAT_64_MAX
    #undef  BOOST_FLOAT_64_MIN
    #undef  BOOST_FLOAT_64_MAX
  #elif(BOOST_CSTDFLOAT_MAXIMUM_AVAILABLE_WIDTH == 80)
    #define BOOST_FLOATMAX_C(x) BOOST_FLOAT80_C(x)
    #define BOOST_FLOATMAX_MIN  BOOST_FLOAT_80_MIN
    #define BOOST_FLOATMAX_MAX  BOOST_FLOAT_80_MAX
    #undef  BOOST_FLOAT_80_MIN
    #undef  BOOST_FLOAT_80_MAX
  #elif(BOOST_CSTDFLOAT_MAXIMUM_AVAILABLE_WIDTH == 128)
    #define BOOST_FLOATMAX_C(x) BOOST_FLOAT128_C(x)
    #define BOOST_FLOATMAX_MIN  BOOST_FLOAT_128_MIN
    #define BOOST_FLOATMAX_MAX  BOOST_FLOAT_128_MAX
    #undef  BOOST_FLOAT_128_MIN
    #undef  BOOST_FLOAT_128_MAX
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

      #undef BOOST_CSTDFLOAT_HAS_FLOAT16_NATIVE_TYPE
    #endif

    #if(BOOST_CSTDFLOAT_HAS_FLOAT32_NATIVE_TYPE == 1)
      typedef BOOST_CSTDFLOAT_FLOAT32_NATIVE_TYPE float32_t;
      typedef boost::float32_t float_fast32_t;
      typedef boost::float32_t float_least32_t;

      BOOST_STATIC_ASSERT(std::numeric_limits<boost::float32_t>::is_iec559    == true);
      BOOST_STATIC_ASSERT(std::numeric_limits<boost::float32_t>::radix        ==    2);
      BOOST_STATIC_ASSERT(std::numeric_limits<boost::float32_t>::digits       ==   24);
      BOOST_STATIC_ASSERT(std::numeric_limits<boost::float32_t>::max_exponent ==  128);

      #undef BOOST_CSTDFLOAT_HAS_FLOAT32_NATIVE_TYPE
    #endif

    #if(BOOST_CSTDFLOAT_HAS_FLOAT64_NATIVE_TYPE == 1)
      typedef BOOST_CSTDFLOAT_FLOAT64_NATIVE_TYPE float64_t;
      typedef boost::float64_t float_fast64_t;
      typedef boost::float64_t float_least64_t;

      BOOST_STATIC_ASSERT(std::numeric_limits<boost::float64_t>::is_iec559    == true);
      BOOST_STATIC_ASSERT(std::numeric_limits<boost::float64_t>::radix        ==    2);
      BOOST_STATIC_ASSERT(std::numeric_limits<boost::float64_t>::digits       ==   53);
      BOOST_STATIC_ASSERT(std::numeric_limits<boost::float64_t>::max_exponent == 1024);

      #undef BOOST_CSTDFLOAT_HAS_FLOAT64_NATIVE_TYPE
    #endif

    #if(BOOST_CSTDFLOAT_HAS_FLOAT80_NATIVE_TYPE == 1)
      typedef BOOST_CSTDFLOAT_FLOAT80_NATIVE_TYPE float80_t;
      typedef boost::float80_t float_fast80_t;
      typedef boost::float80_t float_least80_t;

      BOOST_STATIC_ASSERT(std::numeric_limits<boost::float80_t>::is_iec559    ==  true);
      BOOST_STATIC_ASSERT(std::numeric_limits<boost::float80_t>::radix        ==     2);
      BOOST_STATIC_ASSERT(std::numeric_limits<boost::float80_t>::digits       ==    64);
      BOOST_STATIC_ASSERT(std::numeric_limits<boost::float80_t>::max_exponent == 16384);

      #undef BOOST_CSTDFLOAT_HAS_FLOAT80_NATIVE_TYPE
    #endif

    #if(BOOST_CSTDFLOAT_HAS_FLOAT128_NATIVE_TYPE == 1)
      typedef BOOST_CSTDFLOAT_FLOAT128_NATIVE_TYPE float128_t;
      typedef boost::float128_t float_fast128_t;
      typedef boost::float128_t float_least128_t;

      #if defined(BOOST_MATH_USE_FLOAT128) && !defined(BOOST_CSTDFLOAT_NO_LIBQUADMATH_SUPPORT) && defined(BOOST_CSTDFLOAT_NO_LIBQUADMATH_LIMITS)
      // This configuration does not support std::numeric_limits<boost::float128_t>.
      // So we can not query the values of std::numeric_limits<boost::float128_t>
      // with static assertions.
      #else
      BOOST_STATIC_ASSERT(std::numeric_limits<boost::float128_t>::is_iec559    ==  true);
      BOOST_STATIC_ASSERT(std::numeric_limits<boost::float128_t>::radix        ==     2);
      BOOST_STATIC_ASSERT(std::numeric_limits<boost::float128_t>::digits       ==   113);
      BOOST_STATIC_ASSERT(std::numeric_limits<boost::float128_t>::max_exponent == 16384);
      #endif

      #undef BOOST_CSTDFLOAT_HAS_FLOAT128_NATIVE_TYPE
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

    #undef BOOST_CSTDFLOAT_MAXIMUM_AVAILABLE_WIDTH
  }
  // namespace boost

#endif // _BOOST_CSTDFLOAT_2014_01_09_HPP_
