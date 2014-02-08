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
  #include <limits>
  #include <exception>
  #include <boost/throw_exception.hpp>
  #include <boost/static_assert.hpp>
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

  // Check if __float128 from GCC's libquadmath or if (potentially)
  // ICC's /Qlong-double flag is supported.
  // TODO: Should we allow BOOST_MATH_USE_FLOAT128 for ICC?
  // Here, we use the BOOST_MATH_USE_FLOAT128 pre-processor
  // definition from <boost/math/tools/config.hpp>.
  #if (BOOST_CSTDFLOAT_HAS_FLOAT128_NATIVE_TYPE == 0) && defined(BOOST_MATH_USE_FLOAT128) /*&& defined(FLT128_MIN) && defined(FLT128_MAX) && defined(FLT128_EPSILON) && defined(FLT128_MANT_DIG)*/
    #define BOOST_CSTDFLOAT_FLOAT128_NATIVE_TYPE __float128
    #undef  BOOST_CSTDFLOAT_MAXIMUM_AVAILABLE_WIDTH
    #define BOOST_CSTDFLOAT_MAXIMUM_AVAILABLE_WIDTH 128
    #undef  BOOST_CSTDFLOAT_HAS_FLOAT128_NATIVE_TYPE
    #define BOOST_CSTDFLOAT_HAS_FLOAT128_NATIVE_TYPE  1
    #define BOOST_CSTDFLOAT_FLOAT128_MIN  3.36210314311209350626267781732175260e-4932Q
    #define BOOST_CSTDFLOAT_FLOAT128_MAX  1.18973149535723176508575932662800702e+4932Q
    #define BOOST_CSTDFLOAT_FLOAT128_EPS  1.92592994438723585305597794258492732e-0034Q
    #define BOOST_FLOAT128_C(x)  (x ## Q)

    #if !defined(BOOST_CSTDFLOAT_NO_GCC_FLOAT128_LIMITS)

    // For __float128, implement a specialization of std::numeric_limits<>.

    // Forward declaration of quad square root function.
    extern "C" ::__float128 sqrtq(::__float128);

    namespace std
    {
      template<>
      class numeric_limits< ::__float128 >
      {
      public:
        BOOST_STATIC_CONSTEXPR bool                  is_specialized    = true;
        static                 ::__float128 (min) () BOOST_NOEXCEPT    { return BOOST_CSTDFLOAT_FLOAT128_MIN; }
        static                 ::__float128 (max) () BOOST_NOEXCEPT    { return BOOST_CSTDFLOAT_FLOAT128_MAX; }
        static                 ::__float128 lowest() BOOST_NOEXCEPT    { return -(max)(); }
        BOOST_STATIC_CONSTEXPR int                   digits            = 113;
        BOOST_STATIC_CONSTEXPR int                   digits10          = 34;
        BOOST_STATIC_CONSTEXPR int                   max_digits10      = 36;
        BOOST_STATIC_CONSTEXPR bool                  is_signed         = true;
        BOOST_STATIC_CONSTEXPR bool                  is_integer        = false;
        BOOST_STATIC_CONSTEXPR bool                  is_exact          = false;
        BOOST_STATIC_CONSTEXPR int                   radix             = 2;
        static                 ::__float128          epsilon    ()     { return BOOST_CSTDFLOAT_FLOAT128_EPS; }
        static                 ::__float128          round_error()     { return BOOST_FLOAT128_C(0.5); }
        BOOST_STATIC_CONSTEXPR int                   min_exponent      = -16381;
        BOOST_STATIC_CONSTEXPR int                   min_exponent10    = static_cast<int>((min_exponent * 301L) / 1000L);
        BOOST_STATIC_CONSTEXPR int                   max_exponent      = +16384;
        BOOST_STATIC_CONSTEXPR int                   max_exponent10    = static_cast<int>((max_exponent * 301L) / 1000L);
        BOOST_STATIC_CONSTEXPR bool                  has_infinity      = true;
        BOOST_STATIC_CONSTEXPR bool                  has_quiet_NaN     = true;
        BOOST_STATIC_CONSTEXPR bool                  has_signaling_NaN = false;
        BOOST_STATIC_CONSTEXPR float_denorm_style    has_denorm        = denorm_absent;
        BOOST_STATIC_CONSTEXPR bool                  has_denorm_loss   = false;
        static                 ::__float128          infinity     ()   { return BOOST_FLOAT128_C(1.0) / BOOST_FLOAT128_C(0.0); }
        static                 ::__float128          quiet_NaN    ()   { return ::sqrtq(BOOST_FLOAT128_C(-1.0)); }
        static                 ::__float128          signaling_NaN()   { return BOOST_FLOAT128_C(0.0); }
        static                 ::__float128          denorm_min   ()   { return BOOST_FLOAT128_C(0.0); }
        BOOST_STATIC_CONSTEXPR bool                  is_iec559         = true;
        BOOST_STATIC_CONSTEXPR bool                  is_bounded        = false;
        BOOST_STATIC_CONSTEXPR bool                  is_modulo         = false;
        BOOST_STATIC_CONSTEXPR bool                  traps             = false;
        BOOST_STATIC_CONSTEXPR bool                  tinyness_before   = false;
        BOOST_STATIC_CONSTEXPR float_round_style     round_style       = round_to_nearest;
      };
    }
    #endif

    #if !defined(BOOST_CSTDFLOAT_NO_GCC_FLOAT128_CMATH)

      // For __float128, implement <math.h> functions in boost::cstdfloat::detail.
      // Subsequently *use* these in the global namespace.

      // Begin with some forward function declarations.

      // Forward declaration of quad string print function.
      extern "C" int quadmath_snprintf(char *str, size_t size, const char *format, ...);

      // Forward declarations of quad elementary functions.
      extern "C" ::__float128 ldexpq (::__float128, int);
      extern "C" ::__float128 frexpq (::__float128, int*);
      extern "C" ::__float128 fabsq  (::__float128);
      extern "C" ::__float128 floorq (::__float128);
      extern "C" ::__float128 ceilq  (::__float128);
      extern "C" ::__float128 sqrtq  (::__float128);
      extern "C" ::__float128 truncq (::__float128);
      extern "C" ::__float128 expq   (::__float128);
      extern "C" ::__float128 powq   (::__float128, ::__float128);
      extern "C" ::__float128 logq   (::__float128);
      extern "C" ::__float128 log10q (::__float128);
      extern "C" ::__float128 sinq   (::__float128);
      extern "C" ::__float128 cosq   (::__float128);
      extern "C" ::__float128 tanq   (::__float128);
      extern "C" ::__float128 asinq  (::__float128);
      extern "C" ::__float128 acosq  (::__float128);
      extern "C" ::__float128 atanq  (::__float128);
      extern "C" ::__float128 sinhq  (::__float128);
      extern "C" ::__float128 coshq  (::__float128);
      extern "C" ::__float128 tanhq  (::__float128);
      extern "C" ::__float128 fmodq  (::__float128, ::__float128);
      extern "C" ::__float128 atan2q (::__float128, ::__float128);
      extern "C" ::__float128 lgammaq(::__float128);
      extern "C" ::__float128 tgammaq(::__float128);

      // Put the __float128 <math.h> functions in boost::cstdfloat.
      namespace boost { namespace cstdfloat { namespace detail {
      inline   ::__float128 ldexp (::__float128 x, int n)             { return ::ldexpq (x, n); }
      inline   ::__float128 frexp (::__float128 x, int* pn)           { return ::frexpq (x, pn); }
      inline   ::__float128 fabs  (::__float128 x)                    { return ::fabsq  (x); }
      inline   ::__float128 floor (::__float128 x)                    { return ::floorq (x); }
      inline   ::__float128 ceil  (::__float128 x)                    { return ::ceilq  (x); }
      inline   ::__float128 sqrt  (::__float128 x)                    { return ::sqrtq  (x); }
      inline   ::__float128 trunc (::__float128 x)                    { return ::truncq (x); }
      inline   ::__float128 exp   (::__float128 x)                    { return ::expq   (x); }
      inline   ::__float128 pow   (::__float128 x, ::__float128 a )   { return ::powq   (x, a); }
      inline   ::__float128 log   (::__float128 x)                    { return ::logq   (x); }
      inline   ::__float128 log10 (::__float128 x)                    { return ::log10q (x); }
      inline   ::__float128 sin   (::__float128 x)                    { return ::sinq   (x); }
      inline   ::__float128 cos   (::__float128 x)                    { return ::cosq   (x); }
      inline   ::__float128 tan   (::__float128 x)                    { return ::tanq   (x); }
      inline   ::__float128 asin  (::__float128 x)                    { return ::asinq  (x); }
      inline   ::__float128 acos  (::__float128 x)                    { return ::acosq  (x); }
      inline   ::__float128 atan  (::__float128 x)                    { return ::atanq  (x); }
      inline   ::__float128 sinh  (::__float128 x)                    { return ::sinhq  (x); }
      inline   ::__float128 cosh  (::__float128 x)                    { return ::coshq  (x); }
      inline   ::__float128 tanh  (::__float128 x)                    { return ::tanhq  (x); }
      inline   ::__float128 fmod  (::__float128 a, ::__float128 b )   { return ::fmodq  (a, b); }
      inline   ::__float128 atan2 (::__float128 y, ::__float128 x )   { return ::atan2q (y, x); }
      inline   ::__float128 lgamma(::__float128 x)                    { return ::lgammaq(x); }
      inline   ::__float128 tgamma(::__float128 x)                    { return ::tgammaq(x); }
      } } }  // boost::cstdfloat::detail

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

      // For __float128, implement <cmath> functions in the std namespace.
      namespace std
      {
        inline ::__float128 ldexp (::__float128 x, int n)             { return boost::cstdfloat::detail::ldexp (x, n); }
        inline ::__float128 frexp (::__float128 x, int* pn)           { return boost::cstdfloat::detail::frexp (x, pn); }
        inline ::__float128 fabs  (::__float128 x)                    { return boost::cstdfloat::detail::fabs  (x); }
        inline ::__float128 floor (::__float128 x)                    { return boost::cstdfloat::detail::floor (x); }
        inline ::__float128 ceil  (::__float128 x)                    { return boost::cstdfloat::detail::ceil  (x); }
        inline ::__float128 sqrt  (::__float128 x)                    { return boost::cstdfloat::detail::sqrt  (x); }
        inline ::__float128 trunc (::__float128 x)                    { return boost::cstdfloat::detail::trunc (x); }
        inline ::__float128 exp   (::__float128 x)                    { return boost::cstdfloat::detail::exp   (x); }
        inline ::__float128 pow   (::__float128 x, ::__float128 a )   { return boost::cstdfloat::detail::pow   (x, a); }
        inline ::__float128 log   (::__float128 x)                    { return boost::cstdfloat::detail::log   (x); }
        inline ::__float128 log10 (::__float128 x)                    { return boost::cstdfloat::detail::log10 (x); }
        inline ::__float128 sin   (::__float128 x)                    { return boost::cstdfloat::detail::sin   (x); }
        inline ::__float128 cos   (::__float128 x)                    { return boost::cstdfloat::detail::cos   (x); }
        inline ::__float128 tan   (::__float128 x)                    { return boost::cstdfloat::detail::tan   (x); }
        inline ::__float128 asin  (::__float128 x)                    { return boost::cstdfloat::detail::asin  (x); }
        inline ::__float128 acos  (::__float128 x)                    { return boost::cstdfloat::detail::acos  (x); }
        inline ::__float128 atan  (::__float128 x)                    { return boost::cstdfloat::detail::atan  (x); }
        inline ::__float128 sinh  (::__float128 x)                    { return boost::cstdfloat::detail::sinh  (x); }
        inline ::__float128 cosh  (::__float128 x)                    { return boost::cstdfloat::detail::cosh  (x); }
        inline ::__float128 tanh  (::__float128 x)                    { return boost::cstdfloat::detail::tanh  (x); }
        inline ::__float128 fmod  (::__float128 a, ::__float128 b )   { return boost::cstdfloat::detail::fmod  (a, b); }
        inline ::__float128 atan2 (::__float128 y, ::__float128 x )   { return boost::cstdfloat::detail::atan2 (y, x); }
        inline ::__float128 lgamma(::__float128 x )                   { return boost::cstdfloat::detail::lgamma(x); }
        inline ::__float128 tgamma(::__float128 x )                   { return boost::cstdfloat::detail::tgamma(x); }
      }
    #endif

    #if !defined(BOOST_CSTDFLOAT_NO_GCC_FLOAT128_IOSTREAM)

      // For __float128, implement I/O stream operations.

      #include <algorithm>
      #include <cstddef>
      #include <ostream>

      inline std::ostream& operator<< (std::ostream& os, const ::__float128& x)
      {
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
