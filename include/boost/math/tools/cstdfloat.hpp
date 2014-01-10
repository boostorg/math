///////////////////////////////////////////////////////////////////////////////
// Copyright Christopher Kormanyos 2014.
// Copyright John Maddock 2014.
// Copyright Paul Bristow 2014.
// Distributed under the Boost
// Software License, Version 1.0. (See accompanying file
// LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef _BOOST_CSTDFLOAT_2014_01_09_HPP_
  #define _BOOST_CSTDFLOAT_2014_01_09_HPP_

  #include <float.h>
  #include <limits>
  #include <boost/static_assert.hpp>

  // This is the beginning of the preamble.

  // In this preamble, the preprocessor is used to query certain
  // preprocessor definitions in order to detect the presence of
  // built-in floating-point types having specified widths.

  #define BOOST_CSTDFLOAT_MAXIMUM_AVAILABLE_WIDTH  0

  #if defined(FLT_MAX_EXP) && ((FLT_MAX_EXP == 128) && (FLT_RADIX == 2))
  #define BOOST_CSTDFLOAT_FLOAT32_NATIVE_TYPE  float
  #undef  BOOST_CSTDFLOAT_MAXIMUM_AVAILABLE_WIDTH
  #define BOOST_CSTDFLOAT_MAXIMUM_AVAILABLE_WIDTH 32
  #define BOOST_CSTDFLOAT_HAS_FLOAT32_NATIVE_TYPE  1
  #else
  #define BOOST_CSTDFLOAT_HAS_FLOAT32_NATIVE_TYPE  0
  #endif

  #if defined(DBL_MAX_EXP) && ((DBL_MAX_EXP == 1024) && (FLT_RADIX == 2))
  #define BOOST_CSTDFLOAT_FLOAT64_NATIVE_TYPE  double
  #undef  BOOST_CSTDFLOAT_MAXIMUM_AVAILABLE_WIDTH
  #define BOOST_CSTDFLOAT_MAXIMUM_AVAILABLE_WIDTH 64
  #define BOOST_CSTDFLOAT_HAS_FLOAT64_NATIVE_TYPE  1
  #else
  #define BOOST_CSTDFLOAT_HAS_FLOAT64_NATIVE_TYPE  0
  #endif

  #if (0)
  #define BOOST_CSTDFLOAT_FLOAT80_NATIVE_TYPE  long double
  #undef  BOOST_CSTDFLOAT_MAXIMUM_AVAILABLE_WIDTH
  #define BOOST_CSTDFLOAT_MAXIMUM_AVAILABLE_WIDTH 80
  #define BOOST_CSTDFLOAT_HAS_FLOAT80_NATIVE_TYPE  1
  #else
  #define BOOST_CSTDFLOAT_HAS_FLOAT80_NATIVE_TYPE  0
  #endif

  #if (0)
  #define BOOST_CSTDFLOAT_FLOAT128_NATIVE_TYPE  long double
  #undef  BOOST_CSTDFLOAT_MAXIMUM_AVAILABLE_WIDTH
  #define BOOST_CSTDFLOAT_MAXIMUM_AVAILABLE_WIDTH 128
  #define BOOST_CSTDFLOAT_HAS_FLOAT128_NATIVE_TYPE  1
  #else
  #define BOOST_CSTDFLOAT_HAS_FLOAT128_NATIVE_TYPE  0
  #endif

  // This is the end of preamble and the beginning of the typedefs.

  // Here, we define the floating-point typedefs having specified widths
  // based on the proeprocessor analysis in the preamble above.
  // An assessment of IEEE-754 comformance is also carried out.

  #if(BOOST_CSTDFLOAT_HAS_FLOAT32_NATIVE_TYPE == 1)
  typedef BOOST_CSTDFLOAT_FLOAT32_NATIVE_TYPE  float32_t;
  BOOST_STATIC_ASSERT(std::numeric_limits<::float32_t>::is_iec559 == true);
  #endif

  #if(BOOST_CSTDFLOAT_HAS_FLOAT64_NATIVE_TYPE == 1)
  typedef BOOST_CSTDFLOAT_FLOAT64_NATIVE_TYPE  float64_t;
  BOOST_STATIC_ASSERT(std::numeric_limits<::float64_t>::is_iec559 == true);
  #endif

  #if(BOOST_CSTDFLOAT_HAS_FLOAT80_NATIVE_TYPE == 1)
  typedef BOOST_CSTDFLOAT_FLOAT80_NATIVE_TYPE  float80_t;
  BOOST_STATIC_ASSERT(std::numeric_limits<::float80_t>::is_iec559 == true);
  #endif

  #if(BOOST_CSTDFLOAT_HAS_FLOAT128_NATIVE_NATIVE_TYPE == 1)
  typedef BOOST_CSTDFLOAT_FLOAT128_NATIVE_NATIVE_TYPE float128_t;
  BOOST_STATIC_ASSERT(std::numeric_limits<::float128_t>::is_iec559 == true);
  #endif

  #if  (BOOST_CSTDFLOAT_MAXIMUM_AVAILABLE_WIDTH == 32)
  typedef float32_t floatmax_t;
  #elif(BOOST_CSTDFLOAT_MAXIMUM_AVAILABLE_WIDTH == 64)
  typedef float64_t floatmax_t;
  #elif(BOOST_CSTDFLOAT_MAXIMUM_AVAILABLE_WIDTH == 80)
  typedef float80_t floatmax_t;
  #elif(BOOST_CSTDFLOAT_MAXIMUM_AVAILABLE_WIDTH == 128)
  typedef float128_t floatmax_t;
  #else
  #error undefined floating-point width for cstdfloat
  #endif

  // Here, we inject the floating-point typedefs having
  // specified widths into the namespace boost.
  namespace boost
  {
    #if(BOOST_CSTDFLOAT_HAS_FLOAT32_NATIVE_TYPE == 1)
    using ::float32_t;
    #endif

    #if(BOOST_CSTDFLOAT_HAS_FLOAT64_NATIVE_TYPE == 1)
    using ::float64_t;
    #endif

    #if(BOOST_CSTDFLOAT_HAS_FLOAT80_NATIVE_TYPE == 1)
    using ::float80_t;
    #endif

    #if(BOOST_CSTDFLOAT_HAS_FLOAT128_NATIVE_TYPE == 1)
    using ::float128_t;
    #endif

    #if(BOOST_CSTDFLOAT_MAXIMUM_AVAILABLE_WIDTH != 0)
    using ::floatmax_t;
    #endif
  }

  // Here, we define macros that are used for initializing floating-point
  // literal values for the floating-point typedefs having specified widths.

  #if(BOOST_CSTDFLOAT_HAS_FLOAT32_NATIVE_TYPE == 1)
  #define FLOAT32_C(x)  (x ## F)
  #endif

  #if(BOOST_CSTDFLOAT_HAS_FLOAT64_NATIVE_TYPE == 1)
  #define FLOAT64_C(x)  (x)
  #endif

  #if(BOOST_CSTDFLOAT_HAS_FLOAT80_NATIVE_TYPE == 1)
  #define FLOAT80_C(x)  (x ## L)
  #endif

  #if(BOOST_CSTDFLOAT_HAS_FLOAT128_NATIVE_TYPE == 1)
  #define FLOAT128_C(x) (x ## L)
  #endif

#endif // _BOOST_CSTDFLOAT_2014_01_09_HPP_
