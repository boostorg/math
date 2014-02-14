//  Copyright Christopher Kormanyos 2014.
//  Copyright John Maddock 2014.
//  Copyright Paul A. Bristow 2014.

//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cmath>
#include <complex>
#include <limits>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <boost/cstdfloat.hpp>

#ifdef _MSC_VER
#  pragma warning(disable : 4127) // conditional expression is constant.
#  pragma warning(disable : 4512) // assignment operator could not be generated.
#  pragma warning(disable : 4996) // use -D_SCL_SECURE_NO_WARNINGS.
#endif

#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp> // Boost.Test
#include <boost/test/floating_point_comparison.hpp>

//
// DESCRIPTION:
// ~~~~~~~~~~~~
//
// This file tests the implements floating-point typedefs having
// specified widths, as described in N3626 (proposed for C++14).

// Spot tests of boost::float32_t and boost::float64_t, and where available
// boost::float80_t and boost::float128_t. Check the formal behavior of
// the types. Also check selected values of <cmath> and <complex> functions
// for boost::float128_t.

// For more information on <boost/cstdfloat.hpp> and the corresponding
// proposal, see "Floating-Point Typedefs Having Specified Widths".
// as described in N3626 (proposed for C++14).
// See: http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2013/n3626.pdf


#define TEST_CSTDFLOAT_SANITY_CHECK(the_digits)                                                 \
void test_cstdfloat_sanity_check_##the_digits##_func()                                          \
{                                                                                               \
  typedef boost::float##the_digits##_t float_type;                                              \
                                                                                                \
  BOOST_CONSTEXPR_OR_CONST int my_digits10 = std::numeric_limits<float_type>::digits10;         \
                                                                                                \
  {                                                                                             \
    BOOST_CONSTEXPR_OR_CONST float_type x = BOOST_FLOAT##the_digits##_C(1.0) / 3;               \
    std::stringstream ss;                                                                       \
    ss << std::setprecision(my_digits10)                                                        \
       << x;                                                                                    \
    std::string str = "0.";                                                                     \
    str += std::string(std::string::size_type(my_digits10), char('3'));                         \
    BOOST_CHECK_EQUAL( ss.str(), str );                                                         \
  }                                                                                             \
  {                                                                                             \
    BOOST_CONSTEXPR_OR_CONST float_type x = BOOST_FLOAT##the_digits##_C(2.0) / 3;               \
    std::stringstream ss;                                                                       \
    ss << std::setprecision(my_digits10)                                                        \
       << x;                                                                                    \
    std::string str = "0.";                                                                     \
    str += std::string(std::string::size_type(my_digits10 - 1), char('6'));                     \
    str += "7";                                                                                 \
    BOOST_CHECK_EQUAL( ss.str(), str );                                                         \
  }                                                                                             \
  {                                                                                             \
    BOOST_CONSTEXPR_OR_CONST float_type x = BOOST_FLOAT##the_digits##_C(1.0) / 0;               \
    BOOST_CHECK_EQUAL( x, std::numeric_limits<float_type>::infinity );                          \
  }                                                                                             \
  {                                                                                             \
    using std::sqrt;                                                                            \
    BOOST_CONSTEXPR_OR_CONST float_type x = sqrt(BOOST_FLOAT##the_digits##_C(-1.0);             \
    BOOST_CHECK_EQUAL( bool(x == x), false );                                                   \
  }                                                                                             \
}

namespace test_cstdfloat
{
  #if defined(BOOST_FLOAT16_C)
  TEST_CSTDFLOAT_SANITY_CHECK(16)
  #endif

  #if defined(BOOST_FLOAT32_C)
  TEST_CSTDFLOAT_SANITY_CHECK(32)
  #endif

  #if defined(BOOST_FLOAT64_C)
  TEST_CSTDFLOAT_SANITY_CHECK(64)
  #endif

  #if defined(BOOST_FLOAT80_C)
  TEST_CSTDFLOAT_SANITY_CHECK(80)
  #endif

  #if defined(BOOST_FLOAT128_C)
  TEST_CSTDFLOAT_SANITY_CHECK(128)
  #endif

  #if defined(BOOST_FLOAT128_C)
  void extend_check_128_func()
  {
  }
  #endif

  #if defined(BOOST_FLOATMAX_C)
  BOOST_CONSTEXPR_OR_CONST int has_floatmax_t = 1;
  #else
  BOOST_CONSTEXPR_OR_CONST int has_floatmax_t = 0;
  #endif
}

BOOST_AUTO_TEST_CASE(test_main)
{
  // Perform basic sanity checks that verify both the existence of the proper
  // macros as well as the correct digit handling for a given floating-point
  // typedef having specified width.
  BOOST_CHECK_EQUAL( test_cstdfloat::has_floatmax_t, 1 );

  #if defined(BOOST_FLOAT16_C)
  test_cstdfloat::sanity_check_16_func();
  #endif

  #if defined(BOOST_FLOAT32_C)
  test_cstdfloat::sanity_check_32_func();
  #endif

  #if defined(BOOST_FLOAT64_C)
  test_cstdfloat::sanity_check_64_func();
  #endif

  #if defined(BOOST_FLOAT80_C)
  test_cstdfloat::sanity_check_80_func();
  #endif

  #if defined(BOOST_FLOAT128_C)
  test_cstdfloat::sanity_check_128_func();
  #endif

  // Perform an extended check of boost::float128_t including a variety
  // of functions from the C++ standard library.
  #if defined(BOOST_FLOAT128_C)
  test_cstdfloat::extend_check_128_func();
  #endif
}
