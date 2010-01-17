// Copyright John Maddock 2008

// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/math/concepts/real_concept.hpp> // for real_concept
#include <boost/math/special_functions/sign.hpp>

#include <boost/test/test_exec_monitor.hpp> // Boost.Test
#include <boost/test/results_collector.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <iostream>
   using std::cout;
   using std::endl;
   using std::setprecision;

template <class RealType>
void test_spots(RealType /*T*/, const char* /*type_name*/)
{
   // Basic sanity checks.
   RealType a = 0;
   RealType b = 1;
   RealType c = -1;
   BOOST_CHECK_EQUAL((boost::math::signbit)(a), 0);
   BOOST_CHECK_EQUAL((boost::math::sign)(a), 0);
   BOOST_CHECK_EQUAL((boost::math::copysign)(b, a), RealType(1));
   BOOST_CHECK_EQUAL((boost::math::copysign)(c, a), RealType(1));
   a = 1;
   BOOST_CHECK_EQUAL((boost::math::signbit)(a), 0);
   BOOST_CHECK_EQUAL((boost::math::sign)(a), 1);
   BOOST_CHECK_EQUAL((boost::math::copysign)(b, a), RealType(1));
   BOOST_CHECK_EQUAL((boost::math::copysign)(c, a), RealType(1));
   a = -1;
   BOOST_CHECK((boost::math::signbit)(a) != 0);
   BOOST_CHECK_EQUAL((boost::math::sign)(a), -1);
   BOOST_CHECK_EQUAL((boost::math::copysign)(b, a), RealType(-1));
   BOOST_CHECK_EQUAL((boost::math::copysign)(c, a), RealType(-1));
   a = boost::math::tools::max_value<RealType>();
   BOOST_CHECK_EQUAL((boost::math::signbit)(a), 0);
   BOOST_CHECK_EQUAL((boost::math::sign)(a), 1);
   BOOST_CHECK_EQUAL((boost::math::copysign)(b, a), RealType(1));
   BOOST_CHECK_EQUAL((boost::math::copysign)(c, a), RealType(1));
   a = -boost::math::tools::max_value<RealType>();
   BOOST_CHECK((boost::math::signbit)(a) != 0);
   BOOST_CHECK_EQUAL((boost::math::sign)(a), -1);
   BOOST_CHECK_EQUAL((boost::math::copysign)(b, a), RealType(-1));
   BOOST_CHECK_EQUAL((boost::math::copysign)(c, a), RealType(-1));
   if(std::numeric_limits<RealType>::has_infinity)
   {
      a = std::numeric_limits<RealType>::infinity();
      BOOST_CHECK_EQUAL((boost::math::signbit)(a), 0);
      BOOST_CHECK_EQUAL((boost::math::sign)(a), 1);
      BOOST_CHECK_EQUAL((boost::math::copysign)(b, a), RealType(1));
      BOOST_CHECK_EQUAL((boost::math::copysign)(c, a), RealType(1));
      a = -std::numeric_limits<RealType>::infinity();
      BOOST_CHECK((boost::math::signbit)(a) != 0);
      BOOST_CHECK_EQUAL((boost::math::sign)(a), -1);
      BOOST_CHECK_EQUAL((boost::math::copysign)(b, a), RealType(-1));
      BOOST_CHECK_EQUAL((boost::math::copysign)(c, a), RealType(-1));
   }
#if !defined(__SUNPRO_CC) && !defined(BOOST_INTEL)
   if(std::numeric_limits<RealType>::has_quiet_NaN)
   {
      a = std::numeric_limits<RealType>::quiet_NaN();
      BOOST_CHECK_EQUAL((boost::math::signbit)(a), 0);
      BOOST_CHECK_EQUAL((boost::math::sign)(a), 1);
      BOOST_CHECK_EQUAL((boost::math::copysign)(b, a), RealType(1));
      BOOST_CHECK_EQUAL((boost::math::copysign)(c, a), RealType(1));
      a = -std::numeric_limits<RealType>::quiet_NaN();
      BOOST_CHECK((boost::math::signbit)(a) != 0);
      BOOST_CHECK_EQUAL((boost::math::sign)(a), -1);
      BOOST_CHECK_EQUAL((boost::math::copysign)(b, a), RealType(-1));
      BOOST_CHECK_EQUAL((boost::math::copysign)(c, a), RealType(-1));
   }
#endif
}


int test_main(int, char* [])
{
   // Basic sanity-check spot values.
   // (Parameter value, arbitrarily zero, only communicates the floating point type).
   test_spots(0.0F, "float"); // Test float. OK at decdigits = 0 tolerance = 0.0001 %
   test_spots(0.0, "double"); // Test double. OK at decdigits 7, tolerance = 1e07 %
#ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
   test_spots(0.0L, "long double"); // Test long double.
   test_spots(boost::math::concepts::real_concept(0), "real_concept"); // Test real_concept.
#else
   std::cout << "<note>The long double tests have been disabled on this platform "
      "either because the long double overloads of the usual math functions are "
      "not available at all, or because they are too inaccurate for these tests "
      "to pass.</note>" << std::cout;
#endif

   return 0;
} // int test_main(int, char* [])

