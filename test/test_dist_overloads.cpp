// Copyright John Maddock 2006.

// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

// test_students_t.cpp

#define BOOST_MATH_THROW_ON_DOMAIN_ERROR
#define BOOST_MATH_THROW_ON_OVERFLOW_ERROR

#ifdef _MSC_VER
#  pragma warning(disable: 4127) // conditional expression is constant.
#  pragma warning(disable: 4100) // unreferenced formal parameter.
#  pragma warning(disable: 4512) // assignment operator could not be generated.
#  pragma warning(disable: 4510) // default constructor could not be generated.
#  pragma warning(disable: 4610) // can never be instantiated - user defined constructor required.
//#  pragma warning(disable: 4535) // calling _set_se_translator() requires /EHa (in Boost.test)
// Enable C++ Exceptions Yes With SEH Exceptions (/EHa) prevents warning 4535.
#endif

#include <boost/math/concepts/real_concept.hpp> // for real_concept
#include <boost/math/distributions/normal.hpp>
	 using boost::math::normal_distribution;

#include <boost/test/included/test_exec_monitor.hpp> // Boost.Test
#include <boost/test/floating_point_comparison.hpp>


#include <iostream>
	using std::cout;
	using std::endl;
	using std::setprecision;

template <class RealType>
void test_spots(RealType T)
{
   // Basic sanity checks,
   // 2eps as a persentage:
   RealType tolerance = boost::math::tools::epsilon<RealType>() * 2 * 100;

	cout << "Tolerance for type " << typeid(T).name()  << " is " << tolerance << " %" << endl;

   for(int i = -4; i <= 4; ++i)
   {
      BOOST_CHECK_CLOSE(
         ::boost::math::cdf(normal_distribution<RealType>(), i), 
         ::boost::math::cdf(normal_distribution<RealType>(), static_cast<RealType>(i)), 
         tolerance);
      BOOST_CHECK_CLOSE(
         ::boost::math::pdf(normal_distribution<RealType>(), i), 
         ::boost::math::pdf(normal_distribution<RealType>(), static_cast<RealType>(i)), 
         tolerance);
      BOOST_CHECK_CLOSE(
         ::boost::math::cdf(complement(normal_distribution<RealType>(), i)), 
         ::boost::math::cdf(complement(normal_distribution<RealType>(), static_cast<RealType>(i))), 
         tolerance);
      BOOST_CHECK_CLOSE(
         ::boost::math::hazard(normal_distribution<RealType>(), i), 
         ::boost::math::hazard(normal_distribution<RealType>(), static_cast<RealType>(i)), 
         tolerance);
      BOOST_CHECK_CLOSE(
         ::boost::math::chf(normal_distribution<RealType>(), i), 
         ::boost::math::chf(normal_distribution<RealType>(), static_cast<RealType>(i)), 
         tolerance);
   }
   for(float f = 0.01f; f < 1; f += 0.01f)
   {
      BOOST_CHECK_CLOSE(
         ::boost::math::quantile(normal_distribution<RealType>(), f), 
         ::boost::math::quantile(normal_distribution<RealType>(), static_cast<RealType>(f)), 
         tolerance);
      BOOST_CHECK_CLOSE(
         ::boost::math::quantile(complement(normal_distribution<RealType>(), f)), 
         ::boost::math::quantile(complement(normal_distribution<RealType>(), static_cast<RealType>(f))), 
         tolerance);
   }
} // template <class RealType>void test_spots(RealType)

int test_main(int, char* [])
{
	 // Basic sanity-check spot values.
	// (Parameter value, arbitrarily zero, only communicates the floating point type).
  test_spots(0.0F); // Test float. OK at decdigits = 0 tolerance = 0.0001 %
  test_spots(0.0); // Test double. OK at decdigits 7, tolerance = 1e07 %
#ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
  test_spots(0.0L); // Test long double.
#if !BOOST_WORKAROUND(__BORLANDC__, BOOST_TESTED_AT(0x582))
  test_spots(boost::math::concepts::real_concept(0.)); // Test real concept.
#endif
#else
   std::cout << "<note>The long double tests have been disabled on this platform "
      "either because the long double overloads of the usual math functions are "
      "not available at all, or because they are too inaccurate for these tests "
      "to pass.</note>" << std::cout;
#endif

   return 0;
} // int test_main(int, char* [])

