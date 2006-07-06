// test_binomial.cpp

// Copyright Paul A. Bristow 2006.

// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

// Basic sanity test for binomial Cumulative Distribution Function.

#define BOOST_MATH_THROW_ON_DOMAIN_ERROR

#ifdef _MSC_VER
#  pragma warning(disable: 4127) // conditional expression is constant.
#  pragma warning(disable: 4100) // unreferenced formal parameter.
#  pragma warning(disable: 4512) // assignment operator could not be generated.
#  pragma warning(disable: 4996) // 'std::char_traits<char>::copy' was declared deprecated.
#endif

#include <boost/math/special_functions/binomial.hpp> // for binomial
	using ::boost::math::binomial;
	using ::boost::math::binomial_c;
	using ::boost::math::binomial_inv;
#include <boost/math/concepts/real_concept.hpp> // for real_concept
	using ::boost::math::concepts::real_concept;

#include <boost/test/included/test_exec_monitor.hpp> // for test_main
#include <boost/test/floating_point_comparison.hpp> // for BOOST_CHECK_CLOSE

#include <iostream>
  using std::cout;
  using std::endl;
#include <limits>
  using std::numeric_limits;

template <class RealType> // Any floating-point type RealType.
void test_spots(RealType)
{
  // Basic sanity checks, tolerance is about numeric_limits<RealType>::digits10 decimal places,
	// guaranteed for type RealType, eg 6 for float, 15 for double,
	// expressed as a percentage (so -2) for BOOST_CHECK_CLOSE,

	int decdigits = numeric_limits<RealType>::digits10;
	decdigits -= 3; // Perhaps allow some decimal digits margin of numerical error.
	RealType tolerance = static_cast<RealType>(std::pow(10., -(decdigits-2))); // 1e-6 (-2 so as %)
	tolerance *= 1; // Allow some bit(s) small margin (2 means + or - 1 bit) of numerical error.
	// Typically 2e-13% = 2e-15 as fraction for double.

	// Sources of spot test values:

  // MathCAD defines ppois(k, lambda) as k integer, k >=0  (Note allows 0?)
  // P = pbinom(30, 500, 0.05) = 0.869147702104609


  // Test binomial
   BOOST_CHECK_CLOSE(binomial(
         static_cast<RealType>(30),  // k - as floating-point.
         static_cast<RealType>(500),  // n  < k as floating-point.
         static_cast<RealType>(0.05)),  // probability of success.
      static_cast<RealType>(0.869147702104609), // probability.
			tolerance);

   BOOST_CHECK_CLOSE(binomial(
         static_cast<RealType>(250),  // k - as floating-point.
         static_cast<RealType>(500),  // n  < k as floating-point.
         static_cast<RealType>(0.05)),  // probability of success.
      static_cast<RealType>(1), // probability.
			tolerance);

   BOOST_CHECK_CLOSE(binomial(
         static_cast<RealType>(470),  // k - as floating-point.
         static_cast<RealType>(500),  // n  < k as floating-point.
         static_cast<RealType>(0.95)),  // probability of success.
      static_cast<RealType>(0.176470742656766), // probability.
			tolerance);

   BOOST_CHECK_CLOSE(binomial(
         static_cast<RealType>(400),  // k - as floating-point.
         static_cast<RealType>(500),  // n  < k as floating-point.
         static_cast<RealType>(0.05)),  // probability of success.
      static_cast<RealType>(1), // probability.
			tolerance);

   BOOST_CHECK_CLOSE(binomial_c(
         static_cast<RealType>(400),  // k - as floating-point.
         static_cast<RealType>(500),  // n  < k as floating-point.
         static_cast<RealType>(0.95)),  // probability of success.
      static_cast<RealType>(1), // probability.
			tolerance);

   BOOST_CHECK_CLOSE(binomial(
         static_cast<RealType>(5),  // k - as floating-point.
         static_cast<RealType>(500),  // n  < k as floating-point.
         static_cast<RealType>(0.05)),  // probability of success.
      static_cast<RealType>(9.181808267643E-7), // probability.
			tolerance);

   BOOST_CHECK_CLOSE(binomial_inv(
         static_cast<RealType>(470),  // k - as floating-point.
         static_cast<RealType>(500),  // n  < k as floating-point.
         static_cast<RealType>(0.176470742656766)),  // probability of success.
      static_cast<RealType>(0.95), // probability.
			tolerance);

      BOOST_CHECK_CLOSE(binomial(
         static_cast<RealType>(1),  // k - as floating-point.
         static_cast<RealType>(2),  // n  < k as floating-point.
         static_cast<RealType>(0.5)),  // probability of success.
      static_cast<RealType>(0.75), // probability.
			tolerance);


} // template <class RealType>void test_spots(RealType)

int test_main(int, char* [])
{
	// Basic sanity-check spot values.
#ifdef BOOST_MATH_THROW_ON_DOMAIN_ERROR
  cout << "BOOST_MATH_THROW_ON_DOMAIN_ERROR" << " is defined to throw on domain error." << endl;
#else
  cout << "BOOST_MATH_THROW_ON_DOMAIN_ERROR" << " is NOT defined, so NO throw on domain error." << endl;
#endif

	// (Parameter value, arbitrarily zero, only communicates the floating point type).
	test_spots(0.0F); // Test float.
	test_spots(0.0); // Test double.
	test_spots(0.0L); // Test long double.
	test_spots(boost::math::concepts::real_concept(0.)); // Test real concept.

	return 0;
} // int test_main(int, char* [])

/*

Output:

Compiling...
test_binomial.cpp
Linking...
Autorun "i:\boost-06-05-03-1300\libs\math\test\Math_test\debug\test_binomial.exe"
Running 1 test case...
BOOST_MATH_THROW_ON_DOMAIN_ERROR is defined to throw on domain error.
*** No errors detected
Build Time 0:06
Build log was saved at "file://i:\boost-06-05-03-1300\libs\math\test\Math_test\test_binomial\Debug\BuildLog.htm"
test_binomial - 0 error(s), 0 warning(s)
========== Build: 1 succeeded, 0 failed, 0 up-to-date, 0 skipped ==========

*/