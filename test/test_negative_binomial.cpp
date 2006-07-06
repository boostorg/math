// test_negative_binomial.cpp

// Copyright Paul A. Bristow 2006.

// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

// Basic sanity test for Negative Binomial Cumulative Distribution Function.

#define BOOST_MATH_THROW_ON_DOMAIN_ERROR

#ifdef _MSC_VER
#  pragma warning(disable: 4127) // conditional expression is constant.
#  pragma warning(disable: 4100) // unreferenced formal parameter.
#  pragma warning(disable: 4512) // assignment operator could not be generated.
#  pragma warning(disable: 4996) // 'std::char_traits<char>::copy' was declared deprecated.
#  pragma warning(disable: 4244) // conversion from 'double' to 'const float', possible loss of data.
#endif

#include <boost/math/special_functions/negative_binomial.hpp> // for negative_binomial
	using ::boost::math::negative_binomial;
	using ::boost::math::negative_binomial_c;
	using ::boost::math::negative_binomial_inv;
#include <boost/math/concepts/real_concept.hpp> // for real_concept
	using ::boost::math::concepts::real_concept;

#include <boost/test/included/test_exec_monitor.hpp> // for test_main
#include <boost/test/floating_point_comparison.hpp> // for BOOST_CHECK_CLOSE

#include <iostream>
  using std::cout;
  using std::endl;
#include <limits>
  using std::numeric_limits;

template <class RealType> // RealType is any floating-point type.
void test_spots(RealType)
{
  // Basic sanity checks, tolerance is about numeric_limits<RealType>::digits10 decimal places,
	// guaranteed for type RealType, eg 6 for float, 15 for double,
	// expressed as a percentage (so -2) for BOOST_CHECK_CLOSE,

	int decdigits = numeric_limits<RealType>::digits10;
	decdigits -= 1; // Perhaps allow some decimal digits margin of numerical error.
	RealType tolerance = static_cast<RealType>(std::pow(10., -(decdigits-2))); // 1e-6 (-2 so as %)
	tolerance *= 1; // Allow some bit(s) small margin (2 means + or - 1 bit) of numerical error.
	// Typically 2e-13% = 2e-15 as fraction for double.

	// Sources of spot test values:

  // MathCAD
  //   pnbinom(k, n, p)
  // Returns cumulative negative binomial distribution.
  //   qnbinom(p, n, q) 
  // Returns the inverse negative binomial distribution function,
  // the smallest integer k so that pnbinom(k, n, q) >= p

  // Many be some combinations for which the result is 'exact', or at least is to 40 decimal digits.
	// 40 decimal digits includes 128-bit significand User Defined Floating-Point types,
	// but these are as yet undiscovered for negative binomial. TODO
	
  // Test negative binomial.
  // These are from MathCAD

   //BOOST_CHECK_CLOSE(negative_binomial(
   //      static_cast<RealType>(1), // k. 
   //      static_cast<RealType>(2), // n
   //      static_cast<RealType>(0.5)),  // mean
   //   static_cast<RealType>(3) / // probability 3/4.
   //   static_cast<RealType>(4.), // q
			//tolerance);

   //BOOST_CHECK_CLOSE(negative_binomial(
   //      static_cast<RealType>(1), // k
   //      static_cast<RealType>(2), // n
   //      static_cast<RealType>(0.5)), // p
   //      static_cast<RealType>(0.75), // pnbinom(1,2,0.5) = 0.5
			//tolerance);

   //BOOST_CHECK_CLOSE(negative_binomial_c(
   //      static_cast<RealType>(1),  // k
   //      static_cast<RealType>(2), // n
   //      static_cast<RealType>(0.5)), // p
   //      static_cast<RealType>(0.25),  // q
			//tolerance);
   //BOOST_CHECK_CLOSE(negative_binomial_inv(
   //      static_cast<RealType>(1),  // 
   //      static_cast<RealType>(2), // n
   //      static_cast<RealType>(0.75)), // q
   //      static_cast<RealType>(0.5), // p
			//tolerance);

   BOOST_CHECK_CLOSE(negative_binomial(
         static_cast<RealType>(1), // k. 
         static_cast<RealType>(2), // n
         static_cast<RealType>(0.5)),  // mean
      static_cast<RealType>(0.5), // q
			tolerance);

   BOOST_CHECK_CLOSE(negative_binomial(
         static_cast<RealType>(1), // k. 
         static_cast<RealType>(2), // n
         static_cast<RealType>(0.5)),  // mean
      static_cast<RealType>(0.5), // q
			tolerance);
   BOOST_CHECK_CLOSE(negative_binomial(
         static_cast<RealType>(1), // k
         static_cast<RealType>(2), // n
         static_cast<RealType>(0.5)), // p
         static_cast<RealType>(0.5), // pnbinom(1,2,0.5) = 0.5
			tolerance);

   BOOST_CHECK_CLOSE(negative_binomial_c(
         static_cast<RealType>(1),  // k
         static_cast<RealType>(2), // n
         static_cast<RealType>(0.5)), // p
         static_cast<RealType>(0.5),  // q
			tolerance);
   BOOST_CHECK_CLOSE(negative_binomial_inv(
         static_cast<RealType>(1),  // 
         static_cast<RealType>(2), // n
         static_cast<RealType>(0.5)), // q
         static_cast<RealType>(0.5), // p
			1e-6); // tolerance );


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




*/