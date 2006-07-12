// test_negative_binomial.cpp

// Copyright Paul A. Bristow 2006.

// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

// Basic sanity test for Negative Binomial Cumulative Distribution Function.

#define BOOST_MATH_THROW_ON_DOMAIN_ERROR

#ifdef _MSC_VER
//#  pragma warning(disable: 4127) // conditional expression is constant.
//#  pragma warning(disable: 4100) // unreferenced formal parameter.
//#  pragma warning(disable: 4512) // assignment operator could not be generated.
//#  pragma warning(disable: 4996) // 'std::char_traits<char>::copy' was declared deprecated.
//#  pragma warning(disable: 4244) // conversion from 'double' to 'const float', possible loss of data.
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

  // MathWorld definition: negaBin[n, p]
  // number of failures that occur before achieving n successes, where probability of success in a trial is p.

  // Cephes probability that k or fewer failures preceeds the nth success.

  // Wikipedia For k + r trials with success probability p,
  // gives the probability of k failures and r successes, with sucess on the last trial.
  // Or probability of number of failures before the rth success.

  // Many be some combinations for which the result is 'exact', or at least is to 40 decimal digits.
	// 40 decimal digits includes 128-bit significand User Defined Floating-Point types,
	// but these are as yet undiscovered for negative binomial. TODO
	
  // Test negative binomial.
  // negative_binomial(k, n, p)
   BOOST_CHECK_CLOSE(negative_binomial(
         static_cast<RealType>(0), // k. pnbinom(1,2,0.5) = 0.5
         static_cast<RealType>(2), // n
         static_cast<RealType>(0.5)),  // mean
      static_cast<RealType>(0.25),  // probability 1/4
			tolerance);

   BOOST_CHECK_CLOSE(negative_binomial(
         static_cast<RealType>(1), // k
         static_cast<RealType>(4), // n
         static_cast<RealType>(0.5)), // p
         static_cast<RealType>(0.1875), // 
			tolerance);

   BOOST_CHECK_CLOSE(negative_binomial(
         static_cast<RealType>(1),  // k
         static_cast<RealType>(2), // n
         static_cast<RealType>(0.5)), // pnbinom(1,2,0.5) = 0.5
         static_cast<RealType>(0.5),  //
			tolerance);

   BOOST_CHECK_CLOSE(negative_binomial(
         static_cast<RealType>(2),  // k
         static_cast<RealType>(2), // n
         static_cast<RealType>(0.5)), // pnbinom(1,2,0.5) = 0.5
         static_cast<RealType>(0.6875),  //
			tolerance);
    BOOST_CHECK_CLOSE(negative_binomial(
         static_cast<RealType>(1),  // k
         static_cast<RealType>(3), // n
         static_cast<RealType>(0.5)), // pnbinom(1,3,0.5) = 0.5
         static_cast<RealType>(0.3125),  // 0.6875
			tolerance);

     BOOST_CHECK_CLOSE(negative_binomial(
         static_cast<RealType>(2),  // k
         static_cast<RealType>(2), // n
         static_cast<RealType>(0.1)), // 
         static_cast<RealType>(0.0523),  //
			tolerance);

      BOOST_CHECK_CLOSE(negative_binomial_c(
         static_cast<RealType>(2),  // k
         static_cast<RealType>(2), // n
         static_cast<RealType>(0.1)), // 
         static_cast<RealType>(1- 0.0523),  //
			tolerance);

   BOOST_CHECK_CLOSE(negative_binomial(
         static_cast<RealType>(26), // k
         static_cast<RealType>(5), // n
         static_cast<RealType>(0.4)), // p
         static_cast<RealType>(0.998968624661119), 
			tolerance);                                    

   BOOST_CHECK_CLOSE(negative_binomial(
         static_cast<RealType>(25), // k
         static_cast<RealType>(5), // n
         static_cast<RealType>(0.4)), // p
         static_cast<RealType>(0.998489925933617), 
			tolerance);

    BOOST_CHECK_CLOSE(negative_binomial_c( // complement.
         static_cast<RealType>(26), // k
         static_cast<RealType>(5), // n
         static_cast<RealType>(0.4)), // p
         static_cast<RealType>(1 - 0.998968624661119), 
			tolerance * 100);    // fails at tolerance * 10
    // 0.0010313753388809799 
    // 0.0010313753388808430
    // Shows a loss of about 3 decimal digits accuracy.

   // negative_binomial_inv(k, n, p) tests.
    BOOST_CHECK_CLOSE(negative_binomial(
         static_cast<RealType>(5), // k.
         static_cast<RealType>(10), // n
         static_cast<RealType>(0.5)),  // trial probability
      static_cast<RealType>(0.15087890625000011),  // probability of failure.
			tolerance);

   BOOST_CHECK_CLOSE(negative_binomial_inv(
         static_cast<RealType>(5), // k.
         static_cast<RealType>(10), // n
         static_cast<RealType>(0.15087890625000011)), // probability of failure = negative_binomial(5, 10, 0.5).
      static_cast<RealType>(0.5),  // result is probability of success
			tolerance);

    BOOST_CHECK_CLOSE(negative_binomial(
         static_cast<RealType>(5), // k failures.
         static_cast<RealType>(10), // n successes.
         static_cast<RealType>(0.1)),  // trial probability of success.
      static_cast<RealType>(1.8662024800000033e-007),  //
			tolerance);

   BOOST_CHECK_CLOSE(negative_binomial_inv(
         static_cast<RealType>(5), // k failures.
         static_cast<RealType>(10), // n successes.
         static_cast<RealType>(1.8662024800000033e-007)), // probability of failure = negative_binomial(5, 10, 0.1)
      static_cast<RealType>(0.1),  // results is trial probability
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

  //for (int i = 0; i < 10; i++)
  //{
  //  cout << i << ' ' << negative_binomial(i, 2, 0.5) << endl;
  //}
  //cout << endl;

  //cout.precision(17);
  //for (double p = 0.; p < 1; p+= 0.1)
  //{
  //  double y = negative_binomial(5, 10, p);
  //  double z = negative_binomial_inv(5, 10, y);
  //  cout << p << ' '<<  y << ' ' << z << ' ' << p - z << endl;
  //}

	// (Parameter value, arbitrarily zero, only communicates the floating point type).
	test_spots(0.0F); // Test float.
	test_spots(0.0); // Test double.
	test_spots(0.0L); // Test long double.
	test_spots(boost::math::concepts::real_concept(0.)); // Test real concept.

	return 0;
} // int test_main(int, char* [])

/*

Output:

------ Build started: Project: test_negative_binomial, Configuration: Debug Win32 ------
Compiling...
test_negative_binomial.cpp
Linking...
Autorun "i:\boost-06-05-03-1300\libs\math\test\Math_test\debug\test_negative_binomial.exe"
Running 1 test case...
BOOST_MATH_THROW_ON_DOMAIN_ERROR is defined to throw on domain error.

0 0.25
1 0.5
2 0.6875
3 0.8125
4 0.890625
5 0.9375
6 0.964844
7 0.980469
8 0.989258
9 0.994141

0 0 -1.#INF 1.#INF
0.10000000000000001 1.8662024800000033e-007 0.10000000000000002 -1.3877787807814457e-017
0.20000000000000001 0.00011322566246400017 0.20000000000000004 -2.7755575615628914e-017
0.30000000000000004 0.0036525210084360034 0.29999999999999999 5.5511151231257827e-017
0.40000000000000002 0.033833302884352018 0.39999999999999997 5.5511151231257827e-017
0.5 0.15087890625000011 0.5 0
0.59999999999999998 0.40321555041484813 0.60000000000000009 -1.1102230246251565e-016
0.69999999999999996 0.72162144020436414 0.70000000000000007 -1.1102230246251565e-016
0.79999999999999993 0.93894857038233592 0.79999999999999982 1.1102230246251565e-016
0.89999999999999991 0.99775032991495194 0.89999999999999947 4.4408920985006262e-016
0.99999999999999989 1 1.#INF -1.#INF
*** No errors detected
Build Time 0:05
Build log was saved at "file://i:\boost-06-05-03-1300\libs\math\test\Math_test\test_negative_binomial\Debug\BuildLog.htm"
test_negative_binomial - 0 error(s), 0 warning(s)
========== Build: 1 succeeded, 0 failed, 0 up-to-date, 0 skipped ==========



*/