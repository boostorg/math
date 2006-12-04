// Copyright Paul Bristow 2006.
// Copyright John Maddock 2006.

// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

// test_uniform.cpp

#define BOOST_MATH_THROW_ON_DOMAIN_ERROR
#define BOOST_MATH_THROW_ON_OVERFLOW_ERROR
//#define BOOST_MATH_THROW_ON_UNDERFLOW_ERROR 
// Ignore underflow to zero.

#ifdef _MSC_VER
#  pragma warning(disable: 4127) // conditional expression is constant.
#  pragma warning(disable: 4100) // unreferenced formal parameter.
#  pragma warning(disable: 4512) // assignment operator could not be generated.
#  pragma warning(disable: 4510) // default constructor could not be generated.
#  pragma warning(disable: 4610) // can never be instantiated - user defined constructor required.
#  pragma warning(disable: 4180) // qualifier applied to function type has no meaning; ignored.
#  if !(defined _SCL_SECURE_NO_DEPRECATE) || (_SCL_SECURE_NO_DEPRECATE == 0)
#    pragma warning(disable: 4996) // 'std::char_traits<char>::copy' was declared deprecated.
#  endif
#endif

#include <boost/math/concepts/real_concept.hpp> // for real_concept
#include <boost/test/included/test_exec_monitor.hpp> // Boost.Test
#include <boost/test/floating_point_comparison.hpp>

#include <boost/math/distributions/uniform.hpp>
	 using boost::math::uniform_distribution;
#include <boost/math/tools/test.hpp> 

#include <iostream>
	using std::cout;
	using std::endl;
	using std::setprecision;
#include <limits>
  using std::numeric_limits;

template <class RealType>
void check_uniform(RealType lower, RealType upper, RealType x, RealType p, RealType q, RealType tol)
{
   BOOST_CHECK_CLOSE_FRACTION(
      ::boost::math::cdf(
         uniform_distribution<RealType>(lower, upper),   // distribution.
         x),  // random variable.
         p,    // probability.
			tol);   // tolerance.
   BOOST_CHECK_CLOSE_FRACTION(
      ::boost::math::cdf(
         complement(
            uniform_distribution<RealType>(lower, upper), // distribution.
            x)),    // random variable.
         q,    // probability complement.
			tol);  // tolerance.
   BOOST_CHECK_CLOSE_FRACTION(
      ::boost::math::quantile(
         uniform_distribution<RealType>(lower, upper),  // distribution.
         p),   // probability.
         x,  // random variable.
			tol);  // tolerance.
   BOOST_CHECK_CLOSE_FRACTION(
      ::boost::math::quantile(
         complement(
            uniform_distribution<RealType>(lower, upper),  // distribution.
            q)),     // probability complement.
         x,                                             // random variable.
			tol);  // tolerance.
} // void check_uniform

template <class RealType>
void test_spots(RealType T)
{
   // Basic santity checks
   //
   // These test values were generated for the normal distribution
   // using the online calculator at 
   // http://espse.ed.psu.edu/edpsych/faculty/rhale/hale/507Mat/statlets/free/pdist.htm
   //
   // Tolerance is just over 5 decimal digits expressed as a fraction:
   // that's the limit of the test data.
	RealType tolerance = 2e-5f;  
	cout << "Tolerance for type " << typeid(T).name()  << " is " << tolerance << "." << endl;

   using std::exp;

   // Tests for PDF
   //
   BOOST_CHECK_CLOSE_FRACTION( // x == upper
      pdf(uniform_distribution<RealType>(0, 1), static_cast<RealType>(0)), 
      static_cast<RealType>(1), 
      tolerance);
   BOOST_CHECK_CLOSE_FRACTION( // x == lower
      pdf(uniform_distribution<RealType>(0, 1), static_cast<RealType>(1)), 
      static_cast<RealType>(1), 
      tolerance);
   BOOST_CHECK_CLOSE_FRACTION( // x > upper
      pdf(uniform_distribution<RealType>(0, 1), static_cast<RealType>(-1)), 
      static_cast<RealType>(0), 
      tolerance); 
   BOOST_CHECK_CLOSE_FRACTION( // x < lower
      pdf(uniform_distribution<RealType>(0, 1), static_cast<RealType>(2)), 
      static_cast<RealType>(0), 
      tolerance);

   if(std::numeric_limits<RealType>::has_infinity)
    { // BOOST_CHECK tests for infinity using std::numeric_limits<>::infinity()
      // Note that infinity is not implemented for real_concept, so these tests
      // are only done for types, like built-in float, double.. that have infinity.
      // Note that these assume that  BOOST_MATH_THROW_ON_OVERFLOW_ERROR is NOT defined.
      // #define BOOST_MATH_THROW_ON_OVERFLOW_ERROR would give a throw here.
      // #define BOOST_MATH_THROW_ON_DOMAIN_ERROR IS defined, so the throw path
      // of error handling is tested below with BOOST_CHECK_THROW tests.

     BOOST_CHECK_EQUAL( // x == infinity should be OK.
       pdf(uniform_distribution<RealType>(0, 1), static_cast<RealType>(std::numeric_limits<RealType>::infinity())), 
       static_cast<RealType>(0));

     BOOST_CHECK_EQUAL( // x == minus infinity should be OK too.
       pdf(uniform_distribution<RealType>(0, 1), static_cast<RealType>(-std::numeric_limits<RealType>::infinity())), 
       static_cast<RealType>(0));
   }
   if(std::numeric_limits<RealType>::has_quiet_NaN)
   { // BOOST_CHECK tests for NaN using std::numeric_limits<>::has_quiet_NaN() - should throw.
     BOOST_CHECK_THROW(
       pdf(uniform_distribution<RealType>(0, 1), static_cast<RealType>(std::numeric_limits<RealType>::quiet_NaN())), 
       std::domain_error);
     BOOST_CHECK_THROW(
       pdf(uniform_distribution<RealType>(0, 1), static_cast<RealType>(-std::numeric_limits<RealType>::quiet_NaN())), 
       std::domain_error);
   } // test for x = NaN using std::numeric_limits<>::quiet_NaN()

   // cdf
   BOOST_CHECK_EQUAL( // x < lower
      cdf(uniform_distribution<RealType>(0, 1), static_cast<RealType>(-1)), 
      static_cast<RealType>(0) );
   BOOST_CHECK_CLOSE_FRACTION(
      cdf(uniform_distribution<RealType>(0, 1), static_cast<RealType>(0)), 
      static_cast<RealType>(0), 
      tolerance);
   BOOST_CHECK_CLOSE_FRACTION(
      cdf(uniform_distribution<RealType>(0, 1), static_cast<RealType>(0.5)), 
      static_cast<RealType>(0.5), 
      tolerance);
   BOOST_CHECK_CLOSE_FRACTION(
      cdf(uniform_distribution<RealType>(0, 1), static_cast<RealType>(0.1)), 
      static_cast<RealType>(0.1), 
      tolerance);
   BOOST_CHECK_CLOSE_FRACTION(
      cdf(uniform_distribution<RealType>(0, 1), static_cast<RealType>(0.9)), 
      static_cast<RealType>(0.9), 
      tolerance);
   BOOST_CHECK_EQUAL( // x > upper
      cdf(uniform_distribution<RealType>(0, 1), static_cast<RealType>(2)), 
      static_cast<RealType>(1));

  // cdf complement
   BOOST_CHECK_EQUAL( // x < lower
      cdf(complement(uniform_distribution<RealType>(0, 1), static_cast<RealType>(0))), 
      static_cast<RealType>(1));
   BOOST_CHECK_EQUAL( // x == 0
      cdf(complement(uniform_distribution<RealType>(0, 1), static_cast<RealType>(0))), 
      static_cast<RealType>(1));
   BOOST_CHECK_CLOSE_FRACTION( // x = 0.1
      cdf(complement(uniform_distribution<RealType>(0, 1), static_cast<RealType>(0.1))), 
      static_cast<RealType>(0.9), 
      tolerance);
   BOOST_CHECK_CLOSE_FRACTION( // x = 0.5
      cdf(complement(uniform_distribution<RealType>(0, 1), static_cast<RealType>(0.5))), 
      static_cast<RealType>(0.5), 
      tolerance);
   BOOST_CHECK_EQUAL( // x == 1
      cdf(complement(uniform_distribution<RealType>(0, 1), static_cast<RealType>(1))), 
      static_cast<RealType>(0));
   BOOST_CHECK_EQUAL( // x > upper
      cdf(complement(uniform_distribution<RealType>(0, 1), static_cast<RealType>(2))), 
      static_cast<RealType>(1));

   // quantile

   BOOST_CHECK_CLOSE_FRACTION(
      quantile(uniform_distribution<RealType>(0, 1), static_cast<RealType>(0.9)), 
      static_cast<RealType>(0.9), 
      tolerance);
   BOOST_CHECK_CLOSE_FRACTION(
      quantile(uniform_distribution<RealType>(0, 1), static_cast<RealType>(0.1)), 
      static_cast<RealType>(0.1), 
      tolerance);
   BOOST_CHECK_CLOSE_FRACTION(
      quantile(uniform_distribution<RealType>(0, 1), static_cast<RealType>(0.5)), 
      static_cast<RealType>(0.5), 
      tolerance);
   BOOST_CHECK_CLOSE_FRACTION(
      quantile(uniform_distribution<RealType>(0, 1), static_cast<RealType>(0)), 
      static_cast<RealType>(0), 
      tolerance);
   BOOST_CHECK_CLOSE_FRACTION(
      quantile(uniform_distribution<RealType>(0, 1), static_cast<RealType>(1)), 
      static_cast<RealType>(1), 
      tolerance);

   // quantile complement

   BOOST_CHECK_CLOSE_FRACTION(
      quantile(complement(uniform_distribution<RealType>(0, 1), static_cast<RealType>(0.1))), 
      static_cast<RealType>(0.9), 
      tolerance);
   BOOST_CHECK_CLOSE_FRACTION(
      quantile(complement(uniform_distribution<RealType>(0, 1), static_cast<RealType>(0.9))), 
      static_cast<RealType>(0.1), 
      tolerance);
   BOOST_CHECK_CLOSE_FRACTION(
      quantile(complement(uniform_distribution<RealType>(0, 1), static_cast<RealType>(0.5))), 
      static_cast<RealType>(0.5), 
      tolerance);
   BOOST_CHECK_CLOSE_FRACTION(
      quantile(complement(uniform_distribution<RealType>(0, 1), static_cast<RealType>(0))), 
      static_cast<RealType>(1), 
      tolerance);
   BOOST_CHECK_CLOSE_FRACTION(
      quantile(complement(uniform_distribution<RealType>(0, 1), static_cast<RealType>(1))), 
      static_cast<RealType>(1), 
      tolerance);

   // Some tests using a different location & scale, neight zero or unity.
   BOOST_CHECK_CLOSE_FRACTION( // x == mid
      pdf(uniform_distribution<RealType>(-1, 2), static_cast<RealType>(1)), 
      static_cast<RealType>(0.3333333333333333333333333333333333333333333333333333), 
      tolerance);

   BOOST_CHECK_CLOSE_FRACTION( // x == upper
      pdf(uniform_distribution<RealType>(-1, 2), static_cast<RealType>(+2)), 
      static_cast<RealType>(0.3333333333333333333333333333333333333333333333333333),  // 1 / (2 - -1) = 1/3
      tolerance);

   BOOST_CHECK_CLOSE_FRACTION( // x == lower
      cdf(uniform_distribution<RealType>(-1, 2), static_cast<RealType>(-1)), 
      static_cast<RealType>(0), 
      tolerance);
   BOOST_CHECK_CLOSE_FRACTION( // x == upper
      cdf(uniform_distribution<RealType>(-1, 2), static_cast<RealType>(0)), 
      static_cast<RealType>(0.3333333333333333333333333333333333333333333333333333), 
      tolerance);

   BOOST_CHECK_CLOSE_FRACTION( // x == upper
      cdf(uniform_distribution<RealType>(-1, 2), static_cast<RealType>(1)), 
      static_cast<RealType>(0.6666666666666666666666666666666666666666666666666667), 
      tolerance);

   BOOST_CHECK_CLOSE_FRACTION( // x == lower
      cdf(uniform_distribution<RealType>(-1, 2), static_cast<RealType>(2)), 
      static_cast<RealType>(1), 
      tolerance);

   BOOST_CHECK_CLOSE_FRACTION( // x == upper
      quantile(uniform_distribution<RealType>(-1, 2), static_cast<RealType>(0.6666666666666666666666666666666666666666666666666667)), 
      static_cast<RealType>(1),
      tolerance);

      check_uniform(
      static_cast<RealType>(0),       // lower
      static_cast<RealType>(1),       // upper
      static_cast<RealType>(0.5),     // x
      static_cast<RealType>(0.5),     // p
      static_cast<RealType>(1 - 0.5), // q
      tolerance);

      // Some Not-standard uniform tests.
      check_uniform(
      static_cast<RealType>(-1),    // lower
      static_cast<RealType>(1),     // upper
      static_cast<RealType>(0),     // x
      static_cast<RealType>(0.5),   // p
      static_cast<RealType>(1 - 0.5), // q = 1 - p
      tolerance);

      check_uniform(
      static_cast<RealType>(1),    // lower
      static_cast<RealType>(3),     // upper
      static_cast<RealType>(2),     // x
      static_cast<RealType>(0.5),   // p
      static_cast<RealType>(1 - 0.5), // q = 1 - p
      tolerance);

      check_uniform(
      static_cast<RealType>(-1),    // lower
      static_cast<RealType>(2),     // upper
      static_cast<RealType>(1),     // x
      static_cast<RealType>(0.66666666666666666666666666666666666666666667),   // p
      static_cast<RealType>(0.33333333333333333333333333333333333333333333), // q = 1 - p
      tolerance);
   tolerance = (std::max)(
      boost::math::tools::epsilon<RealType>(),
      static_cast<RealType>(boost::math::tools::epsilon<double>())) * 5; // 5 eps as a fraction.
	 cout << "Tolerance (as fraction) for type " << typeid(T).name()  << " is " << tolerance << "." << endl;
   uniform_distribution<RealType> distu01(0, 1);
   RealType x = static_cast<RealType>(0.5);
   using namespace std; // ADL of std names.
   // mean:
   BOOST_CHECK_CLOSE_FRACTION(
      mean(distu01), static_cast<RealType>(0.5), tolerance);
   // variance:
   BOOST_CHECK_CLOSE_FRACTION(
      variance(distu01), static_cast<RealType>(0.0833333333333333333333333333333333333333333), tolerance);
   // std deviation:
   BOOST_CHECK_CLOSE_FRACTION(
    standard_deviation(distu01), sqrt(variance(distu01)), tolerance);
   // hazard:
   BOOST_CHECK_CLOSE_FRACTION(
    hazard(distu01, x), pdf(distu01, x) / cdf(complement(distu01, x)), tolerance);
   // cumulative hazard:
   BOOST_CHECK_CLOSE_FRACTION(
    chf(distu01, x), -log(cdf(complement(distu01, x))), tolerance);
   // coefficient_of_variation:
   BOOST_CHECK_CLOSE_FRACTION(
    coefficient_of_variation(distu01), standard_deviation(distu01) / mean(distu01), tolerance);
   // mode:
   BOOST_CHECK_CLOSE_FRACTION(
    mode(distu01), static_cast<RealType>(0), tolerance);
   // skewness:
   BOOST_CHECK_EQUAL(
    skewness(distu01), static_cast<RealType>(0));
   // kertosis:
   BOOST_CHECK_CLOSE_FRACTION(
    kurtosis(distu01), kurtosis_excess(distu01) + static_cast<RealType>(3), tolerance);
   // kertosis excess:
   BOOST_CHECK_CLOSE_FRACTION(
    kurtosis_excess(distu01), static_cast<RealType>(-1.2), tolerance);

   if(std::numeric_limits<RealType>::has_infinity)
  { // BOOST_CHECK tests for infinity using std::numeric_limits<>::infinity()
    // Note that infinity is not implemented for real_concept, so these tests
    // are only done for types, like built-in float, double.. that have infinity.
    // Note that these assume that  BOOST_MATH_THROW_ON_OVERFLOW_ERROR is NOT defined.
    // #define BOOST_MATH_THROW_ON_OVERFLOW_ERROR would give a throw here.
    // #define BOOST_MATH_THROW_ON_DOMAIN_ERROR IS defined, so the throw path
    // of error handling is tested below with BOOST_CHECK_THROW tests.

    BOOST_CHECK_EQUAL(pdf(distu01, std::numeric_limits<RealType>::infinity()),  0);
    BOOST_CHECK_EQUAL(pdf(distu01, -std::numeric_limits<RealType>::infinity()),  0);
   } // test for infinity using std::numeric_limits<>::infinity()
   else
   { // real_concept case, does has_infinfity == false, so can't check it throws.
     // cout << std::numeric_limits<RealType>::infinity() << ' '
     // << boost::math::fpclassify(std::numeric_limits<RealType>::infinity()) << endl;
     // value of std::numeric_limits<RealType>::infinity() is zero, so FPclassify is zero,
     // so (boost::math::isfinite)(std::numeric_limits<RealType>::infinity()) does not detect infinity.
     // so these tests would never throw.
     //BOOST_CHECK_THROW(pdf(distu01, std::numeric_limits<RealType>::infinity()),  std::domain_error);
     //BOOST_CHECK_THROW(pdf(distu01, std::numeric_limits<RealType>::quiet_NaN()),  std::domain_error);
     // BOOST_CHECK_THROW(pdf(distu01, boost::math::tools::max_value<RealType>() * 2),  std::domain_error); // Doesn't throw.
     BOOST_CHECK_EQUAL(pdf(distu01, boost::math::tools::max_value<RealType>()), 0); 
   }
   // Special cases:
   BOOST_CHECK(pdf(distu01, 0) == 1);
   BOOST_CHECK(cdf(distu01, 0) == 0);
   BOOST_CHECK(pdf(distu01, 1) == 1);
   BOOST_CHECK(cdf(distu01, 1) == 1);
   BOOST_CHECK(cdf(complement(distu01, 0)) == 1);
   BOOST_CHECK(cdf(complement(distu01, 1)) == 0);
   BOOST_CHECK(quantile(distu01, 0) == 0);
   BOOST_CHECK(quantile(complement(distu01, 0)) == 1);
   BOOST_CHECK(quantile(distu01, 1) == 1);
   BOOST_CHECK(quantile(complement(distu01, 1)) == 1);

   // Error checks:
   if(std::numeric_limits<RealType>::has_quiet_NaN)
   { // BOOST_CHECK tests for quiet_NaN (not for real_concept, for example - see notes above).
     BOOST_CHECK_THROW(uniform_distribution<RealType>(0, std::numeric_limits<RealType>::quiet_NaN()), std::domain_error);
     BOOST_CHECK_THROW(uniform_distribution<RealType>(0, -std::numeric_limits<RealType>::quiet_NaN()), std::domain_error);
   }
   BOOST_CHECK_THROW(uniform_distribution<RealType>(1, 0), std::domain_error); // lower > upper!

} // template <class RealType>void test_spots(RealType)

int test_main(int, char* [])
{
  // Check that can construct uniform distribution using the two convenience methods:
  using namespace boost::math;
  uniform myustd; // Using typedef
  // == uniform_distribution<double> myustd;
  BOOST_CHECK_EQUAL(myustd.lower(), 0); // Check defaults.
  BOOST_CHECK_EQUAL(myustd.upper(), 1);
	uniform_distribution<> myu01(0, 1); // Using default RealType double.
  BOOST_CHECK_EQUAL(myu01.lower(), 0); // Check defaults again.
  BOOST_CHECK_EQUAL(myu01.upper(), 1);

	 // Basic sanity-check spot values.
	// (Parameter value, arbitrarily zero, only communicates the floating point type).
  test_spots(0.0F); // Test float. OK at decdigits = 0 tolerance = 0.0001 %
  test_spots(0.0); // Test double. OK at decdigits 7, tolerance = 1e07 %
#ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
  test_spots(0.0L); // Test long double.
#if !BOOST_WORKAROUND(__BORLANDC__, BOOST_TESTED_AT(0x0582))
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

/*

Output:

------ Build started: Project: test_uniform, Configuration: Debug Win32 ------
Compiling...
test_uniform.cpp
Linking...
Autorun "i:\boost-06-05-03-1300\libs\math\test\Math_test\debug\test_uniform.exe"
Running 1 test case...
Tolerance for type float is 2e-005.
Tolerance (as fraction) for type float is 5.96046e-007.
Tolerance for type double is 2e-005.
Tolerance (as fraction) for type double is 1.11022e-015.
Tolerance for type long double is 2e-005.
Tolerance (as fraction) for type long double is 1.11022e-015.
Tolerance for type class boost::math::concepts::real_concept is 2e-005.
Tolerance (as fraction) for type class boost::math::concepts::real_concept is 1.11022e-015.
*** No errors detected
Build Time 0:06
Build log was saved at "file://i:\boost-06-05-03-1300\libs\math\test\Math_test\test_uniform\Debug\BuildLog.htm"
test_uniform - 0 error(s), 0 warning(s)
========== Build: 1 succeeded, 0 failed, 0 up-to-date, 0 skipped ==========

*/


