// Copyright Paul A. Bristow 2006.

// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

// test_students_t.cpp

// http://en.wikipedia.org/wiki/Student%27s_t_distribution

// Basic sanity test for Student's t probability (quintile) (0. < p < 1).
// and Student's t probability Quintile (0. < p < 1).

#define BOOST_MATH_THROW_ON_DOMAIN_ERROR

#ifdef _MSC_VER
#  pragma warning(disable: 4127) // conditional expression is constant.
#  pragma warning(disable: 4100) // unreferenced formal parameter.
#  pragma warning(disable: 4512) // assignment operator could not be generated.
#  if !(defined _SCL_SECURE_NO_DEPRECATE) || (_SCL_SECURE_NO_DEPRECATE == 0)
#     pragma warning(disable: 4996) // 'std::char_traits<char>::copy' was declared deprecated.
      // #define _SCL_SECURE_NO_DEPRECATE = 1 // avoid C4996 warning.
#  endif
//#  pragma warning(disable: 4244) // conversion from 'double' to 'float', possible loss of data.
#endif

#include <boost/test/included/test_exec_monitor.hpp> // Boost.Test
#include <boost/test/floating_point_comparison.hpp>

#include <boost/math/special_functions/students_t.hpp>
	 using boost::math::students_t_inv;
	 using boost::math::students_t;
#include <boost/math/concepts/real_concept.hpp> // for real_concept

#include <iostream>
	using std::cout;
	using std::endl;
	using std::setprecision;
#include <limits>
  using std::numeric_limits;

template <class FPT>
void test_spots(FPT)
{
  // Basic sanity checks, tolerance is about numeric_limits<FPT>::digits10 decimal places,
	// expressed as a percentage (so -2) for BOOST_CHECK_CLOSE,

	int decdigits = numeric_limits<FPT>::digits10;	// guaranteed for type FPT, eg 6 for float, 15 for double,
	decdigits -= 6; // Perhaps allow some decimal digits margin of numerical error.
	FPT tolerance = static_cast<FPT>(std::pow(10., -(decdigits-2))); // (-2 so as %)
	tolerance *= 1; // Allow some bit(s) small margin (2 means + or - 1 bit) of numerical error.

	cout << "tolerance = " << tolerance << " %" << endl;

	// FPT tolerance = static_cast<FPT>(std::pow(10., -(6-2))); // 1e-6 (as %)
	// Some tests only pass at 1e-5 because probability value is less accurate,
	// a digit in 6th decimal place, although calculated using 
	// a t-distribution generator (calimed 6 decimal digits) at
  // http://faculty.vassar.edu/lowry/VassarStats.html
	// http://faculty.vassar.edu/lowry/tsamp.html
	// df = 5, +/-t = 2.0, 1-tailed = 0.050970, 2-tailed = 0.101939

   //BOOST_CHECK_CLOSE(
   //   ::boost::math::students_t(
   //      static_cast<FPT>(5.),  // df
   //      static_cast<FPT>(-2.)),  // t
   //   static_cast<FPT>(0.050970), // probability.
	//tolerance); // need 5-2

	// http://en.wikipedia.org/wiki/Student%27s_t_distribution#Table_of_selected_values
  // Using tabulated value of t = 3.182 for 0.975, 3 df, one-sided.

	// http://www.mth.kcl.ac.uk/~shaww/web_page/papers/Tdistribution06.pdf refers to:

	// A lookup table of quantiles of the FPT distribution
  // for 1 to 25 in steps of 0.1 is provided in CSV form at:
  // www.mth.kcl.ac.uk/~shaww/web_page/papers/Tsupp/tquantiles.csv
	// gives accurate t of -3.1824463052837 and 3 degrees of freedom.
	// Values below are from this source, saved as tquantiles.xls.
	// DF are across the columns, probabilities down the rows
	// and the t- values (quantiles) are shown.
	// These values are probably accurate to nearly 64-bit double
  // (perhaps 14 decimal digits).

   BOOST_CHECK_CLOSE(
      ::boost::math::students_t(
         static_cast<FPT>(2.),  // degrees_of_freedom
         static_cast<FPT>(-6.96455673428326)),  // t
      static_cast<FPT>(0.01), // probability.
			tolerance);

   BOOST_CHECK_CLOSE(
      ::boost::math::students_t(
         static_cast<FPT>(5.),  // degrees_of_freedom
         static_cast<FPT>(-3.36492999890721)),  // t
      static_cast<FPT>(0.01), // probability.
			tolerance);

   BOOST_CHECK_CLOSE(
      ::boost::math::students_t(
         static_cast<FPT>(1.),  // degrees_of_freedom
         static_cast<FPT>(-31830.988607907)),  // t
      static_cast<FPT>(0.00001), // probability.
			tolerance);

   BOOST_CHECK_CLOSE(
      ::boost::math::students_t(
         static_cast<FPT>(25.),  // degrees_of_freedom
         static_cast<FPT>(-5.2410429995425)),  // t
      static_cast<FPT>(0.00001), // probability.
			tolerance);

   BOOST_CHECK_CLOSE(
      ::boost::math::students_t(
         static_cast<FPT>(1.),  // degrees_of_freedom
         static_cast<FPT>(-63661.97723)),  // t
      static_cast<FPT>(0.000005), // probability.
			tolerance);

    BOOST_CHECK_CLOSE(
      ::boost::math::students_t(
         static_cast<FPT>(5.),  // degrees_of_freedom
         static_cast<FPT>(-17.89686614)),  // t
      static_cast<FPT>(0.000005), // probability.
			tolerance);

    BOOST_CHECK_CLOSE(
      ::boost::math::students_t(
         static_cast<FPT>(25.),  // degrees_of_freedom
         static_cast<FPT>(-5.510848412)),  // t
      static_cast<FPT>(0.000005), // probability.
			tolerance);

  BOOST_CHECK_CLOSE(
      ::boost::math::students_t(
         static_cast<FPT>(10.),  // degrees_of_freedom
         static_cast<FPT>(-1.812461123)),  // t
      static_cast<FPT>(0.05), // probability.
			tolerance);

  BOOST_CHECK_CLOSE(
      ::boost::math::students_t(
         static_cast<FPT>(10.),  // degrees_of_freedom
         static_cast<FPT>(1.812461123)),  // t
      static_cast<FPT>(0.95), // probability.
			tolerance);

  BOOST_CHECK_CLOSE(
      ::boost::math::students_t(
         static_cast<FPT>(10.),  // degrees_of_freedom
         static_cast<FPT>(9.751995491)),  // t
      static_cast<FPT>(0.999999), // probability.
			tolerance);

  BOOST_CHECK_CLOSE(
      ::boost::math::students_t(
         static_cast<FPT>(10.),  // degrees_of_freedom - for ALL degrees_of_freedom!
         static_cast<FPT>(0.)),  // t
      static_cast<FPT>(0.5), // probability.
			tolerance);

	// Student's t Inverse function tests.
  // Special cases

  BOOST_CHECK_EQUAL(boost::math::students_t_inv(
         static_cast<FPT>(1.),  // degrees_of_freedom (ignored).
         static_cast<FPT>(0)),  //  probability == half - special case.
         -numeric_limits<FPT>::infinity()); // t == -infinity.
			
  BOOST_CHECK_EQUAL(boost::math::students_t_inv(
         static_cast<FPT>(1.),  // degrees_of_freedom (ignored).
         static_cast<FPT>(1)),  //  probability == half - special case.
         +numeric_limits<FPT>::infinity()); // t == +infinity.

  BOOST_CHECK_EQUAL(boost::math::students_t_inv(
         static_cast<FPT>(1.),  // degrees_of_freedom (ignored).
         static_cast<FPT>(0.5)),  //  probability == half - special case.
         static_cast<FPT>(0)); // t == zero.

  BOOST_CHECK_CLOSE(boost::math::students_t_inv(
         static_cast<FPT>(1.),  // degrees_of_freedom (ignored).
         static_cast<FPT>(0.5)),  //  probability == half - special case.
      static_cast<FPT>(0), // t == zero.
			tolerance);

   BOOST_CHECK_CLOSE( // Tests of p middling.
      ::boost::math::students_t(
         static_cast<FPT>(5.),  // degrees_of_freedom
         static_cast<FPT>(-0.559429644)),  // t
      static_cast<FPT>(0.3), // probability.
			tolerance);
   BOOST_CHECK_CLOSE(
      ::boost::math::students_t_inv(
         static_cast<FPT>(5.),  // degrees_of_freedom
         static_cast<FPT>(0.3)),  // probability.
      static_cast<FPT>(-0.559429644), // t
			tolerance);

   BOOST_CHECK_CLOSE( // Tests of p high.
      ::boost::math::students_t(
         static_cast<FPT>(5.),  // degrees_of_freedom
         static_cast<FPT>(1.475884049)),  // t
      static_cast<FPT>(0.9), // probability.
			tolerance);
   BOOST_CHECK_CLOSE(
      ::boost::math::students_t_inv(
         static_cast<FPT>(5.),  // degrees_of_freedom
         static_cast<FPT>(0.9)),  // probability.
      static_cast<FPT>(1.475884049), // t
			tolerance);

   BOOST_CHECK_CLOSE( // Tests of p low.
      ::boost::math::students_t(
         static_cast<FPT>(5.),  // degrees_of_freedom
         static_cast<FPT>(-1.475884049)),  // t
      static_cast<FPT>(0.1), // probability.
			tolerance);
   BOOST_CHECK_CLOSE(
      ::boost::math::students_t_inv(
         static_cast<FPT>(5.),  // degrees_of_freedom
         static_cast<FPT>(0.1)),  // probability.
      static_cast<FPT>(-1.475884049), // t
			tolerance);

   BOOST_CHECK_CLOSE(
      ::boost::math::students_t(
         static_cast<FPT>(2.),  // degrees_of_freedom
         static_cast<FPT>(-6.96455673428326)),  // t
      static_cast<FPT>(0.01), // probability.
			tolerance);
   BOOST_CHECK_CLOSE(
      ::boost::math::students_t_inv(
         static_cast<FPT>(2.),  // degrees_of_freedom
         static_cast<FPT>(0.01)),  // probability.
      static_cast<FPT>(-6.96455673428326), // t
			tolerance);

} // template <class FPT>void test_spots(FPT)

int test_main(int, char* [])
{
	 // Basic sanity-check spot values.
	// (Parameter value, arbitrarily zero, only communicates the floating point type).
   test_spots(0.0F); // Test float. OK at decdigits = 0 tolerance = 0.0001 %
  // test_spots(0.0); // Test double. OK at decdigits 7, tolerance = 1e07 %
   //test_spots(0.0L); // Test long double.
   //test_spots(boost::math::concepts::real_concept(0.)); // Test real concept.

   return 0;
} // int test_main(int, char* [])

/*

Output:

Compiling...
test_students_t.cpp
Linking...
Embedding manifest...
Autorun "i:\boost-06-05-03-1300\libs\math\test\Math_test\debug\test_students_t.exe"
Running 1 test case...
tolerance = 100 %
*** No errors detected



*/

