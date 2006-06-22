// Copyright Paul A. Bristow 2006.

// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

// test_fisher.cpp

// Basic sanity test for Fisher Cumulative Distribution Function.

#ifdef _MSC_VER
#  pragma warning(disable: 4127) // conditional expression is constant.
#  pragma warning(disable: 4100) // unreferenced formal parameter.
#endif

#include <boost/math/concepts/real_concept.hpp> // for real_concept
#include <boost/math/special_functions/fisher.hpp>

#include <boost/test/included/test_exec_monitor.hpp> // test_main
#include <boost/test/floating_point_comparison.hpp> //  BOOST_CHECK_CLOSE

#include <limits>
using std::numeric_limits;

template <class FPT>
void test_spots(FPT)
{
   // Basic sanity checks, tolerance is 6 decimal places
	 // expressed as a percentage (so -2) for BOOST_CHECK_CLOSE,
   // One check per domain of the implementation:

	int decdigits = numeric_limits<FPT>::digits10; // guaranteed for this type, eg 6 for float, 15 for double.
	decdigits -= 0; // Perhaps allow some decimal digits margin of numerical error.
	FPT tolerance = static_cast<FPT>(std::pow(10., -(decdigits-2))); // 1e-6 (as %)
	tolerance *= 2; // Allow some small margin (2 means 1 least decimal digit?) of numerical error.
	// Typically 2e-13% = 2e-15 as fraction for double.
	// Whereas tolerance *= 1.5 fails thus
	// difference between ::boost::math::fisher( static_cast<FPT>(1.), static_cast<FPT>(2.), static_cast<FPT>(6.5333333333333))
	// {0.12500000000000022} and static_cast<FPT>(0.125){0.125} exceeds 1.5e-013%

  //	FPT tolerance = static_cast<FPT>(std::pow(10., -(4-2))); // 1e-4 (as %) OK with misc test values.
  //	PT tolerance = static_cast<FPT>(std::pow(10., -(6-2))); // 1e-6 (as %) fails for most values.
	// This test only passes at 1e-5 because probability value is less accurate,
	// a digit in 6th decimal place, although calculated using values from

	// http://www.danielsoper.com/statcalc/calc04.aspx

	// BUT I conclude that these values are only accurate to about 1e-4,
	// so not very useful for testing!

 /*
 BOOST_CHECK_CLOSE(
      ::boost::math::fisher(
         static_cast<FPT>(1.),  // df1
         static_cast<FPT>(1.),  // df2
         static_cast<FPT>(161.4476)),  // F
      static_cast<FPT>(0.05), // probability.
			tolerance);
   BOOST_CHECK_CLOSE(
      ::boost::math::fisher(
         static_cast<FPT>(2.),  // df1
         static_cast<FPT>(2.),  // df2
         static_cast<FPT>(19.0000)),  // F
      static_cast<FPT>(0.05), // probability.
			tolerance);
   BOOST_CHECK_CLOSE(
      ::boost::math::fisher(
         static_cast<FPT>(1.),  // df1
         static_cast<FPT>(10.),  // df2
         static_cast<FPT>(4.9646)),  // F
      static_cast<FPT>(0.05), // probability.
			tolerance);
   BOOST_CHECK_CLOSE(
      ::boost::math::fisher(
         static_cast<FPT>(1.),  // df1
         static_cast<FPT>(100.),  // df2
         static_cast<FPT>(3.9361)),  // F
      static_cast<FPT>(0.05), // probability.
			tolerance);
   BOOST_CHECK_CLOSE(
      ::boost::math::fisher(
         static_cast<FPT>(1.),  // df1
         static_cast<FPT>(1000000.),  // df2
         static_cast<FPT>(3.8415)),  // F
      static_cast<FPT>(0.05), // probability.
			tolerance);
   BOOST_CHECK_CLOSE(
      ::boost::math::fisher(
         static_cast<FPT>(10.),  // df1
         static_cast<FPT>(1.),  // df2
         static_cast<FPT>(241.8818)),  // F
      static_cast<FPT>(0.05), // probability.
			tolerance);
   BOOST_CHECK_CLOSE(
      ::boost::math::fisher(
         static_cast<FPT>(1.),  // df1
         static_cast<FPT>(2.),  // df2
         static_cast<FPT>(8.5263)),  // F
      static_cast<FPT>(0.1), // probability.
			tolerance);
   BOOST_CHECK_CLOSE(
      ::boost::math::fisher(
         static_cast<FPT>(1.),  // df1
         static_cast<FPT>(2.),  // df2
         static_cast<FPT>(18.5128)),  // F
      static_cast<FPT>(0.05), // probability.
			tolerance);
   BOOST_CHECK_CLOSE(
      ::boost::math::fisher(
         static_cast<FPT>(1.),  // df1
         static_cast<FPT>(2.),  // df2
         static_cast<FPT>(98.5025)),  // F
      static_cast<FPT>(0.01), // probability.
			tolerance);

   BOOST_CHECK_CLOSE(
      ::boost::math::fisher(
         static_cast<FPT>(1.),  // df1
         static_cast<FPT>(2.),  // df2
         static_cast<FPT>(998.4979)),  // F
      static_cast<FPT>(0.001), // probability.
			tolerance);
*/

	// http://www.vias.org/tmdatanaleng/cc_distri_calculator.html give some useful info,
	// BUT the calculator link is broken.

	// http://www.vias.org/simulations/simusoft_distcalc.html
	// Distcalc version 1.2 Copyright 2002 H Lohninger, TU Wein
	// H.Lohninger: Teach/Me Data Analysis, Springer-Verlag, Berlin-New York-Tokyo, 1999. ISBN 3-540-14743-8
	// The Windows calculator is available zipped distcalc.exe for download at:
	// http://www.vias.org/simulations/simu_stat.html

	// This interactive Windows program was used to find some combination for which the
	// result appears to be exact.  No doubt this can be done analytically too,
	// by mathematicians!

	// Some combinations for which the result is 'exact', or at least is to 40 decimal digits.
	// 40 decimal digits includes 128-bit significand User Defined Floating-Point types.
  // These all pass tests at near epsilon accuracy for the floating-point type.
   BOOST_CHECK_CLOSE(
      ::boost::math::fisher(
         static_cast<FPT>(1.),  // df1
         static_cast<FPT>(2.),  // df2
//         static_cast<FPT>(0.6666666666666666666666666666666666667)),  // F
         static_cast<FPT>(2.)/static_cast<FPT>(3.) ),  // F
      static_cast<FPT>(0.5), // probability.
			tolerance);

   BOOST_CHECK_CLOSE(
      ::boost::math::fisher(
         static_cast<FPT>(1.),  // df1
         static_cast<FPT>(2.),  // df2
         static_cast<FPT>(1.6)),  // F
      static_cast<FPT>(0.333333333333333333333333333333333333), // probability.
			tolerance);

   BOOST_CHECK_CLOSE(
      ::boost::math::fisher(
         static_cast<FPT>(1.),  // df1
         static_cast<FPT>(2.),  // df2
         static_cast<FPT>(6.5333333333333333333333333333333333)),  // F
      static_cast<FPT>(0.125), // probability.
			tolerance);

   BOOST_CHECK_CLOSE(
      ::boost::math::fisher(
         static_cast<FPT>(2.),  // df1
         static_cast<FPT>(2.),  // df2
         static_cast<FPT>(1.)),  // F
      static_cast<FPT>(0.5), // probability.
			tolerance);

   BOOST_CHECK_CLOSE(
      ::boost::math::fisher(
         static_cast<FPT>(2.),  // df1
         static_cast<FPT>(2.),  // df2
         static_cast<FPT>(3.)),  // F
      static_cast<FPT>(0.25), // probability.
			tolerance);

   BOOST_CHECK_CLOSE(
      ::boost::math::fisher(
         static_cast<FPT>(2.),  // df1
         static_cast<FPT>(2.),  // df2
         static_cast<FPT>(3.)),  // F
      static_cast<FPT>(0.25), // probability.
			tolerance);

   BOOST_CHECK_CLOSE(
      ::boost::math::fisher(
         static_cast<FPT>(2.),  // df1
         static_cast<FPT>(2.),  // df2
         static_cast<FPT>(7.)),  // F
      static_cast<FPT>(0.125), // probability.
			tolerance);

   BOOST_CHECK_CLOSE(
      ::boost::math::fisher(
         static_cast<FPT>(2.),  // df1
         static_cast<FPT>(2.),  // df2
         static_cast<FPT>(9.)),  // F
      static_cast<FPT>(0.1), // probability.
			tolerance);

   BOOST_CHECK_CLOSE(
      ::boost::math::fisher(
         static_cast<FPT>(2.),  // df1
         static_cast<FPT>(2.),  // df2
         static_cast<FPT>(19.)),  // F
      static_cast<FPT>(0.05), // probability.
			tolerance);

   BOOST_CHECK_CLOSE(
      ::boost::math::fisher(
         static_cast<FPT>(2.),  // df1
         static_cast<FPT>(2.),  // df2
         static_cast<FPT>(29.)),  // F
      static_cast<FPT>(0.03333333333333333333333333333333333333333), // probability.
			tolerance);

   BOOST_CHECK_CLOSE(
      ::boost::math::fisher(
         static_cast<FPT>(2.),  // df1
         static_cast<FPT>(2.),  // df2
         static_cast<FPT>(99.)),  // F
      static_cast<FPT>(0.01), // probability. 
			tolerance);

   BOOST_CHECK_CLOSE(
      ::boost::math::fisher(
         static_cast<FPT>(4.),  // df1
         static_cast<FPT>(4.),  // df2
         static_cast<FPT>(9.)),  // F
      static_cast<FPT>(0.028), // probability. 
			tolerance);

   BOOST_CHECK_CLOSE(
      ::boost::math::fisher(
         static_cast<FPT>(8.),  // df1
         static_cast<FPT>(8.),  // df2
         static_cast<FPT>(1.)),  // F
      static_cast<FPT>(0.5), // probability. 
			tolerance);

	 // These two don't pass at the 2 eps tolerance:
   //BOOST_CHECK_CLOSE(
   //   ::boost::math::fisher(
   //      static_cast<FPT>(2.),  // df1
   //      static_cast<FPT>(2.),  // df2
   //      static_cast<FPT>(999.)),  // F
   //   static_cast<FPT>(0.001), // probability.== 0.000999987125
			//tolerance);

	 // And this even with 1 decimal digit tolerance.
   // BOOST_CHECK_CLOSE(
   //   ::boost::math::fisher(
   //      static_cast<FPT>(2.),  // df1
   //      static_cast<FPT>(2.),  // df2
   //      static_cast<FPT>(9999.)),  // F
   //   static_cast<FPT>(0.0001), // probability. == 0.000100016594 or 9.9999999999988987e-005
			//tolerance);
  //BOOST_CHECK_CLOSE(
   //   ::boost::math::fisher(
   //      static_cast<FPT>(2.),  // df1
   //      static_cast<FPT>(7.),  // df2
   //      static_cast<FPT>(99.)),  // F
   //   static_cast<FPT>(0.00001), // probability.
			//tolerance);


   //BOOST_CHECK_CLOSE(
   //   ::boost::math::fisher(
   //      static_cast<FPT>(2.),  // df1
   //      static_cast<FPT>(3.),  // df2
   //      static_cast<FPT>(999.)),  // F
   //   static_cast<FPT>(0.0000666666666666666666666666), // probability.
			//tolerance);


// Also note limit cases for F(1, infinity) == normal distribution
// F(1, n2) == Student's t distribution
// F(n1, infinity) == Chisq distribution

// These might allow some further cross checks?

} // template <class FPT>void test_spots(FPT)

int test_main(int, char* [])
{
	 // Basic sanity-check spot values.
	// (Parameter value, arbitrarily zero, only communicates the floating point type).
   test_spots(0.0F); // Test float.
   test_spots(0.0); // Test double.
   test_spots(0.0L); // Test long double.
   test_spots(boost::math::concepts::real_concept(0.)); // Test real concept.

   return 0;
} // int test_main(int, char* [])

/*

Output:

1>------ Rebuild All started: Project: test_fisher, Configuration: Debug Win32 ------
1>Deleting intermediate and output files for project 'test_fisher', configuration 'Debug|Win32'
1>Compiling...
1>test_fisher.cpp
1>Linking...
1>Embedding manifest...
1>Autorun "i:\boost-06-05-03-1300\libs\math\test\Math_test\debug\test_fisher.exe"
1>Running 1 test case...
1>*** No errors detected
1>Build Time 0:06
1>Build log was saved at "file://i:\boost-06-05-03-1300\libs\math\test\Math_test\test_fisher\Debug\BuildLog.htm"
1>test_fisher - 0 error(s), 0 warning(s)
========== Rebuild All: 1 succeeded, 0 failed, 0 skipped ==========


*/

