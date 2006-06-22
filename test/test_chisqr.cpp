// test_chisqr.cpp

// Copyright Paul A. Bristow 2006.

// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)


// Basic sanity test for chisqr Cumulative Distribution Function.

#define BOOST_MATH_THROW_ON_DOMAIN_ERROR

#ifdef _MSC_VER
#  pragma warning(disable: 4127) // conditional expression is constant.
#  pragma warning(disable: 4100) // unreferenced formal parameter.
#endif

#include <boost/math/special_functions/chisqr.hpp>
	using ::boost::math::chisqr;
#include <boost/math/concepts/real_concept.hpp> // for real_concept
	using ::boost::math::concepts::real_concept;

#include <boost/test/included/test_exec_monitor.hpp> // for test_main
#include <boost/test/floating_point_comparison.hpp> // for BOOST_CHECK_CLOSE

#include <limits>
  using std::numeric_limits;

template <class FPT> // Any floating-point type FPT.
void test_spots(FPT)
{
  // Basic sanity checks, tolerance is about numeric_limits<FPT>::digits10 decimal places,
	// guaranteed for type FPT, eg 6 for float, 15 for double,
	// expressed as a percentage (so -2) for BOOST_CHECK_CLOSE,

	int decdigits = numeric_limits<FPT>::digits10;
	decdigits -= 0; // Perhaps allow some decimal digits margin of numerical error.
	FPT tolerance = static_cast<FPT>(std::pow(10., -(decdigits-2))); // 1e-6 (-2 so as %)
	tolerance *= 1; // Allow some bit(s) small margin (2 means + or - 1 bit) of numerical error.
	// Typically 2e-13% = 2e-15 as fraction for double.

	// Sources of spot test values:

	// http://www.vias.org/tmdatanaleng/cc_distri_calculator.html give some useful info,
	// BUT the calculator link is broken.

	// http://www.vias.org/simulations/simusoft_distcalc.html
	// Distcalc version 1.2 Copyright 2002 H Lohninger, TU Wein
	// H.Lohninger: Teach/Me Data Analysis, Springer-Verlag, Berlin-New York-Tokyo, 1999. ISBN 3-540-14743-8
	// The Windows calculator is available zipped distcalc.exe for download at:
	// http://www.vias.org/simulations/simu_stat.html

	// Also usd a Java Chi Square calculator at:
	// http://www.stat.sc.edu/~west/applets/chisqdemo.html
	// but seems even less accurate.

	// Many be some combinations for which the result is 'exact', or at least is to 40 decimal digits.
	// 40 decimal digits includes 128-bit significand User Defined Floating-Point types,
	// but these are as yet undiscovered for chi-sqr.
	
	// Best source of accurate values is
	// Mathworld online calculator (40 decimal digits precision, suitable for up to 128-bit significands)
	// http://functions.wolfram.com/webMathematica/FunctionEvaluation.jsp?name=GammaRegularized
	// GammaRegularized is same as gamma incomplete, gamma or gamma_Q(a, x) or Q(a, z).

   BOOST_CHECK_CLOSE(
      chisqr(
         static_cast<FPT>(5.),  // degrees_of_freedom - as floating-point.
         static_cast<FPT>(11.0705)),  // chisqr
      static_cast<FPT>(0.04999995542804363455469285561840275602943), // probability.
			// Source Mathworld online calculator (40 decimal digits precision)
			// Q(5.0/2, 11.0705/2) = 0.049999955428043634554692855618402756029434129418384
			tolerance);

   BOOST_CHECK_CLOSE(
      chisqr(
         static_cast<FPT>(4.),  // degrees_of_freedom - as floating-point.
         static_cast<FPT>(4.)),  // chisqr
      static_cast<FPT>(0.4060058497098380756819984849174532102229), // probability.
			// Q(4/2, 4/2) = 0.049999955428043634554692855618402756029434129418384
			tolerance);

   BOOST_CHECK_CLOSE(
      chisqr(
         static_cast<FPT>(1),  // degrees_of_freedom - as floating-point.
         static_cast<FPT>(10.)),  // chisqr
      static_cast<FPT>(0.001565402258002549677499803978385902310435), // probability.
			tolerance);

   BOOST_CHECK_CLOSE(
      chisqr(
         static_cast<FPT>(10.),  // degrees_of_freedom - as floating-point.
         static_cast<FPT>(1.)),  // chisqr
      static_cast<FPT>(0.9998278843700441592218882959620240287207), // probability.
			tolerance);

   BOOST_CHECK_CLOSE(
      chisqr(
         static_cast<FPT>(100.),  // degrees_of_freedom - as floating-point.
         static_cast<FPT>(5.)),  // chisqr
      static_cast<FPT>(0.999999999999999999814527311613020069945), // probability.
			tolerance);
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

Running 1 test case...

*** No errors detected
Press any key to continue . . .

*/