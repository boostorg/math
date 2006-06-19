// Copyright Paul A. Bristow 2006.

// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

// test_students_t.cpp

// http://en.wikipedia.org/wiki/Student%27s_t_distribution

// Basic sanity test for Student's t probability (quintile) (0. < p < 1).

#ifdef _MSC_VER
#  pragma warning(disable: 4127) // conditional expression is constant.
#  pragma warning(disable: 4100) // unreferenced formal parameter.
#endif

#include <boost/math/special_functions/students_t.hpp>
#include <boost/math/concepts/real_concept.hpp> // for real_concept

#include <boost/test/included/test_exec_monitor.hpp> // Test
#include <boost/test/floating_point_comparison.hpp>

template <class FPT>
void test_spots(FPT)
{
   //
   // Basic sanity checks, tolerance is 6 decimal places
	 // expressed as a percentage (so -2) for BOOST_CHECK_CLOSE,
   // One check per domain of the implementation:

	FPT tolerance = static_cast<FPT>(std::pow(10., -(6-2))); // 1e-6 (as %)
	 // This test only passes at 1e-5 because probability value is less accurate,
	// a digit in 6th decimla place, although calculated using 
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
	// and the t- values are shown.


   BOOST_CHECK_CLOSE(
      ::boost::math::students_t(
         static_cast<FPT>(2.),  // df
         static_cast<FPT>(-6.96455673428326)),  // t
      static_cast<FPT>(0.01), // probability.
			tolerance);
   BOOST_CHECK_CLOSE(
      ::boost::math::students_t(
         static_cast<FPT>(5.),  // df
         static_cast<FPT>(-3.36492999890721)),  // t
      static_cast<FPT>(0.01), // probability.
			tolerance);

   BOOST_CHECK_CLOSE(
      ::boost::math::students_t(
         static_cast<FPT>(1.),  // df
         static_cast<FPT>(-31830.988607907)),  // t
      static_cast<FPT>(0.00001), // probability.
			tolerance);

   BOOST_CHECK_CLOSE(
      ::boost::math::students_t(
         static_cast<FPT>(25.),  // df
         static_cast<FPT>(-5.2410429995425)),  // t
      static_cast<FPT>(0.00001), // probability.
			tolerance);

   BOOST_CHECK_CLOSE(
      ::boost::math::students_t(
         static_cast<FPT>(1.),  // df
         static_cast<FPT>(-63661.97723)),  // t
      static_cast<FPT>(0.000005), // probability.
			tolerance);

    BOOST_CHECK_CLOSE(
      ::boost::math::students_t(
         static_cast<FPT>(5.),  // df
         static_cast<FPT>(-17.89686614)),  // t
      static_cast<FPT>(0.000005), // probability.
			tolerance);

    BOOST_CHECK_CLOSE(
      ::boost::math::students_t(
         static_cast<FPT>(25.),  // df
         static_cast<FPT>(-5.510848412)),  // t
      static_cast<FPT>(0.000005), // probability.
			tolerance);

  BOOST_CHECK_CLOSE(
      ::boost::math::students_t(
         static_cast<FPT>(10.),  // df
         static_cast<FPT>(-1.812461123)),  // t
      static_cast<FPT>(0.05), // probability.
			tolerance);

  BOOST_CHECK_CLOSE(
      ::boost::math::students_t(
         static_cast<FPT>(10.),  // df
         static_cast<FPT>(1.812461123)),  // t
      static_cast<FPT>(0.95), // probability.
			tolerance);

  BOOST_CHECK_CLOSE(
      ::boost::math::students_t(
         static_cast<FPT>(10.),  // df
         static_cast<FPT>(9.751995491)),  // t
      static_cast<FPT>(0.999999), // probability.
			tolerance);

  BOOST_CHECK_CLOSE(
      ::boost::math::students_t(
         static_cast<FPT>(10.),  // df - for ALL df!
         static_cast<FPT>(0.)),  // t
      static_cast<FPT>(0.5), // probability.
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



1>------ Build started: Project: test_students_t, Configuration: Debug Win32 ------
1>Compiling...
1>test_students_t.cpp
1>MSVC++ compiler Version 8.0
1>Linking...
1>Embedding manifest...
1>Autorun "i:\boost-06-05-03-1300\libs\math\test\Math_test\debug\test_students_t.exe"
1>Running 1 test case...
1>*** No errors detected
1>Build Time 0:05
1>Build log was saved at "file://i:\boost-06-05-03-1300\libs\math\test\Math_test\test_students_t\Debug\BuildLog.htm"
1>test_students_t - 0 error(s), 0 warning(s)
========== Build: 1 succeeded, 0 failed, 0 up-to-date, 0 skipped ==========


*/

