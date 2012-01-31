// test_owens_t.cpp

// Copyright Paul A. Bristow 2012.
// Copyright Benjamin Sobotta 2012.

// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

// Tested using some 30 decimal digit accuracy values from:
// Fast and accurate calculation of Owen's T-function
// Mike Patefield, and David Tandy
// Journal of Statistical Software, 5 (5), 1-25 (2000).
// http://www.jstatsoft.org/v05/a05/paper  Table 3, page 15
// Values of T(h,a) accurate to thirty figures were calculated using 128 bit arithmetic by
// evaluating (9) with m = 48, the summation over k being continued until additional terms did
// not alter the result. The resultant values Tacc(h,a) say, were validated by evaluating (8) with
// m = 48 (i.e. 96 point Gaussian quadrature).

#ifdef _MSC_VER
#  pragma warning (disable : 4127) // conditional expression is constant
#  pragma warning (disable : 4305) // 'initializing' : truncation from 'double' to 'const float'
// ?? TODO get rid of these warnings?
#endif

#include <boost/math/concepts/real_concept.hpp> // for real_concept.
using ::boost::math::concepts::real_concept;

#include <boost/math/special_functions/owens_t.hpp> // for owens_t function.
using boost::math::owens_t;

#include <boost/test/test_exec_monitor.hpp> // for test_main.
#include <boost/test/floating_point_comparison.hpp> // for BOOST_CHECK_CLOSE and  BOOST_CHECK_CLOSE_fraction.
#include <boost/array.hpp>

#include "libs/math/test/handle_test_result.hpp"
#include "libs/math/test/table_type.hpp"
#include "libs/math/test/functor.hpp"

#include <iostream>
using std::cout;
using std::endl;
#include <limits>
using std::numeric_limits;

void expected_results()
{
   //
   // Define the max and mean errors expected for
   // various compilers and platforms.
   //
   const char* largest_type;
#ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
   if(boost::math::policies::digits<double, boost::math::policies::policy<> >() == boost::math::policies::digits<long double, boost::math::policies::policy<> >())
   {
      largest_type = "(long\\s+)?double";
   }
   else
   {
      largest_type = "long double";
   }
#else
   largest_type = "(long\\s+)?double";
#endif

   //
   // Catch all cases come last:
   //
   add_expected_result(
      ".*",                          // compiler
      ".*",                          // stdlib
      ".*",                          // platform
      largest_type,                  // test type(s)
      ".*",      // test data group
      "boost::math::owens_t", 10, 5);  // test function
   add_expected_result(
      ".*",                          // compiler
      ".*",                          // stdlib
      ".*",                          // platform
      "real_concept",                  // test type(s)
      ".*",      // test data group
      "boost::math::owens_t", 10, 5);  // test function
   //
   // Finish off by printing out the compiler/stdlib/platform names,
   // we do this to make it easier to mark up expected error rates.
   //
   std::cout << "Tests run with " << BOOST_COMPILER << ", " 
      << BOOST_STDLIB << ", " << BOOST_PLATFORM << std::endl;
}


template <class RealType>
void test_spot(
     RealType h,    //
     RealType a,    //
     RealType tol)   // Test tolerance
{
  BOOST_CHECK_CLOSE_FRACTION(owens_t(h, a), 3.89119302347013668966224771378e-2L, tol);
}


template <class RealType> // Any floating-point type RealType.
void test_spots(RealType)
{
  // Basic sanity checks, test data is as accurate as long double,
  // so set tolerance to a few epsilon expressed as a fraction.
  RealType tolerance = boost::math::tools::epsilon<RealType>() * 30; // most OK with 3 eps tolerance.
  cout << "Tolerance = " << tolerance << "." << endl;

  using  ::boost::math::owens_t;
  using namespace std; // ADL of std names.

  // Checks of six sub-methods T1 to T6.
  BOOST_CHECK_CLOSE_FRACTION(owens_t(static_cast<RealType>(0.0625L), static_cast<RealType>(0.25L)), static_cast<RealType>(3.89119302347013668966224771378e-2L), tolerance);  // T1
  BOOST_CHECK_CLOSE_FRACTION(owens_t(static_cast<RealType>(6.5L), static_cast<RealType>(0.4375L)), static_cast<RealType>(2.00057730485083154100907167685E-11L), tolerance); // T2
  BOOST_CHECK_CLOSE_FRACTION(owens_t(static_cast<RealType>(7L), static_cast<RealType>( 0.96875L)), static_cast<RealType>(6.39906271938986853083219914429E-13L), tolerance); // T3
  BOOST_CHECK_CLOSE_FRACTION(owens_t(static_cast<RealType>(4.78125L), static_cast<RealType>(0.0625L)), static_cast<RealType>(1.06329748046874638058307112826E-7L), tolerance); // T4
  BOOST_CHECK_CLOSE_FRACTION(owens_t(static_cast<RealType>(2.L), static_cast<RealType>(0.5L)), static_cast<RealType>(8.62507798552150713113488319155E-3L), tolerance); // T5
  BOOST_CHECK_CLOSE_FRACTION(owens_t(static_cast<RealType>(1.L), static_cast<RealType>(0.9999975L)), static_cast<RealType>(6.67418089782285927715589822405E-2L), tolerance); // T6
  //BOOST_CHECK_CLOSE_FRACTION(owens_t(static_cast<RealType>(L), static_cast<RealType>(L)), static_cast<RealType>(L), tolerance);

  // Mathematica spot values.
  BOOST_CHECK_CLOSE_FRACTION(owens_t(static_cast<RealType>(0.4375L), static_cast<RealType>(6.5L)), static_cast<RealType>(0.1654013012544939L), tolerance);
  BOOST_CHECK_CLOSE_FRACTION(owens_t(static_cast<RealType>(0.96875L), static_cast<RealType>(7.L)), static_cast<RealType>(0.0831674847460297L), tolerance);
  BOOST_CHECK_CLOSE_FRACTION(owens_t(static_cast<RealType>(0.0625L), static_cast<RealType>(4.78125L)), static_cast<RealType>(0.2157118581989799L), tolerance);
  BOOST_CHECK_CLOSE_FRACTION(owens_t(static_cast<RealType>(0.5L), static_cast<RealType>(2.L)), static_cast<RealType>(0.1415806036539784L), tolerance);


//   BOOST_CHECK_CLOSE_FRACTION(owens_t(static_cast<RealType>(L), static_cast<RealType>(L)), static_cast<RealType>(L), tolerance);

  // Spots values using Mathematica but only 17 decimal digits.
  BOOST_CHECK_CLOSE_FRACTION(owens_t(static_cast<RealType>(6.5L), static_cast<RealType>(0.4375L)), static_cast<RealType>(2.000577304850829E-11L), tolerance);
  BOOST_CHECK_CLOSE_FRACTION(owens_t(static_cast<RealType>(0.4375L), static_cast<RealType>(6.5L)), static_cast<RealType>(0.1654013012544939L), tolerance);
  BOOST_CHECK_CLOSE_FRACTION(owens_t(static_cast<RealType>(7.L), static_cast<RealType>(0.96875L)), static_cast<RealType>(6.399062719389912e-13L), tolerance);
  BOOST_CHECK_CLOSE_FRACTION(owens_t(static_cast<RealType>(0.96875L), static_cast<RealType>(7.L)), static_cast<RealType>(0.0831674847460297L), tolerance);
  BOOST_CHECK_CLOSE_FRACTION(owens_t(static_cast<RealType>(4.78125L), static_cast<RealType>(0.0625L)), static_cast<RealType>(1.063297480468747e-7L), tolerance);
  BOOST_CHECK_CLOSE_FRACTION(owens_t(static_cast<RealType>(0.0625L), static_cast<RealType>(4.78125L)), static_cast<RealType>(0.2157118581989799L), tolerance);
  BOOST_CHECK_CLOSE_FRACTION(owens_t(static_cast<RealType>(2.L), static_cast<RealType>(0.5L)), static_cast<RealType>(0.00862507798552151L), tolerance);
  BOOST_CHECK_CLOSE_FRACTION(owens_t(static_cast<RealType>(0.5L), static_cast<RealType>(2L)), static_cast<RealType>(0.1415806036539784L), tolerance);


} // template <class RealType>void test_spots(RealType)

template <class Real, class T>
void do_test_owens_t(const T& data, const char* type_name, const char* test_name)
{
   typedef typename T::value_type row_type;
   typedef Real                   value_type;

   typedef value_type (*pg)(value_type, value_type);
#if defined(BOOST_MATH_NO_DEDUCED_FUNCTION_POINTERS)
   pg funcp = boost::math::owens_t<value_type>;
#else
   pg funcp = boost::math::owens_t;
#endif

   boost::math::tools::test_result<value_type> result;

   std::cout << "Testing " << test_name << " with type " << type_name
      << "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";

   //
   // test hermite against data:
   //
   result = boost::math::tools::test_hetero<Real>(
      data, 
      bind_func<Real>(funcp, 0, 1), 
      extract_result<Real>(2));
   handle_test_result(result, data[result.worst()], result.worst(), type_name, "boost::math::owens_t", test_name);

   std::cout << std::endl;
}

template <class T>
void test_owens_t(T, const char* name)
{
   //
   // The actual test data is rather verbose, so it's in a separate file
   //
   // The contents are as follows, each row of data contains
   // three items, input value a, input value b and erf(a, b):
   // 
#  include "owens_t.ipp"

   do_test_owens_t<T>(owens_t, name, "Owens T (medium small values)");

#include "owens_t_large_data.ipp"

   do_test_owens_t<T>(owens_t_large_data, name, "Owens T (large and diverse values)");
}


int test_main(int, char* [])
{
  BOOST_MATH_CONTROL_FP;

  expected_results();

  // Basic sanity-check spot values.

  // (Parameter value, arbitrarily zero, only communicates the floating point type).
  test_spots(0.0F); // Test float.
  test_spots(0.0); // Test double.
#ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
  test_spots(0.0L); // Test long double.
#if !BOOST_WORKAROUND(__BORLANDC__, BOOST_TESTED_AT(0x582))
  test_spots(boost::math::concepts::real_concept(0.)); // Test real concept.
#endif
#endif

  test_owens_t(0.0F, "float"); // Test float.
  test_owens_t(0.0, "double"); // Test double.
#ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
  test_owens_t(0.0L, "long double"); // Test long double.
#if !BOOST_WORKAROUND(__BORLANDC__, BOOST_TESTED_AT(0x582))
  test_owens_t(boost::math::concepts::real_concept(0.), "real_concept"); // Test real concept.
#endif
#endif
  return 0;
} // int test_main(int, char* [])

/*

Output:
------ Rebuild All started: Project: test_owens_t_sandbox, Configuration: Release Win32 ------
  test_owens_t.cpp
  Generating code
  Finished generating code
  test_owens_t_sandbox.vcxproj -> J:\Cpp\MathToolkit\test\Math_test\Release\test_owens_t_sandbox.exe
  Running 1 test case...
  Platform: Win32
  Compiler: Microsoft Visual C++ version 10.0
  STL     : Dinkumware standard library version 520
  Boost   : 1.49.0
  Tolerance = 3.57628e-007.
  Tolerance = 6.66134e-016.
  Tolerance = 6.66134e-016.
  Tolerance = 6.66134e-016.

  *** No errors detected
========== Rebuild All: 1 succeeded, 0 failed, 0 skipped ==========


------ Build started: Project: test_owens_t_sandbox, Configuration: Debug Win32 ------
  test_owens_t.cpp
  test_owens_t_sandbox.vcxproj -> J:\Cpp\MathToolkit\test\Math_test\Debug\test_owens_t_sandbox.exe
  Running 1 test case...
  Platform: Win32
  Compiler: Microsoft Visual C++ version 10.0
  STL     : Dinkumware standard library version 520
  Boost   : 1.49.0
  Tolerance = 3.57628e-006.
  Tolerance = 6.66134e-015.
  Tolerance = 6.66134e-015.
  Tolerance = 6.66134e-015.

  *** No errors detected
========== Build: 1 succeeded, 0 failed, 0 up-to-date, 0 skipped ==========

*/



