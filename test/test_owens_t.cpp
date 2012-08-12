// Copyright Paul A. Bristow 2012.
// Copyright Benjamin Sobotta 2012.
// Copyright John Maddock 2012.

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

// test_owens_t.cpp

#ifdef _MSC_VER
#  pragma warning (disable : 4127) // conditional expression is constant
#endif

#include <boost/math/concepts/real_concept.hpp> // for real_concept.
using ::boost::math::concepts::real_concept;

#include <boost/math/special_functions/owens_t.hpp> // for owens_t function.
using boost::math::owens_t;

#include <boost/test/test_exec_monitor.hpp> // for test_main.
#include <boost/test/floating_point_comparison.hpp> // for BOOST_CHECK_CLOSE and  BOOST_CHECK_CLOSE_fraction.

#include <iostream>
using std::cout;
using std::endl;
#include <limits>
using std::numeric_limits;

template <class RealType>
void test_spot(
     RealType h,    // 
     RealType a,    // Chi Square statistic
     RealType tol)   // Test tolerance
{
  BOOST_CHECK_CLOSE_FRACTION(owens_t(h, a), 3.89119302347013668966224771378e-2, tol);
}


template <class RealType> // Any floating-point type RealType.
void test_spots(RealType)
{
  // Basic sanity checks, test data is as accurate as long double,
  // so set tolerance to a few epsilon expressed as a fraction.
  RealType tolerance = boost::math::tools::epsilon<RealType>() * 3;
  cout << "Tolerance = " << tolerance << "." << endl;

  using  ::boost::math::owens_t;

    using namespace std; // ADL of std names.
 
   BOOST_CHECK_CLOSE_FRACTION(owens_t(static_cast<RealType>(0.0625L), static_cast<RealType>(0.25L)), static_cast<RealType>(3.89119302347013668966224771378e-2L), tolerance);


} // template <class RealType>void test_spots(RealType)

int test_main(int, char* [])
{
  BOOST_MATH_CONTROL_FP;

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
  return 0;
} // int test_main(int, char* [])

/*

Output:

test_owens_t.cpp
  Generating code
  Finished generating code
  test_owens_t.vcxproj -> J:\Cpp\MathToolkit\test\Math_test\Release\test_owens_t.exe
  Running 1 test case...
  Tolerance = 6.66134e-016.
  Tolerance = 6.66134e-016.
  Tolerance = 6.66134e-016.
  
  *** No errors detected

*/



