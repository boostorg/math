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

#include <boost/math/constants/constants.hpp>
#include <boost/math/tools/test.hpp> 

#include <iostream>
	using std::cout;
	using std::endl;
	using std::setprecision;
#include <limits>
  using std::numeric_limits;

template <class RealType>
void test_spots(RealType T)
{
   // Basic santity checks for constants.

	RealType tolerance = static_cast<RealType>(2e-15);  // double
	cout << "Tolerance for type " << typeid(T).name()  << " is " << tolerance << "." << endl;

	 using namespace boost::math::constants;
   using namespace std; // Help ADL of std exp, log...
   using std::exp;

   BOOST_CHECK_CLOSE_FRACTION(static_cast<RealType>(3.14159265358979323846264338327950288419716939937510), pi<RealType>(), tolerance); 
   BOOST_CHECK_CLOSE_FRACTION(static_cast<RealType>(sqrt(3.14159265358979323846264338327950288419716939937510)), root_pi<RealType>(), tolerance); 
   BOOST_CHECK_CLOSE_FRACTION(static_cast<RealType>(sqrt(3.14159265358979323846264338327950288419716939937510/2)), root_half_pi<RealType>(), tolerance); 
   BOOST_CHECK_CLOSE_FRACTION(static_cast<RealType>(sqrt(3.14159265358979323846264338327950288419716939937510 * 2)), root_two_pi<RealType>(), tolerance); 
   BOOST_CHECK_CLOSE_FRACTION(static_cast<RealType>(sqrt(log(4.))), root_ln_four<RealType>(), tolerance); 
   BOOST_CHECK_CLOSE_FRACTION(static_cast<RealType>(2.71828182845904523536028747135266249775724709369995), e<RealType>(), tolerance); 
   BOOST_CHECK_CLOSE_FRACTION(static_cast<RealType>(0.5), half<RealType>(), tolerance); 
   BOOST_CHECK_CLOSE_FRACTION(static_cast<RealType>(0.57721566490153286060651209008240243104259335), euler<RealType>(), tolerance); 
   BOOST_CHECK_CLOSE_FRACTION(static_cast<RealType>(sqrt(2.)), root_two<RealType>(), tolerance); 
   BOOST_CHECK_CLOSE_FRACTION(static_cast<RealType>(log(2.)), ln_two<RealType>(), tolerance); 
   BOOST_CHECK_CLOSE_FRACTION(static_cast<RealType>(log(log(2.))), ln_ln_two<RealType>(), tolerance); 
   BOOST_CHECK_CLOSE_FRACTION(static_cast<RealType>(1)/3, third<RealType>(), tolerance); 
   BOOST_CHECK_CLOSE_FRACTION(static_cast<RealType>(2)/3, twothirds<RealType>(), tolerance); 
   BOOST_CHECK_CLOSE_FRACTION(static_cast<RealType>(0.14159265358979323846264338327950288419716939937510), pi_minus_three<RealType>(), tolerance); 
   BOOST_CHECK_CLOSE_FRACTION(static_cast<RealType>(4. - 3.14159265358979323846264338327950288419716939937510), four_minus_pi<RealType>(), tolerance); 
   BOOST_CHECK_CLOSE_FRACTION(static_cast<RealType>(pow((4. - 3.14159265358979323846264338327950288419716939937510), 1.5)), pow23_four_minus_pi<RealType>(), tolerance); 
   BOOST_CHECK_CLOSE_FRACTION(static_cast<RealType>(exp(-0.5)), exp_minus_half<RealType>(), tolerance); 

} // template <class RealType>void test_spots(RealType)

int test_main(int, char* [])
{

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

test_constants.cpp
Linking...
Autorun "i:\boost-06-05-03-1300\libs\math\test\Math_test\debug\test_constants.exe"
Running 1 test case...
Tolerance for type float is 2e-015.
Tolerance for type double is 2e-015.
Tolerance for type long double is 2e-015.
Tolerance for type class boost::math::concepts::real_concept is 2e-015.
*** No errors detected




*/


