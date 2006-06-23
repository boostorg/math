// test_error_handling.cpp

// Demo how error handling.

// Copyright Paul A. Bristow 2006.

// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_SPECIAL_ERROR_HPP
#define BOOST_MATH_SPECIAL_ERROR_HPP

#define BOOST_MATH_THROW_ON_DOMAIN_ERROR

#ifdef _MSC_VER
#  pragma warning(disable: 4127) // conditional expression is constant.
#  pragma warning(disable: 4512) // assignment operator could not be generated.
#  pragma warning(disable: 4996) // std::char_traits<char>::copy' was declared deprecated.
#endif

// Boost
#include <boost/math/tools/error_handling.hpp> // for domain_error.
	using ::boost::math::tools::domain_error;
	using ::boost::math::tools::pole_error;
	using ::boost::math::tools::overflow_error;
	using ::boost::math::tools::underflow_error;
	using ::boost::math::tools::denorm_error;
#include <boost/math/special_functions/fpclassify.hpp>
	using boost::math::isnan;
// std
#include <iostream>
	using std::cout;
	using std::endl;
#include <limits>
  using std::numeric_limits;
#include <stdexcept>
	using std::exception;

// Forward declaration of error handling functions.
namespace boost
{
	namespace math
	{
		namespace tools
		{
			// Math Error handling:
			// error is always set to a value in errno.h (ANSI 4.5.1),
			// for example to EDOM or ERANGE. (from cerrno)
			// (except for denormalized value).

			// Domain error: An argument is outside it's allowed range.

			// If  BOOST_MATH_THROW_ON_DOMAIN_ERROR defined
			//   raise error to throw,
			// else
			//   try to return a NaN (but throw is there isn't a NaN for the floating-point type).

			template <class T>
			inline T domain_error(const char* function, const char* message);

			template <class T> // Optionally also a value, usually the offending one,
			// uses Boost.Format so include %1% in message for formatted output.
			inline T domain_error(const char* function, const char* message, const T& value);

			// Evaluation at a pole, this is currently treated the same as a domain error:
			template <class T>
			inline T pole_error(const char* function, const char* message);

			// Result too large to be represented in type T:
			// return infinity, or throw if BOOST_MATH_THROW_ON_OVERFLOW_ERROR defined.
			template <class T>
			inline T overflow_error(const char* function, const char* message);

			// Result too small to be represented in type T, called only when we know the result is not actually zero:
			// return zero, or throw if BOOST_MATH_THROW_ON_UNDERFLOW_ERROR defined.

			template <class T>
			inline T underflow_error(const char* function, const char* message);

			// Result is denormalised:
			template <class T>
			inline T denorm_error(T const& t, const char* function, const char* message);
			// return denormalized value or throw if BOOST_MATH_THROW_ON_DENORM_ERROR defined.
		} // namespace tools
	} // namespace math
} // 	namespace boost

#include <boost/test/included/test_exec_monitor.hpp> // for test_main
#include <boost/test/floating_point_comparison.hpp> // for BOOST_CHECK_CLOSE

#include <boost/math/concepts/real_concept.hpp> // for real_concept
using ::boost::math::concepts::real_concept;

template <class FPT> // Any floating-point type FPT.
FPT test_error(FPT)
{
	cout << "Current function is " << BOOST_CURRENT_FUNCTION << endl;
	using boost::math::tools::denorm_error;

	// FPT r = domain_error<FPT>(BOOST_CURRENT_FUNCTION, "Test domain_error message!");
	// cout << "r = " << r << endl; // r = 1.#QNAN
	// BOOST_ASSERT(isnan(r));
	// BOOST_ASSERT(r == numeric_limits<FPT>::quiet_NaN()); // always fails!!!
	
	FPT x = static_cast<FPT>(1.2345678901234567890123456789012345678980);
	FPT r = domain_error<FPT>(BOOST_CURRENT_FUNCTION, "Test domain_error message, argument = %1%", x);
	return r;
	// Output from test_error(0.0F); is:
	// Current function is float __cdecl test_error<float>(float)
	// unknown location(0):
	// fatal error in "test_main_caller( argc, argv )":
	// std::domain_error: Error in function float __cdecl test_error<float>(float):
	// Test domain_error message! argument = 1.23456788
	// Note 9 decimal digits, max_digits10 for float (17 for double, 33 for quad_float...)

} // template <class FPT>void test_error(FPT)

int test_main(int, char* [])
{
	// Test error handling.
	// (Parameter value, arbitrarily zero, only communicates the floating point type FPT).
	// test_error(0.0F); // Test float.
		test_error(0.0F); // Test float.

	try
	{
		test_error(0.0F); // Test float.
	}
	catch (exception ex)
	{
		cout << "Caught exception: " << ex.what() << std::endl;
		// Current function is float __cdecl test_error<float>(float)
		// Caught exception: Unknown exception

	}
	//test_error(0.0); // Test double.
	//test_error(0.0L); // Test long double.
	//test_error(boost::math::concepts::real_concept(0.)); // Test real concept.

	return 0;
} // int test_main(int, char* [])

#endif // BOOST_MATH_SPECIAL_ERROR_HPP


/*

Output:

Running 1 test case...
unknown location(0): fatal error in "test_main_caller( argc, argv )":
std::domain_error: Error in function void __cdecl test_error<float>(float):
Test domain error message!

*** 1 failure detected in test suite "Test Program"
Press any key to continue . . .


Running 1 test case...
Current function is float __cdecl test_error<float>(float)
Caught exception: Unknown exception

*** No errors detected
Press any key to continue . . .


Autorun "i:\boost-06-05-03-1300\libs\math\test\Math_test\debug\test_error_handling.exe"
Running 1 test case...
Current function is float __cdecl test_error<float>(float)
unknown location(0): fatal error in "test_main_caller( argc, argv )": std::domain_error: Error in function float __cdecl test_error<float>(float): Test domain_error message, argument = 1.23456788
*** 1 failure detected in test suite "Test Program"
Project : error PRJ0019: A tool returned an error code from "Autorun "i:\boost-06-05-03-1300\libs\math\test\Math_test\debug\test_error_handling.exe""




*/