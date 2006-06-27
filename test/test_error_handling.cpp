// test_error_handling.cpp

// Test error handling.

// Copyright Paul A. Bristow 2006.
// Copyright John Maddock 2006.

// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

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
	using ::boost::math::tools::logic_error;
#include <boost/math/special_functions/chisqr.hpp> // chisqr
   using ::boost::math::chisqr;
#include <boost/math/special_functions/gamma.hpp> // tgamma
   using ::boost::math::tgamma;
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

#include <boost/test/included/test_exec_monitor.hpp> // for test_main
#include <boost/test/floating_point_comparison.hpp> // for BOOST_CHECK_CLOSE

#include <boost/math/concepts/real_concept.hpp> // for real_concept
   using ::boost::math::concepts::real_concept;

//
// The macro BOOST_CHECK_THROW_MSG is the same as BOOST_CHECK_THROW
// which is to say it verifies that the exception we expect is
// thrown from the function under test, but it also prints the message
// contained in the thrown exception so we can manually inspect the
// quality of the message.
//
// The expanded code is:
//
//   try
//   { 
//      code; 
//   } 
//   catch(const std::exception& e)
//   {
//      std::cout << 
//          "Message from thrown exception was:\n   " << e.what() << std::endl; 
//   }
//   BOOST_CHECK_THROW(code, t);
//
#define BOOST_CHECK_THROW_MSG(code,t)\
   try{ code; } catch(const std::exception& e){\
   std::cout << "Message from thrown exception was:\n   " << e.what() << std::endl; }\
   BOOST_CHECK_THROW(code, t);

template <class FPT> // Any floating-point type FPT.
void test_error(FPT)
{
	cout << "Current function is " << BOOST_CURRENT_FUNCTION << endl;

   BOOST_CHECK_THROW_MSG(domain_error<FPT>(BOOST_CURRENT_FUNCTION, 0), std::domain_error);
   BOOST_CHECK_THROW_MSG(domain_error<FPT>(BOOST_CURRENT_FUNCTION, "Out of range argument"), std::domain_error);
   BOOST_CHECK_THROW_MSG(domain_error<FPT>(BOOST_CURRENT_FUNCTION, 0, static_cast<FPT>(3.124567890123456789012345678901L)), std::domain_error);
   BOOST_CHECK_THROW_MSG(domain_error<FPT>(BOOST_CURRENT_FUNCTION, "Out of range argument %1% in test invocation", static_cast<FPT>(3.124567890123456789012345678901L)), std::domain_error);
   BOOST_CHECK_THROW_MSG(logic_error<FPT>(BOOST_CURRENT_FUNCTION, 0), std::logic_error);
   BOOST_CHECK_THROW_MSG(logic_error<FPT>(BOOST_CURRENT_FUNCTION, "Internal error"), std::logic_error);
   BOOST_CHECK_THROW_MSG(logic_error<FPT>(BOOST_CURRENT_FUNCTION, 0, static_cast<FPT>(3.124567890123456789012345678901L)), std::logic_error);
   BOOST_CHECK_THROW_MSG(logic_error<FPT>(BOOST_CURRENT_FUNCTION, "Internal error, computed result was %1%, but should be in the range [0,1]", static_cast<FPT>(3.124567890123456789012345678901L)), std::logic_error);
   BOOST_CHECK_THROW_MSG(chisqr(-1, static_cast<FPT>(1)), std::domain_error);
   BOOST_CHECK_THROW_MSG(chisqr(static_cast<FPT>(1), static_cast<FPT>(-1)), std::domain_error);
   BOOST_CHECK_THROW_MSG(chisqr(static_cast<FPT>(0), static_cast<FPT>(1)), std::domain_error);
   BOOST_CHECK_THROW_MSG(tgamma(static_cast<FPT>(0)), std::domain_error);
   BOOST_CHECK_THROW_MSG(tgamma(static_cast<FPT>(-10)), std::domain_error);
   BOOST_CHECK_THROW_MSG(tgamma(static_cast<FPT>(-10123457772243.0)), std::domain_error);

   cout << endl;

} // template <class FPT>void test_error(FPT)

int test_main(int, char* [])
{
	// Test error handling.
	// (Parameter value, arbitrarily zero, only communicates the floating point type FPT).
	test_error(0.0F); // Test float.
	test_error(0.0); // Test double.
	test_error(0.0L); // Test long double.
	test_error(real_concept(0.0L)); // Test concepts.
	return 0;
} // int test_main(int, char* [])

