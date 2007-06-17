
// Copyright Paul A. Bristow 2006.
// Copyright John Maddock 2006.

// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

//
// Define some custom error handlers:
//
struct user_defined_error{};

namespace boost{ namespace math{ namespace policy{

template <class T>
T user_domain_error(const char* , const char* , const T& )
{
   throw user_defined_error();
}

template <class T>
T user_pole_error(const char* , const char* , const T& )
{
   throw user_defined_error();
}

template <class T>
T user_overflow_error(const char* , const char* , const T& )
{
   throw user_defined_error();
}

template <class T>
T user_underflow_error(const char* , const char* , const T& )
{
   throw user_defined_error();
}

template <class T>
T user_denorm_error(const char* , const char* , const T& )
{
   throw user_defined_error();
}

template <class T>
T user_evaluation_error(const char* , const char* , const T& )
{
   throw user_defined_error();
}

}}} // namespaces


#include <boost/math/policy/policy.hpp>
#include <boost/math/policy/error_handling.hpp>
#include <boost/math/concepts/real_concept.hpp>
#include <boost/test/included/test_exec_monitor.hpp> // for test_main
//
// Define some policies:
//
using namespace boost::math::policy;
policy<
   domain_error<throw_on_error>,
   pole_error<throw_on_error>,
   overflow_error<throw_on_error>,
   underflow_error<throw_on_error>,
   denorm_error<throw_on_error>,
   evaluation_error<throw_on_error> > throw_policy;
policy<
   domain_error<errno_on_error>,
   pole_error<errno_on_error>,
   overflow_error<errno_on_error>,
   underflow_error<errno_on_error>,
   denorm_error<errno_on_error>,
   evaluation_error<errno_on_error> > errno_policy;
policy<
   domain_error<user_error>,
   pole_error<user_error>,
   overflow_error<user_error>,
   underflow_error<user_error>,
   denorm_error<user_error>,
   evaluation_error<user_error> > user_policy;
policy<> default_policy;

#define TEST_EXCEPTION(expression, exception)\
   BOOST_CHECK_THROW(expression, exception);\
   try{ expression; }catch(const exception& e){ std::cout << e.what() << std::endl; }

template <class T>
void test_error(T)
{
   const char* func = "boost::math::test_function<%1%>(%1%, %1%, %1%)";
   const char* msg1 = "Error while handling value %1%";
   const char* msg2 = "Error message goes here...";

   TEST_EXCEPTION(boost::math::policy::raise_domain_error(func, msg1, T(0.0), throw_policy), std::domain_error);
   TEST_EXCEPTION(boost::math::policy::raise_domain_error(func, 0, T(0.0), throw_policy), std::domain_error);
   TEST_EXCEPTION(boost::math::policy::raise_pole_error(func, msg1, T(0.0), throw_policy), std::domain_error);
   TEST_EXCEPTION(boost::math::policy::raise_pole_error(func, 0, T(0.0), throw_policy), std::domain_error);
   TEST_EXCEPTION(boost::math::policy::raise_overflow_error<T>(func, msg2, throw_policy), std::overflow_error);
   TEST_EXCEPTION(boost::math::policy::raise_overflow_error<T>(func, 0, throw_policy), std::overflow_error);
   TEST_EXCEPTION(boost::math::policy::raise_underflow_error<T>(func, msg2, throw_policy), std::underflow_error);
   TEST_EXCEPTION(boost::math::policy::raise_underflow_error<T>(func, 0, throw_policy), std::underflow_error);
   TEST_EXCEPTION(boost::math::policy::raise_denorm_error<T>(func, msg2, T(0), throw_policy), std::underflow_error);
   TEST_EXCEPTION(boost::math::policy::raise_denorm_error<T>(func, 0, T(0), throw_policy), std::underflow_error);
   TEST_EXCEPTION(boost::math::policy::raise_evaluation_error(func, msg1, T(1.25), throw_policy), boost::math::evaluation_error);
   TEST_EXCEPTION(boost::math::policy::raise_evaluation_error(func, 0, T(1.25), throw_policy), boost::math::evaluation_error);
   //
   // Now try user error handlers: these should all throw user_error():
   //
   BOOST_CHECK_THROW(boost::math::policy::raise_domain_error(func, msg1, T(0.0), user_policy), user_defined_error);
   BOOST_CHECK_THROW(boost::math::policy::raise_pole_error(func, msg1, T(0.0), user_policy), user_defined_error);
   BOOST_CHECK_THROW(boost::math::policy::raise_overflow_error<T>(func, msg2, user_policy), user_defined_error);
   BOOST_CHECK_THROW(boost::math::policy::raise_underflow_error<T>(func, msg2, user_policy), user_defined_error);
   BOOST_CHECK_THROW(boost::math::policy::raise_denorm_error<T>(func, msg2, T(0), user_policy), user_defined_error);
   BOOST_CHECK_THROW(boost::math::policy::raise_evaluation_error(func, msg1, T(0.0), user_policy), user_defined_error);

}

int test_main(int, char* [])
{
	// Test error handling.
	// (Parameter value, arbitrarily zero, only communicates the floating point type FPT).
	test_error(0.0F); // Test float.
	test_error(0.0); // Test double.
	test_error(0.0L); // Test long double.
   test_error(boost::math::concepts::real_concept(0.0L)); // Test concepts.
	return 0;
} // int test_main(int, char* [])


