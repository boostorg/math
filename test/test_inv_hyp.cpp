//  (C) Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/math/concepts/real_concept.hpp>
#include <boost/math/special_functions/acosh.hpp>
#include <boost/math/special_functions/asinh.hpp>
#include <boost/math/special_functions/atanh.hpp>
#include <boost/test/included/test_exec_monitor.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/math/tools/stats.hpp>
#include <boost/math/tools/test.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/type_traits/is_floating_point.hpp>
#include <boost/array.hpp>
#include "functor.hpp"

#include "handle_test_result.hpp"

//
// DESCRIPTION:
// ~~~~~~~~~~~~
//
// This file tests the inverse hyperbolic functions. There are two sets of tests:
// 1) Sanity checks: comparison to test values created with the
// online calculator at functions.wolfram.com
// 2) Accuracy tests use values generated with NTL::RR at 
// 1000-bit precision and our generic versions of these functions.
//
// Note that when this file is first run on a new platform many of
// these tests will fail: the default accuracy is 1 epsilon which
// is too tight for most platforms.  In this situation you will 
// need to cast a human eye over the error rates reported and make
// a judgement as to whether they are acceptable.  Either way please
// report the results to the Boost mailing list.  Acceptable rates of
// error are marked up below as a series of regular expressions that
// identify the compiler/stdlib/platform/data-type/test-data/test-function
// along with the maximum expected peek and RMS mean errors for that
// test.
//

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

   add_expected_result(
      ".*",                          // compiler
      ".*",                          // stdlib
      ".*",                          // platform
      "real_concept",                // test type(s)
      ".*",                          // test data group
      ".*", 4, 2);                   // test function
   add_expected_result(
      ".*",                          // compiler
      ".*",                          // stdlib
      ".*",                          // platform
      ".*",                          // test type(s)
      ".*",                          // test data group
      ".*", 2, 1);                   // test function

   std::cout << "Tests run with " << BOOST_COMPILER << ", " 
      << BOOST_STDLIB << ", " << BOOST_PLATFORM << std::endl;
}

template <class T>
void do_test_asinh(const T& data, const char* type_name, const char* test_name)
{
   //
   // test asinh(T) against data:
   //
   using namespace std;
   typedef typename T::value_type row_type;
   typedef typename row_type::value_type value_type;

   std::cout << test_name << " with type " << type_name << std::endl;

   typedef value_type (*pg)(value_type);
#if defined(BOOST_MATH_NO_DEDUCED_FUNCTION_POINTERS)
   pg funcp = boost::math::asinh<value_type>;
#else
   pg funcp = boost::math::asinh;
#endif

   boost::math::tools::test_result<value_type> result;
   //
   // test asinh against data:
   //
   result = boost::math::tools::test(
      data,
      bind_func(funcp, 0),
      extract_result(1));
   handle_test_result(result, data[result.worst()], result.worst(), type_name, "boost::math::asinh", test_name);
   std::cout << std::endl;
}

template <class T>
void do_test_acosh(const T& data, const char* type_name, const char* test_name)
{
   //
   // test acosh(T) against data:
   //
   using namespace std;
   typedef typename T::value_type row_type;
   typedef typename row_type::value_type value_type;

   std::cout << test_name << " with type " << type_name << std::endl;

   typedef value_type (*pg)(value_type);
#if defined(BOOST_MATH_NO_DEDUCED_FUNCTION_POINTERS)
   pg funcp = boost::math::acosh<value_type>;
#else
   pg funcp = boost::math::acosh;
#endif

   boost::math::tools::test_result<value_type> result;
   //
   // test acosh against data:
   //
   result = boost::math::tools::test(
      data,
      bind_func(funcp, 0),
      extract_result(1));
   handle_test_result(result, data[result.worst()], result.worst(), type_name, "boost::math::acosh", test_name);
   std::cout << std::endl;
}

template <class T>
void do_test_atanh(const T& data, const char* type_name, const char* test_name)
{
   //
   // test atanh(T) against data:
   //
   using namespace std;
   typedef typename T::value_type row_type;
   typedef typename row_type::value_type value_type;

   std::cout << test_name << " with type " << type_name << std::endl;

   typedef value_type (*pg)(value_type);
#if defined(BOOST_MATH_NO_DEDUCED_FUNCTION_POINTERS)
   pg funcp = boost::math::atanh<value_type>;
#else
   pg funcp = boost::math::atanh;
#endif

   boost::math::tools::test_result<value_type> result;
   //
   // test atanh against data:
   //
   result = boost::math::tools::test(
      data,
      bind_func(funcp, 0),
      extract_result(1));
   handle_test_result(result, data[result.worst()], result.worst(), type_name, "boost::math::atanh", test_name);
   std::cout << std::endl;
}

template <class T>
void test_inv_hyperbolics(T, const char* name)
{
   //
   // The actual test data is rather verbose, so it's in a separate file
   //
#include "asinh_data.ipp"
   do_test_asinh(asinh_data, name, "asinh");
#include "acosh_data.ipp"
   do_test_acosh(acosh_data, name, "acosh");
#include "atanh_data.ipp"
   do_test_atanh(atanh_data, name, "atanh");
}

extern "C" double zetac(double);

template <class T>
void test_spots(T, const char* t)
{
   std::cout << "Testing basic sanity checks for type " << t << std::endl;
   //
   // Basic sanity checks, tolerance is either 5 or 10 epsilon 
   // expressed as a percentage:
   //
   T tolerance = boost::math::tools::epsilon<T>() * 100 *
      (boost::is_floating_point<T>::value ? 5 : 10);
   //BOOST_CHECK_CLOSE(::boost::math::zeta(static_cast<T>(0.125)), static_cast<T>(-0.63277562349869525529352526763564627152686379131122L), tolerance);
   (void)tolerance;
}

int test_main(int, char* [])
{
   expected_results();
   BOOST_MATH_CONTROL_FP;

   test_spots(0.0f, "float");
   test_spots(0.0, "double");
#ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
   test_spots(0.0L, "long double");
   test_spots(boost::math::concepts::real_concept(0.1), "real_concept");
#else
   std::cout << "<note>The long double tests have been disabled on this platform "
      "either because the long double overloads of the usual math functions are "
      "not available at all, or because they are too inaccurate for these tests "
      "to pass.</note>" << std::cout;
#endif

   test_inv_hyperbolics(0.1F, "float");
   test_inv_hyperbolics(0.1, "double");
#ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
   test_inv_hyperbolics(0.1L, "long double");
   test_inv_hyperbolics(boost::math::concepts::real_concept(0.1), "real_concept");
#else
   std::cout << "<note>The long double tests have been disabled on this platform "
      "either because the long double overloads of the usual math functions are "
      "not available at all, or because they are too inaccurate for these tests "
      "to pass.</note>" << std::cout;
#endif
   return 0;
}



