//  (C) Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/math/concepts/real_concept.hpp>
#include <boost/test/included/test_exec_monitor.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/math/special_functions/erf.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/type_traits/is_floating_point.hpp>
#include <boost/array.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>

#include "test_erf_hooks.hpp"
#include "handle_test_result.hpp"

//
// DESCRIPTION:
// ~~~~~~~~~~~~
//
// This file tests the functions erf, erfc, and the inverses
// erf_inv and erfc_inv.  There are two sets of tests, spot
// tests which compare our results with selected values computed
// using the online special function calculator at 
// functions.wolfram.com, while the bulk of the accuracy tests
// use values generated with NTL::RR at 1000-bit precision
// and our generic versions of these functions.
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
   add_expected_result(
      ".*",                          // compiler
      ".*",                          // stdlib
      ".*",                          // platform
      "real_concept",                // test type(s)
      "Erf Function:.*",             // test data group
      "boost::math::erfc?", 20, 6);   // test function
   add_expected_result(
      ".*",                           // compiler
      ".*",                           // stdlib
      ".*",                           // platform
      "real_concept",                 // test type(s)
      "Inverse Erfc.*",               // test data group
      "boost::math::erfc_inv", 80, 10);  // test function


   //
   // Catch all cases come last:
   //
   add_expected_result(
      ".*",                          // compiler
      ".*",                          // stdlib
      ".*",                          // platform
      ".*",                          // test type(s)
      "Erf Function:.*",             // test data group
      "boost::math::erfc?", 2, 2);   // test function
   add_expected_result(
      ".*aCC.*",                     // compiler
      ".*",                          // stdlib
      ".*",                          // platform
      ".*",                          // test type(s)
      "Inverse Erfc.*",               // test data group
      "boost::math::erfc_inv", 80, 10);  // test function
   add_expected_result(
      ".*",                          // compiler
      ".*",                          // stdlib
      ".*",                          // platform
      ".*",                          // test type(s)
      "Inverse Erf.*",               // test data group
      "boost::math::erfc?_inv", 18, 4);  // test function

   //
   // Finish off by printing out the compiler/stdlib/platform names,
   // we do this to make it easier to mark up expected error rates.
   //
   std::cout << "Tests run with " << BOOST_COMPILER << ", " 
      << BOOST_STDLIB << ", " << BOOST_PLATFORM << std::endl;
}

template <class T>
void do_test_erf(const T& data, const char* type_name, const char* test_name)
{
   typedef typename T::value_type row_type;
   typedef typename row_type::value_type value_type;

   typedef value_type (*pg)(value_type);
   pg funcp = boost::math::erf;

   boost::math::tools::test_result<value_type> result;

   std::cout << "Testing " << test_name << " with type " << type_name
      << "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";

   //
   // test erf against data:
   //
   result = boost::math::tools::test(
      data, 
      boost::lambda::bind(funcp, 
         boost::lambda::ret<value_type>(boost::lambda::_1[0])), 
      boost::lambda::ret<value_type>(boost::lambda::_1[1]));
   handle_test_result(result, data[result.worst()], result.worst(), type_name, "boost::math::erf", test_name);
#ifdef TEST_OTHER
   if(::boost::is_floating_point<value_type>::value){
      funcp = other::erf;
      result = boost::math::tools::test(
         data, 
         boost::lambda::bind(funcp, 
            boost::lambda::ret<value_type>(boost::lambda::_1[0])), 
         boost::lambda::ret<value_type>(boost::lambda::_1[1]));
      print_test_result(result, data[result.worst()], result.worst(), type_name, "other::erf");
   }
#endif
   //
   // test erfc against data:
   //
   funcp = boost::math::erfc;
   result = boost::math::tools::test(
      data, 
      boost::lambda::bind(funcp, 
         boost::lambda::ret<value_type>(boost::lambda::_1[0])), 
      boost::lambda::ret<value_type>(boost::lambda::_1[2]));
   handle_test_result(result, data[result.worst()], result.worst(), type_name, "boost::math::erfc", test_name);
#ifdef TEST_OTHER
   if(::boost::is_floating_point<value_type>::value){
      funcp = other::erfc;
      result = boost::math::tools::test(
         data, 
         boost::lambda::bind(funcp, 
            boost::lambda::ret<value_type>(boost::lambda::_1[0])), 
         boost::lambda::ret<value_type>(boost::lambda::_1[2]));
      print_test_result(result, data[result.worst()], result.worst(), type_name, "other::erfc");
   }
#endif
   std::cout << std::endl;
}

template <class T>
void do_test_erf_inv(const T& data, const char* type_name, const char* test_name)
{
   typedef typename T::value_type row_type;
   typedef typename row_type::value_type value_type;

   typedef value_type (*pg)(value_type);
   pg funcp = boost::math::erf;

   boost::math::tools::test_result<value_type> result;
   std::cout << "Testing " << test_name << " with type " << type_name
      << "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
   //
   // test erf_inv against data:
   //
   funcp = boost::math::erf_inv;
   result = boost::math::tools::test(
      data, 
      boost::lambda::bind(funcp, 
         boost::lambda::ret<value_type>(boost::lambda::_1[0])), 
      boost::lambda::ret<value_type>(boost::lambda::_1[1]));
   handle_test_result(result, data[result.worst()], result.worst(), type_name, "boost::math::erf_inv", test_name);
   std::cout << std::endl;
}

template <class T>
void do_test_erfc_inv(const T& data, const char* type_name, const char* test_name)
{
   typedef typename T::value_type row_type;
   typedef typename row_type::value_type value_type;

   typedef value_type (*pg)(value_type);
   pg funcp = boost::math::erf;

   boost::math::tools::test_result<value_type> result;
   std::cout << "Testing " << test_name << " with type " << type_name
      << "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
   //
   // test erfc_inv against data:
   //
   funcp = boost::math::erfc_inv;
   result = boost::math::tools::test(
      data, 
      boost::lambda::bind(funcp, 
         boost::lambda::ret<value_type>(boost::lambda::_1[0])), 
      boost::lambda::ret<value_type>(boost::lambda::_1[1]));
   handle_test_result(result, data[result.worst()], result.worst(), type_name, "boost::math::erfc_inv", test_name);
   std::cout << std::endl;
}

template <class T>
void test_erf(T, const char* name)
{
   //
   // The actual test data is rather verbose, so it's in a separate file
   //
   // The contents are as follows, each row of data contains
   // three items, input value a, input value b and erf(a, b):
   // 
#  include "erf_small_data.ipp"

   do_test_erf(erf_small_data, name, "Erf Function: Small Values");

#  include "erf_data.ipp"

   do_test_erf(erf_data, name, "Erf Function: Medium Values");

#  include "erf_large_data.ipp"

   do_test_erf(erf_large_data, name, "Erf Function: Large Values");

#  include "erf_inv_data.ipp"

   do_test_erf_inv(erf_inv_data, name, "Inverse Erf Function");

#  include "erfc_inv_data.ipp"

   do_test_erfc_inv(erfc_inv_data, name, "Inverse Erfc Function");

#  include "erfc_inv_big_data.ipp"

   if(std::numeric_limits<T>::min_exponent <= -4500)
   {
      do_test_erfc_inv(erfc_inv_big_data, name, "Inverse Erfc Function: extreme values");
   }
}

template <class T>
void test_spots(T, const char* t)
{
   std::cout << "Testing basic sanity checks for type " << t << std::endl;
   //
   // basic sanity checks, tolerance is 10 epsilon expressed as a percentage:
   //
   T tolerance = boost::math::tools::epsilon<T>() * 1000;
   BOOST_CHECK_CLOSE(::boost::math::erfc(static_cast<T>(0.125)), static_cast<T>(0.85968379519866618260697055347837660181302041685015L), tolerance);
   BOOST_CHECK_CLOSE(::boost::math::erfc(static_cast<T>(0.5)), static_cast<T>(0.47950012218695346231725334610803547126354842424204L), tolerance);
   BOOST_CHECK_CLOSE(::boost::math::erfc(static_cast<T>(1)), static_cast<T>(0.15729920705028513065877936491739074070393300203370L), tolerance);
   BOOST_CHECK_CLOSE(::boost::math::erfc(static_cast<T>(5)), static_cast<T>(1.5374597944280348501883434853833788901180503147234e-12L), tolerance);
   BOOST_CHECK_CLOSE(::boost::math::erfc(static_cast<T>(-0.125)), static_cast<T>(1.1403162048013338173930294465216233981869795831498L), tolerance);
   BOOST_CHECK_CLOSE(::boost::math::erfc(static_cast<T>(-0.5)), static_cast<T>(1.5204998778130465376827466538919645287364515757580L), tolerance);
   BOOST_CHECK_CLOSE(::boost::math::erfc(static_cast<T>(0)), static_cast<T>(1), tolerance);

   BOOST_CHECK_CLOSE(::boost::math::erf(static_cast<T>(0.125)), static_cast<T>(0.14031620480133381739302944652162339818697958314985L), tolerance);
   BOOST_CHECK_CLOSE(::boost::math::erf(static_cast<T>(0.5)), static_cast<T>(0.52049987781304653768274665389196452873645157575796L), tolerance);
   BOOST_CHECK_CLOSE(::boost::math::erf(static_cast<T>(1)), static_cast<T>(0.84270079294971486934122063508260925929606699796630L), tolerance);
   BOOST_CHECK_CLOSE(::boost::math::erf(static_cast<T>(5)), static_cast<T>(0.9999999999984625402055719651498116565146166211099L), tolerance);
   BOOST_CHECK_CLOSE(::boost::math::erf(static_cast<T>(-0.125)), static_cast<T>(-0.14031620480133381739302944652162339818697958314985L), tolerance);
   BOOST_CHECK_CLOSE(::boost::math::erf(static_cast<T>(-0.5)), static_cast<T>(-0.52049987781304653768274665389196452873645157575796L), tolerance);
   BOOST_CHECK_CLOSE(::boost::math::erf(static_cast<T>(0)), static_cast<T>(0), tolerance);

   tolerance = boost::math::tools::epsilon<T>() * 100 * 200; // 200 eps %.
#if defined(__CYGWIN__)
   // some platforms long double is only reliably accurate to double precision:
   if(sizeof(T) == sizeof(long double))
      tolerance = boost::math::tools::epsilon<double>() * 100 * 200; // 200 eps %.
#endif
  
   for(T i = -0.95f; i < 1; i += 0.125f)
   {
      T inv = boost::math::erf_inv(i);
      T b = boost::math::erf(inv);
      BOOST_CHECK_CLOSE(b, i, tolerance);
   }
   for(T j = 0.125f; j < 2; j += 0.125f)
   {
      T inv = boost::math::erfc_inv(j);
      T b = boost::math::erfc(inv);
      BOOST_CHECK_CLOSE(b, j, tolerance);
   }
}

int test_main(int, char* [])
{
   test_spots(0.0F, "float");
   test_spots(0.0, "double");
#ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
   test_spots(0.0L, "long double");
   test_spots(boost::math::concepts::real_concept(0.1), "real_concept");
#endif

   expected_results();

   test_erf(0.1F, "float");
   test_erf(0.1, "double");
#ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
   test_erf(0.1L, "long double");
   test_erf(boost::math::concepts::real_concept(0.1), "real_concept");
#else
   std::cout << "<note>The long double tests have been disabled on this platform "
      "either because the long double overloads of the usual math functions are "
      "not available at all, or because they are too inaccurate for these tests "
      "to pass.</note>" << std::cout;
#endif
   return 0;
}


