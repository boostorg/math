//  (C) Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/math/concepts/real_concept.hpp>
#include <boost/test/included/test_exec_monitor.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/math/special_functions/beta.hpp>
#include <boost/math/tools/stats.hpp>
#include <boost/math/tools/test.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/type_traits/is_floating_point.hpp>
#include <boost/array.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>

#include "test_beta_hooks.hpp"
#include "handle_test_result.hpp"

//
// DESCRIPTION:
// ~~~~~~~~~~~~
//
// This file tests the incomplete beta function inverses 
// ibeta_inv and ibetac_inv. There are three sets of tests:
// 1) Spot tests which compare our results with selected values 
// computed using the online special function calculator at 
// functions.wolfram.com, 
// 2) TODO!!!! Accuracy tests use values generated with NTL::RR at 
// 1000-bit precision and our generic versions of these functions.
// 3) Round trip sanity checks, use the test data for the forward
// functions, and verify that we can get (approximately) back
// where we started.
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

   //
   // Catch all cases come last:
   //
   // TODO!!!!
   add_expected_result(
      ".*",                          // compiler
      ".*",                          // stdlib
      ".*",                          // platform
      ".*",                          // test type(s)
      "Inverse Erf.*",               // test data group
      "boost::math::erfc?_inv", 14, 4);  // test function

   //
   // Finish off by printing out the compiler/stdlib/platform names,
   // we do this to make it easier to mark up expected error rates.
   //
   std::cout << "Tests run with " << BOOST_COMPILER << ", " 
      << BOOST_STDLIB << ", " << BOOST_PLATFORM << std::endl;
}

template <class T>
void test_inverses(const T& data)
{
   using namespace std;
   typedef typename T::value_type row_type;
   typedef typename row_type::value_type value_type;

   value_type precision = static_cast<value_type>(ldexp(1.0, 1-boost::math::tools::digits<value_type>()/2)) * 100;
   if(boost::math::tools::digits<value_type>() < 50)
      precision = 1;   // 1% or two decimal digits, all we can hope for when the input is truncated

   for(unsigned i = 0; i < data.size(); ++i)
   {
      //
      // These inverse tests are thrown off if the output of the
      // incomplete beta is too close to 1: basically there is insuffient
      // information left in the value we're using as input to the inverse
      // to be able to get back to the original value.
      //
      if(data[i][5] == 0)
         BOOST_CHECK_EQUAL(boost::math::ibeta_inv(data[i][0], data[i][1], data[i][5]), value_type(0));
      else if((1 - data[i][5] > 0.001) && (fabs(data[i][5]) >= boost::math::tools::min_value<value_type>()))
      {
         value_type inv = boost::math::ibeta_inv(data[i][0], data[i][1], data[i][5]);
         BOOST_CHECK_CLOSE(data[i][2], inv, precision);
      }
      else if(1 == data[i][5])
         BOOST_CHECK_EQUAL(boost::math::ibeta_inv(data[i][0], data[i][1], data[i][5]), value_type(1));

      if(data[i][6] == 0)
         BOOST_CHECK_EQUAL(boost::math::ibetac_inv(data[i][0], data[i][1], data[i][6]), value_type(1));
      else if((1 - data[i][6] > 0.001) && (fabs(data[i][6]) >= boost::math::tools::min_value<value_type>()))
      {
         value_type inv = boost::math::ibetac_inv(data[i][0], data[i][1], data[i][6]);
         BOOST_CHECK_CLOSE(data[i][2], inv, precision);
      }
      else if(data[i][6] == 1)
         BOOST_CHECK_EQUAL(boost::math::ibetac_inv(data[i][0], data[i][1], data[i][6]), value_type(0));
   }
}

template <class T>
void test_beta(T, const char* name)
{
   //
   // The actual test data is rather verbose, so it's in a separate file
   //
   // The contents are as follows, each row of data contains
   // five items, input value a, input value b, integration limits x, beta(a, b, x) and ibeta(a, b, x):
   //
#  include "ibeta_small_data.ipp"

   test_inverses(ibeta_small_data);

#  include "ibeta_data.ipp"

   test_inverses(ibeta_data);

#  include "ibeta_large_data.ipp"

   test_inverses(ibeta_large_data);
}

template <class T>
void test_spots(T)
{
   //
   // basic sanity checks, tolerance is 100 epsilon expressed as a percentage:
   //
   T tolerance = boost::math::tools::epsilon<T>() * 10000;
   BOOST_CHECK_CLOSE(
      ::boost::math::ibeta_inv(
         static_cast<T>(1),
         static_cast<T>(2),
         static_cast<T>(0.5)),
      static_cast<T>(0.29289321881345247559915563789515096071516406231153L), tolerance);
   BOOST_CHECK_CLOSE(
      ::boost::math::ibeta_inv(
         static_cast<T>(3),
         static_cast<T>(0.5),
         static_cast<T>(0.5)),
      static_cast<T>(0.92096723292382700385142816696980724853063433975470L), tolerance);
   BOOST_CHECK_CLOSE(
      ::boost::math::ibeta_inv(
         static_cast<T>(20.125),
         static_cast<T>(0.5),
         static_cast<T>(0.5)),
      static_cast<T>(0.98862133312917003480022776106012775747685870929920L), tolerance);
   BOOST_CHECK_CLOSE(
      ::boost::math::ibeta_inv(
         static_cast<T>(40),
         static_cast<T>(80),
         static_cast<T>(0.5)),
      static_cast<T>(0.33240456430025026300937492802591128972548660643778L), tolerance);
}

int test_main(int, char* [])
{
   expected_results();
#ifdef TEST_GSL
   gsl_set_error_handler_off();
#endif
   test_spots(0.0F);
   test_spots(0.0);
#ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
   test_spots(0.0L);
   test_spots(boost::math::concepts::real_concept(0.1));
#endif

   test_beta(0.1F, "float");
   test_beta(0.1, "double");
#ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
   test_beta(0.1L, "long double");
   test_beta(boost::math::concepts::real_concept(0.1), "real_concept");
#else
   std::cout << "<note>The long double tests have been disabled on this platform "
      "either because the long double overloads of the usual math functions are "
      "not available at all, or because they are too inaccurate for these tests "
      "to pass.</note>" << std::cout;
#endif
   return 0;
}



