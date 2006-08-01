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
// This file tests the incomplete beta functions beta, 
// betac, ibeta and ibetac.  There are two sets of tests, spot
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
   const char* largest_type;
#ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
   if(boost::math::tools::digits<double>() == boost::math::tools::digits<long double>())
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
      "[^|]*",                          // compiler
      "[^|]*",                          // stdlib
      "[^|]*",                          // platform
      largest_type,                     // test type(s)
      "(?i).*small.*",                      // test data group
      ".*", 20, 10);  // test function
   add_expected_result(
      "[^|]*",                          // compiler
      "[^|]*",                          // stdlib
      "[^|]*",                          // platform
      largest_type,                     // test type(s)
      "(?i).*medium.*",                     // test data group
      ".*", 150, 50);  // test function
   add_expected_result(
      "[^|]*",                          // compiler
      "[^|]*",                          // stdlib
      "[^|]*",                          // platform
      largest_type,                     // test type(s)
      "(?i).*large.*",                      // test data group
      ".*", 5000, 500);                 // test function

   add_expected_result(
      "[^|]*",                          // compiler
      "[^|]*",                          // stdlib
      "[^|]*",                          // platform
      "real_concept",                   // test type(s)
      "(?i).*small.*",                      // test data group
      ".*", 30, 15);  // test function
   add_expected_result(
      "[^|]*",                          // compiler
      "[^|]*",                          // stdlib
      "[^|]*",                          // platform
      "real_concept",                   // test type(s)
      "(?i).*medium.*",                     // test data group
      ".*", 100, 50);  // test function
   add_expected_result(
      "[^|]*",                          // compiler
      "[^|]*",                          // stdlib
      "[^|]*",                          // platform
      "real_concept",                   // test type(s)
      "(?i).*large.*",                      // test data group
      ".*", 200000, 50000);             // test function
   //
   // Finish off by printing out the compiler/stdlib/platform names,
   // we do this to make it easier to mark up expected error rates.
   //
   std::cout << "Tests run with " << BOOST_COMPILER << ", " 
      << BOOST_STDLIB << ", " << BOOST_PLATFORM << std::endl;
}

template <class T>
void do_test_beta(const T& data, const char* type_name, const char* test_name)
{
   typedef typename T::value_type row_type;
   typedef typename row_type::value_type value_type;

   typedef value_type (*pg)(value_type, value_type, value_type);
   pg funcp = boost::math::beta;

   boost::math::tools::test_result<value_type> result;

   std::cout << "Testing " << test_name << " with type " << type_name
      << "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";

   //
   // test beta against data:
   //
   result = boost::math::tools::test(
      data,
      boost::lambda::bind(funcp,
         boost::lambda::ret<value_type>(boost::lambda::_1[0]),
         boost::lambda::ret<value_type>(boost::lambda::_1[1]),
         boost::lambda::ret<value_type>(boost::lambda::_1[2])),
      boost::lambda::ret<value_type>(boost::lambda::_1[3]));
   handle_test_result(result, data[result.worst()], result.worst(), type_name, "boost::math::beta", test_name);

   funcp = boost::math::betac;
   result = boost::math::tools::test(
      data,
      boost::lambda::bind(funcp,
         boost::lambda::ret<value_type>(boost::lambda::_1[0]),
         boost::lambda::ret<value_type>(boost::lambda::_1[1]),
         boost::lambda::ret<value_type>(boost::lambda::_1[2])),
      boost::lambda::ret<value_type>(boost::lambda::_1[4]));
   handle_test_result(result, data[result.worst()], result.worst(), type_name, "boost::math::betac", test_name);

   funcp = boost::math::ibeta;
   result = boost::math::tools::test(
      data,
      boost::lambda::bind(funcp,
         boost::lambda::ret<value_type>(boost::lambda::_1[0]),
         boost::lambda::ret<value_type>(boost::lambda::_1[1]),
         boost::lambda::ret<value_type>(boost::lambda::_1[2])),
      boost::lambda::ret<value_type>(boost::lambda::_1[5]));
   handle_test_result(result, data[result.worst()], result.worst(), type_name, "boost::math::ibeta", test_name);

   funcp = boost::math::ibetac;
   result = boost::math::tools::test(
      data,
      boost::lambda::bind(funcp,
         boost::lambda::ret<value_type>(boost::lambda::_1[0]),
         boost::lambda::ret<value_type>(boost::lambda::_1[1]),
         boost::lambda::ret<value_type>(boost::lambda::_1[2])),
      boost::lambda::ret<value_type>(boost::lambda::_1[6]));
   handle_test_result(result, data[result.worst()], result.worst(), type_name, "boost::math::ibetac", test_name);
#ifdef TEST_OTHER
   if(::boost::is_floating_point<value_type>::value){
      funcp = other::ibeta;
      result = boost::math::tools::test(
         data,
         boost::lambda::bind(funcp,
            boost::lambda::ret<value_type>(boost::lambda::_1[0]),
            boost::lambda::ret<value_type>(boost::lambda::_1[1]),
            boost::lambda::ret<value_type>(boost::lambda::_1[2])),
         boost::lambda::ret<value_type>(boost::lambda::_1[5]));
      print_test_result(result, data[result.worst()], result.worst(), type_name, "other::ibeta");
   }
#endif
   std::cout << std::endl;
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

   do_test_beta(ibeta_small_data, name, "Incomplete Beta Function: Small Values");

#  include "ibeta_data.ipp"

   do_test_beta(ibeta_data, name, "Incomplete Beta Function: Medium Values");

#  include "ibeta_large_data.ipp"

   do_test_beta(ibeta_large_data, name, "Incomplete Beta Function: Large and Diverse Values");
}

template <class T>
void test_spots(T)
{
   //
   // basic sanity checks, tolerance is 30 epsilon expressed as a percentage:
   //
   T tolerance = boost::math::tools::epsilon<T>() * 3000;
   BOOST_CHECK_CLOSE(
      ::boost::math::ibeta(
         static_cast<T>(0.015964560210704803),
         static_cast<T>(1.1846856068586931e-005),
         static_cast<T>(0.69176378846168518)),
      static_cast<T>(0.0007508604820642986204162462167319506309750L), tolerance);
   BOOST_CHECK_CLOSE(
      ::boost::math::ibeta(
         static_cast<T>(42.434902191162109),
         static_cast<T>(0.30012050271034241),
         static_cast<T>(0.91574394702911377)),
      static_cast<T>(0.002844243156314242058287766324242276991912L), tolerance);
   BOOST_CHECK_CLOSE(
      ::boost::math::ibeta(
         static_cast<T>(9.7131776809692383),
         static_cast<T>(99.406852722167969),
         static_cast<T>(0.083912998437881470)),
      static_cast<T>(0.4612716118626361034813232775095335302198L), tolerance);
   BOOST_CHECK_CLOSE(
      ::boost::math::ibeta(
         static_cast<T>(72.5),
         static_cast<T>(1.125),
         static_cast<T>(0.75)),
      static_cast<T>(1.3423066982487051710597194786268004978931316494920e-9L), tolerance);
   BOOST_CHECK_CLOSE(
      ::boost::math::ibeta(
         static_cast<T>(4.9854421615600586),
         static_cast<T>(1.0665277242660522),
         static_cast<T>(0.75997146964073181)),
      static_cast<T>(0.2755954254731642667260071858810487404614L), tolerance);
   BOOST_CHECK_CLOSE(
      ::boost::math::ibeta(
         static_cast<T>(6.8127136230468750),
         static_cast<T>(1.0562920570373535),
         static_cast<T>(0.17416560649871826)),
      static_cast<T>(7.702362015088558153029455563361002570531e-6L), tolerance);
   BOOST_CHECK_CLOSE(
      ::boost::math::ibeta(
         static_cast<T>(0.48983201384544373),
         static_cast<T>(0.22512593865394592),
         static_cast<T>(0.20032680034637451)),
      static_cast<T>(0.170905142698145967653807992508983970176L), tolerance);
   BOOST_CHECK_CLOSE(
      ::boost::math::ibeta(
         static_cast<T>(4.0498137474060059),
         static_cast<T>(0.15403440594673157),
         static_cast<T>(0.65370121598243713)),
      static_cast<T>(0.0172702040689452906446803217247250156007L), tolerance);
   BOOST_CHECK_CLOSE(
      ::boost::math::ibeta(
         static_cast<T>(7.2695474624633789),
         static_cast<T>(0.11902070045471191),
         static_cast<T>(0.80036874115467072)),
      static_cast<T>(0.013346136714187857821168127038816508028L), tolerance);
   BOOST_CHECK_CLOSE(
      ::boost::math::ibeta(
         static_cast<T>(2.7266697883605957),
         static_cast<T>(0.011510574258863926),
         static_cast<T>(0.086654007434844971)),
      static_cast<T>(5.812020420972734916187451486321162137375e-6L), tolerance);
   BOOST_CHECK_CLOSE(
      ::boost::math::ibeta(
         static_cast<T>(0.34317314624786377),
         static_cast<T>(0.046342257410287857),
         static_cast<T>(0.75823287665843964)),
      static_cast<T>(0.151317265120184850887504097401768195067L), tolerance);

   BOOST_CHECK_CLOSE(
      ::boost::math::ibeta(
         static_cast<T>(0.34317314624786377),
         static_cast<T>(0.046342257410287857),
         static_cast<T>(0)),
      static_cast<T>(0), tolerance);
   BOOST_CHECK_CLOSE(
      ::boost::math::ibetac(
         static_cast<T>(0.34317314624786377),
         static_cast<T>(0.046342257410287857),
         static_cast<T>(0)),
      static_cast<T>(1), tolerance);
   BOOST_CHECK_CLOSE(
      ::boost::math::ibeta(
         static_cast<T>(0.34317314624786377),
         static_cast<T>(0.046342257410287857),
         static_cast<T>(1)),
      static_cast<T>(1), tolerance);
   BOOST_CHECK_CLOSE(
      ::boost::math::ibetac(
         static_cast<T>(0.34317314624786377),
         static_cast<T>(0.046342257410287857),
         static_cast<T>(1)),
      static_cast<T>(0), tolerance);
   BOOST_CHECK_CLOSE(
      ::boost::math::ibeta(
         static_cast<T>(1),
         static_cast<T>(0.046342257410287857),
         static_cast<T>(0.32)),
      static_cast<T>(0.0177137046180187568703202426065033413304L), tolerance);
   BOOST_CHECK_CLOSE(
      ::boost::math::ibeta(
         static_cast<T>(0.046342257410287857),
         static_cast<T>(1),
         static_cast<T>(0.32)),
      static_cast<T>(0.948565954109602496577407403168592262389L), tolerance);

   // very naive check on derivative:
   using namespace std;  // For ADL of std functions
   tolerance = boost::math::tools::epsilon<T>() * 10000; // 100 eps
   BOOST_CHECK_CLOSE(
      ::boost::math::ibeta_derivative(
         static_cast<T>(2),
         static_cast<T>(3),
         static_cast<T>(0.5)),
         pow(static_cast<T>(0.5), static_cast<T>(2)) * pow(static_cast<T>(0.5), static_cast<T>(1)) / boost::math::beta(static_cast<T>(2), static_cast<T>(3)), tolerance);
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




