// Copyright John Maddock 2006.
// Copyright Paul A. Bristow 2007, 2009
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/math/tools/config.hpp>
#ifndef BOOST_MATH_NO_MP_TESTS

#define BOOST_MATH_OVERFLOW_ERROR_POLICY ignore_error

#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/math/tools/stats.hpp>
#include <boost/math/tools/test.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <array>
#include "functor.hpp"

#include "handle_test_result.hpp"
#include "table_type.hpp"

#ifndef SC_
#define SC_(x) static_cast<typename table_type<T>::type>(BOOST_STRINGIZE(x))
#endif

template <class Real, class T>
void do_test_gamma(const T& data, const char* type_name, const char* test_name)
{
#if !(defined(ERROR_REPORTING_MODE) && (!defined(TGAMMA_FUNCTION_TO_TEST) || !defined(LGAMMA_FUNCTION_TO_TEST)))
   typedef Real                   value_type;

   typedef value_type (*pg)(value_type);
#ifdef TGAMMA_FUNCTION_TO_TEST
   pg funcp = TGAMMA_FUNCTION_TO_TEST;
#elif defined(BOOST_MATH_NO_DEDUCED_FUNCTION_POINTERS)
   pg funcp = boost::math::tgamma<value_type>;
#else
   pg funcp = boost::math::tgamma;
#endif

   boost::math::tools::test_result<value_type> result;

   std::cout << "Testing " << test_name << " with type " << type_name
      << "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";

   //
   // test tgamma against data:
   //
   result = boost::math::tools::test_hetero<Real>(
      data,
      bind_func<Real>(funcp, 0),
      extract_result<Real>(1));
   handle_test_result(result, data[result.worst()], result.worst(), type_name, "tgamma", test_name);
   //
   // test lgamma against data:
   //
#ifdef LGAMMA_FUNCTION_TO_TEST
   funcp = LGAMMA_FUNCTION_TO_TEST;
#elif defined(BOOST_MATH_NO_DEDUCED_FUNCTION_POINTERS)
   funcp = boost::math::lgamma<value_type>;
#else
   funcp = boost::math::lgamma;
#endif
   result = boost::math::tools::test_hetero<Real>(
      data,
      bind_func<Real>(funcp, 0),
      extract_result<Real>(2));
   handle_test_result(result, data[result.worst()], result.worst(), type_name, "lgamma", test_name);

   std::cout << std::endl;
#endif
}

template <class T>
void test_gamma(T, const char* name)
{
   //
   // The actual test data is rather verbose, so it's in a separate file
   //
   // The contents are as follows, each row of data contains
   // three items, input value, gamma and lgamma:
   //
   // gamma and lgamma at integer and half integer values:
   // std::array<std::array<T, 3>, N> factorials;
   //
   // gamma and lgamma for z near 0:
   // std::array<std::array<T, 3>, N> near_0;
   //
   // gamma and lgamma for z near 1:
   // std::array<std::array<T, 3>, N> near_1;
   //
   // gamma and lgamma for z near 2:
   // std::array<std::array<T, 3>, N> near_2;
   //
   // gamma and lgamma for z near -10:
   // std::array<std::array<T, 3>, N> near_m10;
   //
   // gamma and lgamma for z near -55:
   // std::array<std::array<T, 3>, N> near_m55;
   //
   // The last two cases are chosen more or less at random,
   // except that one is even and the other odd, and both are
   // at negative poles.  The data near zero also tests near
   // a pole, the data near 1 and 2 are to probe lgamma as
   // the result -> 0.
   //
#  include "tgamma_mp_data.hpp"

   do_test_gamma<T>(factorials, name, "factorials");
   do_test_gamma<T>(near_0, name, "near 0");
   do_test_gamma<T>(near_1, name, "near 1");
   do_test_gamma<T>(near_2, name, "near 2");
   do_test_gamma<T>(near_m10, name, "near -10");
   do_test_gamma<T>(near_m55, name, "near -55");
}

template <class T>
void test_tgamma1pm1(T, const char* name)
{
   std::cout << "Testing tgamma1pm1 with type " << name << std::endl;
   //
   // https://github.com/boostorg/math/issues/1519
   // The lgamma approximations used internally for these types are not accurate enough
   // to avoid losing digits to cancellation, errors were as high as 1500 epsilon.
   // The arguments are dyadic, reference values calculated with Arb:
   //
   static const std::array<std::array<T, 2>, 8> data = {{
      {{ T(-0.375), SC_(0.4345188480905567756360197394564231366322077722066673307706798580950941973020969146309569665256322935994654337852) }},
      {{ T(0.375), SC_(-0.111086430843774659257572435933755308792224698740403129584327399497564425740749328075507124360869447559309799521) }},
      {{ T(0.75), SC_(-0.080937473151116766153176272477832104861570563918947041577410796282632833679901232072120717838813700063369964821) }},
      {{ T(1.25), SC_(0.1330030963193463474783391112086475009359899009000204585729306248655997636058528505057375427260575938629599223101) }},
      {{ T(1.625), SC_(0.4569332050919717252553325478854297481420860186473965078139717308778300441349421789220656691275952981869570811881) }},
      {{ T(1.75), SC_(0.6083594219855456592319415231637938164922515131418426772395311065053925410601728438737887437820760248891025615618) }},
      {{ T(1.875), SC_(0.7877108988969403106484306138406837807894743299311615258427956959629808411549585601993268102781147763595978133741) }},
      {{ T(1.984375), SC_(0.9714651224878432889438635713198502884452632095994134781726907836796579737404575370382611229556842737635526601037) }},
   }};
   T tolerance = 40 * boost::math::tools::epsilon<T>();
   for(unsigned i = 0; i < data.size(); ++i)
   {
      BOOST_CHECK_CLOSE_FRACTION(boost::math::tgamma1pm1(data[i][0]), data[i][1], tolerance);
   }
}

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
      "cpp_bin_float_100|number<cpp_bin_float<85> >",           // test type(s)
      ".*",                          // test data group
      "lgamma", 600000, 300000);      // test function
   add_expected_result(
      ".*",                          // compiler
      ".*",                          // stdlib
      ".*",                          // platform
      "number<cpp_bin_float<[56]5> >",           // test type(s)
      ".*",                          // test data group
      "lgamma", 7000, 3000);      // test function
   add_expected_result(
      ".*",                          // compiler
      ".*",                          // stdlib
      ".*",                          // platform
      "number<cpp_bin_float<75> >",           // test type(s)
      ".*",                          // test data group
      "lgamma", 40000, 15000);      // test function
   add_expected_result(
      ".*",                          // compiler
      ".*",                          // stdlib
      ".*",                          // platform
      ".*",                          // test type(s)
      ".*",                          // test data group
      "lgamma", 600, 200);            // test function
   add_expected_result(
      ".*",                          // compiler
      ".*",                          // stdlib
      ".*",                          // platform
      ".*",                          // test type(s)
      ".*",                          // test data group
      "[tl]gamma", 120, 50);            // test function
   //
   // Finish off by printing out the compiler/stdlib/platform names,
   // we do this to make it easier to mark up expected error rates.
   //
   std::cout << "Tests run with " << BOOST_COMPILER << ", "
      << BOOST_STDLIB << ", " << BOOST_PLATFORM << std::endl;
}

BOOST_AUTO_TEST_CASE(test_main)
{
   expected_results();
   using namespace boost::multiprecision;
#if !defined(TEST) || (TEST == 1)
   test_gamma(number<cpp_bin_float<38> >(0), "number<cpp_bin_float<38> >");
   test_tgamma1pm1(number<cpp_bin_float<38> >(0), "number<cpp_bin_float<38> >");
   test_gamma(number<cpp_bin_float<45> >(0), "number<cpp_bin_float<45> >");
   test_tgamma1pm1(number<cpp_bin_float<45> >(0), "number<cpp_bin_float<45> >");
#endif
#if !defined(TEST) || (TEST == 2)
   test_gamma(cpp_bin_float_50(0), "cpp_bin_float_50");
   test_tgamma1pm1(cpp_bin_float_50(0), "cpp_bin_float_50");
   test_gamma(number<cpp_bin_float<55> >(0), "number<cpp_bin_float<55> >");
   test_tgamma1pm1(number<cpp_bin_float<55> >(0), "number<cpp_bin_float<55> >");
   test_gamma(number<cpp_bin_float<65> >(0), "number<cpp_bin_float<65> >");
   test_tgamma1pm1(number<cpp_bin_float<65> >(0), "number<cpp_bin_float<65> >");
#endif
#if !defined(TEST) || (TEST == 3)
   test_gamma(number<cpp_bin_float<75> >(0), "number<cpp_bin_float<75> >");
   test_tgamma1pm1(number<cpp_bin_float<75> >(0), "number<cpp_bin_float<75> >");
   test_gamma(number<cpp_bin_float<85> >(0), "number<cpp_bin_float<85> >");
   test_tgamma1pm1(number<cpp_bin_float<85> >(0), "number<cpp_bin_float<85> >");
   test_gamma(cpp_bin_float_100(0), "cpp_bin_float_100");
   test_tgamma1pm1(cpp_bin_float_100(0), "cpp_bin_float_100");
#endif
}
#else // No mp tests
int main(void) { return 0; }
#endif
