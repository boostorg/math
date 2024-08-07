// Copyright John Maddock 2006.
// Copyright Paul A. Bristow 2007.

// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

// test_exponential_dist.cpp

#include "../include_private/boost/math/tools/test.hpp"
#include <boost/math/concepts/real_concept.hpp> // for real_concept
#include <boost/math/distributions/exponential.hpp>
    using boost::math::exponential_distribution;

#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp> // Boost.Test
#include <boost/test/tools/floating_point_comparison.hpp>
#include "test_out_of_range.hpp"

#include <cmath>
#include <type_traits>
#include <iostream>
   using std::cout;
   using std::endl;
   using std::setprecision;
   using std::log;

#if __has_include(<stdfloat>)
#  include <stdfloat>
#endif

template <class RealType>
void test_spot(RealType l, RealType x, RealType p, RealType q, RealType logp, RealType logq, RealType tolerance, RealType logtolerance)
{
   BOOST_IF_CONSTEXPR (std::is_same<RealType, long double>::value)
   {
      logtolerance *= 100;
   }
   
   BOOST_CHECK_CLOSE(
      ::boost::math::cdf(
         exponential_distribution<RealType>(l),
         x),
         p,
         tolerance); // %
   BOOST_CHECK_CLOSE(
      ::boost::math::cdf(
         complement(exponential_distribution<RealType>(l),
         x)),
         q,
         tolerance); // %
   BOOST_CHECK_CLOSE(
      ::boost::math::logcdf(
         exponential_distribution<RealType>(l),
         x),
         logp,
         logtolerance); // %
   BOOST_CHECK_CLOSE(
      ::boost::math::logcdf(
         complement(exponential_distribution<RealType>(l),
         x)),
         logq,
         logtolerance); // %
   if(p < 0.999)
   {
      BOOST_CHECK_CLOSE(
         ::boost::math::quantile(
            exponential_distribution<RealType>(l),
            p),
            x,
            tolerance); // %
   }
   if(q < 0.999)
   {
      BOOST_CHECK_CLOSE(
         ::boost::math::quantile(
            complement(exponential_distribution<RealType>(l),
            q)),
            x,
            tolerance); // %
   }
}

template <class RealType>
void test_spots(RealType T)
{
   // Basic sanity checks.
   // 50 eps as a percentage, up to a maximum of double precision
   // (that's the limit of our test data: obtained by punching
   // numbers into a calculator).
   RealType tolerance = (std::max)(
      static_cast<RealType>(boost::math::tools::epsilon<double>()),
      boost::math::tools::epsilon<RealType>());
   tolerance *= 50 * 100;
   // #  pragma warning(disable: 4100) // unreferenced formal parameter.
   // prevent his spurious warning.
   if (T != 0)
   {
     cout << "Expect parameter T == 0!" << endl;
   }
    cout << "Tolerance for type " << typeid(T).name()  << " is " << tolerance << " %" << endl;

   test_spot(
      static_cast<RealType>(0.5), // lambda
      static_cast<RealType>(0.125), // x
      static_cast<RealType>(0.060586937186524213880289175377695L), // p
      static_cast<RealType>(0.93941306281347578611971082462231L), //q
      static_cast<RealType>(log(static_cast<RealType>(0.060586937186524213880289175377695L))),
      static_cast<RealType>(log(static_cast<RealType>(0.93941306281347578611971082462231L))),
      tolerance,
      tolerance);
   test_spot(
      static_cast<RealType>(0.5), // lambda
      static_cast<RealType>(5), // x
      static_cast<RealType>(0.91791500137610120483047132553284L), // p
      static_cast<RealType>(0.08208499862389879516952867446716L), //q
      static_cast<RealType>(log(static_cast<RealType>(0.91791500137610120483047132553284L))),
      static_cast<RealType>(log(static_cast<RealType>(0.08208499862389879516952867446716L))),
      tolerance,
      tolerance);
   test_spot(
      static_cast<RealType>(2), // lambda
      static_cast<RealType>(0.125), // x
      static_cast<RealType>(0.22119921692859513175482973302168L), // p
      static_cast<RealType>(0.77880078307140486824517026697832L), //q
      static_cast<RealType>(log(static_cast<RealType>(0.22119921692859513175482973302168L))), // log(p)
      static_cast<RealType>(log(static_cast<RealType>(0.77880078307140486824517026697832L))), // log(q)
      tolerance,
      tolerance);
   test_spot(
      static_cast<RealType>(2), // lambda
      static_cast<RealType>(5), // x
      static_cast<RealType>(0.99995460007023751514846440848444L), // p
      static_cast<RealType>(4.5399929762484851535591515560551e-5L), //q
      static_cast<RealType>(-0.00004540096037048920950444635987890882815054L),
      static_cast<RealType>(-10.0000000000000000000000000000000000000000000L),
      tolerance,
      tolerance);

   //
   // Some spot tests generated by MathCAD pexp(x,r):
   //
   test_spot(
      static_cast<RealType>(1), // lambda
      static_cast<RealType>(1), // x
      static_cast<RealType>(6.321205588285580E-001L), // p
      static_cast<RealType>(1-6.321205588285580E-001L), //q
      static_cast<RealType>(-0.4586751453870818910216436450673297018770L), // log(p)
      static_cast<RealType>(-1.0000000000000000000000000000000000000000L),
      tolerance,
      tolerance);
   test_spot(
      static_cast<RealType>(2), // lambda
      static_cast<RealType>(1), // x
      static_cast<RealType>(8.646647167633870E-001L), // p
      static_cast<RealType>(1-8.646647167633870E-001L), //q
      static_cast<RealType>(-0.1454134578688590569726481500994740599617L), //log(p)
      static_cast<RealType>(-2.000000000000000000000000000000000000000L),
      tolerance,
      tolerance);
   test_spot(
      static_cast<RealType>(1), // lambda
      static_cast<RealType>(0.5), // x
      static_cast<RealType>(3.934693402873670E-001L), // p
      static_cast<RealType>(1-3.934693402873670E-001L), //q
      static_cast<RealType>(-0.9327521295671885718946410001485004516325L), 
      static_cast<RealType>(-0.5000000000000000000000000000000000000000L),
      tolerance,
      tolerance);
   test_spot(
      static_cast<RealType>(0.1), // lambda
      static_cast<RealType>(1), // x
      static_cast<RealType>(9.516258196404040E-002L), // p
      static_cast<RealType>(1-9.516258196404040E-002L), //q
      static_cast<RealType>(-2.352168461044090808919497187062612584712L),
      static_cast<RealType>(-0.100000000000000000000000000000000000000L),
      tolerance,
      tolerance);
   test_spot(
      static_cast<RealType>(10), // lambda
      static_cast<RealType>(1), // x
      static_cast<RealType>(9.999546000702380E-001L), // p
      static_cast<RealType>(1-9.999546000702380E-001L), //q
      static_cast<RealType>(-0.00004540096037048920950444635987890882815054L),
      static_cast<RealType>(-10.0000000000000000000000000000000000000000000L),
      tolerance*10000, // we loose four digits to cancellation
      tolerance); 
   test_spot(
      static_cast<RealType>(0.1), // lambda
      static_cast<RealType>(10), // x
      static_cast<RealType>(6.321205588285580E-001L), // p
      static_cast<RealType>(1-6.321205588285580E-001L), //q
      static_cast<RealType>(-0.4586751453870818910216436450673297018770L), // log(p),
      static_cast<RealType>(-1.0000000000000000000000000000000000000000L), // log(q)
      tolerance,
      tolerance);
   test_spot(
      static_cast<RealType>(1), // lambda
      static_cast<RealType>(0.01), // x
      static_cast<RealType>(9.950166250831950E-003L), // p
      static_cast<RealType>(1-9.950166250831950E-003L), //q
      static_cast<RealType>(-4.610166019324896918080084954100938122002L), // log(p),
      static_cast<RealType>(-0.010000000000000000000000000000000000000L), // log(q)
      tolerance,
      tolerance);
   test_spot(
      static_cast<RealType>(1), // lambda
      static_cast<RealType>(0.0001), // x
      static_cast<RealType>(9.999500016666250E-005L), // p
      static_cast<RealType>(1-9.999500016666250E-005L), //q
      static_cast<RealType>(-9.210390371559516069440021374287500922116L), // log(p),
      static_cast<RealType>(-0.000100000000000000000000000000000000000L), // log(q)
      tolerance,
      std::is_same<RealType, float>::value
      #ifdef __STDCPP_FLOAT32_T__
      || std::is_same<RealType, std::float32_t>::value
      #endif
       ? tolerance * 10 : tolerance);
   /*
   // This test data appears to be erroneous, MathCad appears
   // to suffer from cancellation error as x -> 0
   test_spot(
      static_cast<RealType>(1), // lambda
      static_cast<RealType>(0.0000001), // x
      static_cast<RealType>(9.999999499998730E-008L), // p
      static_cast<RealType>(1-9.999999499998730E-008L), //q
      tolerance);
   */

   BOOST_CHECK_CLOSE(
      ::boost::math::pdf(
         exponential_distribution<RealType>(0.5),
         static_cast<RealType>(0.125)),              // x
         static_cast<RealType>(0.46970653140673789305985541231115L),                // probability.
         tolerance); // %
   BOOST_CHECK_CLOSE(
      ::boost::math::pdf(
         exponential_distribution<RealType>(0.5),
         static_cast<RealType>(5)),              // x
         static_cast<RealType>(0.04104249931194939758476433723358L),                // probability.
         tolerance); // %
   BOOST_CHECK_CLOSE(
      ::boost::math::pdf(
         exponential_distribution<RealType>(2),
         static_cast<RealType>(0.125)),              // x
         static_cast<RealType>(1.5576015661428097364903405339566L),                // probability.
         tolerance); // %
   BOOST_CHECK_CLOSE(
      ::boost::math::pdf(
         exponential_distribution<RealType>(2),
         static_cast<RealType>(5)),              // x
         static_cast<RealType>(9.0799859524969703071183031121101e-5L),                // probability.
         tolerance); // %

   BOOST_CHECK_CLOSE(
      ::boost::math::logpdf(
         exponential_distribution<RealType>(0.5),
         static_cast<RealType>(0.125)),              // x
         log(static_cast<RealType>(0.46970653140673789305985541231115L)),                // probability.
         tolerance); // %
   BOOST_CHECK_CLOSE(
      ::boost::math::logpdf(
         exponential_distribution<RealType>(0.5),
         static_cast<RealType>(5)),              // x
         log(static_cast<RealType>(0.04104249931194939758476433723358L)),                // probability.
         tolerance); // %
   BOOST_CHECK_CLOSE(
      ::boost::math::logpdf(
         exponential_distribution<RealType>(2),
         static_cast<RealType>(0.125)),              // x
         log(static_cast<RealType>(1.5576015661428097364903405339566L)),                // probability.
         tolerance); // %
   BOOST_CHECK_CLOSE(
      ::boost::math::logpdf(
         exponential_distribution<RealType>(2),
         static_cast<RealType>(5)),              // x
         log(static_cast<RealType>(9.0799859524969703071183031121101e-5L)),                // probability.
         tolerance); // %

   BOOST_CHECK_CLOSE(
      ::boost::math::mean(
         exponential_distribution<RealType>(2)),
         static_cast<RealType>(0.5),
         tolerance); // %
   BOOST_CHECK_CLOSE(
      ::boost::math::standard_deviation(
         exponential_distribution<RealType>(2)),
         static_cast<RealType>(0.5),
         tolerance); // %
   BOOST_CHECK_CLOSE(
      ::boost::math::mode(
         exponential_distribution<RealType>(2)),
         static_cast<RealType>(0),
         tolerance); // %

   BOOST_CHECK_CLOSE(
      ::boost::math::median(
         exponential_distribution<RealType>(4)),
         static_cast<RealType>(0.693147180559945309417232121458176568075500134360255254) / 4,
         tolerance); // %

   BOOST_CHECK_CLOSE(
      ::boost::math::skewness(
         exponential_distribution<RealType>(2)),
         static_cast<RealType>(2),
         tolerance); // %
   BOOST_CHECK_CLOSE(
      ::boost::math::kurtosis(
         exponential_distribution<RealType>(2)),
         static_cast<RealType>(9),
         tolerance); // %
   BOOST_CHECK_CLOSE(
      ::boost::math::kurtosis_excess(
         exponential_distribution<RealType>(2)),
         static_cast<RealType>(6),
         tolerance); // %

   using std::log;
   BOOST_CHECK_CLOSE(
      ::boost::math::entropy(
         exponential_distribution<RealType>(2)),
         static_cast<RealType>(1-log(RealType(2))),
         tolerance); // %

   //
   // Things that are errors:
   //
   exponential_distribution<RealType> dist(0.5);
   BOOST_MATH_CHECK_THROW(
       quantile(dist, RealType(1.0)),
       std::overflow_error);
   BOOST_MATH_CHECK_THROW(
       quantile(complement(dist, RealType(0.0))),
       std::overflow_error);
   BOOST_MATH_CHECK_THROW(
       pdf(dist, RealType(-1)),
       std::domain_error);
   BOOST_MATH_CHECK_THROW(
       cdf(dist, RealType(-1)),
       std::domain_error);
   BOOST_MATH_CHECK_THROW(
       cdf(exponential_distribution<RealType>(-1), RealType(1)),
       std::domain_error);
   BOOST_MATH_CHECK_THROW(
       quantile(dist, RealType(-1)),
       std::domain_error);
   BOOST_MATH_CHECK_THROW(
       quantile(dist, RealType(2)),
       std::domain_error);

   check_out_of_range<exponential_distribution<RealType> >(2);
   BOOST_MATH_CHECK_THROW(exponential_distribution<RealType>(0), std::domain_error);
   BOOST_MATH_CHECK_THROW(exponential_distribution<RealType>(-1), std::domain_error);
   if(std::numeric_limits<RealType>::has_infinity)
   {
      RealType inf = std::numeric_limits<RealType>::infinity();
      BOOST_CHECK_EQUAL(pdf(exponential_distribution<RealType>(2), inf), 0);
      BOOST_CHECK_EQUAL(cdf(exponential_distribution<RealType>(2), inf), 1);
      BOOST_CHECK_EQUAL(cdf(complement(exponential_distribution<RealType>(2), inf)), 0);
   }
} // template <class RealType>void test_spots(RealType)

BOOST_AUTO_TEST_CASE( test_main )
{
  // Check that can generate exponential distribution using the two convenience methods:
   boost::math::exponential mycexp1(1.); // Using typedef
   exponential_distribution<> myexp2(1.); // Using default RealType double.

    // Basic sanity-check spot values.
   // (Parameter value, arbitrarily zero, only communicates the floating point type).
  test_spots(0.0F); // Test float. OK at decdigits = 0 tolerance = 0.0001 %
  test_spots(0.0); // Test double. OK at decdigits 7, tolerance = 1e07 %
#ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
  test_spots(0.0L); // Test long double.
#if !BOOST_WORKAROUND(BOOST_BORLANDC, BOOST_TESTED_AT(0x582)) && !defined(BOOST_MATH_NO_REAL_CONCEPT_TESTS)
  test_spots(boost::math::concepts::real_concept(0.)); // Test real concept.
#endif
#else
   std::cout << "<note>The long double tests have been disabled on this platform "
      "either because the long double overloads of the usual math functions are "
      "not available at all, or because they are too inaccurate for these tests "
      "to pass.</note>" << std::endl;
#endif

   #ifdef __STDCPP_FLOAT32_T__
   test_spots(0.0F32);
   #endif
   #ifdef __STDCPP_FLOAT64_T__
   test_spots(0.0F64);
   #endif 

} // BOOST_AUTO_TEST_CASE( test_main )

/*

Output is:

Running 1 test case...
Tolerance for type float is 0.000596046 %
Tolerance for type double is 1.11022e-012 %
Tolerance for type long double is 1.11022e-012 %
Tolerance for type class boost::math::concepts::real_concept is 1.11022e-012 %
*** No errors detected

*/
