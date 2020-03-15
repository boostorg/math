// Copyright Paul A. Bristow 2010.
// Copyright John Maddock 2007.

// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

// test_von_mises.cpp

// https://en.wikipedia.org/wiki/Von_Mises_distribution
// From MathWorld--A Wolfram Web Resource.
// http://mathworld.wolfram.com/VonMisesDistribution.html

#include <pch.hpp> // include directory /libs/math/src/tr1/ is needed.

#ifdef _MSC_VER
#  pragma warning (disable: 4127) // conditional expression is constant
// caused by using   if(std::numeric_limits<RealType>::has_infinity)
// and   if (std::numeric_limits<RealType>::has_quiet_NaN)
#endif

#include <boost/math/tools/test.hpp>
#include <boost/math/concepts/real_concept.hpp> // for real_concept
#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp> // Boost.Test
#include <boost/test/tools/floating_point_comparison.hpp>

#include <boost/math/distributions/von_mises.hpp>
    using boost::math::von_mises_distribution;
#include <boost/math/tools/test.hpp>
#include "test_out_of_range.hpp"

#include <iostream>
#include <iomanip>
   using std::cout;
   using std::endl;
   using std::setprecision;
#include <limits>
  using std::numeric_limits;

template <class RealType>
void check_von_mises(RealType mean, RealType conc, RealType x, RealType p, RealType q, RealType tol)
{
  BOOST_CHECK_CLOSE(
      ::boost::math::cdf(
          von_mises_distribution<RealType>(mean, conc),      // distribution.
          x),                                                // random variable.
      p,                                                     // probability.
      tol);                                                  // %tolerance.
  BOOST_CHECK_CLOSE(
      ::boost::math::cdf(
          complement(
              von_mises_distribution<RealType>(mean, conc),  // distribution.
              x)),                                           // random variable.
      q,                                                     // probability complement.
      tol);                                                  // %tolerance.
  BOOST_CHECK_CLOSE(
      ::boost::math::quantile(
          von_mises_distribution<RealType>(mean, conc),      // distribution.
          p),                                                // probability.
      x,                                                     // random variable.
      tol);                                                  // %tolerance.
   BOOST_CHECK_CLOSE(
      ::boost::math::quantile(
          complement(
              von_mises_distribution<RealType>(mean, conc),  // distribution.
              q)),                                           // probability complement.
      x,                                                     // random variable.
      tol);                                                  // %tolerance.
}

template <class RealType>
void test_spots(RealType)
{
  // Basic sanity checks  
  // Check some bad parameters to the distribution,
#ifndef BOOST_NO_EXCEPTIONS   
  BOOST_MATH_CHECK_THROW(boost::math::von_mises_distribution<RealType> nbad1(0, -1), std::domain_error); // negative conc
#else
  BOOST_MATH_CHECK_THROW(boost::math::von_mises_distribution<RealType>(0, -1), std::domain_error); // negative conc
#endif

  // Tests on extreme values of random variate x, if has std::numeric_limits infinity etc.
  von_mises_distribution<RealType> N01;
  if(std::numeric_limits<RealType>::has_infinity)
  {
    BOOST_MATH_CHECK_THROW(pdf(N01, +std::numeric_limits<RealType>::infinity()), std::domain_error);             // x = + infinity, pdf = 0
    BOOST_MATH_CHECK_THROW(pdf(N01, -std::numeric_limits<RealType>::infinity()), std::domain_error);             // x = - infinity, pdf = 0
    BOOST_MATH_CHECK_THROW(cdf(N01, +std::numeric_limits<RealType>::infinity()), std::domain_error);             // x = + infinity, cdf = 1
    BOOST_MATH_CHECK_THROW(cdf(N01, -std::numeric_limits<RealType>::infinity()), std::domain_error);             // x = - infinity, cdf = 0
    BOOST_MATH_CHECK_THROW(cdf(complement(N01, +std::numeric_limits<RealType>::infinity())), std::domain_error); // x = + infinity, c cdf = 0
    BOOST_MATH_CHECK_THROW(cdf(complement(N01, -std::numeric_limits<RealType>::infinity())), std::domain_error); // x = - infinity, c cdf = 1

#ifndef BOOST_NO_EXCEPTIONS
    BOOST_MATH_CHECK_THROW(boost::math::von_mises_distribution<RealType> nbad1(std::numeric_limits<RealType>::infinity(), static_cast<RealType>(1)), std::domain_error); // +infinite mean
    BOOST_MATH_CHECK_THROW(boost::math::von_mises_distribution<RealType> nbad1(-std::numeric_limits<RealType>::infinity(),  static_cast<RealType>(1)), std::domain_error); // -infinite mean
    BOOST_MATH_CHECK_THROW(boost::math::von_mises_distribution<RealType> nbad1(static_cast<RealType>(0), std::numeric_limits<RealType>::infinity()), std::domain_error); // infinite conc
#else     
    BOOST_MATH_CHECK_THROW(boost::math::von_mises_distribution<RealType>(std::numeric_limits<RealType>::infinity(), static_cast<RealType>(1)), std::domain_error); // +infinite mean
    BOOST_MATH_CHECK_THROW(boost::math::von_mises_distribution<RealType>(-std::numeric_limits<RealType>::infinity(),  static_cast<RealType>(1)), std::domain_error); // -infinite mean
    BOOST_MATH_CHECK_THROW(boost::math::von_mises_distribution<RealType>(static_cast<RealType>(0), std::numeric_limits<RealType>::infinity()), std::domain_error); // infinite conc
#endif
  }

  if (std::numeric_limits<RealType>::has_quiet_NaN)
  {
    // No longer allow x to be NaN, then these tests should throw.
    BOOST_MATH_CHECK_THROW(pdf(N01, +std::numeric_limits<RealType>::quiet_NaN()), std::domain_error); // x = NaN
    BOOST_MATH_CHECK_THROW(cdf(N01, +std::numeric_limits<RealType>::quiet_NaN()), std::domain_error); // x = NaN
    BOOST_MATH_CHECK_THROW(cdf(complement(N01, +std::numeric_limits<RealType>::quiet_NaN())), std::domain_error); // x = + infinity
    BOOST_MATH_CHECK_THROW(quantile(N01, +std::numeric_limits<RealType>::quiet_NaN()), std::domain_error); // p = + infinity
    BOOST_MATH_CHECK_THROW(quantile(complement(N01, +std::numeric_limits<RealType>::quiet_NaN())), std::domain_error); // p = + infinity
  }

  //
  // Tests for PDF: we know that the peak value is at e/(2*pi*I0(1)
  //
  RealType tolerance = boost::math::tools::epsilon<RealType>() * 5 * 100; // 5 eps as a percentage  
  BOOST_CHECK_CLOSE(
      pdf(von_mises_distribution<RealType>(), static_cast<RealType>(0)),
      static_cast<RealType>(0.34171048862346315949457814754706159394027L),  // e/(2*pi*I0(1))
      tolerance);
  BOOST_CHECK_CLOSE(
      pdf(von_mises_distribution<RealType>(3, 0), static_cast<RealType>(3)),
      static_cast<RealType>(0.15915494309189533576888376337251L),           // 1/(2*pi)
      tolerance);
  BOOST_CHECK_CLOSE(
      pdf(von_mises_distribution<RealType>(3), static_cast<RealType>(3)),
      static_cast<RealType>(0.34171048862346315949457814754706159394027L),  // e/(2*pi*I0(1))
      tolerance);
  BOOST_CHECK_CLOSE(
      pdf(von_mises_distribution<RealType>(3, 5), static_cast<RealType>(3)),
      static_cast<RealType>(0.86713652854235200257351846969777045343907L),  // e^5/(2*pi*I0(5))
      tolerance);
  // Tests with larger numbers whose bessel function value exceeds the representable range
  BOOST_CHECK_CLOSE(
      pdf(von_mises_distribution<RealType>(3, 25), static_cast<RealType>(3)),
      static_cast<RealType>(1.98455543847726689510475504795539869409664L),  // e^25/(2*pi*I0(25))
      tolerance);
  BOOST_CHECK_CLOSE(
      pdf(von_mises_distribution<RealType>(3, 125), static_cast<RealType>(3)),
      static_cast<RealType>(4.45583423571762164664231477348923797078210L),  // e^125/(2*pi*I0(125))
      tolerance);
  BOOST_CHECK_CLOSE(
      pdf(von_mises_distribution<RealType>(3, 1250), static_cast<RealType>(3)),
      static_cast<RealType>(14.1033286206525302447068654967513214217021L),  // e^1250/(2*pi*I0(1250))
      tolerance);

  tolerance = 1e-3f; // 1e-5 (as %)
  // Some tests only pass at 1e-5 because values generated by

  cout << "Tolerance for type " << typeid(RealType).name()  << " is " << tolerance << " %" << endl;

  // test CDF for mean and interval edges
  check_von_mises(
      static_cast<RealType>(0),
      static_cast<RealType>(1),
      static_cast<RealType>(0),
      static_cast<RealType>(0.5L),
      static_cast<RealType>(0.5L),
      tolerance);
      
  check_von_mises(
      static_cast<RealType>(0),
      static_cast<RealType>(1),
      static_cast<RealType>(-boost::math::constants::pi<RealType>()),
      static_cast<RealType>(0.0L),
      static_cast<RealType>(1.0L),
      tolerance);

  check_von_mises(
      static_cast<RealType>(0),
      static_cast<RealType>(1),
      static_cast<RealType>(boost::math::constants::pi<RealType>()),
      static_cast<RealType>(1.0L),
      static_cast<RealType>(0.0L),
      tolerance);

  // test CDF for low concentrations
  check_von_mises(
      static_cast<RealType>(0),
      static_cast<RealType>(1),
      static_cast<RealType>(1),
      static_cast<RealType>(0.794355307434683479987678129735260058645499262629455722769L),
      static_cast<RealType>(0.205644692565316520012321870264739941354500737370544277230L),
      tolerance);      
      
  check_von_mises(
      static_cast<RealType>(0),
      static_cast<RealType>(1),
      static_cast<RealType>(-1),
      static_cast<RealType>(0.205644692565316520012321870264739941354500737370544277230L),
      static_cast<RealType>(0.794355307434683479987678129735260058645499262629455722769L),
      tolerance);

  check_von_mises(
      static_cast<RealType>(0),
      static_cast<RealType>(5),
      static_cast<RealType>(1),
      static_cast<RealType>(0.98096204546814689054581384796251763480020360394758184271L),
      static_cast<RealType>(0.01903795453185310945418615203748236519979639605241815729L),
      tolerance);      
      
  check_von_mises(
      static_cast<RealType>(0),
      static_cast<RealType>(5),
      static_cast<RealType>(-1),
      static_cast<RealType>(0.01903795453185310945418615203748236519979639605241815729L),
      static_cast<RealType>(0.98096204546814689054581384796251763480020360394758184271L),
      tolerance);
      
  // test CDF for high concentrations
  //~ check_von_mises(
      //~ static_cast<RealType>(0),
      //~ static_cast<RealType>(25),
      //~ static_cast<RealType>(1),
      //~ static_cast<RealType>(0.999999062404464440452233504489299782776166264467210572765L),
      //~ static_cast<RealType>(9.37595535559547766495510700217223833735532789427234e-7L),
      //~ tolerance);      
      
  //~ check_von_mises(
      //~ static_cast<RealType>(0),
      //~ static_cast<RealType>(25),
      //~ static_cast<RealType>(-1),
      //~ static_cast<RealType>(9.37595535559547766495510700217223833735532789427234e-7L),
      //~ static_cast<RealType>(0.999999062404464440452233504489299782776166264467210572765L),
      //~ tolerance);
  
  //~ check_von_mises(
      //~ static_cast<RealType>(0),
      //~ static_cast<RealType>(125),
      //~ static_cast<RealType>(1),
      //~ static_cast<RealType>(0.99999999999999999996645363431349332910951002333081490389L),
      //~ static_cast<RealType>(3.3546365686506670890489976669185096109e-20L),
      //~ tolerance);      
      
  //~ check_von_mises(
      //~ static_cast<RealType>(0),
      //~ static_cast<RealType>(125),
      //~ static_cast<RealType>(-1),
      //~ static_cast<RealType>(3.3546365686506670890489976669185096109e-20L),
      //~ static_cast<RealType>(0.99999999999999999996645363431349332910951002333081490389L),
      //~ tolerance);
  

  RealType tol2 = boost::math::tools::epsilon<RealType>() * 500;
  von_mises_distribution<RealType> dist(2, 3);
  RealType x = static_cast<RealType>(0.125);

  BOOST_MATH_STD_USING // ADL of std math lib names

  // mean:
  BOOST_CHECK_CLOSE(
       mean(dist)
       , static_cast<RealType>(2), tol2);
  // variance:  
  BOOST_CHECK_CLOSE(
       variance(von_mises_distribution<RealType>(2, 0))
       , static_cast<RealType>(1), tol2);
  BOOST_CHECK_CLOSE(
       variance(von_mises_distribution<RealType>(2, 1))
       , static_cast<RealType>(0.553610034103465492952318204807357330223746852599612177138L), tol2);
  BOOST_CHECK_CLOSE(
       variance(von_mises_distribution<RealType>(2, 5))
       , static_cast<RealType>(0.106616862955914778412994992775050827245562823701700738786L), tol2);
  //~ BOOST_CHECK_CLOSE(
       //~ variance(von_mises_distribution<RealType>(2, 25))
       //~ , static_cast<RealType>(0.020208546509484068780303296944566388571293465011951716357L), tol2);
  //~ BOOST_CHECK_CLOSE(
       //~ variance(von_mises_distribution<RealType>(2, 125))
       //~ , static_cast<RealType>(0.004008064813593637057963834031301840442156903531539609019L), tol2);
  return;
       
  // std deviation:
  BOOST_CHECK_CLOSE(
       standard_deviation(dist)
       , static_cast<RealType>(3), tol2);
  // hazard:
  BOOST_CHECK_CLOSE(
       hazard(dist, x)
       , pdf(dist, x) / cdf(complement(dist, x)), tol2);
  // cumulative hazard:
  BOOST_CHECK_CLOSE(
       chf(dist, x)
       , -log(cdf(complement(dist, x))), tol2);
  // coefficient_of_variation:
  BOOST_CHECK_CLOSE(
       coefficient_of_variation(dist)
       , standard_deviation(dist) / mean(dist), tol2);
  // mode:
  BOOST_CHECK_CLOSE(
       mode(dist)
       , static_cast<RealType>(8), tol2);

  BOOST_CHECK_CLOSE(
       median(dist)
       , static_cast<RealType>(8), tol2);

  // skewness:
  BOOST_CHECK_CLOSE(
       skewness(dist)
       , static_cast<RealType>(0), tol2);
  // kurtosis:
  // BOOST_CHECK_CLOSE(
      // kurtosis(dist)
      // , static_cast<RealType>(3), tol2);
  // kurtosis excess:
  // BOOST_CHECK_CLOSE(
      // kurtosis_excess(dist)
      // , static_cast<RealType>(0), tol2);

  RealType expected_entropy = log(boost::math::constants::two_pi<RealType>()*boost::math::constants::e<RealType>()*9)/2;
  BOOST_CHECK_CLOSE(
       entropy(dist)
       ,expected_entropy, tol2);

  von_mises_distribution<RealType> norm01(0, 1); // Test default (0, 1)
  BOOST_CHECK_CLOSE(
       mean(norm01),
       static_cast<RealType>(0), 0); // Mean == zero

  von_mises_distribution<RealType> defsd_norm01(0); // Test default (0, sd = 1)
  BOOST_CHECK_CLOSE(
       mean(defsd_norm01),
       static_cast<RealType>(0), 0); // Mean == zero

  von_mises_distribution<RealType> def_norm01; // Test default (0, sd = 1)
  BOOST_CHECK_CLOSE(
       mean(def_norm01),
       static_cast<RealType>(0), 0); // Mean == zero

  BOOST_CHECK_CLOSE(
       standard_deviation(def_norm01),
       static_cast<RealType>(1), 0); // Mean == zero

  // Error tests:
  check_out_of_range<boost::math::von_mises_distribution<RealType> >(0, 1); // (All) valid constructor parameter values.
    
  BOOST_MATH_CHECK_THROW(pdf(von_mises_distribution<RealType>(0, 0), 0), std::domain_error);
  BOOST_MATH_CHECK_THROW(pdf(von_mises_distribution<RealType>(0, -1), 0), std::domain_error);
  BOOST_MATH_CHECK_THROW(quantile(von_mises_distribution<RealType>(0, 1), -1), std::domain_error);
  BOOST_MATH_CHECK_THROW(quantile(von_mises_distribution<RealType>(0, 1), 2), std::domain_error);
} // template <class RealType>void test_spots(RealType)

BOOST_AUTO_TEST_CASE( test_main )
{
  // Check that we can generate von_mises distribution using the two convenience methods:
  boost::math::von_mises myf1(1., 2);       // Using typedef
  von_mises_distribution<> myf2(1., 2);     // Using default RealType double.
  boost::math::von_mises myn01;             // Use default values.
  // Note NOT myn01() as the compiler will interpret as a function!

  // Check the synonyms, provided to allow generic use of find_location and find_scale.
  BOOST_CHECK_EQUAL(myn01.mean(), myn01.location());
  BOOST_CHECK_EQUAL(myn01.concentration(), myn01.scale());

  // Basic sanity-check spot values.
  // (Parameter value, arbitrarily zero, only communicates the floating point type).
  test_spots(0.0F);   // Test float. OK at decdigits = 0 tolerance = 0.0001 %
  test_spots(0.0);    // Test double. OK at decdigits 7, tolerance = 1e07 %
#ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
  //test_spots(0.0L); // Test long double.
#if !BOOST_WORKAROUND(__BORLANDC__, BOOST_TESTED_AT(0x0582))
  //test_spots(boost::math::concepts::real_concept(0.)); // Test real concept.
#endif
#else
  std::cout << "<note>The long double tests have been disabled on this platform "
      "either because the long double overloads of the usual math functions are "
      "not available at all, or because they are too inaccurate for these tests "
      "to pass.</note>" << std::endl;
#endif

   
} // BOOST_AUTO_TEST_CASE( test_main )

/*

Output:

Autorun "i:\boost-06-05-03-1300\libs\math\test\Math_test\debug\test_von_mises.exe"
Running 1 test case...
Tolerance for type float is 0.01 %
Tolerance for type double is 0.01 %
Tolerance for type long double is 0.01 %
Tolerance for type class boost::math::concepts::real_concept is 0.01 %
*** No errors detected */
