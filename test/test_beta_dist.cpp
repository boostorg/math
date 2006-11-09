// test_beta_dist.cpp

// Copyright John Maddock 2006.
// Copyright  Paul A. Bristow 2006.

// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

// Basic sanity tests for the beta Distribution.

// http://members.aol.com/iandjmsmith/BETAEX.HTM  beta distribution calculator
// Appreas to be a 64-bit calculator showing 17 decimal digit (last is noisy).
// Similar to mathCAD?


// http://www.nuhertz.com/statmat/distributions.html#Beta
// Pretty graphs and explanations for most distributions.

// http://functions.wolfram.com/webMathematica/FunctionEvaluation.jsp
// provided 40 decimal digits accuracy incomplete beta aka beta regularized == cdf

// http://www.ausvet.com.au/pprev/content.php?page=PPscript
// mode 0.75 	5/95% 0.9 	alpha 7.39 	beta 3.13
// http://www.epi.ucdavis.edu/diagnostictests/betabuster.html
// Beta Buster also calculates alpha and beta from mode & percentile estimates.
// This is NOT (yet) implemented.

#define BOOST_MATH_THROW_ON_DOMAIN_ERROR

#ifdef _MSC_VER
#  pragma warning(disable: 4127) // conditional expression is constant.
#  pragma warning(disable: 4100) // unreferenced formal parameter.
#  pragma warning(disable: 4512) // assignment operator could not be generated.
#endif

#include <boost/math/distributions/beta.hpp> // for beta_distribution
using boost::math::beta_distribution;
using boost::math::beta;

#include <boost/math/concepts/real_concept.hpp> // for real_concept
using ::boost::math::concepts::real_concept;

#include <boost/test/included/test_exec_monitor.hpp> // for test_main
#include <boost/test/floating_point_comparison.hpp> // for BOOST_CHECK_CLOSE_FRACTION

#include <iostream>
using std::cout;
using std::endl;
#include <limits>
using std::numeric_limits;

template <class RealType>
void test_spot(
     RealType a,    // alpha a
     RealType b,    // beta b
     RealType x,    // Probability 
     RealType P,    // CDF of beta(a, b)
     RealType Q,    // Complement of CDF
     RealType tol)  // Test tolerance.
{
   boost::math::beta_distribution<RealType> abeta(a, b);
   BOOST_CHECK_CLOSE_FRACTION(cdf(abeta, x), P, tol);
   if((P < 0.99) && (Q < 0.99))
   {  // We can only check this if P is not too close to 1,
      // so that we can guarantee that Q is free of error,
      // (and similarly for Q)
      BOOST_CHECK_CLOSE_FRACTION(
         cdf(complement(abeta, x)), Q, tol);
      if(x != 0)
      {
         BOOST_CHECK_CLOSE_FRACTION(
            quantile(abeta, P), x, tol);
      }
      else
      {
         // Just check quantile is very small:
         if((std::numeric_limits<RealType>::max_exponent <= std::numeric_limits<double>::max_exponent)
           && (boost::is_floating_point<RealType>::value))
         {
            // Limit where this is checked: if exponent range is very large we may
            // run out of iterations in our root finding algorithm.
            BOOST_CHECK(quantile(abeta, P) < boost::math::tools::epsilon<RealType>() * 10);
         }
      } // if k
      if(x != 0)
      {
         BOOST_CHECK_CLOSE_FRACTION(quantile(complement(abeta, Q)), x, tol);
      }
      else
      {  // Just check quantile is very small:
         if((std::numeric_limits<RealType>::max_exponent <= std::numeric_limits<double>::max_exponent) && (boost::is_floating_point<RealType>::value))
         {  // Limit where this is checked: if exponent range is very large we may
            // run out of iterations in our root finding algorithm.
            BOOST_CHECK(quantile(complement(abeta, Q)) < boost::math::tools::epsilon<RealType>() * 10);
         }
      } // if x
      // Estimate alpha & beta from mean and variance:

      BOOST_CHECK_CLOSE_FRACTION(
         beta_distribution<RealType>::estimate_alpha(mean(abeta), variance(abeta)),
         abeta.alpha(), tol);
      BOOST_CHECK_CLOSE_FRACTION(
         beta_distribution<RealType>::estimate_beta(mean(abeta), variance(abeta)),
         abeta.beta(), tol);

      // Estimate sample alpha and beta from others:
      BOOST_CHECK_CLOSE_FRACTION(
         beta_distribution<RealType>::estimate_alpha(abeta.beta(), x, P),
         abeta.alpha(), tol);
      BOOST_CHECK_CLOSE_FRACTION(
         beta_distribution<RealType>::estimate_beta(abeta.alpha(), x, P),
         abeta.beta(), tol);
   } // if((P < 0.99) && (Q < 0.99)
  
} // template <class RealType> void test_spot

template <class RealType> // Any floating-point type RealType.
void test_spots(RealType)
{
  // Basic sanity checks with 'known good' values.
  // MathCAD test data is to double precision only,
  // so set tolerance to 100 eps expressed as a fraction, or
  // 100 eps of type double expressed as a fraction,
  // whichever is the larger.

  RealType tolerance = (std::max)
      (boost::math::tools::epsilon<RealType>(),
      static_cast<RealType>(std::numeric_limits<double>::epsilon())); // 0 if real_concept.

   cout << "Boost::math::tools::epsilon = " << boost::math::tools::epsilon<RealType>() <<endl;
   cout << "std::numeric_limits::epsilon = " << std::numeric_limits<RealType>::epsilon() <<endl;
   cout << "epsilon = " << tolerance;

   tolerance *= 1000; // Note: NO * 100 because is fraction, NOT %.
   cout  << ", Tolerance = " << tolerance * 100 << "%." << endl;

  // RealType teneps = boost::math::tools::epsilon<RealType>() * 10;

  // Sources of spot test values:

  // MathCAD defines dbeta(x, s1, s2) pdf, s1 == alpha, s2 = beta, x = x in Wolfram
  // pbeta(x, s1, s2) cdf and qbeta(x, s1, s2) inverse of cdf
  // returns pr(X ,= x) when random variable X
  // has the beta distribution with parameters s1)alpha) and s2(beta).
  // s1 > 0 and s2 >0 and 0 < x < 1 (but allows x == 0! and x == 1!)
  // dbeta(0,1,1) = 0
  // dbeta(0.5,1,1) = 1

  using boost::math::beta_distribution;
  using  ::boost::math::cdf;
  using  ::boost::math::pdf;

  // Tests that should throw:
  BOOST_CHECK_THROW(mode(beta_distribution<RealType>(static_cast<RealType>(1), static_cast<RealType>(1))), std::domain_error);
  // mode is undefined, and throws domain_error!

  BOOST_CHECK_THROW( // For various bad arguments.
       pdf(
          beta_distribution<RealType>(static_cast<RealType>(-1), static_cast<RealType>(1)), // bad alpha < 0.
          static_cast<RealType>(1)), std::domain_error);

  BOOST_CHECK_THROW( 
       pdf(
          beta_distribution<RealType>(static_cast<RealType>(0), static_cast<RealType>(1)), // bad alpha == 0.
          static_cast<RealType>(1)), std::domain_error);

  BOOST_CHECK_THROW(
       pdf(
          beta_distribution<RealType>(static_cast<RealType>(1), static_cast<RealType>(0)), // bad beta == 0.
          static_cast<RealType>(1)), std::domain_error);

  BOOST_CHECK_THROW(
       pdf(
          beta_distribution<RealType>(static_cast<RealType>(1), static_cast<RealType>(-1)), // bad beta < 0.
          static_cast<RealType>(1)), std::domain_error);

  BOOST_CHECK_THROW(
       pdf(
          beta_distribution<RealType>(static_cast<RealType>(1), static_cast<RealType>(1)), // bad x < 0.
          static_cast<RealType>(-1)), std::domain_error);

  BOOST_CHECK_THROW(
       pdf(
          beta_distribution<RealType>(static_cast<RealType>(1), static_cast<RealType>(1)), // bad x > 1.
          static_cast<RealType>(999)), std::domain_error);

  // Some exact pdf values.

  BOOST_CHECK_EQUAL( // a = b = 1 is uniform distribution.
     pdf(beta_distribution<RealType>(static_cast<RealType>(1), static_cast<RealType>(1)),
     static_cast<RealType>(1)),  // x
     static_cast<RealType>(1));
  BOOST_CHECK_EQUAL(
     pdf(beta_distribution<RealType>(static_cast<RealType>(1), static_cast<RealType>(1)),
     static_cast<RealType>(0)),  // x
     static_cast<RealType>(1));
  BOOST_CHECK_EQUAL(
     pdf(beta_distribution<RealType>(static_cast<RealType>(1), static_cast<RealType>(1)),
     static_cast<RealType>(0.5)),  // x
     static_cast<RealType>(1));

  BOOST_CHECK_EQUAL(
     beta_distribution<RealType>(static_cast<RealType>(1), static_cast<RealType>(1)).alpha(),
     static_cast<RealType>(1) ); // 

  BOOST_CHECK_EQUAL(
     mean(beta_distribution<RealType>(static_cast<RealType>(1), static_cast<RealType>(1))),
     static_cast<RealType>(0.5) ); // Exact one half.

  BOOST_CHECK_CLOSE_FRACTION(
     pdf(beta_distribution<RealType>(static_cast<RealType>(2), static_cast<RealType>(2)),
     static_cast<RealType>(0.5)),  // x
     static_cast<RealType>(1.5), // Exactly 3/2
      tolerance);

  BOOST_CHECK_CLOSE_FRACTION(
     pdf(beta_distribution<RealType>(static_cast<RealType>(2), static_cast<RealType>(2)),
     static_cast<RealType>(0.5)),  // x
     static_cast<RealType>(1.5), // Exactly 3/2
      tolerance);

  // CDF
  BOOST_CHECK_CLOSE_FRACTION(
     cdf(beta_distribution<RealType>(static_cast<RealType>(2), static_cast<RealType>(2)),
     static_cast<RealType>(0.1)),  // x
     static_cast<RealType>(0.02800000000000000000000000000000000000000), // Seems exact.
     // http://functions.wolfram.com/webMathematica/FunctionEvaluation.jsp?name=BetaRegularized&ptype=0&z=0.1&a=2&b=2&digits=40
      tolerance);

  BOOST_CHECK_CLOSE_FRACTION(
     cdf(beta_distribution<RealType>(static_cast<RealType>(2), static_cast<RealType>(2)),
     static_cast<RealType>(0.0001)),  // x
     static_cast<RealType>(2.999800000000000000000000000000000000000e-8),
     // http://members.aol.com/iandjmsmith/BETAEX.HTM 2.9998000000004
     // http://functions.wolfram.com/webMathematica/FunctionEvaluation.jsp?name=BetaRegularized&ptype=0&z=0.0001&a=2&b=2&digits=40
      tolerance);


  BOOST_CHECK_CLOSE_FRACTION(
     pdf(beta_distribution<RealType>(static_cast<RealType>(2), static_cast<RealType>(2)),
     static_cast<RealType>(0.0001)),  // x
     static_cast<RealType>(0.0005999400000000004), // http://members.aol.com/iandjmsmith/BETAEX.HTM
      tolerance);


  BOOST_CHECK_CLOSE_FRACTION(
     cdf(beta_distribution<RealType>(static_cast<RealType>(2), static_cast<RealType>(2)),
     static_cast<RealType>(0.9999)),  // x
     static_cast<RealType>(0.999999970002), // http://members.aol.com/iandjmsmith/BETAEX.HTM
     // Wolfram 0.9999999700020000000000000000000000000000
      tolerance);

  BOOST_CHECK_CLOSE_FRACTION(
     cdf(beta_distribution<RealType>(static_cast<RealType>(0.5), static_cast<RealType>(2)),
     static_cast<RealType>(0.9)),  // x
     static_cast<RealType>(0.9961174629530394895796514664963063381217), 
     // Wolfram
      tolerance);

  BOOST_CHECK_CLOSE_FRACTION(
     cdf(beta_distribution<RealType>(static_cast<RealType>(0.5), static_cast<RealType>(0.5)),
     static_cast<RealType>(0.1)),  // x
     static_cast<RealType>(0.2048327646991334516491978475505189480977), 
     // Wolfram
      tolerance);

  BOOST_CHECK_CLOSE_FRACTION(
     cdf(beta_distribution<RealType>(static_cast<RealType>(0.5), static_cast<RealType>(0.5)),
     static_cast<RealType>(0.9)),  // x
     static_cast<RealType>(0.7951672353008665483508021524494810519023), 
     // Wolfram
      tolerance);

  BOOST_CHECK_CLOSE_FRACTION(
     quantile(beta_distribution<RealType>(static_cast<RealType>(0.5), static_cast<RealType>(0.5)),
     static_cast<RealType>(0.7951672353008665483508021524494810519023)),  // x
     static_cast<RealType>(0.9), 
     // Wolfram
     tolerance);

  BOOST_CHECK_CLOSE_FRACTION(
     cdf(beta_distribution<RealType>(static_cast<RealType>(0.5), static_cast<RealType>(0.5)),
     static_cast<RealType>(0.6)),  // x
     static_cast<RealType>(0.5640942168489749316118742861695149357858), 
     // Wolfram
      tolerance);

  BOOST_CHECK_CLOSE_FRACTION(
     quantile(beta_distribution<RealType>(static_cast<RealType>(0.5), static_cast<RealType>(0.5)),
     static_cast<RealType>(0.5640942168489749316118742861695149357858)),  // x
     static_cast<RealType>(0.6), 
     // Wolfram
      tolerance); 


  BOOST_CHECK_CLOSE_FRACTION(
     cdf(beta_distribution<RealType>(static_cast<RealType>(2), static_cast<RealType>(0.5)),
     static_cast<RealType>(0.6)),  // x
     static_cast<RealType>(0.1778078083562213736802876784474931812329), 
     // Wolfram
      tolerance);

  BOOST_CHECK_CLOSE_FRACTION(
     quantile(beta_distribution<RealType>(static_cast<RealType>(2), static_cast<RealType>(0.5)),
     static_cast<RealType>(0.1778078083562213736802876784474931812329)),  // x
     static_cast<RealType>(0.6), 
     // Wolfram
      tolerance); // gives 

  BOOST_CHECK_CLOSE_FRACTION(
     cdf(beta_distribution<RealType>(static_cast<RealType>(1), static_cast<RealType>(1)),
     static_cast<RealType>(0.1)),  // x 
     static_cast<RealType>(0.1),  // 0.1000000000000000000000000000000000000000
     // Wolfram
      tolerance);

  BOOST_CHECK_CLOSE_FRACTION(
     quantile(beta_distribution<RealType>(static_cast<RealType>(1), static_cast<RealType>(1)),
     static_cast<RealType>(0.1)),  // x 
     static_cast<RealType>(0.1),  // 0.1000000000000000000000000000000000000000
     // Wolfram
      tolerance);

  BOOST_CHECK_CLOSE_FRACTION(
     cdf(complement(beta_distribution<RealType>(static_cast<RealType>(0.5), static_cast<RealType>(0.5)),
     static_cast<RealType>(0.1))),  // complement of x
     static_cast<RealType>(0.7951672353008665483508021524494810519023), 
     // Wolfram
      tolerance);

    BOOST_CHECK_CLOSE_FRACTION(
     quantile(beta_distribution<RealType>(static_cast<RealType>(2), static_cast<RealType>(2)),
     static_cast<RealType>(0.0280000000000000000000000000000000000)),  // x
     static_cast<RealType>(0.1), 
     // Wolfram
      tolerance);


  BOOST_CHECK_CLOSE_FRACTION(
     cdf(complement(beta_distribution<RealType>(static_cast<RealType>(2), static_cast<RealType>(2)),
     static_cast<RealType>(0.1))),  // x
     static_cast<RealType>(0.9720000000000000000000000000000000000000), // Exact.
     // Wolfram
      tolerance);

  BOOST_CHECK_CLOSE_FRACTION(
     pdf(beta_distribution<RealType>(static_cast<RealType>(2), static_cast<RealType>(2)),
     static_cast<RealType>(0.9999)),  // x
     static_cast<RealType>(0.0005999399999999344), // http://members.aol.com/iandjmsmith/BETAEX.HTM
      tolerance*10); // Note less accurate.

  //void test_spot(
  //   RealType a,    // alpha a
  //   RealType b,    // beta b
  //   RealType x,    // Probability 
  //   RealType P,    // CDF of beta(a, b)
  //   RealType Q,    // Complement of CDF
  //   RealType tol)  // Test tolerance.

   // These test quantiles and complements, and parameter estimates as well.
  // Spot values using, for example:
  // http://functions.wolfram.com/webMathematica/FunctionEvaluation.jsp?name=BetaRegularized&ptype=0&z=0.1&a=0.5&b=3&digits=40

  test_spot(
     static_cast<RealType>(1),   // alpha a
     static_cast<RealType>(1),   // beta b
     static_cast<RealType>(0.1), // Probability  p
     static_cast<RealType>(0.1), // Probability of result (CDF of beta), P
     static_cast<RealType>(0.9),  // Complement of CDF Q = 1 - P
     tolerance); // Test tolerance.
  test_spot(
     static_cast<RealType>(2),   // alpha a
     static_cast<RealType>(2),   // beta b
     static_cast<RealType>(0.1), // Probability  p
     static_cast<RealType>(0.0280000000000000000000000000000000000), // Probability of result (CDF of beta), P
     static_cast<RealType>(1 - 0.0280000000000000000000000000000000000),  // Complement of CDF Q = 1 - P
     tolerance); // Test tolerance.


  test_spot(
     static_cast<RealType>(2),   // alpha a
     static_cast<RealType>(2),   // beta b
     static_cast<RealType>(0.5), // Probability  p
     static_cast<RealType>(0.5), // Probability of result (CDF of beta), P
     static_cast<RealType>(0.5),  // Complement of CDF Q = 1 - P
     tolerance); // Test tolerance.

  test_spot(
     static_cast<RealType>(2),   // alpha a
     static_cast<RealType>(2),   // beta b
     static_cast<RealType>(0.9), // Probability  p
     static_cast<RealType>(0.972000000000000), // Probability of result (CDF of beta), P
     static_cast<RealType>(1-0.972000000000000),  // Complement of CDF Q = 1 - P
     tolerance); // Test tolerance.

  test_spot(
     static_cast<RealType>(2),   // alpha a
     static_cast<RealType>(2),   // beta b
     static_cast<RealType>(0.01), // Probability  p
     static_cast<RealType>(0.0002980000000000000000000000000000000000000), // Probability of result (CDF of beta), P
     static_cast<RealType>(1-0.0002980000000000000000000000000000000000000),  // Complement of CDF Q = 1 - P
     tolerance); // Test tolerance.

  test_spot(
     static_cast<RealType>(2),   // alpha a
     static_cast<RealType>(2),   // beta b
     static_cast<RealType>(0.001), // Probability  p
     static_cast<RealType>(2.998000000000000000000000000000000000000E-6), // Probability of result (CDF of beta), P
     static_cast<RealType>(1-2.998000000000000000000000000000000000000E-6),  // Complement of CDF Q = 1 - P
     tolerance); // Test tolerance.

  test_spot(
     static_cast<RealType>(2),   // alpha a
     static_cast<RealType>(2),   // beta b
     static_cast<RealType>(0.0001), // Probability  p
     static_cast<RealType>(2.999800000000000000000000000000000000000E-8), // Probability of result (CDF of beta), P
     static_cast<RealType>(1-2.999800000000000000000000000000000000000E-8),  // Complement of CDF Q = 1 - P
     tolerance); // Test tolerance.

  test_spot(
     static_cast<RealType>(2),   // alpha a
     static_cast<RealType>(2),   // beta b
     static_cast<RealType>(0.99), // Probability  p
     static_cast<RealType>(0.9997020000000000000000000000000000000000), // Probability of result (CDF of beta), P
     static_cast<RealType>(1-0.9997020000000000000000000000000000000000),  // Complement of CDF Q = 1 - P
     tolerance); // Test tolerance.

  test_spot(
     static_cast<RealType>(0.5),   // alpha a
     static_cast<RealType>(2),   // beta b
     static_cast<RealType>(0.5), // Probability  p
     static_cast<RealType>(0.8838834764831844055010554526310612991060), // Probability of result (CDF of beta), P
     static_cast<RealType>(1-0.8838834764831844055010554526310612991060),  // Complement of CDF Q = 1 - P
     tolerance); // Test tolerance.

  test_spot(
     static_cast<RealType>(0.5),   // alpha a
     static_cast<RealType>(3.),   // beta b
     static_cast<RealType>(0.7), // Probability  p
     static_cast<RealType>(0.9903963064097119299191611355232156905687), // Probability of result (CDF of beta), P
     static_cast<RealType>(1-0.9903963064097119299191611355232156905687),  // Complement of CDF Q = 1 - P
     tolerance); // Test tolerance.

  test_spot(
     static_cast<RealType>(0.5),   // alpha a
     static_cast<RealType>(3.),   // beta b
     static_cast<RealType>(0.1), // Probability  p
     static_cast<RealType>(0.5545844446520295253493059553548880128511), // Probability of result (CDF of beta), P
     static_cast<RealType>(1-0.5545844446520295253493059553548880128511),  // Complement of CDF Q = 1 - P
     tolerance); // Test tolerance.

} // template <class RealType>void test_spots(RealType)

int test_main(int, char* [])
{
#ifdef BOOST_MATH_THROW_ON_DOMAIN_ERROR
  cout << "BOOST_MATH_THROW_ON_DOMAIN_ERROR" << " is defined to throw on domain error." << endl;
#else
  cout << "BOOST_MATH_THROW_ON_DOMAIN_ERROR" << " is NOT defined, so NO throw on domain error." << endl;
#endif


	// Check that can generate beta distribution using one convenience methods:
	beta_distribution<> mybeta11(1., 1.); // Using default RealType double.
  // but that
	//boost::math::beta mybeta1(1., 1.); // Using typedef fails.
  // error C2039: 'beta' : is not a member of 'boost::math'

  // Basic sanity-check spot values.

  // Some simple checks using double only.
  BOOST_CHECK_EQUAL(mybeta11.alpha(), 1); // 
  BOOST_CHECK_EQUAL(mybeta11.beta(), 1);
  BOOST_CHECK_EQUAL(mean(mybeta11), 0.5); // 1 / (1 + 1) = 1/2 exactly
  BOOST_CHECK_THROW(mode(mybeta11), std::domain_error);
	beta_distribution<> mybeta22(2., 2.); // pdf is dome shape.
  BOOST_CHECK_EQUAL(mode(mybeta22), 0.5); // 2-1 / (2+2-2) = 1/2 exactly.
	beta_distribution<> mybetaH2(0.5, 2.); // 
	beta_distribution<> mybetaH3(0.5, 3.); // 


  // Check a few values using double.
  BOOST_CHECK_EQUAL(pdf(mybeta11, 1), 1); // is uniform unity over 0 to 1,
  BOOST_CHECK_EQUAL(pdf(mybeta11, 0), 1); // including zero and unity.
  BOOST_CHECK_EQUAL(pdf(mybeta11, 0.5), 1);
  BOOST_CHECK_EQUAL(pdf(mybeta11, 0.0001), 1);
  BOOST_CHECK_EQUAL(pdf(mybeta11, 0.9999), 1);
  BOOST_CHECK_CLOSE_FRACTION(cdf(mybeta11, 0.1), 0.1, std::numeric_limits<double>::epsilon());
  BOOST_CHECK_CLOSE_FRACTION(cdf(mybeta11, 0.5), 0.5, std::numeric_limits<double>::epsilon());
  BOOST_CHECK_CLOSE_FRACTION(cdf(mybeta11, 0.9), 0.9, std::numeric_limits<double>::epsilon());
  BOOST_CHECK_EQUAL(cdf(mybeta11, 1), 1.); // Exact unity expected.

  double tol = std::numeric_limits<double>::epsilon() * 10;
  BOOST_CHECK_EQUAL(pdf(mybeta22, 1), 0); // is dome shape.
  BOOST_CHECK_EQUAL(pdf(mybeta22, 0), 0);
  BOOST_CHECK_CLOSE_FRACTION(pdf(mybeta22, 0.5), 1.5, tol); // top of dome, expect exactly 3/2.
  BOOST_CHECK_CLOSE_FRACTION(pdf(mybeta22, 0.0001), 5.9994000000000E-4, tol);
  BOOST_CHECK_CLOSE_FRACTION(pdf(mybeta22, 0.9999), 5.9994000000000E-4, tol*50);

  BOOST_CHECK_EQUAL(cdf(mybeta22, 0.), 0); // cdf is a curved line from 0 to 1.
  BOOST_CHECK_CLOSE_FRACTION(cdf(mybeta22, 0.1), 0.028000000000000, tol);
  BOOST_CHECK_CLOSE_FRACTION(cdf(mybeta22, 0.5), 0.5, tol);
  BOOST_CHECK_CLOSE_FRACTION(cdf(mybeta22, 0.9), 0.972000000000000, tol);
  BOOST_CHECK_CLOSE_FRACTION(cdf(mybeta22, 0.0001), 2.999800000000000000000000000000000000000E-8, tol);
  BOOST_CHECK_CLOSE_FRACTION(cdf(mybeta22, 0.001), 2.998000000000000000000000000000000000000E-6, tol);
  BOOST_CHECK_CLOSE_FRACTION(cdf(mybeta22, 0.01), 0.0002980000000000000000000000000000000000000, tol);
  BOOST_CHECK_CLOSE_FRACTION(cdf(mybeta22, 0.1), 0.02800000000000000000000000000000000000000, tol); // exact
  BOOST_CHECK_CLOSE_FRACTION(cdf(mybeta22, 0.99), 0.9997020000000000000000000000000000000000, tol);

  BOOST_CHECK_EQUAL(cdf(mybeta22, 1), 1.); // Exact unity expected.

  // Complement

  BOOST_CHECK_CLOSE_FRACTION(cdf(complement(mybeta22, 0.9)), 0.028000000000000, tol);

  // quantile.
  BOOST_CHECK_CLOSE_FRACTION(quantile(mybeta22, 0.028), 0.1, tol); 
  BOOST_CHECK_CLOSE_FRACTION(quantile(complement(mybeta22, 1 - 0.028)), 0.1, tol);
  BOOST_CHECK_EQUAL(kurtosis(mybeta11), 3+ kurtosis_excess(mybeta11)); // Check kurtosis_excess = kurtosis - 3;
  BOOST_CHECK_CLOSE_FRACTION(variance(mybeta22), 0.05, tol);
  BOOST_CHECK_CLOSE_FRACTION(mode(mybeta22), 0.5, tol);
  BOOST_CHECK_CLOSE_FRACTION(mean(mybeta22), 0.5, tol);


   //cout << beta_distribution<double>::estimate_alpha(mean(mybeta22), variance(mybeta22)) << endl; // 2
  //cout << beta_distribution<double>::estimate_beta(mean(mybeta22), variance(mybeta22)) << endl; // 2
  BOOST_CHECK_EQUAL(beta_distribution<double>::estimate_alpha(mean(mybeta22), variance(mybeta22)), mybeta22.alpha()); // mean, variance, probability. 
  BOOST_CHECK_EQUAL(beta_distribution<double>::estimate_beta(mean(mybeta22), variance(mybeta22)), mybeta22.beta());// mean, variance, probability. 

 // BOOST_CHECK_CLOSE_FRACTION(ibeta_inva(mybeta22.beta(), 0.8, cdf(mybeta22, 0.8)), mybeta22.alpha(), tol);
  //BOOST_CHECK_CLOSE_FRACTION(ibeta_inva(mybeta22.beta(), 0.8, cdf(mybeta22, 0.8)), mybeta22.alpha(), tol);
  BOOST_CHECK_CLOSE_FRACTION(mybeta22.estimate_alpha(mybeta22.beta(), 0.8, cdf(mybeta22, 0.8)), mybeta22.alpha(), tol);
  BOOST_CHECK_CLOSE_FRACTION(mybeta22.estimate_beta(mybeta22.alpha(), 0.8, cdf(mybeta22, 0.8)), mybeta22.beta(), tol);


	beta_distribution<real_concept> rcbeta22(2, 2); // Using RealType real_concept.
  cout << "numeric_limits<real_concept>::is_specialized " << numeric_limits<real_concept>::is_specialized << endl;
  cout << "numeric_limits<real_concept>::digits " << numeric_limits<real_concept>::digits << endl;
  cout << "numeric_limits<real_concept>::digits10 " << numeric_limits<real_concept>::digits10 << endl;
  cout << "numeric_limits<real_concept>::epsilon " << numeric_limits<real_concept>::epsilon() << endl;

  // Tests for improvements to Boost.test display of errors.
  //BOOST_CHECK_CLOSE_FRACTION(rcbeta22.alpha(), static_cast<real_concept>(2 + 3 * std::numeric_limits<double>::epsilon()), 0);
  //BOOST_CHECK_CLOSE_FRACTION(cdf(rcbeta22, 0.1), static_cast<real_concept>(0.028000000000009), 0);

  // (Parameter value, arbitrarily zero, only communicates the floating point type).
  test_spots(0.0F); // Test float.
  test_spots(0.0); // Test double.
  test_spots(0.0L); // Test long double.
#if !BOOST_WORKAROUND(__BORLANDC__, BOOST_TESTED_AT(0x582))
  test_spots(boost::math::concepts::real_concept(0.)); // Test real concept.
#endif

  return 0;
} // int test_main(int, char* [])

/*

Output is:

------ Build started: Project: test_beta_dist, Configuration: Debug Win32 ------
Compiling...
test_beta_dist.cpp
Linking...
Autorun "i:\boost-06-05-03-1300\libs\math\test\Math_test\debug\test_beta_dist.exe"
Running 1 test case...
BOOST_MATH_THROW_ON_DOMAIN_ERROR is defined to throw on domain error.
numeric_limits<real_concept>::is_specialized 0
numeric_limits<real_concept>::digits 0
numeric_limits<real_concept>::digits10 0
numeric_limits<real_concept>::epsilon 0
Boost::math::tools::epsilon = 1.19209e-007
std::numeric_limits::epsilon = 1.19209e-007
epsilon = 1.19209e-007, Tolerance = 0.0119209%.
Boost::math::tools::epsilon = 2.22045e-016
std::numeric_limits::epsilon = 2.22045e-016
epsilon = 2.22045e-016, Tolerance = 2.22045e-011%.
Boost::math::tools::epsilon = 2.22045e-016
std::numeric_limits::epsilon = 2.22045e-016
epsilon = 2.22045e-016, Tolerance = 2.22045e-011%.
Boost::math::tools::epsilon = 2.22045e-016
std::numeric_limits::epsilon = 0
epsilon = 2.22045e-016, Tolerance = 2.22045e-011%.
*** No errors detected
Build Time 0:07
Build log was saved at "file://i:\boost-06-05-03-1300\libs\math\test\Math_test\test_beta_dist\Debug\BuildLog.htm"
test_beta_dist - 0 error(s), 0 warning(s)
========== Build: 1 succeeded, 0 failed, 0 up-to-date, 0 skipped ==========


*/

