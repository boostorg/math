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

// http://www.ausvet.com.au/pprev/content.php?page=PPscript
// mode 0.75 	5/95% 0.9 	alpha 7.39 	beta 3.13

// http://www.nuhertz.com/statmat/distributions.html#Beta
// Pretty graphs and explanations for most distributions.

// http://functions.wolfram.com/webMathematica/FunctionEvaluation.jsp
// provided 40 decimal digits accuracy incomplete beta aka beta regularized == cdf



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
      if(k != 0)
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
      // Estimate alpha:
      BOOST_CHECK_CLOSE_FRACTION(
         beta_distribution<RealType>::estimate_alpha(N, x, Q),
         p, tol);
      BOOST_CHECK_CLOSE_FRACTION(
         beta_distribution<RealType>::estimate_beta(N, x, P),
         p, tol);

      // Estimate sample alpha and beta:
      BOOST_CHECK_CLOSE_FRACTION(
         beta_distribution<RealType>::estimate_alpha(
            k, p, P),
         N, tol);
      BOOST_CHECK_CLOSE_FRACTION(
         beta_distribution<RealType>::estimate_alpha(
            boost::math::complement(k, p, Q)),
         N, tol);
      BOOST_CHECK_CLOSE_FRACTION(
         beta_distribution<RealType>::estimate_beta(
            k, p, P),
         N, tol);
      BOOST_CHECK_CLOSE_FRACTION(
         beta_distribution<RealType>::estimate_beta(
            boost::math::complement(k, p, Q)),
         N, tol);
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
      static_cast<RealType>(std::numeric_limits<double>::epsilon()));
   tolerance *= 1000; // Note: NO * 100 because is fraction, NOT %.

  cout << "Tolerance = " << tolerance * 100 << "%." << endl;

  //RealType teneps = boost::math::tools::epsilon<RealType>() * 10;

  // Sources of spot test values:

  // MathCAD defines dbeta(x, s1, s2) pdf, s1 == alpha, s2 = beta, x = x in Wolfram
  // pbeta(x, s1, s2) cdf and qbeta(x, s1, s2) inverse of cdf
  // returns pr(X ,= x) when random variable X
  // has the beta distribution with parameters s1)alpha) and s2(beta).

  // s1 > 0 and s2 >0 and 0 < x < 1 but allows x == 0! and x == 1!)
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
     tolerance); // gives 0.5 ??????????????

  BOOST_CHECK_CLOSE_FRACTION(
     cdf(beta_distribution<RealType>(static_cast<RealType>(0.5), static_cast<RealType>(0.5)),
     static_cast<RealType>(0.6)),  // x
     static_cast<RealType>(0.5640942168489749316118742861695149357858), 
     // Wolfram
      tolerance);

  BOOST_CHECK_CLOSE_FRACTION(
     quantile(beta_distribution<RealType>(static_cast<RealType>(2), static_cast<RealType>(0.5)),
     static_cast<RealType>(0.5640942168489749316118742861695149357858)),  // x
     static_cast<RealType>(0.6), 
     // Wolfram
      tolerance); // gives 


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

} // template <class RealType>void test_spots(RealType)

int test_main(int, char* [])
{
	// Check that can generate beta distribution using one convenience methods:
	beta_distribution<> mybeta11(1., 1.); // Using default RealType double.
  // but that
	//boost::math::beta mybeta1(1., 1.); // Using typedef fails.
  // error C2039: 'beta' : is not a member of 'boost::math'

  // Some simple checks using double only.
  BOOST_CHECK_EQUAL(mybeta11.alpha(), 1); // 
  BOOST_CHECK_EQUAL(mybeta11.beta(), 1);
  BOOST_CHECK_EQUAL(mean(mybeta11), 0.5); // 1 / (1 + 1) = 1/2 exactly
  BOOST_CHECK_THROW(mode(mybeta11), std::domain_error);
	beta_distribution<> mybeta22(2., 2.); // pdf is dome shape.
  BOOST_CHECK_EQUAL(mode(mybeta22), 0.5); // 2-1 / (2+2-2) = 1/2 exactly.
	beta_distribution<> mybetaH2(0.5, 2.); // 

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
  //BOOST_CHECK_CLOSE_FRACTION(cdf(mybeta22, 0.0001), 4.999666666666666666666666666666666666667E-9, tol);
  BOOST_CHECK_CLOSE_FRACTION(cdf(mybeta22, 0.1), 0.02800000000000000000000000000000000000000, tol); // exact
  //BOOST_CHECK_CLOSE_FRACTION(cdf(mybeta22, 0.1), 0.004666666666666666666666666666666666666667, tol);

  // 4.999666666666666666666666666666666666667E-9  Wolfram

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
  BOOST_CHECK_EQUAL(beta_distribution<double>::estimate_alpha(mean(mybeta22), variance(mybeta22), 0.5), mybeta22.alpha()); // mean, variance, probability. 
  BOOST_CHECK_EQUAL(beta_distribution<double>::estimate_beta(mean(mybeta22), variance(mybeta22), 0.5), mybeta22.beta());// mean, variance, probability. 

  

  // Basic sanity-check spot values.
#ifdef BOOST_MATH_THROW_ON_DOMAIN_ERROR
  cout << "BOOST_MATH_THROW_ON_DOMAIN_ERROR" << " is defined to throw on domain error." << endl;
#else
  cout << "BOOST_MATH_THROW_ON_DOMAIN_ERROR" << " is NOT defined, so NO throw on domain error." << endl;
#endif

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
*** No errors detected
Build Time 0:05
Build log was saved at "file://i:\boost-06-05-03-1300\libs\math\test\Math_test\test_beta_dist\Debug\BuildLog.htm"
test_beta_dist - 0 error(s), 0 warning(s)
========== Build: 1 succeeded, 0 failed, 0 up-to-date, 0 skipped ==========


*/

/*

Park

 // These test quantiles and complements as well.
  test_spot(
     static_cast<RealType>(500),                     // Sample size, N
     static_cast<RealType>(30),                      // Number of successes, k
     static_cast<RealType>(0.05),                    // Probability of success, p
     static_cast<RealType>(0.869147702104609),       // Probability of result (CDF), P
     static_cast<RealType>(1 - 0.869147702104609),   // Q = 1 - P
     tolerance);

  test_spot(
     static_cast<RealType>(500),                     // Sample size, N
     static_cast<RealType>(250),                     // Number of successes, k
     static_cast<RealType>(0.05),                    // Probability of success, p
     static_cast<RealType>(1),                       // Probability of result (CDF), P
     static_cast<RealType>(0),   // Q = 1 - P
     tolerance);

  test_spot(
     static_cast<RealType>(500),                     // Sample size, N
     static_cast<RealType>(470),                     // Number of successes, k
     static_cast<RealType>(0.95),                    // Probability of success, p
     static_cast<RealType>(0.176470742656766),       // Probability of result (CDF), P
     static_cast<RealType>(1 - 0.176470742656766),   // Q = 1 - P
     tolerance * 10);                                // Note higher tolerance on this test!

  test_spot(
     static_cast<RealType>(500),                       // Sample size, N
     static_cast<RealType>(400),                       // Number of successes, k
     static_cast<RealType>(0.05),                      // Probability of success, p
     static_cast<RealType>(1),                         // Probability of result (CDF), P
     static_cast<RealType>(0),                         // Q = 1 - P
     tolerance);

  test_spot(
     static_cast<RealType>(500),                       // Sample size, N
     static_cast<RealType>(400),                       // Number of successes, k
     static_cast<RealType>(0.9),                       // Probability of success, p
     static_cast<RealType>(1.80180425681923E-11),      // Probability of result (CDF), P
     static_cast<RealType>(1 - 1.80180425681923E-11),  // Q = 1 - P
     tolerance);

  test_spot(
     static_cast<RealType>(500),                       // Sample size, N
     static_cast<RealType>(5),                         // Number of successes, k
     static_cast<RealType>(0.05),                      // Probability of success, p
     static_cast<RealType>(9.181808267643E-7),         // Probability of result (CDF), P
     static_cast<RealType>(1 - 9.181808267643E-7),     // Q = 1 - P
     tolerance);

  test_spot(
     static_cast<RealType>(2),                       // Sample size, N
     static_cast<RealType>(1),                       // Number of successes, k
     static_cast<RealType>(0.5),                     // Probability of success, p
     static_cast<RealType>(0.75),                    // Probability of result (CDF), P
     static_cast<RealType>(0.25),                    // Q = 1 - P
     tolerance);

  test_spot(
     static_cast<RealType>(8),                       // Sample size, N
     static_cast<RealType>(3),                       // Number of successes, k
     static_cast<RealType>(0.25),                    // Probability of success, p
     static_cast<RealType>(0.8861846923828125),      // Probability of result (CDF), P
     static_cast<RealType>(1 - 0.8861846923828125),  // Q = 1 - P
     tolerance);

  test_spot(
     static_cast<RealType>(8),                       // Sample size, N
     static_cast<RealType>(0),                       // Number of successes, k
     static_cast<RealType>(0.25),                    // Probability of success, p
     static_cast<RealType>(0.1001129150390625),      // Probability of result (CDF), P
     static_cast<RealType>(1 - 0.1001129150390625),  // Q = 1 - P
     tolerance);

  test_spot(
     static_cast<RealType>(8),                       // Sample size, N
     static_cast<RealType>(1),                       // Number of successes, k
     static_cast<RealType>(0.25),                    // Probability of success, p
     static_cast<RealType>(0.36708068847656244),     // Probability of result (CDF), P
     static_cast<RealType>(1 - 0.36708068847656244), // Q = 1 - P
     tolerance);

  test_spot(
     static_cast<RealType>(8),                       // Sample size, N
     static_cast<RealType>(4),                       // Number of successes, k
     static_cast<RealType>(0.25),                    // Probability of success, p
     static_cast<RealType>(0.9727020263671875),      // Probability of result (CDF), P
     static_cast<RealType>(1 - 0.9727020263671875),  // Q = 1 - P
     tolerance);

  test_spot(
     static_cast<RealType>(8),                       // Sample size, N
     static_cast<RealType>(7),                       // Number of successes, k
     static_cast<RealType>(0.25),                    // Probability of success, p
     static_cast<RealType>(0.9999847412109375),      // Probability of result (CDF), P
     static_cast<RealType>(1 - 0.9999847412109375),  // Q = 1 - P
     tolerance);

  // Tests on PDF follow:
  BOOST_CHECK_CLOSE_FRACTION(
     pdf(beta_distribution<RealType>(static_cast<RealType>(20), static_cast<RealType>(0.75)),
     static_cast<RealType>(10)),  // k.
     static_cast<RealType>(0.00992227527967770583927631378173), // 0.00992227527967770583927631378173
     tolerance);

  BOOST_CHECK_CLOSE_FRACTION(
    pdf(beta_distribution<RealType>(static_cast<RealType>(20), static_cast<RealType>(0.5)),
    static_cast<RealType>(10)),  // k.
    static_cast<RealType>(0.17619705200195312500000000000000000000), // get k=10 0.049611376398388612 p = 0.25
    tolerance);

     BOOST_CHECK_CLOSE_FRACTION(
    pdf(beta_distribution<RealType>(static_cast<RealType>(20), static_cast<RealType>(0.25)),
    static_cast<RealType>(10)),  // k.
    static_cast<RealType>(0.00992227527967770583927631378173), // k=10  p = 0.25
    tolerance);

    BOOST_CHECK_CLOSE_FRACTION( // k = 0 use different formula - only exp so more accurate.
    pdf(beta_distribution<RealType>(static_cast<RealType>(20), static_cast<RealType>(0.25)),
    static_cast<RealType>(0)),  // k.
    static_cast<RealType>(0.00317121193893399322405457496643), // k=0  p = 0.25
    tolerance);

    BOOST_CHECK_CLOSE_FRACTION( // k = 20 use different formula - only exp so more accurate.
    pdf(beta_distribution<RealType>(static_cast<RealType>(20), static_cast<RealType>(0.25)),
    static_cast<RealType>(20)),  // k == n.
    static_cast<RealType>(0.00000000000090949470177292823791), // k=20  p = 0.25
    tolerance);

    BOOST_CHECK_CLOSE_FRACTION( // k = 1.
    pdf(beta_distribution<RealType>(static_cast<RealType>(20), static_cast<RealType>(0.25)),
    static_cast<RealType>(1)),  // k.
    static_cast<RealType>(0.02114141292622662149369716644287), // k=1  p = 0.25
    tolerance);

    // Some exact (probably) values.
    BOOST_CHECK_CLOSE_FRACTION(
    pdf(beta_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(0.25)),
    static_cast<RealType>(0)),  // k.
    static_cast<RealType>(0.10011291503906250000000000000000), // k=0  p = 0.25
    tolerance);

    BOOST_CHECK_CLOSE_FRACTION( // k = 1.
    pdf(beta_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(0.25)),
    static_cast<RealType>(1)),  // k.
    static_cast<RealType>(0.26696777343750000000000000000000), // k=1  p = 0.25
    tolerance);

    BOOST_CHECK_CLOSE_FRACTION( // k = 2.
    pdf(beta_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(0.25)),
    static_cast<RealType>(2)),  // k.
    static_cast<RealType>(0.31146240234375000000000000000000), // k=2  p = 0.25
    tolerance);

    BOOST_CHECK_CLOSE_FRACTION( // k = 3.
    pdf(beta_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(0.25)),
    static_cast<RealType>(3)),  // k.
    static_cast<RealType>(0.20764160156250000000000000000000), // k=3  p = 0.25
    tolerance);

    BOOST_CHECK_CLOSE_FRACTION( // k = 7.
    pdf(beta_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(0.25)),
    static_cast<RealType>(7)),  // k.
    static_cast<RealType>(0.00036621093750000000000000000000), // k=7  p = 0.25
    tolerance);

    BOOST_CHECK_CLOSE_FRACTION( // k = 8.
    pdf(beta_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(0.25)),
    static_cast<RealType>(8)),  // k = n.
    static_cast<RealType>(0.00001525878906250000000000000000), // k=8  p = 0.25
    tolerance);

    RealType tol2 = boost::math::tools::epsilon<RealType>() * 5 * 100;  // 5 eps as a persent
    beta_distribution<RealType> dist(static_cast<RealType>(8), static_cast<RealType>(0.25));
    RealType x = static_cast<RealType>(0.125);
    using namespace std; // ADL of std names.
    // mean:
    BOOST_CHECK_CLOSE_FRACTION(
       mean(dist)
       , static_cast<RealType>(8 * 0.25), tol2);
    // variance:
    BOOST_CHECK_CLOSE_FRACTION(
       variance(dist)
       , static_cast<RealType>(8 * 0.25 * 0.75), tol2);
    // std deviation:
    BOOST_CHECK_CLOSE_FRACTION(
       standard_deviation(dist)
       , static_cast<RealType>(sqrt(8 * 0.25L * 0.75L)), tol2);
    // hazard:
    BOOST_CHECK_CLOSE_FRACTION(
       hazard(dist, x)
       , pdf(dist, x) / cdf(complement(dist, x)), tol2);
    // cumulative hazard:
    BOOST_CHECK_CLOSE_FRACTION(
       chf(dist, x)
       , -log(cdf(complement(dist, x))), tol2);
    // coefficient_of_variation:
    BOOST_CHECK_CLOSE_FRACTION(
       coefficient_of_variation(dist)
       , standard_deviation(dist) / mean(dist), tol2);
    // mode:
    BOOST_CHECK_CLOSE_FRACTION(
       mode(dist)
       , static_cast<RealType>(std::floor(9 * 0.25)), tol2);
    // skewness:
    BOOST_CHECK_CLOSE_FRACTION(
       skewness(dist)
       , static_cast<RealType>(0.40824829046386301636621401245098L), tol2);
    // kurtosis:
    BOOST_CHECK_CLOSE_FRACTION(
       kurtosis(dist)
       , static_cast<RealType>(2.9166666666666666666666666666667L), tol2);
    // kurtosis excess:
    BOOST_CHECK_CLOSE_FRACTION(
       kurtosis_excess(dist)
       , static_cast<RealType>(-0.083333333333333333333333333333333L), tol2);

    // special cases for PDF:
    BOOST_CHECK_EQUAL(
       pdf(
          beta_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(0)),
          static_cast<RealType>(0)), static_cast<RealType>(1)
       );
    BOOST_CHECK_EQUAL(
       pdf(
          beta_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(0)),
          static_cast<RealType>(0.0001)), static_cast<RealType>(0)
       );
    BOOST_CHECK_EQUAL(
       pdf(
          beta_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(1)),
          static_cast<RealType>(0.001)), static_cast<RealType>(0)
       );
    BOOST_CHECK_EQUAL(
       pdf(
          beta_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(1)),
          static_cast<RealType>(8)), static_cast<RealType>(1)
       );
    BOOST_CHECK_EQUAL(
       pdf(
          beta_distribution<RealType>(static_cast<RealType>(0), static_cast<RealType>(0.25)),
          static_cast<RealType>(0)), static_cast<RealType>(1)
       );
    BOOST_CHECK_THROW(
       pdf(
          beta_distribution<RealType>(static_cast<RealType>(-1), static_cast<RealType>(0.25)),
          static_cast<RealType>(0)), std::domain_error
       );
    BOOST_CHECK_THROW(
       pdf(
          beta_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(-0.25)),
          static_cast<RealType>(0)), std::domain_error
       );
    BOOST_CHECK_THROW(
       pdf(
          beta_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(1.25)),
          static_cast<RealType>(0)), std::domain_error
       );
    BOOST_CHECK_THROW(
       pdf(
          beta_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(0.25)),
          static_cast<RealType>(-1)), std::domain_error
       );
    BOOST_CHECK_THROW(
       pdf(
          beta_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(0.25)),
          static_cast<RealType>(9)), std::domain_error
       );
    BOOST_CHECK_THROW(
       cdf(
          beta_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(0.25)),
          static_cast<RealType>(-1)), std::domain_error
       );
    BOOST_CHECK_THROW(
       cdf(
          beta_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(0.25)),
          static_cast<RealType>(9)), std::domain_error
       );
    BOOST_CHECK_THROW(
       cdf(
          beta_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(-0.25)),
          static_cast<RealType>(0)), std::domain_error
       );
    BOOST_CHECK_THROW(
       cdf(
          beta_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(1.25)),
          static_cast<RealType>(0)), std::domain_error
       );
    BOOST_CHECK_THROW(
       quantile(
          beta_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(-0.25)),
          static_cast<RealType>(0)), std::domain_error
       );
    BOOST_CHECK_THROW(
       quantile(
          beta_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(1.25)),
          static_cast<RealType>(0)), std::domain_error
       );
    BOOST_CHECK_EQUAL(
       cdf(
          beta_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(0.25)),
          static_cast<RealType>(8)), static_cast<RealType>(1)
       );
    BOOST_CHECK_EQUAL(
       cdf(
          beta_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(0)),
          static_cast<RealType>(7)), static_cast<RealType>(1)
       );
    BOOST_CHECK_EQUAL(
       cdf(
          beta_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(1)),
          static_cast<RealType>(7)), static_cast<RealType>(0)
       );


       */

