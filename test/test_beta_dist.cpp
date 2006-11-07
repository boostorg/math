// test_beta_dist.cpp

// Copyright John Maddock 2006.
// Copyright  Paul A. Bristow 2006.

// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

// Basic sanity test for the beta Distribution.

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
#include <boost/test/floating_point_comparison.hpp> // for BOOST_CHECK_CLOSE

#include <iostream>
using std::cout;
using std::endl;
#include <limits>
using std::numeric_limits;

template <class RealType>
void test_spot(
     RealType N,    // Number of trials
     RealType k,    // Number of successes
     RealType p,    // Probability of success
     RealType P,    // CDF
     RealType Q,    // Complement of CDF
     RealType tol)  // Test tolerance
{
   boost::math::beta_distribution<RealType> bn(N, p);
   BOOST_CHECK_CLOSE(
      cdf(bn, k), P, tol);
   if((P < 0.99) && (Q < 0.99))
   {
      //
      // We can only check this if P is not too close to 1,
      // so that we can guarentee Q is free of error:
      //
      BOOST_CHECK_CLOSE(
         cdf(complement(bn, k)), Q, tol);
      if(k != 0)
      {
         BOOST_CHECK_CLOSE(
            quantile(bn, P), k, tol);
      }
      else
      {
         // Just check quantile is very small:
         if((std::numeric_limits<RealType>::max_exponent <= std::numeric_limits<double>::max_exponent) && (boost::is_floating_point<RealType>::value))
         {
            // Limit where this is checked: if exponent range is very large we may
            // run out of iterations in our root finding algorithm.
            BOOST_CHECK(quantile(bn, P) < boost::math::tools::epsilon<RealType>() * 10);
         }
      }
      if(k != 0)
      {
         BOOST_CHECK_CLOSE(
            quantile(complement(bn, Q)), k, tol);
      }
      else
      {
         // Just check quantile is very small:
         if((std::numeric_limits<RealType>::max_exponent <= std::numeric_limits<double>::max_exponent) && (boost::is_floating_point<RealType>::value))
         {
            // Limit where this is checked: if exponent range is very large we may
            // run out of iterations in our root finding algorithm.
            BOOST_CHECK(quantile(complement(bn, Q)) < boost::math::tools::epsilon<RealType>() * 10);
         }
      }
      // estimate success ratio:
      BOOST_CHECK_CLOSE(
         beta_distribution<RealType>::estimate_lower_bound_on_p(
            N, k, Q),
         p, tol);
      BOOST_CHECK_CLOSE(
         beta_distribution<RealType>::estimate_upper_bound_on_p(
            N, k, P),
         p, tol);

      if(Q < P)
      {
         BOOST_CHECK(
            beta_distribution<RealType>::estimate_lower_bound_on_p(
               N, k, Q)
               <
            beta_distribution<RealType>::estimate_upper_bound_on_p(
               N, k, Q)
               );
      }
      else
      {
         BOOST_CHECK(
            beta_distribution<RealType>::estimate_lower_bound_on_p(
               N, k, P)
               <
            beta_distribution<RealType>::estimate_upper_bound_on_p(
               N, k, P)
               );
      }
      //
      // estimate sample size:
      //
      BOOST_CHECK_CLOSE(
         beta_distribution<RealType>::estimate_number_of_trials(
            k, p, P),
         N, tol);
      BOOST_CHECK_CLOSE(
         beta_distribution<RealType>::estimate_number_of_trials(
            boost::math::complement(k, p, Q)),
         N, tol);
   }

   // Double check consistency of CDF and PDF by computing
   // the finite sum:
   RealType sum = 0;
   for(unsigned i = 0; i <= k; ++i)
      sum += pdf(bn, RealType(i));
   BOOST_CHECK_CLOSE(
      sum, P, tol);
   // And complement as well:
   sum = 0;
   for(RealType i = N; i > k; i -= 1)
      sum += pdf(bn, i);
   if(P < 0.99)
   {
      BOOST_CHECK_CLOSE(
         sum, Q, tol);
   }
   else
   {
      // Not enough information content in P for Q to be meaningful
      RealType tol = (std::max)(2 * Q, boost::math::tools::epsilon<RealType>());
      BOOST_CHECK(sum < tol);
   }
} // template <class RealType> void test_spot

template <class RealType> // Any floating-point type RealType.
void test_spots(RealType)
{
  // Basic sanity checks, MathCAD test data is to double precision only
  // so set tolerance to 100 eps expressed as a fraction, or
  // 100 eps of type double expressed as a fraction,
  // whichever is the larger.

  RealType tolerance = (std::max)
      (boost::math::tools::epsilon<RealType>(),
      static_cast<RealType>(std::numeric_limits<double>::epsilon()));
   tolerance *= 1000; // NO * 100 because is fraction, NOT %.

  cout << "Tolerance = " << tolerance * 100 << "%." << endl;

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
  BOOST_CHECK_EQUAL(cdf(mybeta22, 1), 1.); // Exact unity expected.
  // Complement

  BOOST_CHECK_CLOSE_FRACTION(cdf(complement(mybeta22, 0.9)), 0.028000000000000, tol);

  // quantile.
  BOOST_CHECK_CLOSE_FRACTION(quantile(mybeta22, 0.028), 0.1, tol); 
  BOOST_CHECK_CLOSE_FRACTION(quantile(complement(mybeta22, 1 - 0.028)), 0.1, tol);
  BOOST_CHECK_EQUAL(kurtosis(mybeta11), 3+ kurtosis_excess(mybeta11)); // Check kurtosis_excess = kurtosis - 3;
  BOOST_CHECK_CLOSE(variance(mybeta22), 0.05, tol);
  BOOST_CHECK_CLOSE(mode(mybeta22), 0.5, tol);
  BOOST_CHECK_CLOSE(mean(mybeta22), 0.5, tol);
  BOOST_CHECK_EQUAL(beta_distribution<double>::estimate_alpha(mean(mybeta22), variance(mybeta22), 0.5), mybeta22.alpha()); // mean, variance, probability. 
  BOOST_CHECK_EQUAL(beta_distribution<double>::estimate_beta(mean(mybeta22), variance(mybeta22), 0.5), mybeta22.beta());// mean, variance, probability. 


  // Basic sanity-check spot values.
#ifdef BOOST_MATH_THROW_ON_DOMAIN_ERROR
  cout << "BOOST_MATH_THROW_ON_DOMAIN_ERROR" << " is defined to throw on domain error." << endl;
#else
  cout << "BOOST_MATH_THROW_ON_DOMAIN_ERROR" << " is NOT defined, so NO throw on domain error." << endl;
#endif

  // (Parameter value, arbitrarily zero, only communicates the floating point type).
//  test_spots(0.0F); // Test float.
//  test_spots(0.0); // Test double.
//  test_spots(0.0L); // Test long double.
//#if !BOOST_WORKAROUND(__BORLANDC__, BOOST_TESTED_AT(0x582))
//  test_spots(boost::math::concepts::real_concept(0.)); // Test real concept.
//#endif

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
  BOOST_CHECK_CLOSE(
     pdf(beta_distribution<RealType>(static_cast<RealType>(20), static_cast<RealType>(0.75)),
     static_cast<RealType>(10)),  // k.
     static_cast<RealType>(0.00992227527967770583927631378173), // 0.00992227527967770583927631378173
     tolerance);

  BOOST_CHECK_CLOSE(
    pdf(beta_distribution<RealType>(static_cast<RealType>(20), static_cast<RealType>(0.5)),
    static_cast<RealType>(10)),  // k.
    static_cast<RealType>(0.17619705200195312500000000000000000000), // get k=10 0.049611376398388612 p = 0.25
    tolerance);

     BOOST_CHECK_CLOSE(
    pdf(beta_distribution<RealType>(static_cast<RealType>(20), static_cast<RealType>(0.25)),
    static_cast<RealType>(10)),  // k.
    static_cast<RealType>(0.00992227527967770583927631378173), // k=10  p = 0.25
    tolerance);

    BOOST_CHECK_CLOSE( // k = 0 use different formula - only exp so more accurate.
    pdf(beta_distribution<RealType>(static_cast<RealType>(20), static_cast<RealType>(0.25)),
    static_cast<RealType>(0)),  // k.
    static_cast<RealType>(0.00317121193893399322405457496643), // k=0  p = 0.25
    tolerance);

    BOOST_CHECK_CLOSE( // k = 20 use different formula - only exp so more accurate.
    pdf(beta_distribution<RealType>(static_cast<RealType>(20), static_cast<RealType>(0.25)),
    static_cast<RealType>(20)),  // k == n.
    static_cast<RealType>(0.00000000000090949470177292823791), // k=20  p = 0.25
    tolerance);

    BOOST_CHECK_CLOSE( // k = 1.
    pdf(beta_distribution<RealType>(static_cast<RealType>(20), static_cast<RealType>(0.25)),
    static_cast<RealType>(1)),  // k.
    static_cast<RealType>(0.02114141292622662149369716644287), // k=1  p = 0.25
    tolerance);

    // Some exact (probably) values.
    BOOST_CHECK_CLOSE(
    pdf(beta_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(0.25)),
    static_cast<RealType>(0)),  // k.
    static_cast<RealType>(0.10011291503906250000000000000000), // k=0  p = 0.25
    tolerance);

    BOOST_CHECK_CLOSE( // k = 1.
    pdf(beta_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(0.25)),
    static_cast<RealType>(1)),  // k.
    static_cast<RealType>(0.26696777343750000000000000000000), // k=1  p = 0.25
    tolerance);

    BOOST_CHECK_CLOSE( // k = 2.
    pdf(beta_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(0.25)),
    static_cast<RealType>(2)),  // k.
    static_cast<RealType>(0.31146240234375000000000000000000), // k=2  p = 0.25
    tolerance);

    BOOST_CHECK_CLOSE( // k = 3.
    pdf(beta_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(0.25)),
    static_cast<RealType>(3)),  // k.
    static_cast<RealType>(0.20764160156250000000000000000000), // k=3  p = 0.25
    tolerance);

    BOOST_CHECK_CLOSE( // k = 7.
    pdf(beta_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(0.25)),
    static_cast<RealType>(7)),  // k.
    static_cast<RealType>(0.00036621093750000000000000000000), // k=7  p = 0.25
    tolerance);

    BOOST_CHECK_CLOSE( // k = 8.
    pdf(beta_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(0.25)),
    static_cast<RealType>(8)),  // k = n.
    static_cast<RealType>(0.00001525878906250000000000000000), // k=8  p = 0.25
    tolerance);

    RealType tol2 = boost::math::tools::epsilon<RealType>() * 5 * 100;  // 5 eps as a persent
    beta_distribution<RealType> dist(static_cast<RealType>(8), static_cast<RealType>(0.25));
    RealType x = static_cast<RealType>(0.125);
    using namespace std; // ADL of std names.
    // mean:
    BOOST_CHECK_CLOSE(
       mean(dist)
       , static_cast<RealType>(8 * 0.25), tol2);
    // variance:
    BOOST_CHECK_CLOSE(
       variance(dist)
       , static_cast<RealType>(8 * 0.25 * 0.75), tol2);
    // std deviation:
    BOOST_CHECK_CLOSE(
       standard_deviation(dist)
       , static_cast<RealType>(sqrt(8 * 0.25L * 0.75L)), tol2);
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
       , static_cast<RealType>(std::floor(9 * 0.25)), tol2);
    // skewness:
    BOOST_CHECK_CLOSE(
       skewness(dist)
       , static_cast<RealType>(0.40824829046386301636621401245098L), tol2);
    // kurtosis:
    BOOST_CHECK_CLOSE(
       kurtosis(dist)
       , static_cast<RealType>(2.9166666666666666666666666666667L), tol2);
    // kurtosis excess:
    BOOST_CHECK_CLOSE(
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

