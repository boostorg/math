// test_negative_binomial.cpp

// Copyright Paul A. Bristow 2006.

// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

// Tests for Negative Binomial Distribution.

#define BOOST_MATH_THROW_ON_DOMAIN_ERROR

#ifdef _MSC_VER
#  pragma warning(disable: 4127) // conditional expression is constant.
#  pragma warning(disable: 4100) // unreferenced formal parameter.
#  pragma warning(disable: 4512) // assignment operator could not be generated.
#endif

#include <boost/math/distributions/negative_binomial.hpp> // for negative_binomial_distribution
using boost::math::negative_binomial_distribution;

#include <boost/math/special_functions/gamma.hpp>
  using boost::math::lgamma;  // log gamma

#include <boost/math/concepts/real_concept.hpp> // for real_concept
using ::boost::math::concepts::real_concept;

#include <boost/test/included/test_exec_monitor.hpp> // for test_main
#include <boost/test/floating_point_comparison.hpp> // for BOOST_CHECK_CLOSE

#include <iostream>
using std::cout;
using std::endl;
using std::setprecision;
using std::showpoint;
#include <limits>
using std::numeric_limits;

// Need all functions working to use this test_spot!

template <class RealType>
void test_spot( // Test a single spot value against 'known good' values.
               RealType N,    // Number of successes.
               RealType k,    // Number of failures.
               RealType p,    // Probability of success_fraction.
               RealType P,    // CDF probability.
               RealType Q,    // Complement of CDF.
               RealType tol)  // Test tolerance.
{
  boost::math::negative_binomial_distribution<RealType> bn(N, p);
   BOOST_CHECK_EQUAL(N, bn.successes());
   BOOST_CHECK_EQUAL(p, bn.success_fraction());
   BOOST_CHECK_CLOSE(
     cdf(bn, k), P, tol);

  if((P < 0.99) && (Q < 0.99))
  {
    // We can only check this if P is not too close to 1,
    // so that we can guarantee that Q is free of error:
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
#if 0
    // estimate success ratio:
    BOOST_CHECK_CLOSE(
      negative_binomial_distribution<RealType>::estimate_lower_bound_on_p(
      N, k, Q), 
      p, tol);
    BOOST_CHECK_CLOSE(
      negative_binomial_distribution<RealType>::estimate_upper_bound_on_p(
      N, k, P), 
      p, tol);

    if(Q < P)
    {
      BOOST_CHECK(
        negative_binomial_distribution<RealType>::estimate_lower_bound_on_p(
        N, k, Q)
        < 
        negative_binomial_distribution<RealType>::estimate_upper_bound_on_p(
        N, k, Q)
        );
    }
    else
    {
      BOOST_CHECK(
        negative_binomial_distribution<RealType>::estimate_lower_bound_on_p(
        N, k, P)
        < 
        negative_binomial_distribution<RealType>::estimate_upper_bound_on_p(
        N, k, P)
        );
    }
#endif
    //
    // Estimate sample size:
    //
    //BOOST_CHECK_CLOSE(
    //  negative_binomial_distribution<RealType>::estimate_number_of_successes(
    //  k, p, P), 
    //  N, tol);
    //BOOST_CHECK_CLOSE(
    //  negative_binomial_distribution<RealType>::estimate_number_of_successes(
    //  boost::math::complement(k, p, Q)), 
    //  N, tol);

    // Double check consistency of CDF and PDF by computing the finite sum:
    RealType sum = 0;
    for(unsigned i = 0; i <= k; ++i)
    {
      sum += pdf(bn, RealType(i));
    }
    BOOST_CHECK_CLOSE(sum, P, tol);

    // Complement is not possible since sum is to infinity.
  } // 
} // test_spot

template <class RealType> // Any floating-point type RealType.
void test_spots(RealType)
{
  // Basic sanity checks, test data is to double precision only
  // so set tolerance to 100eps expressed as a percent, or
  // 100eps of type double expressed as a percent, whichever
  // is the larger.

  RealType tolerance = (std::max)
    (boost::math::tools::epsilon<RealType>(),
    static_cast<RealType>(std::numeric_limits<double>::epsilon()));
  tolerance *= 100 * 1000;

  cout << "Tolerance = " << tolerance << "%." << endl;

  RealType tol1eps = boost::math::tools::epsilon<RealType>() * 2; // Very tight, suit exact values.
  RealType tol2eps = boost::math::tools::epsilon<RealType>() * 2; // Tight, suit exact values.
  RealType tol5eps = boost::math::tools::epsilon<RealType>() * 5; // Wider 5 epsilon.

  // Sources of spot test values:

  // MathCAD defines pbinom(k, r, p)
  // returns pr(X , k) when random variable X has the binomial distribution with parameters r and p.
  // 0 <= k 
  // r > 0
  // 0 <= p <= 1
  // P = pbinom(30, 500, 0.05) = 0.869147702104609

  using boost::math::negative_binomial_distribution;
  using  ::boost::math::negative_binomial;
  using  ::boost::math::cdf;
  using  ::boost::math::pdf;

  // Test negative binomial using cdf spot values from MathCAD cdf = pnbinom(k, r, p).
  // These test quantiles and complements as well.

  test_spot(  // pnbinom(1,2,0.5) = 0.5
  static_cast<RealType>(2),   // successes r
  static_cast<RealType>(1),   // Number of failures, k
  static_cast<RealType>(0.5), // Probability of success as fraction, p
  static_cast<RealType>(0.5), // Probability of result (CDF), P
  static_cast<RealType>(0.5),  // complement CCDF Q = 1 - P
  tolerance);

  test_spot( // pbinom(0, 2, 0.25)
  static_cast<RealType>(2),    // successes r
  static_cast<RealType>(0),    // Number of failures, k
  static_cast<RealType>(0.25),   
  static_cast<RealType>(0.0625),                    // Probability of result (CDF), P
  static_cast<RealType>(0.9375),                    // Q = 1 - P
  tolerance);

  test_spot(  // pbinom(48,8,0.25)
  static_cast<RealType>(8),     // successes r
  static_cast<RealType>(48),    // Number of failures, k
  static_cast<RealType>(0.25),                    // Probability of success, p
  static_cast<RealType>(9.826582228110670E-1),     // Probability of result (CDF), P
  static_cast<RealType>(1 - 9.826582228110670E-1),   // Q = 1 - P
  tolerance);

  test_spot(  // pbinom(2,5,0.4)
  static_cast<RealType>(5),     // successes r
  static_cast<RealType>(2),     // Number of failures, k
  static_cast<RealType>(0.4),                    // Probability of success, p
  static_cast<RealType>(9.625600000000020E-2),     // Probability of result (CDF), P
  static_cast<RealType>(1 - 9.625600000000020E-2),   // Q = 1 - P
  tolerance);

  test_spot(  // pbinom(10,100,0.9)
  static_cast<RealType>(100),     // successes r
  static_cast<RealType>(10),     // Number of failures, k
  static_cast<RealType>(0.9),                    // Probability of success, p
  static_cast<RealType>(4.535522887695670E-1),     // Probability of result (CDF), P
  static_cast<RealType>(1 - 4.535522887695670E-1),   // Q = 1 - P
  tolerance);

  test_spot(  // pbinom(1,100,0.991)
  static_cast<RealType>(100),     // successes r
  static_cast<RealType>(1),     // Number of failures, k
  static_cast<RealType>(0.991),                    // Probability of success, p
  static_cast<RealType>(7.693413044217000E-1),     // Probability of result (CDF), P
  static_cast<RealType>(1 - 7.693413044217000E-1),   // Q = 1 - P
  tolerance);
  
  test_spot(  // pbinom(10,100,0.991)
  static_cast<RealType>(100),     // successes r
  static_cast<RealType>(10),     // Number of failures, k
  static_cast<RealType>(0.991),                    // Probability of success, p
  static_cast<RealType>(9.999999940939000E-1),     // Probability of result (CDF), P
  static_cast<RealType>(1 - 9.999999940939000E-1),   // Q = 1 - P
  tolerance);

  test_spot(  // pbinom(100000,100,0.001)
  static_cast<RealType>(100),     // successes r
  static_cast<RealType>(100000),     // Number of failures, k
  static_cast<RealType>(0.001),                    // Probability of success, p
  static_cast<RealType>(5.173047534260320E-1),     // Probability of result (CDF), P
  static_cast<RealType>(1 - 5.173047534260320E-1),   // Q = 1 - P
  tolerance*1000); // *1000 is OK 0.51730475350664229  versus
  //                       MathCAD 0.51730475342603199 differs at 10th decimal digit.

 // End of single spot tests using RealType

  // Tests on cdf:
// MathCAD pbinom k, r, p) == failures, successes, 

  BOOST_CHECK_CLOSE(cdf( 
    negative_binomial_distribution<RealType>(static_cast<RealType>(2), static_cast<RealType>(0.5)), // successes = 2,prob 0.25
    static_cast<RealType>(0) ), // k = 0
    static_cast<RealType>(0.25), // probability 1/4
    tolerance);

  BOOST_CHECK_CLOSE(cdf(complement(
    negative_binomial_distribution<RealType>(static_cast<RealType>(2), static_cast<RealType>(0.5)), // successes = 2,prob 0.25
    static_cast<RealType>(0) )), // k = 0
    static_cast<RealType>(0.75), // probability 3/4
    tolerance);

  // Tests on PDF:
  BOOST_CHECK_CLOSE(
  pdf(negative_binomial_distribution<RealType>(static_cast<RealType>(2), static_cast<RealType>(0.5)),
  static_cast<RealType>(0) ),  // k = 0.
  static_cast<RealType>(0.25), // 0
  tolerance);
 
  BOOST_CHECK_CLOSE(
  pdf(negative_binomial_distribution<RealType>(static_cast<RealType>(4), static_cast<RealType>(0.5)), 
  static_cast<RealType>(0)),  // k = 0.
  static_cast<RealType>(0.0625), // exact 1/16
  tolerance);

  BOOST_CHECK_CLOSE(
  pdf(negative_binomial_distribution<RealType>(static_cast<RealType>(20), static_cast<RealType>(0.25)), 
  static_cast<RealType>(0)),  // k = 0
  static_cast<RealType>(9.094947017729270E-13), // pbinom(0,20,0.25) = 9.094947017729270E-13
  tolerance); 

  BOOST_CHECK_CLOSE(
  pdf(negative_binomial_distribution<RealType>(static_cast<RealType>(20), static_cast<RealType>(0.2)), 
  static_cast<RealType>(0)),  // k = 0
  static_cast<RealType>(1.0485760000000003e-014), // MathCAD 1.048576000000000E-14
  tolerance);

  BOOST_CHECK_CLOSE( 
  pdf(negative_binomial_distribution<RealType>(static_cast<RealType>(10), static_cast<RealType>(0.1)), 
  static_cast<RealType>(0)),  // k = 0.
  static_cast<RealType>(1e-10), // MathCAD says zero, but suffers cancellation error?
  tolerance); 

  BOOST_CHECK_CLOSE( 
  pdf(negative_binomial_distribution<RealType>(static_cast<RealType>(20), static_cast<RealType>(0.1)), 
  static_cast<RealType>(0)),  // k = 0.
  static_cast<RealType>(1e-20), // MathCAD says zero, but suffers cancellation error?
  tolerance); 


  BOOST_CHECK_CLOSE( // .
  pdf(negative_binomial_distribution<RealType>(static_cast<RealType>(20), static_cast<RealType>(0.9)), 
  static_cast<RealType>(0)),  // k.
  static_cast<RealType>(1.215766545905690E-1), // k=20  p = 0.9
  tolerance); 

  // Test on cdf

  BOOST_CHECK_CLOSE( // k = 1.
  cdf(negative_binomial_distribution<RealType>(static_cast<RealType>(20), static_cast<RealType>(0.25)), 
  static_cast<RealType>(1)),  // k =1.
  static_cast<RealType>(1.455191522836700E-11),
  tolerance); 

  BOOST_CHECK_SMALL( // Check within an epsilon with CHECK_SMALL
  cdf(negative_binomial_distribution<RealType>(static_cast<RealType>(20), static_cast<RealType>(0.25)), 
  static_cast<RealType>(1)) -
  static_cast<RealType>(1.455191522836700E-11),
  static_cast<RealType>(tol1eps) );
  
  // Some exact (probably - judging by trailing zeros) values.
  BOOST_CHECK_CLOSE( 
  cdf(negative_binomial_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(0.25)), 
  static_cast<RealType>(0)),  // k.
  static_cast<RealType>(1.525878906250000E-5), 
  tolerance); 

  BOOST_CHECK_CLOSE( 
  cdf(negative_binomial_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(0.25)), 
  static_cast<RealType>(0)),  // k.
  static_cast<RealType>(1.525878906250000E-5), 
  tolerance); 

  BOOST_CHECK_SMALL( 
  cdf(negative_binomial_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(0.25)), 
  static_cast<RealType>(0)) - 
  static_cast<RealType>(1.525878906250000E-5),
  tol2eps );

  BOOST_CHECK_CLOSE( // k = 1.
  cdf(negative_binomial_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(0.25)), 
  static_cast<RealType>(1)),  // k.
  static_cast<RealType>(1.068115234375010E-4), 
  tolerance); 

  BOOST_CHECK_CLOSE( // k = 2.
  cdf(negative_binomial_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(0.25)), 
  static_cast<RealType>(2)),  // k.
  static_cast<RealType>(4.158020019531300E-4),
  tolerance); 

  BOOST_CHECK_CLOSE( // k = 3.
  cdf(negative_binomial_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(0.25)), 
  static_cast<RealType>(3)),  // k.bristow
  static_cast<RealType>(1.188278198242200E-3),
  tolerance); 
  
  BOOST_CHECK_CLOSE( // k = 4.
  cdf(negative_binomial_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(0.25)), 
  static_cast<RealType>(4)),  // k.
  static_cast<RealType>(2.781510353088410E-3),
  tolerance); 

  BOOST_CHECK_CLOSE( // k = 5.
  cdf(negative_binomial_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(0.25)), 
  static_cast<RealType>(5)),  // k.
  static_cast<RealType>(5.649328231811500E-3),
  tolerance); 

  BOOST_CHECK_CLOSE( // k = 6.
  cdf(negative_binomial_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(0.25)), 
  static_cast<RealType>(6)),  // k.
  static_cast<RealType>(1.030953228473680E-2),
  tolerance); 

  BOOST_CHECK_CLOSE( // k = 7.
  cdf(negative_binomial_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(0.25)), 
  static_cast<RealType>(7)),  // k.
  static_cast<RealType>(1.729983836412430E-2),
  tolerance); 

  BOOST_CHECK_CLOSE( // k = 8.
  cdf(negative_binomial_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(0.25)), 
  static_cast<RealType>(8)),  // k = n.
  static_cast<RealType>(2.712995628826370E-2), 
  tolerance); 

  BOOST_CHECK_CLOSE( //
  cdf(negative_binomial_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(0.25)), 
  static_cast<RealType>(48)),  // k 
  static_cast<RealType>(9.826582228110670E-1), 
  tolerance); 

  BOOST_CHECK_CLOSE( //
  cdf(negative_binomial_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(0.25)), 
  static_cast<RealType>(64)),  // k 
  static_cast<RealType>(9.990295004935590E-1), 
  tolerance); 

  BOOST_CHECK_CLOSE( //
  cdf(negative_binomial_distribution<RealType>(static_cast<RealType>(5), static_cast<RealType>(0.4)), 
  static_cast<RealType>(26)),  // k 
  static_cast<RealType>(9.989686246611190E-1), 
  tolerance); 

  BOOST_CHECK_CLOSE( //
  cdf(negative_binomial_distribution<RealType>(static_cast<RealType>(5), static_cast<RealType>(0.4)), 
  static_cast<RealType>(2)),  // k failures
  static_cast<RealType>(9.625600000000020E-2), 
  tolerance); 

  BOOST_CHECK_CLOSE( //
  cdf(negative_binomial_distribution<RealType>(static_cast<RealType>(50), static_cast<RealType>(0.9)), 
  static_cast<RealType>(20)),  // k 
  static_cast<RealType>(9.999970854144170E-1), 
  tolerance); 

  BOOST_CHECK_CLOSE( //
  cdf(negative_binomial_distribution<RealType>(static_cast<RealType>(500), static_cast<RealType>(0.7)), 
  static_cast<RealType>(200)),  // k 
  static_cast<RealType>(2.172846379930550E-1), 
  tolerance* 2); 

  BOOST_CHECK_CLOSE( //
  cdf(negative_binomial_distribution<RealType>(static_cast<RealType>(50), static_cast<RealType>(0.7)), 
  static_cast<RealType>(20)),  // k 
  static_cast<RealType>(4.550203671301790E-1),
  tolerance); 


  negative_binomial_distribution<RealType> dist(static_cast<RealType>(8), static_cast<RealType>(0.25));
  using namespace std; // ADL of std names.
  // mean:
  BOOST_CHECK_CLOSE(
    mean(dist)
    , static_cast<RealType>(8 * ( 1- 0.25) /0.25), tol5eps);
  // variance:
  BOOST_CHECK_CLOSE(
    variance(dist)
    , static_cast<RealType>(8 * (1-0.25) / (0.25 * 0.25)), tol5eps);
  // std deviation:
  BOOST_CHECK_CLOSE(
    standard_deviation(dist)
    , static_cast<RealType>(sqrt(8 * (1-0.25) / (0.25 * 0.25))), tol5eps);
  // hazard:
  RealType x = static_cast<RealType>(0.125);
  BOOST_CHECK_CLOSE(
  hazard(dist, x)
  , pdf(dist, x) / cdf(complement(dist, x)), tol5eps);
  // cumulative hazard:
  BOOST_CHECK_CLOSE(
  chf(dist, x)
  , -log(cdf(complement(dist, x))), tol5eps);
  // coefficient_of_variation:
  BOOST_CHECK_CLOSE(
  coefficient_of_variation(dist)
  , standard_deviation(dist) / mean(dist), tol5eps);

  // Special cases for PDF:
  BOOST_CHECK_EQUAL(
  pdf(
  negative_binomial_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(0)), // 
  static_cast<RealType>(0)),
  static_cast<RealType>(0) );

  BOOST_CHECK_EQUAL(
  pdf(
  negative_binomial_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(0)),
  static_cast<RealType>(0.0001)), 
  static_cast<RealType>(0) );

  BOOST_CHECK_EQUAL(
  pdf(
  negative_binomial_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(1)),
  static_cast<RealType>(0.001)),
  static_cast<RealType>(0) );

  BOOST_CHECK_EQUAL(
  pdf(
  negative_binomial_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(1)),
  static_cast<RealType>(8)),
  static_cast<RealType>(0) );
  BOOST_CHECK_SMALL(
  pdf(
   negative_binomial_distribution<RealType>(static_cast<RealType>(2), static_cast<RealType>(0.25)),
  static_cast<RealType>(0))-
  static_cast<RealType>(0.0625),
  numeric_limits<RealType>::epsilon()); // Expect exact, but not quite.


  // Check that duff arguments throw domain_error:
  BOOST_CHECK_THROW(
  pdf( // Negative successes!
  negative_binomial_distribution<RealType>(static_cast<RealType>(-1), static_cast<RealType>(0.25)),
  static_cast<RealType>(0)), std::domain_error
  );
  BOOST_CHECK_THROW(
  pdf( // Negative success_fraction!
  negative_binomial_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(-0.25)),
  static_cast<RealType>(0)), std::domain_error
  );
  BOOST_CHECK_THROW(
  pdf( // Success_fraction > 1!
  negative_binomial_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(1.25)),
  static_cast<RealType>(0)),
  std::domain_error
  );
  BOOST_CHECK_THROW(
  pdf( // Negative k argument !
  negative_binomial_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(0.25)),
  static_cast<RealType>(-1)),
  std::domain_error
  );
  //BOOST_CHECK_THROW(
  //pdf( // Unlike binomial there is NO limit on k (failures)
  //negative_binomial_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(0.25)),
  //static_cast<RealType>(9)), std::domain_error
  //);
  BOOST_CHECK_THROW(
  cdf(  // Negative k argument !
  negative_binomial_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(0.25)),
  static_cast<RealType>(-1)),
  std::domain_error
  );
  BOOST_CHECK_THROW(
  cdf( // Negative success_fraction!
  negative_binomial_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(-0.25)),
  static_cast<RealType>(0)), std::domain_error
  );
  BOOST_CHECK_THROW(
  cdf( // Success_fraction > 1!
  negative_binomial_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(1.25)),
  static_cast<RealType>(0)), std::domain_error
  );
  BOOST_CHECK_THROW(
  quantile(  // Negative success_fraction!
  negative_binomial_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(-0.25)),
  static_cast<RealType>(0)), std::domain_error
  );
  BOOST_CHECK_THROW(
  quantile( // Success_fraction > 1! 
  negative_binomial_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(1.25)),
  static_cast<RealType>(0)), std::domain_error
  );
  // End of check throwing 'duff' out-of-domain values.

  return;

} // template <class RealType> void test_spots(RealType) // Any floating-point type RealType.

int test_main(int, char* [])
{
  // Check that can generate negative_binomial distribution using the two convenience methods:
  using namespace boost::math;
	negative_binomial mynb1(2., 0.5); // Using typedef - default type is double.
	negative_binomial_distribution<> myf2(2., 0.5); // Using default RealType double.

  // Basic sanity-check spot values.
#ifdef BOOST_MATH_THROW_ON_DOMAIN_ERROR
  cout << "BOOST_MATH_THROW_ON_DOMAIN_ERROR" << " is defined to throw on domain error." << endl;
#else
  cout << "BOOST_MATH_THROW_ON_DOMAIN_ERROR" << " is NOT defined, so NO throw on domain error." << endl;
#endif

  // This is a visual sanity check that everything is OK:
  cout << setprecision(17) << showpoint << endl;
  // Show double max_digits10 precision, including trailing zeros.

  // Test some simple double only examples.
  negative_binomial_distribution<double> my8dist(8., 0.25);
  // 8 successes (r), 0.25 success fraction = 35% or 1 in 4 successes.
  // Note: double values (matching the distribution definition) avoid the need for any casting.

  BOOST_CHECK_EQUAL(my8dist.successes(), static_cast<double>(8));
  BOOST_CHECK_EQUAL(my8dist.success_fraction(), static_cast<double>(1./4.)); // Exact.

  double tol = boost::math::tools::epsilon<double>() * 100 * 10;
  // * 100 for %, so tol is 10 epsilon.
  BOOST_CHECK_SMALL(cdf(my8dist, 2.), 4.1580200195313E-4); 
  BOOST_CHECK_SMALL(cdf(my8dist, 8.), 0.027129956288264);
  BOOST_CHECK_CLOSE(cdf(my8dist, 16.), 0.233795830683125, tol);

  BOOST_CHECK_CLOSE(pdf(
  negative_binomial_distribution<double>(2., 0.5),
  1.),
  static_cast<double>(0.25),
  tol);

  BOOST_CHECK_CLOSE(cdf(complement(
    negative_binomial_distribution<double>(2., 0.5), // 
    1.)), // k 
    static_cast<double>(0.5), // half
    tol);

  BOOST_CHECK_CLOSE(cdf(
    negative_binomial_distribution<double>(2., 0.5),
    1.),
    static_cast<double>(0.5),
    tol);

  BOOST_CHECK_CLOSE(
  quantile(
  negative_binomial_distribution<double>(2., 0.5),
  0.5),
  1., // 
  tol);

  BOOST_CHECK_CLOSE(
  cdf(
  negative_binomial_distribution<double>(8., 0.25),
  16.),
  0.233795830683125, // 
  tol);

  BOOST_CHECK_CLOSE(
  quantile(
  negative_binomial_distribution<double>(8., 0.25),
  0.233795830683125),
  16., // 
  tol);

    // Compare pdf using the simple formulae:
  //   exp(lgamma(r + k) - lgamma(r) - lgamma(k+1)) * pow(p, r) * pow((1-p), k)
  // with
  //   (p/(r+k) * ibeta_derivative(r, k+1, p) (as used in pdf)
  {
   typedef double RealType;
   RealType r = my8dist.successes();
    RealType p = my8dist.success_fraction();
    for (int i = 0; i <= r * 4; i++) // 
    {
      RealType k = static_cast<RealType>(i);
      RealType pmf = exp(lgamma(r + k) - lgamma(r) - lgamma(k+1)) * pow(p, r) * pow((1-p), k);
      BOOST_CHECK_CLOSE(pdf(my8dist, static_cast<RealType>(k)), pmf, tol * 10);
      // 0.0015932321548461931
    // 0.0015932321548461866
  }

  // Double-check consistency of CDF and PDF by computing the finite sum of pdfs:
  RealType sum = 0;
  for(unsigned i = 0; i <= 20; ++i)
  {
    sum += pdf(my8dist, RealType(i));
  }

  cout << setprecision(17) << showpoint << sum <<' '  // 0.40025683281803714 
  << cdf(my8dist, static_cast<RealType>(20)) << endl; // 0.40025683281803681
  BOOST_CHECK_CLOSE(sum, cdf(my8dist, static_cast<RealType>(20)), tol);
  }

  // (Parameter value, arbitrarily zero, only communicates the floating point type).
  test_spots(0.0F); // Test float.
  test_spots(0.0); // Test double.
  test_spots(0.0L); // Test long double.
  //#if !BOOST_WORKAROUND(__BORLANDC__, BOOST_TESTED_AT(0x582))
  //  test_spots(boost::math::concepts::real_concept(0.)); // Test real concept.
  //#endif

  return 0;
} // int test_main(int, char* [])

/*

------ Build started: Project: test_negative_binomial, Configuration: Debug Win32 ------
Compiling...
test_negative_binomial.cpp
Linking...
Autorun "i:\boost-06-05-03-1300\libs\math\test\Math_test\debug\test_negative_binomial.exe"
Running 1 test case...
BOOST_MATH_THROW_ON_DOMAIN_ERROR is defined to throw on domain error.
0.40025683281803714 0.40025683281803681
Tolerance = 0.0119209%.
Tolerance = 2.22045e-011%.
Tolerance = 2.22045e-011%.
*** No errors detected
Build Time 0:07
Build log was saved at "file://i:\boost-06-05-03-1300\libs\math\test\Math_test\test_negative_binomial\Debug\BuildLog.htm"
test_negative_binomial - 0 error(s), 0 warning(s)
========== Build: 1 succeeded, 0 failed, 0 up-to-date, 0 skipped ==========

BUT hangs doing real concept :-(


*/
