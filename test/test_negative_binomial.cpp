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
   BOOST_CHECK_EQUAL(N, bn.success());
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
      // Not enough information content in P for Q to be meaningful!
      RealType tol = (std::max)(2 * Q, boost::math::tools::epsilon<RealType>());
      BOOST_CHECK(sum < tol);
    }
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

  // Test negative binomial using cdf spot values from MathCAD.
  // These test quantiles and complements as well.
  //test_spot(
  //static_cast<RealType>(2),      // Trial n
  //static_cast<RealType>(1),      // Number of failures, k
  //static_cast<RealType>(0.5),    // Probability of success, p
  //static_cast<RealType>(0.25),   // Probability of result (CDF), P
  //static_cast<RealType>(0.75),   // Q = 1 - P
  //tolerance);


  //test_spot(
  //static_cast<RealType>(500),                     // Sample size, N
  //static_cast<RealType>(250),                     // Number of successes, k
  //static_cast<RealType>(0.05),                    // Probability of success, p
  //static_cast<RealType>(1),                       // Probability of result (CDF), P
  //static_cast<RealType>(0),   // Q = 1 - P
  //tolerance);

  //test_spot(
  //static_cast<RealType>(500),                     // Sample size, N
  //static_cast<RealType>(470),                     // Number of successes, k
  //static_cast<RealType>(0.95),                    // Probability of success, p
  //static_cast<RealType>(0.176470742656766),       // Probability of result (CDF), P
  //static_cast<RealType>(1 - 0.176470742656766),   // Q = 1 - P
  //tolerance * 10);                                // Note higher tolerance on this test!

  //test_spot(
  //static_cast<RealType>(500),                       // Sample size, N
  //static_cast<RealType>(400),                       // Number of successes, k
  //static_cast<RealType>(0.05),                      // Probability of success, p
  //static_cast<RealType>(1),                         // Probability of result (CDF), P
  //static_cast<RealType>(0),                         // Q = 1 - P
  //tolerance);

  //test_spot(
  //static_cast<RealType>(500),                       // Sample size, N
  //static_cast<RealType>(400),                       // Number of successes, k
  //static_cast<RealType>(0.9),                       // Probability of success, p
  //static_cast<RealType>(1.80180425681923E-11),      // Probability of result (CDF), P
  //static_cast<RealType>(1 - 1.80180425681923E-11),  // Q = 1 - P
  //tolerance);

  //test_spot(
  //static_cast<RealType>(500),                       // Sample size, N
  //static_cast<RealType>(5),                         // Number of successes, k
  //static_cast<RealType>(0.05),                      // Probability of success, p
  //static_cast<RealType>(9.181808267643E-7),         // Probability of result (CDF), P
  //static_cast<RealType>(1 - 9.181808267643E-7),     // Q = 1 - P
  //tolerance);

  //test_spot(
  //static_cast<RealType>(2),                       // Sample size, N
  //static_cast<RealType>(1),                       // Number of successes, k
  //static_cast<RealType>(0.5),                     // Probability of success, p
  //static_cast<RealType>(0.75),                    // Probability of result (CDF), P
  //static_cast<RealType>(0.25),                    // Q = 1 - P
  //tolerance);

  //test_spot(
  //static_cast<RealType>(8),                       // Sample size, N
  //static_cast<RealType>(3),                       // Number of successes, k
  //static_cast<RealType>(0.25),                    // Probability of success, p
  //static_cast<RealType>(0.8861846923828125),      // Probability of result (CDF), P
  //static_cast<RealType>(1 - 0.8861846923828125),  // Q = 1 - P
  //tolerance);

  //test_spot(
  //static_cast<RealType>(8),                       // Sample size, N
  //static_cast<RealType>(0),                       // Number of successes, k
  //static_cast<RealType>(0.25),                    // Probability of success, p
  //static_cast<RealType>(0.1001129150390625),      // Probability of result (CDF), P
  //static_cast<RealType>(1 - 0.1001129150390625),  // Q = 1 - P
  //tolerance);

  //test_spot(
  //static_cast<RealType>(8),                       // Sample size, N
  //static_cast<RealType>(1),                       // Number of successes, k
  //static_cast<RealType>(0.25),                    // Probability of success, p
  //static_cast<RealType>(0.36708068847656244),     // Probability of result (CDF), P
  //static_cast<RealType>(1 - 0.36708068847656244), // Q = 1 - P
  //tolerance);

  //test_spot(
  //static_cast<RealType>(8),                       // Sample size, N
  //static_cast<RealType>(4),                       // Number of successes, k
  //static_cast<RealType>(0.25),                    // Probability of success, p
  //static_cast<RealType>(0.9727020263671875),      // Probability of result (CDF), P
  //static_cast<RealType>(1 - 0.9727020263671875),  // Q = 1 - P
  //tolerance);

  //test_spot(
  //static_cast<RealType>(8),                       // Sample size, N
  //static_cast<RealType>(7),                       // Number of successes, k
  //static_cast<RealType>(0.25),                    // Probability of success, p
  //static_cast<RealType>(0.9999847412109375),      // Probability of result (CDF), P
  //static_cast<RealType>(1 - 0.9999847412109375),  // Q = 1 - P
  //tolerance);

  // End of single spot tests using RealType

  // Tests on cdf:
// MathCAD pbinom k, r, p) == failures, successes, 

  BOOST_CHECK_CLOSE(cdf( // pnbinom(1,2,0.5) = 0.5
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
  static_cast<RealType>(0) ),  // k.
  static_cast<RealType>(0.25), // 0
  tolerance);
 

  BOOST_CHECK_CLOSE(
  pdf(negative_binomial_distribution<RealType>(static_cast<RealType>(2), static_cast<RealType>(0.5)), 
  static_cast<RealType>(1)),  // k.
  static_cast<RealType>(0.25), // 
  tolerance);

  BOOST_CHECK_CLOSE(
  cdf(negative_binomial_distribution<RealType>(static_cast<RealType>(20), static_cast<RealType>(0.25)), 
  static_cast<RealType>(10)),  // k = 10
  static_cast<RealType>(1.820643361830440E-6), // pbinom(10,20,0.25) = 1.820643361830440E-6
  tolerance); 

  BOOST_CHECK_CLOSE( // k = 0 use different formula - only exp so more accurate.
  cdf(negative_binomial_distribution<RealType>(static_cast<RealType>(20), static_cast<RealType>(0.25)), 
  static_cast<RealType>(0)),  // k.
  static_cast<RealType>(9.094947017729270E-13), // k=0  p = 0.25
  tolerance); 

  BOOST_CHECK_CLOSE( // k = 20 use different formula - only exp so more accurate.
  pdf(negative_binomial_distribution<RealType>(static_cast<RealType>(20), static_cast<RealType>(0.25)), 
  static_cast<RealType>(20)),  // k == n.
  static_cast<RealType>(0.00000000000090949470177292823791), // k=20  p = 0.25
  tolerance); 

  BOOST_CHECK_CLOSE( // k = 1.
  pdf(negative_binomial_distribution<RealType>(static_cast<RealType>(20), static_cast<RealType>(0.25)), 
  static_cast<RealType>(1)),  // k.
  static_cast<RealType>(0.02114141292622662149369716644287), // k=1  p = 0.25
  tolerance); 

  // Some exact (probably - judging by trailing zeros) values.
  BOOST_CHECK_CLOSE( 
  pdf(negative_binomial_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(0.25)), 
  static_cast<RealType>(0)),  // k.
  static_cast<RealType>(0.10011291503906250000000000000000), // k=0  p = 0.25
  tolerance); 

  BOOST_CHECK_CLOSE( // k = 1.
  pdf(negative_binomial_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(0.25)), 
  static_cast<RealType>(1)),  // k.
  static_cast<RealType>(0.26696777343750000000000000000000), // k=1  p = 0.25
  tolerance); 

  BOOST_CHECK_CLOSE( // k = 2.
  pdf(negative_binomial_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(0.25)), 
  static_cast<RealType>(2)),  // k.
  static_cast<RealType>(0.31146240234375000000000000000000), // k=2  p = 0.25
  tolerance); 

  BOOST_CHECK_CLOSE( // k = 3.
  pdf(negative_binomial_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(0.25)), 
  static_cast<RealType>(3)),  // k.
  static_cast<RealType>(0.20764160156250000000000000000000), // k=3  p = 0.25
  tolerance); 

  BOOST_CHECK_CLOSE( // k = 7.
  pdf(negative_binomial_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(0.25)), 
  static_cast<RealType>(7)),  // k.
  static_cast<RealType>(0.00036621093750000000000000000000), // k=7  p = 0.25
  tolerance); 

  BOOST_CHECK_CLOSE( // k = 8.
  pdf(negative_binomial_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(0.25)), 
  static_cast<RealType>(8)),  // k = n.
  static_cast<RealType>(0.00001525878906250000000000000000), // k=8  p = 0.25
  tolerance); 

  RealType tol2 = boost::math::tools::epsilon<RealType>() * 5;
  negative_binomial_distribution<RealType> dist(static_cast<RealType>(8), static_cast<RealType>(0.25));
  using namespace std; // ADL of std names.
  // mean:
  BOOST_CHECK_CLOSE(
    mean(dist)
    , static_cast<RealType>(8 * ( 1- 0.25) /0.25), tol2);
  // variance:
  BOOST_CHECK_CLOSE(
    variance(dist)
    , static_cast<RealType>(8 * (1-0.25) / (0.25 * 0.25)), tol2);
  // std deviation:
  BOOST_CHECK_CLOSE(
    standard_deviation(dist)
    , static_cast<RealType>(sqrt(8 * (1-0.25) / (0.25 * 0.25))), tol2);
  // hazard:
  RealType x = static_cast<RealType>(0.125);
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

  // Special cases for PDF:
  BOOST_CHECK_EQUAL(
  pdf(
  negative_binomial_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(0)), // 
  static_cast<RealType>(0)),
  static_cast<RealType>(0) // Expected pdf.
  );
  BOOST_CHECK_EQUAL(
  pdf(
  negative_binomial_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(0)),
  static_cast<RealType>(0.0001)), // Expected pdf.
  static_cast<RealType>(0)
  );
  BOOST_CHECK_EQUAL(
  pdf(
  negative_binomial_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(1)),
  static_cast<RealType>(0.001)), // Expected pdf.
  static_cast<RealType>(0)
  );
  BOOST_CHECK_EQUAL(
  pdf(
  negative_binomial_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(1)),
  static_cast<RealType>(8)),
  static_cast<RealType>(0) // Expected pdf.
  );
  BOOST_CHECK_EQUAL(
  pdf(
   negative_binomial_distribution<RealType>(static_cast<RealType>(2), static_cast<RealType>(0.25)),
  static_cast<RealType>(0)),
  static_cast<RealType>(1) // Expected pdf.
  );

  // Check that duff arguments throw domain_error:
  BOOST_CHECK_THROW(
  pdf(
  negative_binomial_distribution<RealType>(static_cast<RealType>(-1), static_cast<RealType>(0.25)),
  static_cast<RealType>(0)), std::domain_error
  );
  BOOST_CHECK_THROW(
  pdf(
  negative_binomial_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(-0.25)),
  static_cast<RealType>(0)), std::domain_error
  );
  BOOST_CHECK_THROW(
  pdf(
  negative_binomial_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(1.25)),
  static_cast<RealType>(0)), std::domain_error
  );
  BOOST_CHECK_THROW(
  pdf(
  negative_binomial_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(0.25)),
  static_cast<RealType>(-1)), std::domain_error
  );
  BOOST_CHECK_THROW(
  pdf(
  negative_binomial_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(0.25)),
  static_cast<RealType>(9)), std::domain_error
  );
  BOOST_CHECK_THROW(
  cdf(
  negative_binomial_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(0.25)),
  static_cast<RealType>(-1)), std::domain_error
  );
  BOOST_CHECK_THROW(
  cdf(
  negative_binomial_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(0.25)),
  static_cast<RealType>(9)), std::domain_error
  );
  BOOST_CHECK_THROW(
  cdf(
  negative_binomial_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(-0.25)),
  static_cast<RealType>(0)), std::domain_error
  );
  BOOST_CHECK_THROW(
  cdf(
  negative_binomial_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(1.25)),
  static_cast<RealType>(0)), std::domain_error
  );
  BOOST_CHECK_THROW(
  quantile(
  negative_binomial_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(-0.25)),
  static_cast<RealType>(0)), std::domain_error
  );
  BOOST_CHECK_THROW(
  quantile(
  negative_binomial_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(1.25)),
  static_cast<RealType>(0)), std::domain_error
  );
  BOOST_CHECK_EQUAL(
  cdf(
  negative_binomial_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(0.25)),
  static_cast<RealType>(8)), static_cast<RealType>(1)
  );
  BOOST_CHECK_EQUAL(
  cdf(
  negative_binomial_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(0)),
  static_cast<RealType>(7)), static_cast<RealType>(1)
  );
  BOOST_CHECK_EQUAL(
  cdf(
  negative_binomial_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(1)),
  static_cast<RealType>(7)), static_cast<RealType>(0)
  );
  // End of single checks on exact and check throwing 'duff' out-of-domain values.
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
  //  exp(lgamma(r + k) - lgamma(r) - lgamma(k+1)) * pow(p, r) * pow((1-p), k)
  // with
  // (p/(r+k) * ibeta_derivative(r, k+1, p) (as used in pdf)
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
  // test_spots(0.0F); // Test float.
   test_spots(0.0); // Test double.
  // test_spots(0.0L); // Test long double.
  //#if !BOOST_WORKAROUND(__BORLANDC__, BOOST_TESTED_AT(0x582))
  //  test_spots(boost::math::concepts::real_concept(0.)); // Test real concept.
  //#endif

  return 0;
} // int test_main(int, char* [])

/*

    // Test negative binomial cdf using MathCAD computed values.
    // negative_binomial(k, r, p)
  
    BOOST_CHECK_CLOSE(negative_binomial_distribution<RealType>(dist(2, 0.25)),
    static_cast<RealType>(0), // k. 
    static_cast<RealType>(2), // r
    static_cast<RealType>(0.5)),  // mean
    static_cast<RealType>(0.25),  // probability 1/4
    tol);

    BOOST_CHECK_CLOSE(negative_binomial(
    static_cast<RealType>(1), // k
    static_cast<RealType>(4), // r
    static_cast<RealType>(0.5)), // p
    static_cast<RealType>(0.1875), // 
    tolerance);

    BOOST_CHECK_CLOSE(negative_binomial(
    static_cast<RealType>(1),  // k
    static_cast<RealType>(2), // r
    static_cast<RealType>(0.5)), // pnbinom(1,2,0.5) = 0.5
    static_cast<RealType>(0.5),  //
    tolerance);

    BOOST_CHECK_CLOSE(negative_binomial(
    static_cast<RealType>(2),  // k
    static_cast<RealType>(2), // n
    static_cast<RealType>(0.5)), // pnbinom(1,2,0.5) = 0.5
    static_cast<RealType>(0.6875),  //
    tolerance);
    BOOST_CHECK_CLOSE(negative_binomial(
    static_cast<RealType>(1),  // k
    static_cast<RealType>(3), // n
    static_cast<RealType>(0.5)), // pnbinom(1,3,0.5) = 0.5
    static_cast<RealType>(0.3125),  // 0.6875
    tolerance);

    BOOST_CHECK_CLOSE(negative_binomial(
    static_cast<RealType>(2),  // k
    static_cast<RealType>(2), // n
    static_cast<RealType>(0.1)), // 
    static_cast<RealType>(0.0523),  //
    tolerance);

    BOOST_CHECK_CLOSE(negative_binomial_c(
    static_cast<RealType>(2),  // k
    static_cast<RealType>(2), // n
    static_cast<RealType>(0.1)), // 
    static_cast<RealType>(1- 0.0523),  //
    tolerance);

    BOOST_CHECK_CLOSE(negative_binomial(
    static_cast<RealType>(26), // k
    static_cast<RealType>(5), // n
    static_cast<RealType>(0.4)), // p
    static_cast<RealType>(0.998968624661119), 
    tolerance);                                    

    BOOST_CHECK_CLOSE(negative_binomial(
    static_cast<RealType>(25), // k
    static_cast<RealType>(5), // n
    static_cast<RealType>(0.4)), // p
    static_cast<RealType>(0.998489925933617), 
    tolerance);

    BOOST_CHECK_CLOSE(negative_binomial_c( // complement.
    static_cast<RealType>(26), // k
    static_cast<RealType>(5), // n
    static_cast<RealType>(0.4)), // p
    static_cast<RealType>(1 - 0.998968624661119), 
    tolerance * 100);    // fails at tolerance * 10
    // 0.0010313753388809799 
    */

