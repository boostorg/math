// test_negative_binomial.cpp

// Copyright Paul A. Bristow 2007.
// Copyright John Maddock 2006.

// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

// Tests for Negative Binomial Distribution.

// Note that these defines must be placed BEFORE #includes.
#define BOOST_MATH_OVERFLOW_ERROR_POLICY ignore_error
// because several tests overflow & underflow by design.
#define BOOST_MATH_DISCRETE_QUANTILE_POLICY real

#ifdef _MSC_VER
#  pragma warning(disable: 4127) // conditional expression is constant.
#endif

#if !defined(TEST_FLOAT) && !defined(TEST_DOUBLE) && !defined(TEST_LDOUBLE) && !defined(TEST_REAL_CONCEPT)
#  define TEST_FLOAT
#  define TEST_DOUBLE
#  define TEST_LDOUBLE
#  define TEST_REAL_CONCEPT
#endif

#ifdef BOOST_MATH_ENABLE_SYCL
#include "sycl/sycl.hpp"
#endif

#include <boost/math/tools/config.hpp>
#include "math_unit_test.hpp"

#ifndef BOOST_MATH_NO_REAL_CONCEPT_TESTS
#include <boost/math/concepts/real_concept.hpp> // for real_concept
#endif

#include <boost/math/distributions/negative_binomial.hpp>
#include <boost/math/tools/precision.hpp>
#include "table_type.hpp"
#include "test_out_of_range.hpp"

#include <array>
#include <limits>
#include <type_traits>

using boost::math::negative_binomial_distribution;

// The Boost.Test checks these replaced took tolerances in percent; the tolerances below are
// fractions, the same values divided by 100.

template <class RealType>
void test_spot( // Test a single spot value against 'known good' values.
               RealType N,    // Number of successes.
               RealType k,    // Number of failures.
               RealType p,    // Probability of success_fraction.
               RealType P,    // CDF probability.
               RealType Q,    // Complement of CDF.
               RealType tol)  // Test tolerance.
{
  negative_binomial_distribution<RealType> bn(N, p);
  CHECK_EQUAL(N, bn.successes());
  CHECK_EQUAL(p, bn.success_fraction());
  CHECK_RELATIVE_ERROR(P, cdf(bn, k), tol);

  if((P < 0.99) && (Q < 0.99))
  {
    // We can only check this if P is not too close to 1,
    // so that we can guarantee that Q is free of error:
    CHECK_RELATIVE_ERROR(Q, cdf(complement(bn, k)), tol);
    // Limit where a zero quantile is checked: if exponent range is very large we may
    // run out of iterations in our root finding algorithm.
    bool check_zero = (std::numeric_limits<RealType>::max_exponent <= std::numeric_limits<double>::max_exponent)
        && std::is_floating_point<RealType>::value;
    if(k != 0)
    {
      CHECK_RELATIVE_ERROR(k, quantile(bn, P), tol);
      CHECK_RELATIVE_ERROR(k, quantile(complement(bn, Q)), tol);
    }
    else if(check_zero)
    {
      // Just check quantile is very small:
      CHECK_LE(quantile(bn, P), boost::math::tools::epsilon<RealType>() * 10);
      CHECK_LE(quantile(complement(bn, Q)), boost::math::tools::epsilon<RealType>() * 10);
    }
    // estimate success ratio:
    CHECK_RELATIVE_ERROR(p, negative_binomial_distribution<RealType>::find_lower_bound_on_p(N+k, N, P), tol);
    // Note we bump up the sample size here, purely for the sake of the test,
    // internally the function has to adjust the sample size so that we get
    // the right upper bound, our test undoes this, so we can verify the result.
    CHECK_RELATIVE_ERROR(p, negative_binomial_distribution<RealType>::find_upper_bound_on_p(N+k+1, N, Q), tol);

    // Check that the upper and lower bounds are the right way around, and that they
    // bracket the naive estimate of p = successes / (sample size); with Q when Q < P, else with P:
    RealType alpha = (Q < P) ? Q : P;
    RealType lower = negative_binomial_distribution<RealType>::find_lower_bound_on_p(N+k, N, alpha);
    RealType upper = negative_binomial_distribution<RealType>::find_upper_bound_on_p(N+k, N, alpha);
    CHECK_LE(lower, upper);
    CHECK_LE(lower, RealType(N / (N+k)));
    CHECK_LE(RealType(N / (N+k)), upper);

    // Estimate sample size:
    CHECK_RELATIVE_ERROR(N+k, negative_binomial_distribution<RealType>::find_minimum_number_of_trials(k, p, P), tol);
    CHECK_RELATIVE_ERROR(N+k, negative_binomial_distribution<RealType>::find_maximum_number_of_trials(k, p, Q), tol);

    // Double check consistency of CDF and PDF by computing the finite sum:
    RealType sum = 0;
    for(unsigned i = 0; i <= k; ++i)
    {
      sum += pdf(bn, RealType(i));
    }
    CHECK_RELATIVE_ERROR(P, sum, tol);
    // Complement is not possible since sum is to infinity.
  }
} // test_spot

// The relative condition number |P / (x dP/dx)| of the quantile x of P, from a central difference
// of the cdf: the quantile's relative error is the cdf's times this, which is large near x = 0.
template <class Dist, class RealType>
RealType quantile_condition_number(Dist const& d, RealType x, RealType P, bool complemented)
{
  using std::abs;
  using std::pow;
  if(x == 0)
  {
    return 1;
  }
  // Not tools::cbrt_epsilon, whose function-local static a SYCL kernel can't use:
  RealType h = x * pow(boost::math::tools::epsilon<RealType>(), RealType(1) / 3);
  RealType up = complemented ? cdf(complement(d, RealType(x + h))) : cdf(d, RealType(x + h));
  RealType down = complemented ? cdf(complement(d, RealType(x - h))) : cdf(d, RealType(x - h));
  RealType kappa = abs(P * 2 * h / (x * (up - down)));
  return kappa > 1 ? kappa : RealType(1);
}

template <class RealType> // Any floating-point type RealType.
void test_spots(RealType)
{
  // Basic sanity checks, test data is to double precision only
  // so set tolerance to 100000 eps, or 100000 eps of type double,
  // whichever is the larger.
  RealType tolerance = (std::max)
    (boost::math::tools::epsilon<RealType>(),
    static_cast<RealType>(std::numeric_limits<double>::epsilon()));
  tolerance *= 100000.0f;
  // BOOST_CHECK_SMALL took the percentage tolerance above as an absolute bound:
  RealType small_bound = tolerance * 100;

  RealType tol1eps = boost::math::tools::epsilon<RealType>() * 2 / 100; // Very tight, suit exact values.
  RealType tol5eps = boost::math::tools::epsilon<RealType>() * 5 / 100; // Wider 5 epsilon.

  // Sources of spot test values:

  // MathCAD defines pbinom(k, r, p) (at about 64-bit double precision, about 16 decimal digits)
  // returns pr(X , k) when random variable X has the binomial distribution with parameters r and p.
  // 0 <= k
  // r > 0
  // 0 <= p <= 1
  // P = pbinom(30, 500, 0.05) = 0.869147702104609

  // And functions.wolfram.com

  using ::boost::math::cdf;
  using ::boost::math::pdf;

  // Test negative binomial using cdf spot values from MathCAD cdf = pnbinom(k, r, p).
  // These test quantiles and complements as well.
  // Arguments: successes r, failures k, success fraction p, CDF P, complement Q = 1 - P.

  test_spot(  // pnbinom(1,2,0.5) = 0.5
  static_cast<RealType>(2), static_cast<RealType>(1), static_cast<RealType>(0.5),
  static_cast<RealType>(0.5), static_cast<RealType>(0.5), tolerance);

  test_spot( // pbinom(0, 2, 0.25)
  static_cast<RealType>(2), static_cast<RealType>(0), static_cast<RealType>(0.25),
  static_cast<RealType>(0.0625), static_cast<RealType>(0.9375), tolerance);

  test_spot(  // pbinom(48,8,0.25)
  static_cast<RealType>(8), static_cast<RealType>(48), static_cast<RealType>(0.25),
  static_cast<RealType>(9.826582228110670E-1), static_cast<RealType>(1 - 9.826582228110670E-1), tolerance);

  test_spot(  // pbinom(2,5,0.4)
  static_cast<RealType>(5), static_cast<RealType>(2), static_cast<RealType>(0.4),
  static_cast<RealType>(9.625600000000020E-2), static_cast<RealType>(1 - 9.625600000000020E-2), tolerance);

  test_spot(  // pbinom(10,100,0.9)
  static_cast<RealType>(100), static_cast<RealType>(10), static_cast<RealType>(0.9),
  static_cast<RealType>(4.535522887695670E-1), static_cast<RealType>(1 - 4.535522887695670E-1), tolerance);

  test_spot(  // pbinom(1,100,0.991)
  static_cast<RealType>(100), static_cast<RealType>(1), static_cast<RealType>(0.991),
  static_cast<RealType>(7.693413044217000E-1), static_cast<RealType>(1 - 7.693413044217000E-1), tolerance);

  test_spot(  // pbinom(10,100,0.991)
  static_cast<RealType>(100), static_cast<RealType>(10), static_cast<RealType>(0.991),
  static_cast<RealType>(9.999999940939000E-1), static_cast<RealType>(1 - 9.999999940939000E-1), tolerance);

  if(std::numeric_limits<RealType>::is_specialized)
  { // An extreme value test that takes 3 minutes using the real concept type
    // for which numeric_limits<RealType>::is_specialized == false, deliberately
    // and for which there is no Lanczos approximation defined (also deliberately)
    // giving a very slow computation, but with acceptable accuracy.
    // A possible enhancement might be to use a normal approximation for
    // extreme values, but this is not implemented.
    test_spot(  // pbinom(100000,100,0.001)
    static_cast<RealType>(100), static_cast<RealType>(100000), static_cast<RealType>(0.001),
    static_cast<RealType>(5.173047534260320E-1), static_cast<RealType>(1 - 5.173047534260320E-1),
    tolerance*1000); // *1000 is OK 0.51730475350664229  versus

    // functions.wolfram.com
    //   for I[0.001](100, 100000+1) gives:
    // Wolfram       0.517304753506834882009032744488738352004003696396461766326713
    // JM nonLanczos 0.51730475350664229 differs at the 13th decimal digit.
    // MathCAD       0.51730475342603199 differs at 10th decimal digit.

    // Error tests:
    check_out_of_range<negative_binomial_distribution<RealType> >(20, 0.5);
    CHECK_THROW(negative_binomial_distribution<RealType>(0, 0.5), std::domain_error);
    CHECK_THROW(negative_binomial_distribution<RealType>(-2, 0.5), std::domain_error);
    CHECK_THROW(negative_binomial_distribution<RealType>(20, -0.5), std::domain_error);
    CHECK_THROW(negative_binomial_distribution<RealType>(20, 1.5), std::domain_error);
  }
  // End of single spot tests using RealType

  auto nb = [](double r, double p) { return negative_binomial_distribution<RealType>(static_cast<RealType>(r), static_cast<RealType>(p)); };

  // Tests on PDF at k = 0:
  CHECK_RELATIVE_ERROR(static_cast<RealType>(0.25), pdf(nb(2, 0.5), static_cast<RealType>(0)), tolerance);
  CHECK_RELATIVE_ERROR(static_cast<RealType>(0.0625), pdf(nb(4, 0.5), static_cast<RealType>(0)), tolerance); // exact 1/16
  CHECK_RELATIVE_ERROR(static_cast<RealType>(9.094947017729270E-13), pdf(nb(20, 0.25), static_cast<RealType>(0)), tolerance); // pbinom(0,20,0.25)
  CHECK_RELATIVE_ERROR(static_cast<RealType>(1.0485760000000003e-014), pdf(nb(20, 0.2), static_cast<RealType>(0)), tolerance); // MathCAD 1.048576000000000E-14
  CHECK_RELATIVE_ERROR(static_cast<RealType>(1e-10), pdf(nb(10, 0.1), static_cast<RealType>(0)), tolerance); // MathCAD says zero, but suffers cancellation error?
  CHECK_RELATIVE_ERROR(static_cast<RealType>(1e-20), pdf(nb(20, 0.1), static_cast<RealType>(0)), tolerance); // MathCAD says zero, but suffers cancellation error?
  CHECK_RELATIVE_ERROR(static_cast<RealType>(1.215766545905690E-1), pdf(nb(20, 0.9), static_cast<RealType>(0)), tolerance);

  // Tests on cdf:
  // MathCAD pbinom k, r, p) == failures, successes, probability.
  CHECK_RELATIVE_ERROR(static_cast<RealType>(0.25), cdf(nb(2, 0.5), static_cast<RealType>(0)), tolerance);
  CHECK_RELATIVE_ERROR(static_cast<RealType>(0.75), cdf(complement(nb(2, 0.5), static_cast<RealType>(0))), tolerance);
  CHECK_RELATIVE_ERROR(static_cast<RealType>(1.455191522836700E-11), cdf(nb(20, 0.25), static_cast<RealType>(1)), tolerance);
  CHECK_ABSOLUTE_ERROR(static_cast<RealType>(1.455191522836700E-11), cdf(nb(20, 0.25), static_cast<RealType>(1)), small_bound);

  // Some exact (probably - judging by trailing zeros) values.
  CHECK_RELATIVE_ERROR(static_cast<RealType>(1.525878906250000E-5), cdf(nb(8, 0.25), static_cast<RealType>(0)), tolerance);
  CHECK_ABSOLUTE_ERROR(static_cast<RealType>(1.525878906250000E-5), cdf(nb(8, 0.25), static_cast<RealType>(0)), small_bound);
  const double eight_quarter[][2] = {
    {1, 1.068115234375010E-4}, {2, 4.158020019531300E-4}, {3, 1.188278198242200E-3}, {4, 2.781510353088410E-3},
    {5, 5.649328231811500E-3}, {6, 1.030953228473680E-2}, {7, 1.729983836412430E-2}, {8, 2.712995628826370E-2},
    {48, 9.826582228110670E-1}, {64, 9.990295004935590E-1} };
  for(auto const& row : eight_quarter)
  {
    CHECK_RELATIVE_ERROR(static_cast<RealType>(row[1]), cdf(nb(8, 0.25), static_cast<RealType>(row[0])), tolerance);
  }
  CHECK_RELATIVE_ERROR(static_cast<RealType>(9.989686246611190E-1), cdf(nb(5, 0.4), static_cast<RealType>(26)), tolerance);
  CHECK_RELATIVE_ERROR(static_cast<RealType>(9.625600000000020E-2), cdf(nb(5, 0.4), static_cast<RealType>(2)), tolerance);
  CHECK_RELATIVE_ERROR(static_cast<RealType>(9.999970854144170E-1), cdf(nb(50, 0.9), static_cast<RealType>(20)), tolerance);
  CHECK_RELATIVE_ERROR(static_cast<RealType>(2.172846379930550E-1), cdf(nb(500, 0.7), static_cast<RealType>(200)), tolerance * 2);
  CHECK_RELATIVE_ERROR(static_cast<RealType>(4.550203671301790E-1), cdf(nb(50, 0.7), static_cast<RealType>(20)), tolerance);

  // Tests of other functions, mean and other moments ...

  negative_binomial_distribution<RealType> dist = nb(8, 0.25);
  using namespace std; // ADL of std names.
  CHECK_RELATIVE_ERROR(static_cast<RealType>(8 * (1 - 0.25) /0.25), mean(dist), tol5eps);
  CHECK_RELATIVE_ERROR(static_cast<RealType>(21), mode(dist), tol1eps);
  CHECK_RELATIVE_ERROR(static_cast<RealType>(8 * (1 - 0.25) / (0.25 * 0.25)), variance(dist), tol5eps);
  CHECK_RELATIVE_ERROR(static_cast<RealType>(9.797958971132712392789136298823565567864L), standard_deviation(dist), tol5eps * 100); // using functions.wolfram.com
  CHECK_RELATIVE_ERROR(static_cast<RealType>(0.71443450831176036), skewness(dist), tolerance); // using http://mathworld.wolfram.com/skewness.html
  CHECK_RELATIVE_ERROR(static_cast<RealType>(0.7604166666666666666666666666666666666666L), kurtosis_excess(dist), tol5eps * 100); // using Wikipedia Kurtosis(excess) formula
  CHECK_RELATIVE_ERROR(static_cast<RealType>(3.76041666666666666666666666666666666666666L), kurtosis(dist), tol5eps * 100);
  RealType x = static_cast<RealType>(0.125);
  CHECK_RELATIVE_ERROR(pdf(dist, x) / cdf(complement(dist, x)), hazard(dist, x), tol5eps);
  CHECK_RELATIVE_ERROR(-log(cdf(complement(dist, x))), chf(dist, x), tol5eps);
  CHECK_RELATIVE_ERROR(standard_deviation(dist) / mean(dist), coefficient_of_variation(dist), tol5eps);

  // Special cases for PDF:
  CHECK_EQUAL(pdf(nb(8, 0), static_cast<RealType>(0)), static_cast<RealType>(0));
  CHECK_EQUAL(pdf(nb(8, 0), static_cast<RealType>(0.0001)), static_cast<RealType>(0));
  CHECK_EQUAL(pdf(nb(8, 1), static_cast<RealType>(0.001)), static_cast<RealType>(0));
  CHECK_EQUAL(pdf(nb(8, 1), static_cast<RealType>(8)), static_cast<RealType>(0));
  // Expect exact, but not quite:
  CHECK_ABSOLUTE_ERROR(static_cast<RealType>(0.0625), pdf(nb(2, 0.25), static_cast<RealType>(0)), RealType(2 * boost::math::tools::epsilon<RealType>()));

  // Quantile boundary cases checks:
  CHECK_EQUAL(quantile(dist, static_cast<RealType>(0)), static_cast<RealType>(0)); // zero P < cdf(0) so should be exactly zero.
  CHECK_EQUAL(quantile(dist, boost::math::tools::min_value<RealType>()), static_cast<RealType>(0)); // min P < cdf(0)
  // Small P < cdf(0) so should be near zero; a relative check against zero requires zero:
  CHECK_RELATIVE_ERROR(static_cast<RealType>(0), quantile(dist, boost::math::tools::epsilon<RealType>()), tol5eps);
  CHECK_RELATIVE_ERROR(static_cast<RealType>(0.95854156929288470), quantile(dist, static_cast<RealType>(0.0001)), tolerance);

  if(std::numeric_limits<RealType>::has_infinity)
  { // Infinity is not implemented for real_concept, so these tests are only done for built-in types.
    // These assume that BOOST_MATH_OVERFLOW_ERROR_POLICY is NOT throw_on_error.
    CHECK_EQUAL(quantile(dist, static_cast<RealType>(1)), std::numeric_limits<RealType>::infinity()); // At P == 1 so k failures should be infinite.
    CHECK_EQUAL(quantile(complement(dist, static_cast<RealType>(0))), std::numeric_limits<RealType>::infinity()); // Q zero so P == 1.
  }
  else
  { // real_concept case, which returns max_value rather than infinity.
    CHECK_EQUAL(quantile(dist, static_cast<RealType>(1)), boost::math::tools::max_value<RealType>());
    CHECK_EQUAL(quantile(complement(dist, static_cast<RealType>(0))), boost::math::tools::max_value<RealType>());
  }
  // Q very near to 1 so P nearly 1 so should be large > 384:
  CHECK_GE(quantile(complement(dist, boost::math::tools::min_value<RealType>())), static_cast<RealType>(384));
  CHECK_EQUAL(quantile(dist, static_cast<RealType>(0)), static_cast<RealType>(0)); // P == 0 < cdf(0) so should be zero.

  // Quantile Complement boundary cases:
  CHECK_EQUAL(quantile(complement(dist, static_cast<RealType>(1))), static_cast<RealType>(0)); // Q = 1 so P = 0 < cdf(0)
  CHECK_EQUAL(quantile(complement(dist, static_cast<RealType>(1 - boost::math::tools::epsilon<RealType>()))), static_cast<RealType>(0)); // P == epsilon < cdf(0)

  // Check that duff arguments throw domain_error:
  CHECK_THROW(pdf(nb(-1, 0.25), static_cast<RealType>(0)), std::domain_error); // Negative successes!
  CHECK_THROW(pdf(nb(8, -0.25), static_cast<RealType>(0)), std::domain_error); // Negative success_fraction!
  CHECK_THROW(pdf(nb(8, 1.25), static_cast<RealType>(0)), std::domain_error); // Success_fraction > 1!
  CHECK_THROW(pdf(nb(8, 0.25), static_cast<RealType>(-1)), std::domain_error); // Negative k argument!
  // Unlike binomial there is NO limit on k (failures).
  CHECK_THROW(cdf(nb(8, 0.25), static_cast<RealType>(-1)), std::domain_error); // Negative k argument!
  CHECK_THROW(cdf(nb(8, -0.25), static_cast<RealType>(0)), std::domain_error); // Negative success_fraction!
  CHECK_THROW(cdf(nb(8, 1.25), static_cast<RealType>(0)), std::domain_error); // Success_fraction > 1!
  CHECK_THROW(quantile(nb(8, -0.25), static_cast<RealType>(0)), std::domain_error); // Negative success_fraction!
  CHECK_THROW(quantile(nb(8, 1.25), static_cast<RealType>(0)), std::domain_error); // Success_fraction > 1!
  // End of check throwing 'duff' out-of-domain values.

// The table's SC_ uses BOOST_JOIN, from Boost.Config, which a standalone build doesn't have:
#ifndef BOOST_JOIN
#define BOOST_JOIN(X, Y) BOOST_MATH_JOIN(X, Y)
#endif
#define T RealType
#include "negative_binomial_quantile.ipp"

  for(auto const& row : negative_binomial_quantile_data)
  {
     using namespace boost::math::policies;
     typedef policy<discrete_quantile<boost::math::policies::real> > P1;
     typedef policy<discrete_quantile<integer_round_down> > P2;
     typedef policy<discrete_quantile<integer_round_up> > P3;
     typedef policy<discrete_quantile<integer_round_outwards> > P4;
     typedef policy<discrete_quantile<integer_round_inwards> > P5;
     typedef policy<discrete_quantile<integer_round_nearest> > P6;
     RealType r = static_cast<RealType>(row[0]);
     RealType sf = static_cast<RealType>(row[1]);
     RealType p = static_cast<RealType>(row[2]);
     RealType lower = static_cast<RealType>(row[3]);
     RealType upper = static_cast<RealType>(row[4]);
     RealType tol = boost::math::tools::epsilon<RealType>() * 700;
     if(!std::is_floating_point<RealType>::value)
        tol *= 10;  // no lanczos approximation implies less accuracy
     //
     // Check full real value first, to the cdf's tolerance times the quantile's condition number
     // (548 for r = 18.3, p = 0.906, where the complement's quantile is 0.0061):
     //
     negative_binomial_distribution<RealType, P1> p1(r, sf);
     CHECK_RELATIVE_ERROR(lower, quantile(p1, p), tol * quantile_condition_number(p1, lower, p, false));
     CHECK_RELATIVE_ERROR(upper, quantile(complement(p1, p)), tol * quantile_condition_number(p1, upper, p, true));
     //
     // Now with round down to integer:
     //
     negative_binomial_distribution<RealType, P2> p2(r, sf);
     CHECK_EQUAL(quantile(p2, p), RealType(floor(lower)));
     CHECK_EQUAL(quantile(complement(p2, p)), RealType(floor(upper)));
     //
     // Now with round up to integer:
     //
     negative_binomial_distribution<RealType, P3> p3(r, sf);
     CHECK_EQUAL(quantile(p3, p), RealType(ceil(lower)));
     CHECK_EQUAL(quantile(complement(p3, p)), RealType(ceil(upper)));
     //
     // Now with round to integer "outside":
     //
     negative_binomial_distribution<RealType, P4> p4(r, sf);
     CHECK_EQUAL(quantile(p4, p), RealType(p < 0.5f ? floor(lower) : ceil(lower)));
     CHECK_EQUAL(quantile(complement(p4, p)), RealType(p < 0.5f ? ceil(upper) : floor(upper)));
     //
     // Now with round to integer "inside":
     //
     negative_binomial_distribution<RealType, P5> p5(r, sf);
     CHECK_EQUAL(quantile(p5, p), RealType(p < 0.5f ? ceil(lower) : floor(lower)));
     CHECK_EQUAL(quantile(complement(p5, p)), RealType(p < 0.5f ? floor(upper) : ceil(upper)));
     //
     // Now with round to nearest integer:
     //
     negative_binomial_distribution<RealType, P6> p6(r, sf);
     CHECK_EQUAL(quantile(p6, p), RealType(floor(lower + 0.5f)));
     CHECK_EQUAL(quantile(complement(p6, p)), RealType(floor(upper + 0.5f)));
  }
#undef T
} // template <class RealType> void test_spots(RealType) // Any floating-point type RealType.

template <class RealType>
void test_entropy(RealType)
{
   // Reference values from mpmath at 60 digits: closed forms, or -sum p log(p) for the discrete
   // distributions, which have none. Parameters are exact in every type.
   RealType tol = 16 * boost::math::tools::epsilon<RealType>();
   auto check = [&](RealType a, RealType b, RealType h, RealType tolerance)
   {
      CHECK_RELATIVE_ERROR(h, entropy(negative_binomial_distribution<RealType>(a, b)), tolerance);
   };
   check(static_cast<RealType>(5.0L), static_cast<RealType>(0.25L), static_cast<RealType>(3.39359238247812518456836017160705370313434987L), tol);
   check(static_cast<RealType>(0.5L), static_cast<RealType>(0.5L), static_cast<RealType>(0.943370759269334533806787413559674682998205284L), tol);
   check(static_cast<RealType>(1000.5L), static_cast<RealType>(0.0625L), static_cast<RealType>(7.61305197413188991442493027236853895416590232L), tol);
   check(static_cast<RealType>(1000000.0L), static_cast<RealType>(0.5L), static_cast<RealType>(8.67326702746660536495437769998638212050453356L), tol);
   // About 1e4 terms, summed with compensation:
   check(static_cast<RealType>(0.125L), static_cast<RealType>(0.00390625L), static_cast<RealType>(3.05357381884703876385483729149700634949203959L), tol);
   check(static_cast<RealType>(2.5L), static_cast<RealType>(0.001953125L), static_cast<RealType>(7.96729439209919695998955366089632316751290691L), tol);
   check(static_cast<RealType>(100000000.0L), static_cast<RealType>(0.9990234375L), static_cast<RealType>(7.16451918702377261816611719177425949753810124L), tol);
   CHECK_EQUAL(entropy(negative_binomial_distribution<RealType>(3, 1)), RealType(0));
}

int main()
{
  // Check that can generate negative_binomial distribution using the two convenience methods:
  using namespace boost::math;
  negative_binomial mynb1(2., 0.5); // Using typedef - default type is double.
  negative_binomial_distribution<> myf2(2., 0.5); // Using default RealType double.
  static_cast<void>(mynb1);
  static_cast<void>(myf2);

  // Check accessor functions return exact values for double at least.
  negative_binomial_distribution<double> my8dist(8., 0.25);
  CHECK_EQUAL(my8dist.successes(), static_cast<double>(8));
  CHECK_EQUAL(my8dist.success_fraction(), static_cast<double>(1./4.));

  // (Parameter value, arbitrarily zero, only communicates the floating point type).
#ifdef TEST_FLOAT
  test_spots(0.0F);
  test_entropy(0.0F);
#endif
#ifdef TEST_DOUBLE
  test_spots(0.0);
  test_entropy(0.0);
#endif
#ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
#ifdef TEST_LDOUBLE
  test_spots(0.0L);
  test_entropy(0.0L);
#endif
#ifndef BOOST_MATH_NO_REAL_CONCEPT_TESTS
#ifdef TEST_REAL_CONCEPT
  test_spots(boost::math::concepts::real_concept(0.));
  test_entropy(boost::math::concepts::real_concept(0.));
#endif
#endif
#endif
  return boost::math::test::report_errors();
}
