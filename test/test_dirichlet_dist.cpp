// test_dirichlet_dist.cpp

// Copyright Mrityunjay Tripathi 2020.

// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

// Basic sanity tests for the Dirichlet Distribution.

#ifdef _MSC_VER
#pragma warning(disable : 4127) // conditional expression is constant.
#pragma warning(disable : 4996) // POSIX name for this item is deprecated.
#pragma warning(disable : 4224) // nonstandard extension used : formal parameter 'arg' was previously defined as a type.
#endif

#include <limits>
#include <iostream>
#include <numeric>
#include <boost/test/unit_test.hpp>                       // for test_main
#include <boost/multiprecision/float128.hpp>
#include <boost/math/concepts/real_concept.hpp>           // for real_concept
#include <boost/math/distributions/dirichlet.hpp>         // for dirichlet_distribution
#include <boost/test/tools/floating_point_comparison.hpp> // for BOOST_CHECK_CLOSE_FRACTION
#include "test_out_of_range.hpp"
#include "math_unit_test.hpp"

using boost::math::dirichlet_distribution;
using boost::math::concepts::real_concept;
using std::numeric_limits;
using std::domain_error;

#define BOOST_TEST_MAIN
#define BOOST_MATH_CHECK_THROW

template <class RandomAccessContainer>
void test_spot(
    RandomAccessContainer &&alpha,                      // concentration parameters 'a'
    RandomAccessContainer &&x,                          // quantiles 'x'
    RandomAccessContainer &&mean,                       // mean
    RandomAccessContainer &&mode,                       // mode
    RandomAccessContainer &&var,                        // variance
    RandomAccessContainer &&skewness,                   // skewness
    RandomAccessContainer &&kurtosis,                   // kurtosis
    typename RandomAccessContainer::value_type entropy, // entropy
    typename RandomAccessContainer::value_type pdf,     // pdf
    typename RandomAccessContainer::value_type cdf,     // cdf
    typename RandomAccessContainer::value_type tol)     // Test tolerance.
{
  // using RealType = typename RandomAccessContainer::value_type;
  typedef RandomAccessContainer V;
  boost::math::dirichlet_distribution<V> diri(std::move(alpha));

  V calc_mean = boost::math::mean(diri);
  V calc_variance = boost::math::variance(diri);
  V calc_mode = boost::math::mode(diri);
  V calc_kurtosis = boost::math::kurtosis(diri);
  V calc_skewness = boost::math::skewness(diri);

  BOOST_CHECK_CLOSE_FRACTION(boost::math::pdf(diri, x), pdf, tol);
  BOOST_CHECK_CLOSE_FRACTION(boost::math::cdf(diri, x), cdf, tol);
  BOOST_CHECK_CLOSE_FRACTION(boost::math::entropy(diri), entropy, tol);

  for (decltype(alpha.size()) i = 0; i < alpha.size(); ++i)
  {
    BOOST_CHECK_CLOSE_FRACTION(calc_mean[i], mean[i], tol);
    BOOST_CHECK_CLOSE_FRACTION(calc_variance[i], var[i], tol);
    BOOST_CHECK_CLOSE_FRACTION(calc_kurtosis[i], kurtosis[i], tol);
    BOOST_CHECK_CLOSE_FRACTION(calc_skewness[i], skewness[i], tol);
  }
} // template <class RandomAccessContainer> void test_spot

template <class RandomAccessContainer>
void test_spots(RandomAccessContainer)
{
  typedef RandomAccessContainer V;
  using RealType = typename V::value_type;
  RealType tolerance = std::max(boost::math::tools::epsilon<RealType>(),
                                static_cast<RealType>(std::numeric_limits<double>::epsilon()));
  V alpha(2);
  V x(2);

  // Error checks:
  // Necessary conditions for instantiation:
  // 1. alpha[i] > 0.
  alpha[0] = 0.35;
  alpha[1] = -1.72; // alpha[1] < 0.
  BOOST_MATH_CHECK_THROW(dirichlet_distribution<V>(std::move(alpha)), std::domain_error);
  BOOST_MATH_CHECK_THROW(dirichlet_distribution<V>(std::move(alpha)), std::domain_error);

  // Domain test for pdf. Necessary conditions for pdf:
  // 1. alpha[i] > 0
  // 2. 0 <= x[i] <=1
  // 3. sum(x) <= 1.
  alpha[0] = -0.2;
  alpha[1] = 1.7; // alpha[0] < 0.
  x[0] = 0.5;
  x[1] = 0.5;
  BOOST_MATH_CHECK_THROW(boost::math::pdf(dirichlet_distribution<V>(std::move(alpha)), x), std::domain_error);

  alpha[0] = 1.36;
  alpha[1] = 0.0; // alpha[1] = 0.
  x[0] = 0.47;
  x[1] = 0.53;
  BOOST_MATH_CHECK_THROW(boost::math::pdf(dirichlet_distribution<V>(std::move(alpha)), x), std::domain_error);

  alpha[0] = 1.26;
  alpha[1] = 2.99;
  x[0] = 0.5;
  x[1] = 0.75; // sum(x) > 1.0
  BOOST_MATH_CHECK_THROW(boost::math::pdf(dirichlet_distribution<V>(std::move(alpha)), x), std::domain_error);

  alpha[0] = 1.56;
  alpha[1] = 4.00;
  x[0] = 0.31;
  x[1] = -0.03; // x[1] < 0.
  BOOST_MATH_CHECK_THROW(boost::math::pdf(dirichlet_distribution<V>(std::move(alpha)), x), std::domain_error);

  alpha[0] = 1.56;
  alpha[1] = 4.00;
  x[0] = 0.31;
  x[1] = 1.06; // x[1] > 1.
  BOOST_MATH_CHECK_THROW(boost::math::pdf(dirichlet_distribution<V>(std::move(alpha)), x), std::domain_error);

  // Domain test for cdf. Necessary conditions for cdf:
  // 1. alpha[i] > 0
  // 2. 0 <= x[i] <= 1
  // 3. sum(x) <= 1.
  alpha[0] = 1.56;
  alpha[1] = 4.00;
  x[0] = 0.31;
  x[1] = 1.06; // x[1] > 1.
  BOOST_MATH_CHECK_THROW(boost::math::cdf(dirichlet_distribution<V>(std::move(alpha)), x), std::domain_error);

  alpha[0] = 3.756;
  alpha[1] = 4.91;
  x[0] = 0.31;
  x[1] = -1.06; // x[1] < 0.
  BOOST_MATH_CHECK_THROW(boost::math::cdf(dirichlet_distribution<V>(std::move(alpha)), x), std::domain_error);

  alpha[0] = 1.56;
  alpha[1] = -4.00; // alpha[1] < 0
  x[0] = 0.31;
  x[1] = 0.69;
  BOOST_MATH_CHECK_THROW(boost::math::cdf(dirichlet_distribution<V>(std::move(alpha)), x), std::domain_error);

  alpha[0] = 0.0;
  alpha[1] = 4.00; // alpha[0] = 0.
  x[0] = 0.25;
  x[1] = 0.75;
  BOOST_MATH_CHECK_THROW(boost::math::cdf(dirichlet_distribution<V>(std::move(alpha)), x), std::domain_error);

  alpha[0] = 1.56;
  alpha[1] = 4.00;
  x[0] = 0.31;
  x[1] = 0.71; // sum(x) > 1.
  BOOST_MATH_CHECK_THROW(boost::math::cdf(dirichlet_distribution<V>(std::move(alpha)), x), std::domain_error);

  // Domain test for mode. Necessary conditions for mode:
  // 1. alpha[i] > 1.
  alpha[0] = 1.0;
  alpha[1] = 1.4; // alpha[0] = 1.
  BOOST_MATH_CHECK_THROW(boost::math::mode(dirichlet_distribution<V>(std::move(alpha))), std::domain_error);

  alpha[0] = 1.56;
  alpha[1] = 0.92; // alpha[1] < 1.
  BOOST_MATH_CHECK_THROW(boost::math::mode(dirichlet_distribution<V>(std::move(alpha))), std::domain_error);

  // Some exact values of pdf.
  alpha[0] = 1.0, alpha[1] = 1.0;
  x[0] = 0.5, x[1] = 0.5;
  BOOST_CHECK_EQUAL(boost::math::pdf(dirichlet_distribution<V>(std::move(alpha)), x), static_cast<RealType>(1.0));

  alpha[0] = 2.0, alpha[1] = 2.0;
  x[0] = 0.5, x[1] = 0.5;
  BOOST_CHECK_EQUAL(boost::math::pdf(dirichlet_distribution<V>(std::move(alpha)), x), static_cast<RealType>(1.5));

  // Checking precalculated values on scipy.
  // https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.dirichlet.html
  alpha[0] = static_cast<RealType>(5.778238829L); alpha[1] = static_cast<RealType>(2.55821892973L);
  x[0] = static_cast<RealType>(0.23667289213L); x[1] = static_cast<RealType>(1.0L - 0.23667289213L);
  V mean = {static_cast<RealType>(0.69312878L), static_cast<RealType>(0.30687122L)};
  V mode = {static_cast<RealType>(0.75408675L), static_cast<RealType>(0.24591325L)};
  V var = {static_cast<RealType>(0.0227818L), static_cast<RealType>(0.0227818L)};
  V skewness = {static_cast<RealType>(-0.49515513L), static_cast<RealType>(0.49515513L)};
  V kurtosis = {static_cast<RealType>(-139231.64518864L), static_cast<RealType>(-6993.41057616L)};
  RealType entropy = static_cast<RealType>(17.6747L);
  RealType pdf = static_cast<RealType>(0.05866154L);
  RealType cdf = static_cast<RealType>(0.00071693L);
  tolerance *= 1E+07;

  test_spot(std::move(alpha),
    std::move(x),
    std::move(mean),
    std::move(mode),
    std::move(var),
    std::move(skewness),
    std::move(kurtosis),
    entropy, pdf, cdf, tolerance);

  // No longer allow any parameter to be NaN or inf.
  if (std::numeric_limits<RealType>::has_quiet_NaN)
  {
    RealType not_a_num = std::numeric_limits<RealType>::quiet_NaN();
    alpha[0] = not_a_num; alpha[1] = 0.37;
#ifndef BOOST_NO_EXCEPTIONS
    BOOST_MATH_CHECK_THROW(dirichlet_distribution<V>(std::move(alpha)), std::domain_error);
#else
    BOOST_MATH_CHECK_THROW(dirichlet_distribution<V>(std::move(alpha)), std::domain_error);
#endif

    // Non-finite parameters should throw.
    alpha[0] = 1.67; alpha[1] = 3.8;
    x[0] = not_a_num; x[1] = 0.5;
    dirichlet_distribution<V> w(std::move(alpha));
    BOOST_MATH_CHECK_THROW(boost::math::pdf(w, x), std::domain_error); // x = NaN
    BOOST_MATH_CHECK_THROW(boost::math::cdf(w, x), std::domain_error); // x = NaN
  }                                                          // has_quiet_NaN

  if (std::numeric_limits<RealType>::has_infinity)
  {
    // Attempt to construct from non-finite should throw.
    RealType infinite = std::numeric_limits<RealType>::infinity();
    alpha[0] = infinite;
    alpha[1] = 7.2;
#ifndef BOOST_NO_EXCEPTIONS
    BOOST_MATH_CHECK_THROW(dirichlet_distribution<V>(std::move(alpha)), std::domain_error);
#else
    BOOST_MATH_CHECK_THROW(dirichlet_distribution<V>(std::move(alpha)), std::domain_error);
#endif
    alpha[0] = 1.42; alpha[1] = 7.91;
    x[0] = 0.25; x[1] = infinite;
    dirichlet_distribution<V> w(std::move(alpha));
    BOOST_MATH_CHECK_THROW(boost::math::pdf(w, x), std::domain_error); // x = inf
    BOOST_MATH_CHECK_THROW(boost::math::cdf(w, x), std::domain_error); // x = inf
    x[1] = -infinite;
    BOOST_MATH_CHECK_THROW(boost::math::pdf(w, x), std::domain_error); // x = -inf
    BOOST_MATH_CHECK_THROW(boost::math::cdf(w, x), std::domain_error); // x = -inf
  }
} // test_spots()


BOOST_AUTO_TEST_CASE(test_main)
{
  BOOST_MATH_CONTROL_FP;
  test_spots(std::vector<long double>(0.0L));

  test_spots(std::vector<double>(0.0));

  test_spots(std::vector<float>(0.0F));

// #ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
//   test_spots(); // Test long double.
// #if !BOOST_WORKAROUND(__BORLANDC__, BOOST_TESTED_AT(0x582))
//   test_spots(boost::math::concepts::real_concept(0.)); // Test real concept.
// #endif
} // BOOST_AUTO_TEST_CASE( test_main )