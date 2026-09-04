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

#define BOOST_TEST_MAIN
#define BOOST_MATH_CHECK_THROW
#define BOOST_LIB_DIAGNOSTIC
#define BOOST_TEST_MODULE

#include <limits>
#include <iostream>
#include <boost/test/unit_test.hpp> // for test_main
#include <boost/multiprecision/float128.hpp>
#include <boost/math/concepts/real_concept.hpp>           // for real_concept
#include <boost/math/distributions/dirichlet.hpp>         // for dirichlet_distribution
#include <boost/test/tools/floating_point_comparison.hpp> // for BOOST_CHECK_CLOSE_FRACTION
#include "test_out_of_range.hpp"
#include "math_unit_test.hpp"

using boost::math::dirichlet_distribution;
using boost::math::concepts::real_concept;
using std::domain_error;
using std::numeric_limits;

template <class RandomAccessContainer>
void test_spot(
    RandomAccessContainer &&alpha,                      // concentration parameters 'a'
    RandomAccessContainer &&x,                          // quantiles 'x'
    RandomAccessContainer &&mean,                       // mean
    RandomAccessContainer &&var,                        // variance
    typename RandomAccessContainer::value_type entropy, // entropy
    typename RandomAccessContainer::value_type pdf,     // pdf
    typename RandomAccessContainer::value_type tol)     // Test tolerance.
{
  // using RealType = typename RandomAccessContainer::value_type;
  typedef RandomAccessContainer V;
  boost::math::dirichlet_distribution<V> diri(std::move(alpha));

  V calc_mean = boost::math::mean(diri);
  V calc_variance = boost::math::variance(diri);

  BOOST_CHECK_CLOSE_FRACTION(boost::math::pdf(diri, x), pdf, tol);
  BOOST_CHECK_CLOSE_FRACTION(boost::math::entropy(diri), entropy, tol);

  for (decltype(alpha.size()) i = 0; i < alpha.size(); ++i)
  {
    BOOST_CHECK_CLOSE_FRACTION(calc_mean[i], mean[i], tol);
    BOOST_CHECK_CLOSE_FRACTION(calc_variance[i], var[i], tol);
  }
} // template <class RandomAccessContainer> void test_spot

template <class RandomAccessContainer>
void test_spots()
{
  typedef RandomAccessContainer V;
  using RealType = typename V::value_type;
  RealType tolerance = 1e-8;

  // Error checks:
  // Necessary conditions for instantiation:
  // 1. alpha.size() > 0.
  // 2. alpha[i] > 0.

  V alpha; // alpha.size() == 0.
  V x;     // x.size() == 0.
  BOOST_MATH_CHECK_THROW(dirichlet_distribution<V>(std::move(alpha)), std::domain_error);
  alpha.resize(2);
  alpha[0] = static_cast<RealType>(0.35);
  alpha[1] = static_cast<RealType>(-1.72); // alpha[1] < 0.
  BOOST_MATH_CHECK_THROW(dirichlet_distribution<V>(std::move(alpha)), std::domain_error);

  // Domain test for pdf. Necessary conditions for pdf:
  // 1. alpha[i] > 0.
  // 2. x.size() > 0.
  // 3. 0 <= x[i] <=1.
  // 4. sum(x) <= 1.
  alpha[0] = static_cast<RealType>(0.2);
  alpha[1] = static_cast<RealType>(1.7);
  BOOST_MATH_CHECK_THROW(boost::math::pdf(dirichlet_distribution<V>(std::move(alpha)), x), std::domain_error);
  x[0] = static_cast<RealType>(0.5);
  x[1] = static_cast<RealType>(0.5);
  BOOST_MATH_CHECK_THROW(boost::math::pdf(dirichlet_distribution<V>(std::move(alpha)), x), std::domain_error);

  alpha[0] = static_cast<RealType>(1.36);
  alpha[1] = static_cast<RealType>(0.0); // alpha[1] = 0.
  x[0] = static_cast<RealType>(0.47);
  x[1] = static_cast<RealType>(0.53);
  BOOST_MATH_CHECK_THROW(boost::math::pdf(dirichlet_distribution<V>(std::move(alpha)), x), std::domain_error);

  alpha[0] = static_cast<RealType>(1.26);
  alpha[1] = static_cast<RealType>(2.99);
  x[0] = static_cast<RealType>(0.5);
  x[1] = static_cast<RealType>(0.75); // sum(x) > 1.0
  BOOST_MATH_CHECK_THROW(boost::math::pdf(dirichlet_distribution<V>(std::move(alpha)), x), std::domain_error);

  alpha[0] = static_cast<RealType>(1.56);
  alpha[1] = static_cast<RealType>(4.00);
  x[0] = static_cast<RealType>(0.31);
  x[1] = static_cast<RealType>(-0.03); // x[1] < 0.
  BOOST_MATH_CHECK_THROW(boost::math::pdf(dirichlet_distribution<V>(std::move(alpha)), x), std::domain_error);

  alpha[0] = static_cast<RealType>(1.56);
  alpha[1] = static_cast<RealType>(4.00);
  x[0] = static_cast<RealType>(0.31);
  x[1] = static_cast<RealType>(1.06); // x[1] > 1.
  BOOST_MATH_CHECK_THROW(boost::math::pdf(dirichlet_distribution<V>(std::move(alpha)), x), std::domain_error);

  // Domain test for cdf. Necessary conditions for cdf:
  // 1. alpha[i] > 0
  // 2. 0 <= x[i] <= 1
  // 3. sum(x) <= 1.
  alpha[0] = static_cast<RealType>(1.56);
  alpha[1] = static_cast<RealType>(4.00);
  x[0] = static_cast<RealType>(0.31);
  x[1] = static_cast<RealType>(1.06); // x[1] > 1.
  BOOST_MATH_CHECK_THROW(boost::math::cdf(dirichlet_distribution<V>(std::move(alpha)), x), std::domain_error);

  alpha[0] = static_cast<RealType>(3.756);
  alpha[1] = static_cast<RealType>(4.91);
  x[0] = static_cast<RealType>(0.31);
  x[1] = static_cast<RealType>(-1.06); // x[1] < 0.
  BOOST_MATH_CHECK_THROW(boost::math::cdf(dirichlet_distribution<V>(std::move(alpha)), x), std::domain_error);

  alpha[0] = static_cast<RealType>(1.56);
  alpha[1] = static_cast<RealType>(-4.00); // alpha[1] < 0
  x[0] = static_cast<RealType>(0.31);
  x[1] = static_cast<RealType>(0.69);
  BOOST_MATH_CHECK_THROW(boost::math::cdf(dirichlet_distribution<V>(std::move(alpha)), x), std::domain_error);

  alpha[0] = static_cast<RealType>(0.0);
  alpha[1] = static_cast<RealType>(4.00); // alpha[0] = 0.
  x[0] = static_cast<RealType>(0.25);
  x[1] = static_cast<RealType>(0.75);
  BOOST_MATH_CHECK_THROW(boost::math::cdf(dirichlet_distribution<V>(std::move(alpha)), x), std::domain_error);

  alpha[0] = static_cast<RealType>(1.56);
  alpha[1] = static_cast<RealType>(4.00);
  x[0] = static_cast<RealType>(0.31);
  x[1] = static_cast<RealType>(0.71); // sum(x) > 1.
  BOOST_MATH_CHECK_THROW(boost::math::cdf(dirichlet_distribution<V>(std::move(alpha)), x), std::domain_error);

  // Domain test for mode. Necessary conditions for mode:
  // 1. alpha[i] > 1.
  alpha[0] = static_cast<RealType>(1.0);
  alpha[1] = static_cast<RealType>(1.4); // alpha[0] = 1.
  BOOST_MATH_CHECK_THROW(boost::math::mode(dirichlet_distribution<V>(std::move(alpha))), std::domain_error);

  alpha[0] = static_cast<RealType>(1.56);
  alpha[1] = static_cast<RealType>(0.92); // alpha[1] < 1.
  BOOST_MATH_CHECK_THROW(boost::math::mode(dirichlet_distribution<V>(std::move(alpha))), std::domain_error);

  // Some exact values of pdf.
  alpha[0] = static_cast<RealType>(1.0), alpha[1] = static_cast<RealType>(1.0);
  x[0] = static_cast<RealType>(0.5), x[1] = static_cast<RealType>(0.5);
  BOOST_CHECK_EQUAL(boost::math::pdf(dirichlet_distribution<V>(std::move(alpha)), x), static_cast<RealType>(1.0));

  alpha[0] = static_cast<RealType>(2.0), alpha[1] = static_cast<RealType>(2.0);
  x[0] = static_cast<RealType>(0.5), x[1] = static_cast<RealType>(0.5);
  BOOST_CHECK_EQUAL(boost::math::pdf(dirichlet_distribution<V>(std::move(alpha)), x), static_cast<RealType>(1.5));

  // Checking precalculated values on scipy.
  // https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.dirichlet.html
  alpha[0] = static_cast<RealType>(5.778238829);
  alpha[1] = static_cast<RealType>(2.55821892973);
  x[0] = static_cast<RealType>(0.23667289213);
  x[1] = static_cast<RealType>(0.76332710787);
  V mean = {static_cast<RealType>(0.693128783978901), static_cast<RealType>(0.3068712160210989)};
  V var = {static_cast<RealType>(0.022781795654775592), static_cast<RealType>(0.022781795654775592)};
  RealType entropy = static_cast<RealType>(-0.516646371355904);
  RealType pdf = static_cast<RealType>(0.05866153821852176);

  test_spot(std::move(alpha),
            std::move(x),
            std::move(mean),
            std::move(var),
            entropy, pdf, tolerance);

  alpha[0] = static_cast<RealType>(5.310948003052013);
  alpha[1] = static_cast<RealType>(8.003963132298916);
  x[0] = static_cast<RealType>(0.35042614416132284);
  x[1] = static_cast<RealType>(0.64957385583867716);
  mean[0] = static_cast<RealType>(0.398872207937724);
  mean[1] = static_cast<RealType>(0.601127792062276);
  var[0] = static_cast<RealType>(0.016749888798155716);
  var[1] = static_cast<RealType>(0.016749888798155716);
  pdf = static_cast<RealType>(2.870121181949622);
  entropy = static_cast<RealType>(-0.6347509574442718);

  test_spot(std::move(alpha),
            std::move(x),
            std::move(mean),
            std::move(var),
            entropy, pdf, tolerance);

  alpha[0] = static_cast<RealType>(8.764102220201394);
  alpha[1] = static_cast<RealType>(4.348446856921846);
  x[0] = static_cast<RealType>(0.6037585982123262);
  x[1] = static_cast<RealType>(0.39624140178767375);
  mean[0] = static_cast<RealType>(0.6683751701255137);
  mean[1] = static_cast<RealType>(0.33162482987448627);
  var[0] = static_cast<RealType>(0.015705865813037533);
  var[1] = static_cast<RealType>(0.015705865813037533);
  pdf = static_cast<RealType>(2.473329499915834);
  entropy = static_cast<RealType>(-0.6769547381491741);

  test_spot(std::move(alpha),
            std::move(x),
            std::move(mean),
            std::move(var),
            entropy, pdf, tolerance);

  alpha.resize(3);
  x.resize(3);
  mean.resize(3);
  var.resize(3);
  alpha[0] = static_cast<RealType>(5.622313698848736);
  alpha[1] = static_cast<RealType>(0.3516907178071482);
  alpha[2] = static_cast<RealType>(9.15629985496498);
  x[0] = static_cast<RealType>(0.6571425803855344);
  x[1] = static_cast<RealType>(0.2972004956337586);
  x[2] = static_cast<RealType>(0.04565692398070697);
  mean[0] = static_cast<RealType>(0.37159290374577736);
  mean[1] = static_cast<RealType>(0.023244127249099442);
  mean[2] = static_cast<RealType>(0.6051629690051231);
  var[0] = static_cast<RealType>(0.014476578600094457);
  var[1] = static_cast<RealType>(0.0014075269390591417);
  var[2] = static_cast<RealType>(0.014813158259538361);
  pdf = static_cast<RealType>(4.97846312846897e-08);
  entropy = static_cast<RealType>(-4.047215462643532);

  test_spot(std::move(alpha),
            std::move(x),
            std::move(mean),
            std::move(var),
            entropy, pdf, tolerance);

  alpha.resize(4);
  x.resize(4);
  mean.resize(4);
  var.resize(4);
  alpha[0] = static_cast<RealType>(5.958168192947443);
  alpha[1] = static_cast<RealType>(6.823198187239482);
  alpha[2] = static_cast<RealType>(6.297779996686504);
  alpha[3] = static_cast<RealType>(4.396226676824867);
  x[0] = static_cast<RealType>(0.15589020332495018);
  x[1] = static_cast<RealType>(0.3893497609653562);
  x[2] = static_cast<RealType>(0.060839680922786556);
  x[3] = static_cast<RealType>(0.393920354786907);
  mean[0] = static_cast<RealType>(0.2538050483508204);
  mean[1] = static_cast<RealType>(0.2906534508155371);
  mean[2] = static_cast<RealType>(0.26827177494818794);
  mean[3] = static_cast<RealType>(0.1872697258854546);
  var[0] = static_cast<RealType>(0.007737902313764369);
  var[1] = static_cast<RealType>(0.00842373359916587);
  var[2] = static_cast<RealType>(0.008020389690635378);
  var[3] = static_cast<RealType>(0.006218486448329886);
  pdf = static_cast<RealType>(0.2649374226055107);
  entropy = static_cast<RealType>(-3.4416182654031537);

  test_spot(std::move(alpha),
            std::move(x),
            std::move(mean),
            std::move(var),
            entropy, pdf, tolerance);

  alpha[0] = static_cast<RealType>(3.1779256968768976);
  alpha[1] = static_cast<RealType>(1.355989101047721);
  alpha[2] = static_cast<RealType>(5.594207813755373);
  alpha[3] = static_cast<RealType>(5.9897453525066355);
  x[0] = static_cast<RealType>(0.3388848203529338);
  x[1] = static_cast<RealType>(0.36731530174264704);
  x[2] = static_cast<RealType>(0.11166014002460622);
  x[3] = static_cast<RealType>(0.1821397378798129);
  mean[0] = static_cast<RealType>(0.19716787008915473);
  mean[1] = static_cast<RealType>(0.08412955758545015);
  mean[2] = static_cast<RealType>(0.34708112922785567);
  mean[3] = static_cast<RealType>(0.37162144309753953);
  var[0] = static_cast<RealType>(0.009247220589902615);
  var[1] = static_cast<RealType>(0.004501248361485874);
  var[2] = static_cast<RealType>(0.013238553973887957);
  var[3] = static_cast<RealType>(0.013641824239806118);
  pdf = static_cast<RealType>(0.06803159432725718);
  entropy = static_cast<RealType>(-3.398201562087422);

  test_spot(std::move(alpha),
            std::move(x),
            std::move(mean),
            std::move(var),
            entropy, pdf, tolerance);

  // No longer allow any parameter to be NaN or inf.
  if (std::numeric_limits<RealType>::has_quiet_NaN)
  {
    RealType not_a_num = std::numeric_limits<RealType>::quiet_NaN();
    alpha[0] = not_a_num;
    alpha[1] = static_cast<RealType>(0.37);
#ifndef BOOST_NO_EXCEPTIONS
    BOOST_MATH_CHECK_THROW(dirichlet_distribution<V>(std::move(alpha)), std::domain_error);
#else
    BOOST_MATH_CHECK_THROW(dirichlet_distribution<V>(std::move(alpha)), std::domain_error);
#endif

    // Non-finite parameters should throw.
    alpha[0] = static_cast<RealType>(1.67);
    alpha[1] = static_cast<RealType>(3.8);
    x[0] = not_a_num;
    x[1] = static_cast<RealType>(0.5);
    dirichlet_distribution<V> w(std::move(alpha));
    BOOST_MATH_CHECK_THROW(boost::math::pdf(w, x), std::domain_error); // x = NaN
    BOOST_MATH_CHECK_THROW(boost::math::cdf(w, x), std::domain_error); // x = NaN
  }                                                                    // has_quiet_NaN

  if (std::numeric_limits<RealType>::has_infinity)
  {
    // Attempt to construct from non-finite should throw.
    RealType infinite = std::numeric_limits<RealType>::infinity();
    alpha[0] = infinite;
    alpha[1] = static_cast<RealType>(7.2);
#ifndef BOOST_NO_EXCEPTIONS
    BOOST_MATH_CHECK_THROW(dirichlet_distribution<V>(std::move(alpha)), std::domain_error);
#else
    BOOST_MATH_CHECK_THROW(dirichlet_distribution<V>(std::move(alpha)), std::domain_error);
#endif
    alpha[0] = static_cast<RealType>(1.42);
    alpha[1] = static_cast<RealType>(7.91);
    x[0] = static_cast<RealType>(0.25);
    x[1] = infinite;
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
  test_spots<std::vector<long double>>();

  test_spots<std::vector<double>>();

  test_spots<std::vector<double>>();

  // #ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
  //   test_spots(); // Test long double.
  // #if !BOOST_WORKAROUND(__BORLANDC__, BOOST_TESTED_AT(0x582))
  //   test_spots(boost::math::concepts::real_concept(0.)); // Test real concept.
  // #endif
} // BOOST_AUTO_TEST_CASE( test_main )