//  Copyright Anton Leontev 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0.
//  (See accompanying file LICENSE_1_0.txt or copy at
//   http://www.boost.org/LICENSE_1_0.txt)
//
//  Regression test for spurious intermediate overflow of x*x in the
//  Student's t pdf() and cdf(): for |x| > sqrt(max_value) the functions
//  must return the exact tail values (0 for the pdf, 0/1 for the cdf and
//  its complement) without raising FE_OVERFLOW and without invoking any
//  error handler.

#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/math/tools/test.hpp>

#include <boost/math/distributions/students_t.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/math/tools/precision.hpp>

#ifndef BOOST_MATH_STANDALONE
#include <boost/math/concepts/real_concept.hpp>  // for real_concept
#include <boost/multiprecision/cpp_bin_float.hpp>  // for cpp_bin_float_50
#endif

#include <cfenv>
#include <cmath>
#include <limits>

#pragma STDC FENV_ACCESS ON

template <class RealType>
RealType overflowing_deviate()
{
   BOOST_MATH_STD_USING  // ADL sqrt for built-in and UDT types
   using boost::math::tools::max_value;
   return 2 * sqrt(max_value<RealType>());
}

// For |x| > sqrt(max_value) the mathematically exact results have already
// saturated: pdf == 0, cdf == 0 or 1.  These must be returned directly.
// A throwing policy is used so that any error-handler invocation would
// fail the test with an unexpected exception.
template <class RealType>
void test_tail_values(RealType)
{
   using namespace boost::math;
   using namespace boost::math::policies;

   typedef policy<overflow_error<throw_on_error> > strict_policy;
   typedef students_t_distribution<RealType, strict_policy> dist_t;

   dist_t dist(static_cast<RealType>(5));
   const RealType big = overflowing_deviate<RealType>();

   BOOST_CHECK((boost::math::isfinite)(big));

   BOOST_CHECK_EQUAL(pdf(dist,  big), static_cast<RealType>(0));
   BOOST_CHECK_EQUAL(pdf(dist, -big), static_cast<RealType>(0));

   BOOST_CHECK_EQUAL(cdf(dist,  big), static_cast<RealType>(1));
   BOOST_CHECK_EQUAL(cdf(dist, -big), static_cast<RealType>(0));

   BOOST_CHECK_EQUAL(cdf(complement(dist,  big)), static_cast<RealType>(0));
   BOOST_CHECK_EQUAL(cdf(complement(dist, -big)), static_cast<RealType>(1));
}

// Built-in types: the calls must not leave FE_OVERFLOW set -- avoiding the
// spurious intermediate overflow is the whole point of the guard.
template <class RealType>
void test_no_fp_exceptions(RealType)
{
   using namespace boost::math;

   students_t_distribution<RealType> dist(static_cast<RealType>(5));
   const RealType big = overflowing_deviate<RealType>();

   std::feclearexcept(FE_ALL_EXCEPT);
   RealType p = pdf(dist, big);
   RealType c1 = cdf(dist, big);
   RealType c2 = cdf(dist, -big);
   BOOST_CHECK_EQUAL(std::fetestexcept(FE_OVERFLOW), 0);
   BOOST_CHECK_EQUAL(p, static_cast<RealType>(0));
   BOOST_CHECK_EQUAL(c1, static_cast<RealType>(1));
   BOOST_CHECK_EQUAL(c2, static_cast<RealType>(0));
}

// Ordinary arguments must be unaffected by the guard.
template <class RealType>
void test_no_regression(RealType)
{
   using namespace boost::math;

   students_t_distribution<RealType> dist(static_cast<RealType>(10));
   const RealType tol = boost::math::tools::epsilon<RealType>() * 100;

   BOOST_CHECK_CLOSE_FRACTION(cdf(dist, static_cast<RealType>(0)),
                              static_cast<RealType>(0.5), tol);

   RealType p = pdf(dist, static_cast<RealType>(1.5));
   BOOST_CHECK((boost::math::isfinite)(p));
   BOOST_CHECK(p > 0);

   RealType big_but_safe = static_cast<RealType>(1e6);
   BOOST_CHECK_CLOSE_FRACTION(cdf(dist, big_but_safe),
                              static_cast<RealType>(1), tol);
}

BOOST_AUTO_TEST_CASE(students_t_overflow_test)
{
   test_tail_values(0.0F);
   test_tail_values(0.0);
   test_tail_values(0.0L);

   test_no_fp_exceptions(0.0F);
   test_no_fp_exceptions(0.0);
   test_no_fp_exceptions(0.0L);

   test_no_regression(0.0F);
   test_no_regression(0.0);
   test_no_regression(0.0L);

   // real_concept and multiprecision compute distribution constants via
   // lexical_cast, which standalone mode disables, so both are skipped there.
#ifndef BOOST_MATH_STANDALONE
   test_tail_values(boost::math::concepts::real_concept(0));
   test_no_regression(boost::math::concepts::real_concept(0));

   test_tail_values(boost::multiprecision::cpp_bin_float_50(0));
   test_no_regression(boost::multiprecision::cpp_bin_float_50(0));
#endif
}
