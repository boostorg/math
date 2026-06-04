//  Copyright Anton Leontev 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0.
//  (See accompanying file LICENSE_1_0.txt or copy at
//   http://www.boost.org/LICENSE_1_0.txt)
//
//  Regression test for intermediate overflow of x*x in the Student's t
//  distribution pdf() and cdf().

#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/math/tools/test.hpp>

#include <boost/math/distributions/students_t.hpp>
#include <boost/math/policies/policy.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/math/tools/precision.hpp>

#ifndef BOOST_MATH_STANDALONE
#include <boost/math/concepts/real_concept.hpp>  // for real_concept
#include <boost/multiprecision/cpp_bin_float.hpp>  // for cpp_bin_float_50
#endif

#include <cerrno>
#include <cmath>
#include <limits>

template <class RealType>
RealType overflowing_deviate()
{
   BOOST_MATH_STD_USING  // pull in std::sqrt for built-in types via ADL
   using boost::math::tools::max_value;
   return 2 * sqrt(max_value<RealType>());
}

template <class RealType>
void test_throws(RealType)
{
   using namespace boost::math;
   using namespace boost::math::policies;

   typedef policy<overflow_error<throw_on_error> > throw_policy;
   typedef students_t_distribution<RealType, throw_policy> dist_t;

   dist_t dist(static_cast<RealType>(5));
   const RealType big = overflowing_deviate<RealType>();

   BOOST_CHECK((boost::math::isfinite)(big));

   BOOST_CHECK_THROW(pdf(dist,  big), std::overflow_error);
   BOOST_CHECK_THROW(pdf(dist, -big), std::overflow_error);
   BOOST_CHECK_THROW(cdf(dist,  big), std::overflow_error);
   BOOST_CHECK_THROW(cdf(dist, -big), std::overflow_error);
}

template <class RealType>
void test_errno(RealType)
{
   using namespace boost::math;
   using namespace boost::math::policies;

   typedef policy<overflow_error<errno_on_error> > errno_policy;
   typedef students_t_distribution<RealType, errno_policy> dist_t;

   dist_t dist(static_cast<RealType>(5));
   const RealType big = overflowing_deviate<RealType>();

   errno = 0;
   RealType r = pdf(dist, big);
   BOOST_CHECK_EQUAL(errno, ERANGE);
   BOOST_CHECK(!(boost::math::isfinite)(r));

   errno = 0;
   r = cdf(dist, big);
   BOOST_CHECK_EQUAL(errno, ERANGE);
   BOOST_CHECK(!(boost::math::isfinite)(r));
}

template <class RealType>
void test_no_regression(RealType)
{
   using namespace boost::math;
   using namespace boost::math::policies;

   typedef policy<overflow_error<errno_on_error> > errno_policy;
   typedef students_t_distribution<RealType, errno_policy> dist_t;

   dist_t dist(static_cast<RealType>(10));
   const RealType tol = boost::math::tools::epsilon<RealType>() * 100;

   // No errno assertions on these success paths: C11 7.5/3 allows any
   // library call to set errno even on success, and under UBSan the
   // diagnostic printer itself clobbers errno (isatty() on a redirected
   // stderr yields ENOTTY).  Policy-driven errno semantics are verified
   // separately in test_errno().
   BOOST_CHECK_CLOSE_FRACTION(cdf(dist, static_cast<RealType>(0)),
                              static_cast<RealType>(0.5), tol);

   RealType p = pdf(dist, static_cast<RealType>(1.5));
   BOOST_CHECK((boost::math::isfinite)(p));
   BOOST_CHECK(p > 0);

   RealType big_but_safe = static_cast<RealType>(1e6);
   RealType c = cdf(dist, big_but_safe);
   BOOST_CHECK_CLOSE_FRACTION(c, static_cast<RealType>(1), tol);
}

BOOST_AUTO_TEST_CASE(students_t_overflow_test)
{
   test_throws(0.0F);
   test_throws(0.0);
   test_throws(0.0L);

   test_errno(0.0F);
   test_errno(0.0);
   test_errno(0.0L);

   test_no_regression(0.0F);
   test_no_regression(0.0);
   test_no_regression(0.0L);

   // Non-standard types: these must compile (exercises ADL fabs) and behave.
   // errno semantics are only meaningful for built-in types, so only the
   // policy-driven throw and the no-regression cases are exercised here.
   // real_concept and multiprecision compute distribution constants via
   // lexical_cast, which standalone mode disables, so both are skipped there.
   #ifndef BOOST_MATH_STANDALONE
   test_throws(boost::math::concepts::real_concept(0));
   test_no_regression(boost::math::concepts::real_concept(0));

   test_throws(boost::multiprecision::cpp_bin_float_50(0));
   test_no_regression(boost::multiprecision::cpp_bin_float_50(0));
   #endif
}
