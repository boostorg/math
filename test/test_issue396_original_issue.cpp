///////////////////////////////////////////////////////////////////
//  Copyright Christopher Kormanyos 2020.
//  Copyright John Maddock 2020.
//  Distributed under the Boost Software License,
//  Version 1.0. (See accompanying file LICENSE_1_0.txt
//  or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/math/distributions/beta.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>

#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp> // Boost.Test

// In issue 396, a bug regarding overflow was traced back to the tgamma function.
// This test file tests the fix in 433 and its corresponding original bug report code.

namespace local {

bool test_cdf____for_issue396_bug_report___()
{
  typedef boost::math::beta_distribution<boost::multiprecision::cpp_dec_float_100> b_dist_type;

  b_dist_type b_dist(1.0, 2.0);

  boost::multiprecision::cpp_dec_float_100 test = -boost::math::cdf(boost::math::complement(b_dist, 0.5));

  using std::fabs;

  const boost::multiprecision::cpp_dec_float_100 control = boost::multiprecision::cpp_dec_float_100(-0.25F);

  const boost::multiprecision::cpp_dec_float_100 closeness = fabs(1 - fabs(test / control));

  const bool result_is_ok = (closeness < std::numeric_limits<boost::multiprecision::cpp_dec_float_100>::epsilon() * 10U);

  return result_is_ok;
}

}

BOOST_AUTO_TEST_CASE(test_cdf____for_issue396_bug_report____tag)
{
  const bool cpp_bug_check_is_ok = local::test_cdf____for_issue396_bug_report___();

  BOOST_CHECK(cpp_bug_check_is_ok);
}
