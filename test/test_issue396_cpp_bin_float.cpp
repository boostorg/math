///////////////////////////////////////////////////////////////////
//  Copyright Christopher Kormanyos 2020.
//  Copyright John Maddock 2020.
//  Distributed under the Boost Software License,
//  Version 1.0. (See accompanying file LICENSE_1_0.txt
//  or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include "test_issue396_data_and_checker.ipp"
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>

#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp> // Boost.Test

// In issue 396, a bug regarding overflow was traced back to the tgamma function.
// This test file tests the fix in 433 and its corresponding original bug report code.

namespace local
{

bool test_tgamma_for_issue396_cpp_bin_float()
{
  typedef boost::multiprecision::number<boost::multiprecision::cpp_bin_float< 20U>, boost::multiprecision::et_off> cpp_bin_float_type_020;
  typedef boost::multiprecision::number<boost::multiprecision::cpp_bin_float< 30U>, boost::multiprecision::et_off> cpp_bin_float_type_030;
  typedef boost::multiprecision::number<boost::multiprecision::cpp_bin_float< 40U>, boost::multiprecision::et_off> cpp_bin_float_type_040;
  typedef boost::multiprecision::number<boost::multiprecision::cpp_bin_float< 49U>, boost::multiprecision::et_off> cpp_bin_float_type_049;
  typedef boost::multiprecision::number<boost::multiprecision::cpp_bin_float< 50U>, boost::multiprecision::et_off> cpp_bin_float_type_050;
  typedef boost::multiprecision::number<boost::multiprecision::cpp_bin_float< 51U>, boost::multiprecision::et_off> cpp_bin_float_type_051;
  typedef boost::multiprecision::number<boost::multiprecision::cpp_bin_float<101U>, boost::multiprecision::et_off> cpp_bin_float_type_101;
  typedef boost::multiprecision::number<boost::multiprecision::cpp_bin_float<501U>, boost::multiprecision::et_off> cpp_bin_float_type_501;

  const bool cpp_bin_float_020_is_ok = boost::math::detail::test_issue396_value_checker<cpp_bin_float_type_020>();
  const bool cpp_bin_float_030_is_ok = boost::math::detail::test_issue396_value_checker<cpp_bin_float_type_030>();
  const bool cpp_bin_float_040_is_ok = boost::math::detail::test_issue396_value_checker<cpp_bin_float_type_040>();
  const bool cpp_bin_float_049_is_ok = boost::math::detail::test_issue396_value_checker<cpp_bin_float_type_049>();
  const bool cpp_bin_float_050_is_ok = boost::math::detail::test_issue396_value_checker<cpp_bin_float_type_050>();
  const bool cpp_bin_float_051_is_ok = boost::math::detail::test_issue396_value_checker<cpp_bin_float_type_051>();
  const bool cpp_bin_float_101_is_ok = boost::math::detail::test_issue396_value_checker<cpp_bin_float_type_101>();
  const bool cpp_bin_float_501_is_ok = boost::math::detail::test_issue396_value_checker<cpp_bin_float_type_501>();

  const bool cpp_bin_float_is_ok = (   cpp_bin_float_020_is_ok
                                    && cpp_bin_float_030_is_ok
                                    && cpp_bin_float_040_is_ok
                                    && cpp_bin_float_049_is_ok
                                    && cpp_bin_float_050_is_ok
                                    && cpp_bin_float_051_is_ok
                                    && cpp_bin_float_101_is_ok
                                    && cpp_bin_float_501_is_ok);

  return cpp_bin_float_is_ok;
}

}

BOOST_AUTO_TEST_CASE(test_tgamma_for_issue396_cpp_bin_float_tag)
{
  const bool cpp_bin_float_is_ok = local::test_tgamma_for_issue396_cpp_bin_float();

  BOOST_CHECK(cpp_bin_float_is_ok);
}
