///////////////////////////////////////////////////////////////////
//  Copyright Christopher Kormanyos 2020.
//  Copyright John Maddock 2020.
//  Distributed under the Boost Software License,
//  Version 1.0. (See accompanying file LICENSE_1_0.txt
//  or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#define BOOST_MATH_DISABLE_DEPRECATED_03_WARNING

#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>

#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp> // Boost.Test

// In issue 396, a bug regarding overflow was traced back to the tgamma function.
// This test file tests the fix in 433 and its corresponding original bug report code.

namespace local {

template<class BigFloatType>
bool test_issue396_value_checker()
{
  typedef BigFloatType floating_point_type;

  // Table[N[Gamma[(1/2) + (10^n)], 103], {n, 0, 3, 1}]

  const boost::array<floating_point_type, 4U> control_values =
  {{
    floating_point_type("0.8862269254527580136490837416705725913987747280611935641069038949264556422955160906874753283692723327081"),
    floating_point_type("1.13327838894878556733457416558889247556029830827515977660872341452948339005600415371763053872760729065835E6"),
    floating_point_type("9.320963104082716608349109809141910437906497038162361154016117519412076597761162355221807605383606022360999E156"),
    floating_point_type("1.27230119569505546418224418037744456950663470986552782839399298388048086183891436363933143173336221543437E2566")
  }};

  boost::uint32_t ten_pow_n = UINT32_C(1);

  const floating_point_type tol = std::numeric_limits<floating_point_type>::epsilon() * UINT32_C(5000);

  bool result_is_ok = true;

  for(std::size_t i = 0U; i < control_values.size(); ++i)
  {
    const floating_point_type g = boost::math::tgamma(boost::math::constants::half<floating_point_type>() + ten_pow_n);

    ten_pow_n *= UINT32_C(10);

    using std::fabs;

    const floating_point_type closeness = fabs(1 - (g / control_values[i]));

    result_is_ok &= (closeness < tol);
  }

  return result_is_ok;
}

}

bool test_tgamma_for_issue396_cpp_dec_float020()
{
  typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float<20U>, boost::multiprecision::et_off> cpp_dec_float_type_020;

  const bool cpp_dec_float_020_is_ok = local::test_issue396_value_checker<cpp_dec_float_type_020>();

  return cpp_dec_float_020_is_ok;
}

bool test_tgamma_for_issue396_cpp_dec_float030()
{
  typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float<30U>, boost::multiprecision::et_off> cpp_dec_float_type_030;

  const bool cpp_dec_float_030_is_ok = local::test_issue396_value_checker<cpp_dec_float_type_030>();

  return cpp_dec_float_030_is_ok;
}

bool test_tgamma_for_issue396_cpp_dec_float040()
{
  typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float<40U>, boost::multiprecision::et_off> cpp_dec_float_type_040;

  const bool cpp_dec_float_040_is_ok = local::test_issue396_value_checker<cpp_dec_float_type_040>();

  return cpp_dec_float_040_is_ok;
}

bool test_tgamma_for_issue396_cpp_dec_float049()
{
  typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float<49U>, boost::multiprecision::et_off> cpp_dec_float_type_049;

  const bool cpp_dec_float_049_is_ok = local::test_issue396_value_checker<cpp_dec_float_type_049>();

  return cpp_dec_float_049_is_ok;
}

bool test_tgamma_for_issue396_cpp_dec_float050()
{
  typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float<50U>, boost::multiprecision::et_off> cpp_dec_float_type_050;

  const bool cpp_dec_float_050_is_ok = local::test_issue396_value_checker<cpp_dec_float_type_050>();

  return cpp_dec_float_050_is_ok;
}

bool test_tgamma_for_issue396_cpp_dec_float101()
{
  typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float<101U>, boost::multiprecision::et_off> cpp_dec_float_type_101;

  const bool cpp_dec_float_101_is_ok = local::test_issue396_value_checker<cpp_dec_float_type_101>();

  return cpp_dec_float_101_is_ok;
}

BOOST_AUTO_TEST_CASE(test_tgamma_for_issue396_cpp_dec_float020_tag)
{
  const bool cpp_dec_float_020_is_ok = test_tgamma_for_issue396_cpp_dec_float020();

  BOOST_CHECK(cpp_dec_float_020_is_ok);
}

BOOST_AUTO_TEST_CASE(test_tgamma_for_issue396_cpp_dec_float030_tag)
{
  const bool cpp_dec_float_030_is_ok = test_tgamma_for_issue396_cpp_dec_float030();

  BOOST_CHECK(cpp_dec_float_030_is_ok);
}

BOOST_AUTO_TEST_CASE(test_tgamma_for_issue396_cpp_dec_float040_tag)
{
  const bool cpp_dec_float_040_is_ok = test_tgamma_for_issue396_cpp_dec_float040();

  BOOST_CHECK(cpp_dec_float_040_is_ok);
}

BOOST_AUTO_TEST_CASE(test_tgamma_for_issue396_cpp_dec_float049_tag)
{
  const bool cpp_dec_float_049_is_ok = test_tgamma_for_issue396_cpp_dec_float049();

  BOOST_CHECK(cpp_dec_float_049_is_ok);
}

BOOST_AUTO_TEST_CASE(test_tgamma_for_issue396_cpp_dec_float050_tag)
{
  const bool cpp_dec_float_050_is_ok = test_tgamma_for_issue396_cpp_dec_float050();

  BOOST_CHECK(cpp_dec_float_050_is_ok);
}

BOOST_AUTO_TEST_CASE(test_tgamma_for_issue396_cpp_dec_float101_tag)
{
  const bool cpp_dec_float_101_is_ok = test_tgamma_for_issue396_cpp_dec_float101();

  BOOST_CHECK(cpp_dec_float_101_is_ok);
}
