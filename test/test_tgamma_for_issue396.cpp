///////////////////////////////////////////////////////////////////
//  Copyright John Maddock 2020.
//  Copyright Christopher Kormanyos 2020.
//  Distributed under the Boost Software License,
//  Version 1.0. (See accompanying file LICENSE_1_0.txt
//  or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/array.hpp>
#include <boost/math/constants/constants.hpp>
#include "boost/math/distributions/beta.hpp"
#include <boost/math/special_functions/gamma.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>

#define BOOST_TEST_MODULE test_tgamma_for_issue396
#include <boost/test/included/unit_test.hpp>

// In issue 396, a bug regarding overflow was traced back to the tgamma function.
// This test file tests the fix in 433 and its corresponding original bug report code.

namespace local {

bool test_cdf____for_issue396_bug_report___()
{
  using b_dist_type = boost::math::beta_distribution<boost::multiprecision::cpp_dec_float_100>;

  b_dist_type b_dist(1.0, 2.0);

  auto test = -boost::math::cdf(boost::math::complement(b_dist, 0.5));

  using std::fabs;

  const boost::multiprecision::cpp_dec_float_100 control = boost::multiprecision::cpp_dec_float_100(-0.25F);

  const boost::multiprecision::cpp_dec_float_100 closeness = fabs(1 - fabs(test / control));

  const bool result_is_ok = (closeness < std::numeric_limits<boost::multiprecision::cpp_dec_float_100>::epsilon() * 10U);

  return result_is_ok;
}

template<class BigFloatType>
bool test_tgamma_for_issue396()
{
  using floating_point_type = BigFloatType;

  // Table[N[Gamma[(1/2) + (10^n)], 103], {n, 0, 3, 1}]

  const boost::array<floating_point_type, 4U> control =
  {{
    floating_point_type("0.8862269254527580136490837416705725913987747280611935641069038949264556422955160906874753283692723327081"),
    floating_point_type("1.133278388948785567334574165588892475560298308275159776608723414529483390056004153717630538727607290658E6"),
    floating_point_type("9.320963104082716608349109809141910437906497038162361154016117519412076597761162355221807605383606022361E156"),
    floating_point_type("1.272301195695055464182244180377444569506634709865527828393992983880480861838914363639331431733362215434E2566")
  }};

  boost::uint32_t ten_pow_n = UINT32_C(1);

  const floating_point_type tol = std::numeric_limits<floating_point_type>::epsilon() * UINT32_C(100000);

  bool result_is_ok = true;

  for(typename boost::array<floating_point_type, 4U>::size_type i = 0U; i < control.size(); ++i)
  {
    const floating_point_type g = boost::math::tgamma(boost::math::constants::half<floating_point_type>() + ten_pow_n);

    ten_pow_n *= UINT32_C(10);

    using std::fabs;

    const floating_point_type closeness = fabs(1 - (g / control[i]));

    result_is_ok &= (closeness < tol);
  }

  return result_is_ok;
}

bool test_tgamma_for_issue396_cpp_dec_float()
{
  using cpp_dec_float_type_37 = boost::multiprecision::number<boost::multiprecision::cpp_dec_float<37U>, boost::multiprecision::et_off>;
  using cpp_dec_float_type_38 = boost::multiprecision::number<boost::multiprecision::cpp_dec_float<38U>, boost::multiprecision::et_off>;
  using cpp_dec_float_type_39 = boost::multiprecision::number<boost::multiprecision::cpp_dec_float<39U>, boost::multiprecision::et_off>;
  using cpp_dec_float_type_40 = boost::multiprecision::number<boost::multiprecision::cpp_dec_float<40U>, boost::multiprecision::et_off>;
  using cpp_dec_float_type_41 = boost::multiprecision::number<boost::multiprecision::cpp_dec_float<41U>, boost::multiprecision::et_off>;
  using cpp_dec_float_type_42 = boost::multiprecision::number<boost::multiprecision::cpp_dec_float<42U>, boost::multiprecision::et_off>;
  using cpp_dec_float_type_43 = boost::multiprecision::number<boost::multiprecision::cpp_dec_float<43U>, boost::multiprecision::et_off>;

  const bool cpp_dec_float_37_is_ok = test_tgamma_for_issue396<cpp_dec_float_type_37>();
  const bool cpp_dec_float_38_is_ok = test_tgamma_for_issue396<cpp_dec_float_type_38>();
  const bool cpp_dec_float_39_is_ok = test_tgamma_for_issue396<cpp_dec_float_type_39>();
  const bool cpp_dec_float_40_is_ok = test_tgamma_for_issue396<cpp_dec_float_type_40>();
  const bool cpp_dec_float_41_is_ok = test_tgamma_for_issue396<cpp_dec_float_type_41>();
  const bool cpp_dec_float_42_is_ok = test_tgamma_for_issue396<cpp_dec_float_type_42>();
  const bool cpp_dec_float_43_is_ok = test_tgamma_for_issue396<cpp_dec_float_type_43>();

  const bool cpp_dec_float_is_ok =    cpp_dec_float_37_is_ok
                                   && cpp_dec_float_38_is_ok
                                   && cpp_dec_float_39_is_ok
                                   && cpp_dec_float_40_is_ok
                                   && cpp_dec_float_41_is_ok
                                   && cpp_dec_float_42_is_ok
                                   && cpp_dec_float_43_is_ok;

  return cpp_dec_float_is_ok;
}

bool test_tgamma_for_issue396_cpp_bin_float()
{
  using cpp_bin_float_type_37 = boost::multiprecision::number<boost::multiprecision::cpp_bin_float<37U>, boost::multiprecision::et_off>;
  using cpp_bin_float_type_38 = boost::multiprecision::number<boost::multiprecision::cpp_bin_float<38U>, boost::multiprecision::et_off>;
  using cpp_bin_float_type_39 = boost::multiprecision::number<boost::multiprecision::cpp_bin_float<39U>, boost::multiprecision::et_off>;
  using cpp_bin_float_type_40 = boost::multiprecision::number<boost::multiprecision::cpp_bin_float<40U>, boost::multiprecision::et_off>;
  using cpp_bin_float_type_41 = boost::multiprecision::number<boost::multiprecision::cpp_bin_float<41U>, boost::multiprecision::et_off>;
  using cpp_bin_float_type_42 = boost::multiprecision::number<boost::multiprecision::cpp_bin_float<42U>, boost::multiprecision::et_off>;
  using cpp_bin_float_type_43 = boost::multiprecision::number<boost::multiprecision::cpp_bin_float<43U>, boost::multiprecision::et_off>;

  const bool cpp_bin_float_37_is_ok = test_tgamma_for_issue396<cpp_bin_float_type_37>();
  const bool cpp_bin_float_38_is_ok = test_tgamma_for_issue396<cpp_bin_float_type_38>();
  const bool cpp_bin_float_39_is_ok = test_tgamma_for_issue396<cpp_bin_float_type_39>();
  const bool cpp_bin_float_40_is_ok = test_tgamma_for_issue396<cpp_bin_float_type_40>();
  const bool cpp_bin_float_41_is_ok = test_tgamma_for_issue396<cpp_bin_float_type_41>();
  const bool cpp_bin_float_42_is_ok = test_tgamma_for_issue396<cpp_bin_float_type_42>();
  const bool cpp_bin_float_43_is_ok = test_tgamma_for_issue396<cpp_bin_float_type_43>();

  const bool cpp_bin_float_is_ok =    cpp_bin_float_37_is_ok
                                   && cpp_bin_float_38_is_ok
                                   && cpp_bin_float_39_is_ok
                                   && cpp_bin_float_40_is_ok
                                   && cpp_bin_float_41_is_ok
                                   && cpp_bin_float_42_is_ok
                                   && cpp_bin_float_43_is_ok;

  return cpp_bin_float_is_ok;
}

} // namespace local

BOOST_AUTO_TEST_CASE(test_cdf____for_issue396_bug_report____tag)
{
  const bool cpp_bug_report_is_ok = local::test_cdf____for_issue396_bug_report___();
  BOOST_CHECK(cpp_bug_report_is_ok);
}

BOOST_AUTO_TEST_CASE(test_tgamma_for_issue396_cpp_dec_float_tag)
{
  const bool cpp_dec_float_is_ok = local::test_tgamma_for_issue396_cpp_dec_float();
  BOOST_CHECK(cpp_dec_float_is_ok);
}

BOOST_AUTO_TEST_CASE(test_tgamma_for_issue396_cpp_bin_float_tag)
{
  const bool cpp_bin_float_is_ok = local::test_tgamma_for_issue396_cpp_bin_float();
  BOOST_CHECK(cpp_bin_float_is_ok);
}
