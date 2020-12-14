///////////////////////////////////////////////////////////////////
//  Copyright Christopher Kormanyos 2020.
//  Copyright John Maddock 2020.
//  Distributed under the Boost Software License,
//  Version 1.0. (See accompanying file LICENSE_1_0.txt
//  or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp> // Boost.Test

#include <boost/array.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/distributions/beta.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>

// In issue 396, a bug regarding overflow was traced back to the tgamma function.
// This test file tests the fix in 433 and its corresponding original bug report code.

namespace local {

namespace detail {

template<class BigFloatType>
bool test_tgamma_for_issue396_value_checker()
{
  using floating_point_type = BigFloatType;

  // Table[N[Gamma[(1/2) + (10^n)], 503], {n, 0, 3, 1}]

  const boost::array<floating_point_type, 4U> control =
  {{
    floating_point_type("0.88622692545275801364908374167057259139877472806119356410690389492645564229551609068747532836927233270811341181214128533311807643286221130126254685480139353423101884932655256142496258651447541311446604768963398140008731950767573986025835009509261700929272348724745632015696088776295310820270966625045319920380686673873757671683399489468292591820439772558258086938002953369671589566640492742312409245102732742609780662578082373375752136938052805399806355360503018602224183618264830685404716174941583421211"),
    floating_point_type("1.1332783889487855673345741655888924755602983082751597766087234145294833900560041537176305387276072906583502717008932373348895801731780765775979953796646009714415152490764416630481375706606053932396039541459764525989187023837695167161085523804417015113740063535865261183579508922972990386756543208549178543857406373798865630303794109491220205170302558277398183764099268751365861892723863412249690833216320407918186480305202146014474770321625907339955121137559264239090240758401696425720048012081453338360E6"),
    floating_point_type("9.3209631040827166083491098091419104379064970381623611540161175194120765977611623552218076053836060223609993676387199220631835256331102029826429784793420637988460945604451237342972023988743201341318701614328454618664952897316247603329530308777063116667275003586843755354841307657702809317290363831151480295446074722690100652644579131609996151999119113967501099655433566352849645431012667388627160383486515144610582794470005796689975604764040892168183647321540427819244511610500074895473959438490375652158E156"),
    floating_point_type("1.2723011956950554641822441803774445695066347098655278283939929838804808618389143636393314317333622154343715992535881414698586440455330620652019981627229614973177953241634213768203151670660953863412381880742653187501307209325406338924004280546485392703623101051957976224599412003938216329590158926122017907280168159527761842471509358725974702333390709735919152262756462872191402491961250987725812831155116532550035967994387094267848607390288008530653715254376729558412833771092612838971719786622446726968E2566")
  }};

  boost::uint32_t ten_pow_n = UINT32_C(1);

  const floating_point_type tol = std::numeric_limits<floating_point_type>::epsilon() * UINT32_C(5000);

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

} // namespace detail

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

bool test_tgamma_for_issue396_cpp_dec_float()
{
  using cpp_dec_float_type_020 = boost::multiprecision::number<boost::multiprecision::cpp_dec_float< 20U>, boost::multiprecision::et_off>;
  using cpp_dec_float_type_030 = boost::multiprecision::number<boost::multiprecision::cpp_dec_float< 30U>, boost::multiprecision::et_off>;
  using cpp_dec_float_type_040 = boost::multiprecision::number<boost::multiprecision::cpp_dec_float< 40U>, boost::multiprecision::et_off>;
  using cpp_dec_float_type_049 = boost::multiprecision::number<boost::multiprecision::cpp_dec_float< 49U>, boost::multiprecision::et_off>;
  using cpp_dec_float_type_050 = boost::multiprecision::number<boost::multiprecision::cpp_dec_float< 50U>, boost::multiprecision::et_off>;
  using cpp_dec_float_type_051 = boost::multiprecision::number<boost::multiprecision::cpp_dec_float< 51U>, boost::multiprecision::et_off>;
  using cpp_dec_float_type_101 = boost::multiprecision::number<boost::multiprecision::cpp_dec_float<101U>, boost::multiprecision::et_off>;
  using cpp_dec_float_type_501 = boost::multiprecision::number<boost::multiprecision::cpp_dec_float<501U>, boost::multiprecision::et_off>;

  const bool cpp_dec_float_020_is_ok = detail::test_tgamma_for_issue396_value_checker<cpp_dec_float_type_020>();
  const bool cpp_dec_float_030_is_ok = detail::test_tgamma_for_issue396_value_checker<cpp_dec_float_type_030>();
  const bool cpp_dec_float_040_is_ok = detail::test_tgamma_for_issue396_value_checker<cpp_dec_float_type_040>();
  const bool cpp_dec_float_049_is_ok = detail::test_tgamma_for_issue396_value_checker<cpp_dec_float_type_049>();
  const bool cpp_dec_float_050_is_ok = detail::test_tgamma_for_issue396_value_checker<cpp_dec_float_type_050>();
  const bool cpp_dec_float_051_is_ok = detail::test_tgamma_for_issue396_value_checker<cpp_dec_float_type_051>();
  const bool cpp_dec_float_101_is_ok = detail::test_tgamma_for_issue396_value_checker<cpp_dec_float_type_101>();
  const bool cpp_dec_float_501_is_ok = detail::test_tgamma_for_issue396_value_checker<cpp_dec_float_type_501>();

  const bool cpp_dec_float_is_ok = (   cpp_dec_float_020_is_ok
                                    && cpp_dec_float_030_is_ok
                                    && cpp_dec_float_040_is_ok
                                    && cpp_dec_float_049_is_ok
                                    && cpp_dec_float_050_is_ok
                                    && cpp_dec_float_051_is_ok
                                    && cpp_dec_float_101_is_ok
                                    && cpp_dec_float_501_is_ok);

  return cpp_dec_float_is_ok;
}

bool test_tgamma_for_issue396_cpp_bin_float()
{
  using cpp_bin_float_type_020 = boost::multiprecision::number<boost::multiprecision::cpp_bin_float< 20U>, boost::multiprecision::et_off>;
  using cpp_bin_float_type_030 = boost::multiprecision::number<boost::multiprecision::cpp_bin_float< 30U>, boost::multiprecision::et_off>;
  using cpp_bin_float_type_040 = boost::multiprecision::number<boost::multiprecision::cpp_bin_float< 40U>, boost::multiprecision::et_off>;
  using cpp_bin_float_type_049 = boost::multiprecision::number<boost::multiprecision::cpp_bin_float< 49U>, boost::multiprecision::et_off>;
  using cpp_bin_float_type_050 = boost::multiprecision::number<boost::multiprecision::cpp_bin_float< 50U>, boost::multiprecision::et_off>;
  using cpp_bin_float_type_051 = boost::multiprecision::number<boost::multiprecision::cpp_bin_float< 51U>, boost::multiprecision::et_off>;
  using cpp_bin_float_type_101 = boost::multiprecision::number<boost::multiprecision::cpp_bin_float<101U>, boost::multiprecision::et_off>;
  using cpp_bin_float_type_501 = boost::multiprecision::number<boost::multiprecision::cpp_bin_float<501U>, boost::multiprecision::et_off>;

  const bool cpp_bin_float_020_is_ok = detail::test_tgamma_for_issue396_value_checker<cpp_bin_float_type_020>();
  const bool cpp_bin_float_030_is_ok = detail::test_tgamma_for_issue396_value_checker<cpp_bin_float_type_030>();
  const bool cpp_bin_float_040_is_ok = detail::test_tgamma_for_issue396_value_checker<cpp_bin_float_type_040>();
  const bool cpp_bin_float_049_is_ok = detail::test_tgamma_for_issue396_value_checker<cpp_bin_float_type_049>();
  const bool cpp_bin_float_050_is_ok = detail::test_tgamma_for_issue396_value_checker<cpp_bin_float_type_050>();
  const bool cpp_bin_float_051_is_ok = detail::test_tgamma_for_issue396_value_checker<cpp_bin_float_type_051>();
  const bool cpp_bin_float_101_is_ok = detail::test_tgamma_for_issue396_value_checker<cpp_bin_float_type_101>();
  const bool cpp_bin_float_501_is_ok = detail::test_tgamma_for_issue396_value_checker<cpp_bin_float_type_501>();

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

} // namespace local

BOOST_AUTO_TEST_CASE(test_cdf____for_issue396_bug_report____tag)
{
  const bool cpp_bug_check_is_ok = local::test_cdf____for_issue396_bug_report___();

  BOOST_CHECK(cpp_bug_check_is_ok);
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
