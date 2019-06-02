//           Copyright Matthew Pulver 2018 - 2019.
// Distributed under the Boost Software License, Version 1.0.
//      (See accompanying file LICENSE_1_0.txt or copy at
//           https://www.boost.org/LICENSE_1_0.txt)

#include "test_autodiff.hpp"
#include <boost/utility/identity_type.hpp>

BOOST_AUTO_TEST_SUITE(test_autodiff_3)

BOOST_AUTO_TEST_CASE_TEMPLATE(atanh_test, T, all_float_types) {
  const T eps = 3000 * test_constants_t<T>::pct_epsilon(); // percent
  constexpr unsigned m = 5;
  const T cx = 0.5;
  auto x = make_fvar<T, m>(cx);
  auto y = atanh(x);
  // BOOST_CHECK_EQUAL(y.derivative(0) , atanh(cx)); // fails due to overload
  BOOST_CHECK_CLOSE(y.derivative(0u), atanh(static_cast<T>(x)), eps);
  BOOST_CHECK_CLOSE(y.derivative(1u), static_cast<T>(4) / 3, eps);
  BOOST_CHECK_CLOSE(y.derivative(2u), static_cast<T>(16) / 9, eps);
  BOOST_CHECK_CLOSE(y.derivative(3u), static_cast<T>(224) / 27, eps);
  BOOST_CHECK_CLOSE(y.derivative(4u), static_cast<T>(1280) / 27, eps);
  BOOST_CHECK_CLOSE(y.derivative(5u), static_cast<T>(31232) / 81, eps);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(atan_test, T, all_float_types) {
  BOOST_MATH_STD_USING
  using namespace boost;

  const T cx = 1.0;
  constexpr unsigned m = 5;
  const auto x = make_fvar<T, m>(cx);
  auto y = atan(x);
  const auto eps = boost::math::tools::epsilon<T>();
  BOOST_CHECK_CLOSE(y.derivative(0u), boost::math::constants::pi<T>() / 4, eps);
  BOOST_CHECK_CLOSE(y.derivative(1u), T(0.5), eps);
  BOOST_CHECK_CLOSE(y.derivative(2u), T(-0.5), eps);
  BOOST_CHECK_CLOSE(y.derivative(3u), T(0.5), eps);
  BOOST_CHECK_CLOSE(y.derivative(4u), T(0), eps);
  BOOST_CHECK_CLOSE(y.derivative(5u), T(-3), eps);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(erf_test, T, all_float_types) {
  BOOST_MATH_STD_USING
  using namespace boost;

  const T eps = 300 * 100 * boost::math::tools::epsilon<T>(); // percent
  const T cx = 1.0;
  constexpr unsigned m = 5;
  const auto x = make_fvar<T, m>(cx);
  auto y = erf(x);
  BOOST_CHECK_CLOSE(y.derivative(0u), erf(static_cast<T>(x)), eps);
  BOOST_CHECK_CLOSE(
      y.derivative(1u),
      T(2) / (math::constants::e<T>() * math::constants::root_pi<T>()), eps);
  BOOST_CHECK_CLOSE(
      y.derivative(2u),
      T(-4) / (math::constants::e<T>() * math::constants::root_pi<T>()), eps);
  BOOST_CHECK_CLOSE(
      y.derivative(3u),
      T(4) / (math::constants::e<T>() * math::constants::root_pi<T>()), eps);
  BOOST_CHECK_CLOSE(
      y.derivative(4u),
      T(8) / (math::constants::e<T>() * math::constants::root_pi<T>()), eps);
  BOOST_CHECK_CLOSE(
      y.derivative(5u),
      T(-40) / (math::constants::e<T>() * math::constants::root_pi<T>()), eps);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(sinc_test, T, bin_float_types) {
  BOOST_MATH_STD_USING
  const T eps = 20000 * boost::math::tools::epsilon<T>(); // percent
  const T cx = 1;
  constexpr unsigned m = 5;
  auto x = make_fvar<T, m>(cx);
  auto y = sinc(x);
  BOOST_CHECK_CLOSE(y.derivative(0u), sin(cx), eps);
  BOOST_CHECK_CLOSE(y.derivative(1u), cos(cx) - sin(cx), eps);
  BOOST_CHECK_CLOSE(y.derivative(2u), sin(cx) - 2 * cos(cx), eps);
  BOOST_CHECK_CLOSE(y.derivative(3u), T(5) * cos(cx) - T(3) * sin(cx), eps);
  BOOST_CHECK_CLOSE(y.derivative(4u), T(13) * sin(cx) - T(20) * cos(cx), eps);
  BOOST_CHECK_CLOSE(y.derivative(5u), T(101) * cos(cx) - T(65) * sin(cx), eps);
  // Test at x = 0
  auto y2 = sinc(make_fvar<T, 10>(0));
  BOOST_CHECK_CLOSE(y2.derivative(0u), T(1), eps);
  BOOST_CHECK_CLOSE(y2.derivative(1u), T(0), eps);
  BOOST_CHECK_CLOSE(y2.derivative(2u), -cx / T(3), eps);
  BOOST_CHECK_CLOSE(y2.derivative(3u), T(0), eps);
  BOOST_CHECK_CLOSE(y2.derivative(4u), cx / T(5), eps);
  BOOST_CHECK_CLOSE(y2.derivative(5u), T(0), eps);
  BOOST_CHECK_CLOSE(y2.derivative(6u), -cx / T(7), eps);
  BOOST_CHECK_CLOSE(y2.derivative(7u), T(0), eps);
  BOOST_CHECK_CLOSE(y2.derivative(8u), cx / T(9), eps);
  BOOST_CHECK_CLOSE(y2.derivative(9u), T(0), eps);
  BOOST_CHECK_CLOSE(y2.derivative(10u), -cx / T(11), eps);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(sinh_and_cosh, T, bin_float_types) {
  BOOST_MATH_STD_USING
  const T eps = 300 * boost::math::tools::epsilon<T>(); // percent
  const T cx = 1;
  constexpr unsigned m = 5;
  auto x = make_fvar<T, m>(cx);
  auto s = sinh(x);
  auto c = cosh(x);
  BOOST_CHECK_CLOSE(s.derivative(0u), sinh(static_cast<T>(x)), eps);
  BOOST_CHECK_CLOSE(c.derivative(0u), cosh(static_cast<T>(x)), eps);
  for (auto i : boost::irange(m + 1)) {
    BOOST_CHECK_CLOSE(s.derivative(i), static_cast<T>(i % 2 == 1 ? c : s), eps);
    BOOST_CHECK_CLOSE(c.derivative(i), static_cast<T>(i % 2 == 1 ? s : c), eps);
  }
}

#ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
BOOST_AUTO_TEST_CASE_TEMPLATE(tanh_test, T, all_float_types) {
  using bmp::fabs;
  using bmp::tanh;
  using detail::fabs;
  using detail::tanh;
  using std::fabs;
  using std::tanh;
  constexpr std::array<const char *, 6> tanh_derivatives{
      {"0."
       "76159415595576488811945828260479359041276859725793655159681050012195324"
       "457663848345894752167367671442190275970155",
       "0."
       "41997434161402606939449673904170144491718672823077095471331144024458989"
       "95240483056156940088623187260",
       "-0."
       "63970000844922450018849176930384395321921136306079914494299856318702069"
       "34885434644440069533372017992",
       "0."
       "62162668077129626310653042872222339967572411755445418563968706335816206"
       "22188951465548376863495698762",
       "0."
       "66509104475050167773507148092106234992757132833203125448814929383096463"
       "47626843278089998045994094537",
       "-5."
       "55689355847371979760458290231697200987383372116293456019531342394708989"
       "7942786231796317250984197038"}};
  const T cx = 1;
  constexpr std::size_t m = 5;
  auto x = make_fvar<T, m>(cx);
  auto t = tanh(x);
  for (auto i : boost::irange(tanh_derivatives.size())) {
    BOOST_TEST_WARN(isNearZero(t.derivative(i) -
                               boost::lexical_cast<T>(tanh_derivatives[i])));
  }
}
#endif

BOOST_AUTO_TEST_CASE_TEMPLATE(tan_test, T, bin_float_types) {
  BOOST_MATH_STD_USING
  const T eps = 800 * boost::math::tools::epsilon<T>(); // percent
  const T cx = boost::math::constants::third_pi<T>();
  const T root_three = boost::math::constants::root_three<T>();
  constexpr unsigned m = 5;
  const auto x = make_fvar<T, m>(cx);
  auto y = tan(x);
  BOOST_CHECK_CLOSE(y.derivative(0u), root_three, eps);
  BOOST_CHECK_CLOSE(y.derivative(1u), T(4), eps);
  BOOST_CHECK_CLOSE(y.derivative(2u), T(8) * root_three, eps);
  BOOST_CHECK_CLOSE(y.derivative(3u), T(80), eps);
  BOOST_CHECK_CLOSE(y.derivative(4u), T(352) * root_three, eps);
  BOOST_CHECK_CLOSE(y.derivative(5u), T(5824), eps);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(fmod_test, T, bin_float_types) {
  BOOST_MATH_STD_USING
  constexpr unsigned m = 3;
  const T cx = 3.25;
  const T cy = 0.5;
  auto x = make_fvar<T, m>(cx);
  auto y = fmod(x, autodiff_fvar<T, m>(cy));
  BOOST_CHECK_EQUAL(y.derivative(0u), T(0.25));
  BOOST_CHECK_EQUAL(y.derivative(1u), T(1));
  BOOST_CHECK_EQUAL(y.derivative(2u), T(0));
  BOOST_CHECK_EQUAL(y.derivative(3u), T(0));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(round_and_trunc, T, all_float_types) {
  BOOST_MATH_STD_USING
  constexpr unsigned m = 3;
  const T cx = 3.25;
  auto x = make_fvar<T, m>(cx);
  auto y = round(x);
  BOOST_CHECK_EQUAL(y.derivative(0u), round(cx));
  BOOST_CHECK_EQUAL(y.derivative(1u), T(0));
  BOOST_CHECK_EQUAL(y.derivative(2u), T(0));
  BOOST_CHECK_EQUAL(y.derivative(3u), T(0));
  y = trunc(x);
  BOOST_CHECK_EQUAL(y.derivative(0u), trunc(cx));
  BOOST_CHECK_EQUAL(y.derivative(1u), T(0));
  BOOST_CHECK_EQUAL(y.derivative(2u), T(0));
  BOOST_CHECK_EQUAL(y.derivative(3u), T(0));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(iround_and_itrunc, T, all_float_types) {
  BOOST_MATH_STD_USING
  using namespace boost::math;
  constexpr unsigned m = 3;
  const T cx = 3.25;
  auto x = make_fvar<T, m>(cx);
  int y = iround(x);
  BOOST_CHECK_EQUAL(y, iround(cx));
  y = itrunc(x);
  BOOST_CHECK_EQUAL(y, itrunc(cx));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(lambert_w0_test, T, all_float_types) {
  const T eps = 1000 * boost::math::tools::epsilon<T>(); // percent
  constexpr unsigned m = 10;
  const T cx = 3;
  // Mathematica: N[Table[D[ProductLog[x], {x, n}], {n, 0, 10}] /. x -> 3, 52]
  constexpr std::array<const char *, m + 1> answers{
      {"1.049908894964039959988697070552897904589466943706341",
       "0.1707244807388472968312949774415522047470762509741737",
       "-0.04336545501146252734105411312976167858858970875797718",
       "0.02321456264324789334313200360870492961288748451791104",
       "-0.01909049778427783072663170526188353869136655225133878",
       "0.02122935002563637629500975949987796094687564718834156",
       "-0.02979093848448877259041971538394953658978044986784643",
       "0.05051290266216717699803334605370337985567016837482099",
       "-0.1004503154972645060971099914384090562800544486549660",
       "0.2292464437392250211967939182075930820454464472006425",
       "-0.5905839053125614593682763387470654123192290838719517"}};
  auto x = make_fvar<T, m>(cx);
  auto y = lambert_w0(x);
  for (auto i : boost::irange(m + 1)) {
    const T answer = boost::lexical_cast<T>(answers[i]);
    BOOST_CHECK_CLOSE(y.derivative(i), answer, eps);
  }
  // const T cx0 = -1 / boost::math::constants::e<T>();
  // auto edge = lambert_w0(make_fvar<T,m>(cx0));
  // std::cout << "edge = " << edge << std::endl;
  // edge = depth(1)(-1,inf,-inf,inf,-inf,inf,-inf,inf,-inf,inf,-inf)
  // edge = depth(1)(-1,inf,-inf,inf,-inf,inf,-inf,inf,-inf,inf,-inf)
  // edge =
  // depth(1)(-1,3.68935e+19,-9.23687e+57,4.62519e+96,-2.89497e+135,2.02945e+174,-1.52431e+213,1.19943e+252,-9.75959e+290,8.14489e+329,-6.93329e+368)
}

BOOST_AUTO_TEST_SUITE_END()
