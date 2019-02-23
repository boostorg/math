//           Copyright Matthew Pulver 2018 - 2019.
// Distributed under the Boost Software License, Version 1.0.
//      (See accompanying file LICENSE_1_0.txt or copy at
//           https://www.boost.org/LICENSE_1_0.txt)

#include "test_autodiff.hpp"

BOOST_AUTO_TEST_SUITE(test_autodiff_3)

struct atanh_test_test
{
  template<typename T>
  void operator()(const T&) const
  {
    const T eps = 300*std::numeric_limits<T>::epsilon(); // percent
    using std::atanh;
    constexpr int m = 5;
    const T cx = 0.5;
    auto x = make_fvar<T,m>(cx);
    auto y = atanh(x);
    // BOOST_REQUIRE(y.derivative(0) == atanh(cx)); // fails due to overload
    BOOST_REQUIRE(y.derivative(0) == atanh(static_cast<T>(x)));
    BOOST_REQUIRE_CLOSE(y.derivative(1), static_cast<T>(4)/3, eps);
    BOOST_REQUIRE_CLOSE(y.derivative(2), static_cast<T>(16)/9, eps);
    BOOST_REQUIRE_CLOSE(y.derivative(3), static_cast<T>(224)/27, eps);
    BOOST_REQUIRE_CLOSE(y.derivative(4), static_cast<T>(1280)/27, eps);
    BOOST_REQUIRE_CLOSE(y.derivative(5), static_cast<T>(31232)/81, eps);
  }
};

BOOST_AUTO_TEST_CASE(atanh_test)
{
    boost::fusion::for_each(bin_float_types, atanh_test_test());
}

struct atan_test_test
{
  template<typename T>
  void operator()(const T&) const
  {
    constexpr int m = 5;
    constexpr float cx = 1.0;
    const auto x = make_fvar<T,m>(cx);
    auto y = atan(x);
    BOOST_REQUIRE(y.derivative(0) == boost::math::constants::pi<T>()/4);
    BOOST_REQUIRE(y.derivative(1) == 0.5);
    BOOST_REQUIRE(y.derivative(2) == -0.5);
    BOOST_REQUIRE(y.derivative(3) == 0.5);
    BOOST_REQUIRE(y.derivative(4) == 0.0);
    BOOST_REQUIRE(y.derivative(5) == -3.0);
  }
};

BOOST_AUTO_TEST_CASE(atan_test)
{
    boost::fusion::for_each(bin_float_types, atan_test_test());
    boost::fusion::for_each(multiprecision_float_types, atan_test_test());
}

struct erf_test_test
{
  template<typename T>
  void operator()(const T&) const
  {
    const T eps = 300*std::numeric_limits<T>::epsilon(); // percent
    using std::erf;
    using namespace boost;
    constexpr int m = 5;
    constexpr float cx = 1.0;
    const auto x = make_fvar<T,m>(cx);
    auto y = erf(x);
    BOOST_REQUIRE(y.derivative(0) == erf(static_cast<T>(x)));
    BOOST_REQUIRE_CLOSE(y.derivative(1), 2/(math::constants::e<T>()*math::constants::root_pi<T>()), eps);
    BOOST_REQUIRE_CLOSE(y.derivative(2), -4/(math::constants::e<T>()*math::constants::root_pi<T>()), eps);
    BOOST_REQUIRE_CLOSE(y.derivative(3), 4/(math::constants::e<T>()*math::constants::root_pi<T>()), eps);
    BOOST_REQUIRE_CLOSE(y.derivative(4), 8/(math::constants::e<T>()*math::constants::root_pi<T>()), eps);
    BOOST_REQUIRE_CLOSE(y.derivative(5), -40/(math::constants::e<T>()*math::constants::root_pi<T>()), eps);
  }
};

BOOST_AUTO_TEST_CASE(erf_test)
{
    boost::fusion::for_each(bin_float_types, erf_test_test());
    boost::fusion::for_each(multiprecision_float_types, erf_test_test());
}

struct sinc_test_test
{
  template<typename T>
  void operator()(const T&) const
  {
    const T eps = 20000*std::numeric_limits<T>::epsilon(); // percent
    using std::sin;
    using std::cos;
    constexpr int m = 5;
    const T cx = 1;
    auto x = make_fvar<T,m>(cx);
    auto y = sinc(x);
    BOOST_REQUIRE_CLOSE(y.derivative(0), sin(cx), eps);
    BOOST_REQUIRE_CLOSE(y.derivative(1), cos(cx)-sin(cx), eps);
    BOOST_REQUIRE_CLOSE(y.derivative(2), sin(cx)-2*cos(cx), eps);
    BOOST_REQUIRE_CLOSE(y.derivative(3), 5*cos(cx)-3*sin(cx), eps);
    BOOST_REQUIRE_CLOSE(y.derivative(4), 13*sin(cx)-20*cos(cx), eps);
    BOOST_REQUIRE_CLOSE(y.derivative(5), 101*cos(cx)-65*sin(cx), eps);
    // Test at x = 0
    auto y2 = sinc(make_fvar<T,10>(0));
    BOOST_REQUIRE_CLOSE(y2.derivative(0), 1, eps);
    BOOST_REQUIRE_CLOSE(y2.derivative(1), 0, eps);
    BOOST_REQUIRE_CLOSE(y2.derivative(2), -cx/3, eps);
    BOOST_REQUIRE_CLOSE(y2.derivative(3), 0, eps);
    BOOST_REQUIRE_CLOSE(y2.derivative(4), cx/5, eps);
    BOOST_REQUIRE_CLOSE(y2.derivative(5), 0, eps);
    BOOST_REQUIRE_CLOSE(y2.derivative(6), -cx/7, eps);
    BOOST_REQUIRE_CLOSE(y2.derivative(7), 0, eps);
    BOOST_REQUIRE_CLOSE(y2.derivative(8), cx/9, eps);
    BOOST_REQUIRE_CLOSE(y2.derivative(9), 0, eps);
    BOOST_REQUIRE_CLOSE(y2.derivative(10), -cx/11, eps);
  }
};

BOOST_AUTO_TEST_CASE(sinc_test)
{
    boost::fusion::for_each(bin_float_types, sinc_test_test());
}

struct sinh_and_cosh_test
{
  template<typename T>
  void operator()(const T&) const
  {
    const T eps = 300*std::numeric_limits<T>::epsilon(); // percent
    using std::sinh;
    using std::cosh;
    constexpr int m = 5;
    const T cx = 1;
    auto x = make_fvar<T,m>(cx);
    auto s = sinh(x);
    auto c = cosh(x);
    BOOST_REQUIRE_CLOSE(s.derivative(0), sinh(static_cast<T>(x)), eps);
    BOOST_REQUIRE_CLOSE(c.derivative(0), cosh(static_cast<T>(x)), eps);
    for (size_t i=0 ; i<=m ; ++i)
    {
        BOOST_REQUIRE_CLOSE(s.derivative(i), static_cast<T>(i&1?c:s), eps);
        BOOST_REQUIRE_CLOSE(c.derivative(i), static_cast<T>(i&1?s:c), eps);
    }
  }
};

BOOST_AUTO_TEST_CASE(sinh_and_cosh)
{
    boost::fusion::for_each(bin_float_types, sinh_and_cosh_test());
}

struct tan_test_test
{
  template<typename T>
  void operator()(const T&) const
  {
    const T eps = 800*std::numeric_limits<T>::epsilon(); // percent
    using std::sqrt;
    constexpr int m = 5;
    const T cx = boost::math::constants::third_pi<T>();
    const T root_three = boost::math::constants::root_three<T>();
    const auto x = make_fvar<T,m>(cx);
    auto y = tan(x);
    BOOST_REQUIRE_CLOSE(y.derivative(0), root_three, eps);
    BOOST_REQUIRE_CLOSE(y.derivative(1), 4.0, eps);
    BOOST_REQUIRE_CLOSE(y.derivative(2), 8*root_three, eps);
    BOOST_REQUIRE_CLOSE(y.derivative(3), 80.0, eps);
    BOOST_REQUIRE_CLOSE(y.derivative(4), 352*root_three, eps);
    BOOST_REQUIRE_CLOSE(y.derivative(5), 5824.0, eps);
  }
};

BOOST_AUTO_TEST_CASE(tan_test)
{
    boost::fusion::for_each(bin_float_types, tan_test_test());
}

struct fmod_test_test
{
  template<typename T>
  void operator()(const T&) const
  {
    constexpr int m = 3;
    constexpr float cx = 3.25;
    const T cy = 0.5;
    auto x = make_fvar<T,m>(cx);
    auto y = fmod(x,autodiff_fvar<T,m>(cy));
    BOOST_REQUIRE(y.derivative(0) == 0.25);
    BOOST_REQUIRE(y.derivative(1) == 1.0);
    BOOST_REQUIRE(y.derivative(2) == 0.0);
    BOOST_REQUIRE(y.derivative(3) == 0.0);
  }
};

BOOST_AUTO_TEST_CASE(fmod_test)
{
    boost::fusion::for_each(bin_float_types, fmod_test_test());
    boost::fusion::for_each(multiprecision_float_types, fmod_test_test());
}

struct round_and_trunc_test
{
  template<typename T>
  void operator()(const T&) const
  {
    using std::round;
    using std::trunc;
    constexpr int m = 3;
    constexpr float cx = 3.25;
    auto x = make_fvar<T,m>(cx);
    auto y = round(x);
    BOOST_REQUIRE(y.derivative(0) == round(cx));
    BOOST_REQUIRE(y.derivative(1) == 0.0);
    BOOST_REQUIRE(y.derivative(2) == 0.0);
    BOOST_REQUIRE(y.derivative(3) == 0.0);
    y = trunc(x);
    BOOST_REQUIRE(y.derivative(0) == trunc(cx));
    BOOST_REQUIRE(y.derivative(1) == 0.0);
    BOOST_REQUIRE(y.derivative(2) == 0.0);
    BOOST_REQUIRE(y.derivative(3) == 0.0);
  }
};

BOOST_AUTO_TEST_CASE(round_and_trunc)
{
    boost::fusion::for_each(bin_float_types, round_and_trunc_test());
    boost::fusion::for_each(multiprecision_float_types, round_and_trunc_test());
}

struct iround_and_itrunc_test
{
  template<typename T>
  void operator()(const T&) const
  {
    using namespace boost::math;
    constexpr int m = 3;
    constexpr float cx = 3.25;
    auto x = make_fvar<T,m>(cx);
    int y = iround(x);
    BOOST_REQUIRE(y == iround(cx));
    y = itrunc(x);
    BOOST_REQUIRE(y == itrunc(cx));
  }
};

BOOST_AUTO_TEST_CASE(iround_and_itrunc)
{
    boost::fusion::for_each(bin_float_types, iround_and_itrunc_test());
    boost::fusion::for_each(multiprecision_float_types, iround_and_itrunc_test());
}

BOOST_AUTO_TEST_SUITE_END()
