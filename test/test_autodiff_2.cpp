//           Copyright Matthew Pulver 2018 - 2019.
// Distributed under the Boost Software License, Version 1.0.
//      (See accompanying file LICENSE_1_0.txt or copy at
//           https://www.boost.org/LICENSE_1_0.txt)

#include "test_autodiff.hpp"

BOOST_AUTO_TEST_SUITE(test_autodiff_2)

struct one_over_one_plus_x_squared_test
{
  template<typename T>
  void operator()(const T&) const
  {
    constexpr int m = 4;
    constexpr float cx = 1.0;
    auto f = make_fvar<T,m>(cx);
    //f = 1 / ((f *= f) += 1);
    f = ((f *= f) += 1).inverse();
    BOOST_REQUIRE(f.derivative(0) == 0.5);
    BOOST_REQUIRE(f.derivative(1) == -0.5);
    BOOST_REQUIRE(f.derivative(2) == 0.5);
    BOOST_REQUIRE(f.derivative(3) == 0.0);
    BOOST_REQUIRE(f.derivative(4) == -3.0);
  }
};

BOOST_AUTO_TEST_CASE(one_over_one_plus_x_squared)
{
    boost::fusion::for_each(bin_float_types, one_over_one_plus_x_squared_test());
    boost::fusion::for_each(multiprecision_float_types, one_over_one_plus_x_squared_test());
}

struct exp_test_test
{
  template<typename T>
  void operator()(const T&) const
  {
    using std::exp;
    constexpr int m = 4;
    const T cx = 2.0;
    const auto x = make_fvar<T,m>(cx);
    auto y = exp(x);
    for (int i=0 ; i<=m ; ++i)
    {
        //std::cout.precision(100);
        //std::cout << "y.derivative("<<i<<") = " << y.derivative(i) << ", std::exp(cx) = " << std::exp(cx) << std::endl;
        BOOST_REQUIRE(y.derivative(i) == exp(cx));
    }
  }
};

BOOST_AUTO_TEST_CASE(exp_test)
{
    boost::fusion::for_each(bin_float_types, exp_test_test());
    boost::fusion::for_each(multiprecision_float_types, exp_test_test());
}

struct pow_test_test
{
  template<typename T>
  void operator()(const T&) const
  {
    const T eps = 201*std::numeric_limits<T>::epsilon(); // percent
    using std::exp;
    using std::log;
    using std::pow;
    constexpr int m = 5;
    constexpr int n = 4;
    const T cx = 2.0;
    const T cy = 3.0;
    const auto x = make_fvar<T,m>(cx);
    const auto y = make_fvar<T,m,n>(cy);
    auto z0 = pow(x,cy);
    BOOST_REQUIRE(z0.derivative(0) == pow(cx,cy));
    BOOST_REQUIRE(z0.derivative(1) == cy*pow(cx,cy-1));
    BOOST_REQUIRE(z0.derivative(2) == cy*(cy-1)*pow(cx,cy-2));
    BOOST_REQUIRE(z0.derivative(3) == cy*(cy-1)*(cy-2)*pow(cx,cy-3));
    BOOST_REQUIRE(z0.derivative(4) == 0.0);
    BOOST_REQUIRE(z0.derivative(5) == 0.0);
    auto z1 = pow(cx,y);
    BOOST_REQUIRE_CLOSE(z1.derivative(0,0), pow(cx,cy), eps);
    for (int j=1 ; j<=n ; ++j)
        BOOST_REQUIRE_CLOSE(z1.derivative(0,j), pow(log(cx),j)*exp(cy*log(cx)), eps);
    for (int i=1 ; i<=m ; ++i)
        for (int j=0 ; j<=n ; ++j)
            BOOST_REQUIRE(z1.derivative(i,j) == 0.0);
    auto z2 = pow(x,y);
    for (int j=0 ; j<=n ; ++j)
        BOOST_REQUIRE_CLOSE(z2.derivative(0,j), pow(cx,cy)*pow(log(cx),j), eps);
    for (int j=0 ; j<=n ; ++j)
        BOOST_REQUIRE_CLOSE(z2.derivative(1,j), pow(cx,cy-1)*pow(log(cx),j-1)*(cy*log(cx)+j), eps);
    BOOST_REQUIRE_CLOSE(z2.derivative(2,0), pow(cx,cy-2)*cy*(cy-1), eps);
    BOOST_REQUIRE_CLOSE(z2.derivative(2,1), pow(cx,cy-2)*(cy*(cy-1)*log(cx)+2*cy-1), eps);
    for (int j=2 ; j<=n ; ++j)
        BOOST_REQUIRE_CLOSE(z2.derivative(2,j), pow(cx,cy-2)*pow(log(cx),j-2)*(j*(2*cy-1)*log(cx)+(j-1)*j+(cy-1)*cy*pow(log(cx),2)), eps);
    BOOST_REQUIRE_CLOSE(z2.derivative(2,4), pow(cx,cy-2)*pow(log(cx),2)*(4*(2*cy-1)*log(cx)+(4-1)*4+(cy-1)*cy*pow(log(cx),2)), eps);
  }
};

BOOST_AUTO_TEST_CASE(pow_test)
{
    boost::fusion::for_each(bin_float_types, pow_test_test());
}

struct sqrt_test_test
{
  template<typename T>
  void operator()(const T&) const
  {
    using std::sqrt;
    using std::pow;
    constexpr int m = 5;
    constexpr float cx = 4.0;
    auto x = make_fvar<T,m>(cx);
    auto y = sqrt(x);
    BOOST_REQUIRE(y.derivative(0) == sqrt(cx));
    BOOST_REQUIRE(y.derivative(1) == 0.5*pow(cx,-0.5));
    BOOST_REQUIRE(y.derivative(2) == -0.5*0.5*pow(cx,-1.5));
    BOOST_REQUIRE(y.derivative(3) == 0.5*0.5*1.5*pow(cx,-2.5));
    BOOST_REQUIRE(y.derivative(4) == -0.5*0.5*1.5*2.5*pow(cx,-3.5));
    BOOST_REQUIRE(y.derivative(5) == 0.5*0.5*1.5*2.5*3.5*pow(cx,-4.5));
    x = make_fvar<T,m>(0);
    y = sqrt(x);
    //std::cout << "sqrt(0) = " << y << std::endl; // (0,inf,-inf,inf,-inf,inf)
    BOOST_REQUIRE(y.derivative(0) == 0.0);
    for (int i=1; i<=m ; ++i)
        BOOST_REQUIRE(y.derivative(i) == (i&1?1:-1)*std::numeric_limits<T>::infinity());
  }
};

BOOST_AUTO_TEST_CASE(sqrt_test)
{
    boost::fusion::for_each(bin_float_types, sqrt_test_test());
    boost::fusion::for_each(multiprecision_float_types, sqrt_test_test());
}

struct log_test_test
{
  template<typename T>
  void operator()(const T&) const
  {
    using std::log;
    using std::pow;
    constexpr int m = 5;
    const T cx = 2.0;
    auto x = make_fvar<T,m>(cx);
    auto y = log(x);
    BOOST_REQUIRE(y.derivative(0) == log(cx));
    BOOST_REQUIRE(y.derivative(1) == 1/cx);
    BOOST_REQUIRE(y.derivative(2) == -1/pow(cx,2));
    BOOST_REQUIRE(y.derivative(3) == 2/pow(cx,3));
    BOOST_REQUIRE(y.derivative(4) == -6/pow(cx,4));
    BOOST_REQUIRE(y.derivative(5) == 24/pow(cx,5));
    x = make_fvar<T,m>(0);
    y = log(x);
    //std::cout << "log(0) = " << y << std::endl; // log(0) = depth(1)(-inf,inf,-inf,inf,-inf,inf)
    for (int i=0; i<=m ; ++i)
        BOOST_REQUIRE(y.derivative(i) == (i&1?1:-1)*std::numeric_limits<T>::infinity());
  }
};

BOOST_AUTO_TEST_CASE(log_test)
{
    boost::fusion::for_each(bin_float_types, log_test_test());
    boost::fusion::for_each(multiprecision_float_types, log_test_test());
}

struct ylogx_test
{
  template<typename T>
  void operator()(const T&) const
  {
    using std::log;
    using std::pow;
    const T eps = 100*std::numeric_limits<T>::epsilon(); // percent
    constexpr int m = 5;
    constexpr int n = 4;
    const T cx = 2.0;
    const T cy = 3.0;
    const auto x = make_fvar<T,m>(cx);
    const auto y = make_fvar<T,m,n>(cy);
    auto z = y*log(x);
    BOOST_REQUIRE(z.derivative(0,0) == cy*log(cx));
    BOOST_REQUIRE(z.derivative(0,1) == log(cx));
    BOOST_REQUIRE(z.derivative(0,2) == 0.0);
    BOOST_REQUIRE(z.derivative(0,3) == 0.0);
    BOOST_REQUIRE(z.derivative(0,4) == 0.0);
    for (size_t i=1 ; i<=m ; ++i)
        BOOST_REQUIRE_CLOSE(z.derivative(i,0), pow(-1,i-1)*boost::math::factorial<T>(i-1)*cy/pow(cx,i), eps);
    for (size_t i=1 ; i<=m ; ++i)
        BOOST_REQUIRE_CLOSE(z.derivative(i,1), pow(-1,i-1)*boost::math::factorial<T>(i-1)/pow(cx,i), eps);
    for (size_t i=1 ; i<=m ; ++i)
        for (size_t j=2 ; j<=n ; ++j)
            BOOST_REQUIRE(z.derivative(i,j) == 0.0);
    auto z1 = exp(z);
    // RHS is confirmed by
    // https://www.wolframalpha.com/input/?i=D%5Bx%5Ey,%7Bx,2%7D,%7By,4%7D%5D+%2F.+%7Bx-%3E2.0,+y-%3E3.0%7D
    BOOST_REQUIRE_CLOSE(z1.derivative(2,4),
        pow(cx,cy-2)*pow(log(cx),2)*(4*(2*cy-1)*log(cx)+(4-1)*4+(cy-1)*cy*pow(log(cx),2)), eps);
  }
};

BOOST_AUTO_TEST_CASE(ylogx)
{
    boost::fusion::for_each(bin_float_types, ylogx_test());
}

struct frexp_test_test
{
  template<typename T>
  void operator()(const T&) const
  {
    using std::frexp;
    using std::exp2;
    constexpr int m = 3;
    const T cx = 3.5;
    const auto x = make_fvar<T,m>(cx);
    int exp, testexp;
    auto y = frexp(x,&exp);
    BOOST_REQUIRE(y.derivative(0) == frexp(cx,&testexp));
    BOOST_REQUIRE(exp == testexp);
    BOOST_REQUIRE(y.derivative(1) == exp2(-exp));
    BOOST_REQUIRE(y.derivative(2) == 0.0);
    BOOST_REQUIRE(y.derivative(3) == 0.0);
  }
};

BOOST_AUTO_TEST_CASE(frexp_test)
{
    boost::fusion::for_each(bin_float_types, frexp_test_test());
    boost::fusion::for_each(multiprecision_float_types, frexp_test_test());
}

struct ldexp_test_test
{
  template<typename T>
  void operator()(const T&) const
  {
    using std::ldexp;
    using std::exp2;
    constexpr int m = 3;
    const T cx = 3.5;
    const auto x = make_fvar<T,m>(cx);
    constexpr int exp = 3;
    auto y = ldexp(x,exp);
    BOOST_REQUIRE(y.derivative(0) == ldexp(cx,exp));
    BOOST_REQUIRE(y.derivative(1) == exp2(exp));
    BOOST_REQUIRE(y.derivative(2) == 0.0);
    BOOST_REQUIRE(y.derivative(3) == 0.0);
  }
};

BOOST_AUTO_TEST_CASE(ldexp_test)
{
    boost::fusion::for_each(bin_float_types, ldexp_test_test());
    boost::fusion::for_each(multiprecision_float_types, ldexp_test_test());
}

struct cos_and_sin_test
{
  template<typename T>
  void operator()(const T&) const
  {
    using std::cos;
    using std::sin;
    const T eps = 200*std::numeric_limits<T>::epsilon(); // percent
    constexpr int m = 5;
    const T cx = boost::math::constants::third_pi<T>();
    const auto x = make_fvar<T,m>(cx);
    auto cos5 = cos(x);
    BOOST_REQUIRE_CLOSE(cos5.derivative(0), cos(cx), eps);
    BOOST_REQUIRE_CLOSE(cos5.derivative(1), -sin(cx), eps);
    BOOST_REQUIRE_CLOSE(cos5.derivative(2), -cos(cx), eps);
    BOOST_REQUIRE_CLOSE(cos5.derivative(3), sin(cx), eps);
    BOOST_REQUIRE_CLOSE(cos5.derivative(4), cos(cx), eps);
    BOOST_REQUIRE_CLOSE(cos5.derivative(5), -sin(cx), eps);
    auto sin5 = sin(x);
    BOOST_REQUIRE_CLOSE(sin5.derivative(0), sin(cx), eps);
    BOOST_REQUIRE_CLOSE(sin5.derivative(1), cos(cx), eps);
    BOOST_REQUIRE_CLOSE(sin5.derivative(2), -sin(cx), eps);
    BOOST_REQUIRE_CLOSE(sin5.derivative(3), -cos(cx), eps);
    BOOST_REQUIRE_CLOSE(sin5.derivative(4), sin(cx), eps);
    BOOST_REQUIRE_CLOSE(sin5.derivative(5), cos(cx), eps);
    // Test Order = 0 for codecov
    auto cos0 = cos(make_fvar<T,0>(cx));
    BOOST_REQUIRE_CLOSE(cos0.derivative(0), cos(cx), eps);
    auto sin0 = sin(make_fvar<T,0>(cx));
    BOOST_REQUIRE_CLOSE(sin0.derivative(0), sin(cx), eps);
  }
};

BOOST_AUTO_TEST_CASE(cos_and_sin)
{
    boost::fusion::for_each(bin_float_types, cos_and_sin_test());
}

struct acos_test_test
{
  template<typename T>
  void operator()(const T&) const
  {
    const T eps = 300*std::numeric_limits<T>::epsilon(); // percent
    using std::acos;
    using std::pow;
    using std::sqrt;
    constexpr int m = 5;
    const T cx = 0.5;
    auto x = make_fvar<T,m>(cx);
    auto y = acos(x);
    BOOST_REQUIRE_CLOSE(y.derivative(0), acos(cx), eps);
    BOOST_REQUIRE_CLOSE(y.derivative(1), -1/sqrt(1-cx*cx), eps);
    BOOST_REQUIRE_CLOSE(y.derivative(2), -cx/pow(1-cx*cx,1.5), eps);
    BOOST_REQUIRE_CLOSE(y.derivative(3), -(2*cx*cx+1)/pow(1-cx*cx,2.5), eps);
    BOOST_REQUIRE_CLOSE(y.derivative(4), -3*cx*(2*cx*cx+3)/pow(1-cx*cx,3.5), eps);
    BOOST_REQUIRE_CLOSE(y.derivative(5), -(24*(cx*cx+3)*cx*cx+9)/pow(1-cx*cx,4.5), eps);
  }
};

BOOST_AUTO_TEST_CASE(acos_test)
{
    boost::fusion::for_each(bin_float_types, acos_test_test());
}

struct acosh_test_test
{
  template<typename T>
  void operator()(const T&) const
  {
    const T eps = 300*std::numeric_limits<T>::epsilon(); // percent
    using std::acosh;
    constexpr int m = 5;
    const T cx = 2;
    auto x = make_fvar<T,m>(cx);
    auto y = acosh(x);
    //BOOST_REQUIRE(y.derivative(0) == acosh(cx)); // FAILS! acosh(2) is overloaded for integral types
    BOOST_REQUIRE(y.derivative(0) == acosh(static_cast<T>(x)));
    BOOST_REQUIRE_CLOSE(y.derivative(1), 1/boost::math::constants::root_three<T>(), eps);
    BOOST_REQUIRE_CLOSE(y.derivative(2), -2/(3*boost::math::constants::root_three<T>()), eps);
    BOOST_REQUIRE_CLOSE(y.derivative(3), 1/boost::math::constants::root_three<T>(), eps);
    BOOST_REQUIRE_CLOSE(y.derivative(4), -22/(9*boost::math::constants::root_three<T>()), eps);
    BOOST_REQUIRE_CLOSE(y.derivative(5), 227/(27*boost::math::constants::root_three<T>()), eps);
  }
};

BOOST_AUTO_TEST_CASE(acosh_test)
{
    boost::fusion::for_each(bin_float_types, acosh_test_test());
}

struct asin_test_test
{
  template<typename T>
  void operator()(const T&) const
  {
    const T eps = 300*std::numeric_limits<T>::epsilon(); // percent
    using std::asin;
    using std::pow;
    using std::sqrt;
    constexpr int m = 5;
    const T cx = 0.5;
    auto x = make_fvar<T,m>(cx);
    auto y = asin(x);
    BOOST_REQUIRE_CLOSE(y.derivative(0), asin(cx), eps);
    BOOST_REQUIRE_CLOSE(y.derivative(1), 1/sqrt(1-cx*cx), eps);
    BOOST_REQUIRE_CLOSE(y.derivative(2), cx/pow(1-cx*cx,1.5), eps);
    BOOST_REQUIRE_CLOSE(y.derivative(3), (2*cx*cx+1)/pow(1-cx*cx,2.5), eps);
    BOOST_REQUIRE_CLOSE(y.derivative(4), 3*cx*(2*cx*cx+3)/pow(1-cx*cx,3.5), eps);
    BOOST_REQUIRE_CLOSE(y.derivative(5), (24*(cx*cx+3)*cx*cx+9)/pow(1-cx*cx,4.5), eps);
  }
};

BOOST_AUTO_TEST_CASE(asin_test)
{
    boost::fusion::for_each(bin_float_types, asin_test_test());
}

struct asin_infinity_test
{
  template<typename T>
  void operator()(const T&) const
  {
    const T eps = 100*std::numeric_limits<T>::epsilon(); // percent
    constexpr int m = 5;
    auto x = make_fvar<T,m>(1);
    auto y = asin(x);
    //std::cout << "asin(1) = " << y << std::endl; // depth(1)(1.5707963267949,inf,inf,-nan,-nan,-nan)
    BOOST_REQUIRE_CLOSE(y.derivative(0), boost::math::constants::half_pi<T>(), eps); // MacOS is not exact
    BOOST_REQUIRE(y.derivative(1) == std::numeric_limits<T>::infinity());
  }
};

BOOST_AUTO_TEST_CASE(asin_infinity)
{
    boost::fusion::for_each(bin_float_types, asin_infinity_test());
    boost::fusion::for_each(multiprecision_float_types, asin_infinity_test());
}

struct asin_derivative_test
{
  template<typename T>
  void operator()(const T&) const
  {
    const T eps = 300*std::numeric_limits<T>::epsilon(); // percent
    using std::pow;
    using std::sqrt;
    constexpr int m = 4;
    const T cx = 0.5;
    auto x = make_fvar<T,m>(cx);
    auto y = 1-x*x;
    BOOST_REQUIRE(y.derivative(0) == 1-cx*cx);
    BOOST_REQUIRE(y.derivative(1) == -2*cx);
    BOOST_REQUIRE(y.derivative(2) == -2);
    BOOST_REQUIRE(y.derivative(3) == 0);
    BOOST_REQUIRE(y.derivative(4) == 0);
    y = sqrt(y);
    BOOST_REQUIRE(y.derivative(0) == sqrt(1-cx*cx));
    BOOST_REQUIRE_CLOSE(y.derivative(1), -cx/sqrt(1-cx*cx), eps);
    BOOST_REQUIRE_CLOSE(y.derivative(2), -1/pow(1-cx*cx,1.5), eps);
    BOOST_REQUIRE_CLOSE(y.derivative(3), -3*cx/pow(1-cx*cx,2.5), eps);
    BOOST_REQUIRE_CLOSE(y.derivative(4), -(12*cx*cx+3)/pow(1-cx*cx,3.5), eps);
    y = y.inverse(); // asin'(x) = 1 / sqrt(1-x*x).
    BOOST_REQUIRE_CLOSE(y.derivative(0), 1/sqrt(1-cx*cx), eps);
    BOOST_REQUIRE_CLOSE(y.derivative(1), cx/pow(1-cx*cx,1.5), eps);
    BOOST_REQUIRE_CLOSE(y.derivative(2), (2*cx*cx+1)/pow(1-cx*cx,2.5), eps);
    BOOST_REQUIRE_CLOSE(y.derivative(3), 3*cx*(2*cx*cx+3)/pow(1-cx*cx,3.5), eps);
    BOOST_REQUIRE_CLOSE(y.derivative(4), (24*(cx*cx+3)*cx*cx+9)/pow(1-cx*cx,4.5), eps);
  }
};

BOOST_AUTO_TEST_CASE(asin_derivative)
{
    boost::fusion::for_each(bin_float_types, asin_derivative_test());
}

struct asinh_test_test
{
  template<typename T>
  void operator()(const T&) const
  {
    const T eps = 300*std::numeric_limits<T>::epsilon(); // percent
    using std::asinh;
    constexpr int m = 5;
    const T cx = 1;
    auto x = make_fvar<T,m>(cx);
    auto y = asinh(x);
    //BOOST_REQUIRE(y.derivative(0) == asinh(cx)); // Fails for gcc-mingw - similar to acosh()?
    BOOST_REQUIRE(y.derivative(0) == asinh(static_cast<T>(x)));
    BOOST_REQUIRE_CLOSE(y.derivative(1), 1/boost::math::constants::root_two<T>(), eps);
    BOOST_REQUIRE_CLOSE(y.derivative(2), -1/(2*boost::math::constants::root_two<T>()), eps);
    BOOST_REQUIRE_CLOSE(y.derivative(3), 1/(4*boost::math::constants::root_two<T>()), eps);
    BOOST_REQUIRE_CLOSE(y.derivative(4), 3/(8*boost::math::constants::root_two<T>()), eps);
    BOOST_REQUIRE_CLOSE(y.derivative(5), -39/(16*boost::math::constants::root_two<T>()), eps);
  }
};

BOOST_AUTO_TEST_CASE(asinh_test)
{
    boost::fusion::for_each(bin_float_types, asinh_test_test());
}

BOOST_AUTO_TEST_SUITE_END()
