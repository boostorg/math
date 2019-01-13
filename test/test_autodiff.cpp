//           Copyright Matthew Pulver 2018 - 2019.
// Distributed under the Boost Software License, Version 1.0.
//      (See accompanying file LICENSE_1_0.txt or copy at
//           https://www.boost.org/LICENSE_1_0.txt)

#include <boost/math/differentiation/autodiff.hpp>

#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/special_functions/fpclassify.hpp> // isnan
#include <boost/math/special_functions/round.hpp> // iround
#include <boost/math/special_functions/trunc.hpp> // itrunc
#include <boost/multiprecision/cpp_dec_float.hpp>

#define BOOST_TEST_MODULE test_autodiff
#include <boost/test/included/unit_test.hpp>

#include <iostream>

template<typename W,typename X,typename Y,typename Z>
boost::math::differentiation::autodiff::promote<W,X,Y,Z>
    mixed_partials_f(const W& w, const X& x, const Y& y, const Z& z)
{
    using namespace std;
    return exp(w*sin(x*log(y)/z) + sqrt(w*z/(x*y))) + w*w/tan(z);
}

// Equations and function/variable names are from
// https://en.wikipedia.org/wiki/Greeks_(finance)#Formulas_for_European_option_Greeks
//
// Standard normal probability density function
template<typename X>
X phi(const X& x)
{
  return boost::math::constants::one_div_root_two_pi<double>()*exp(-0.5*x*x);
}

// Standard normal cumulative distribution function
template<typename X>
X Phi(const X& x)
{
  return 0.5*erfc(-boost::math::constants::one_div_root_two<double>()*x);
}

enum CP { call, put };

// Assume zero annual dividend yield (q=0).
template<typename Price,typename Sigma,typename Tau,typename Rate>
boost::math::differentiation::autodiff::promote<Price,Sigma,Tau,Rate>
    black_scholes_option_price(CP cp, double K, const Price& S, const Sigma& sigma, const Tau& tau, const Rate& r)
{
  using namespace std;
  const auto d1 = (log(S/K) + (r+sigma*sigma/2)*tau) / (sigma*sqrt(tau));
  const auto d2 = (log(S/K) + (r-sigma*sigma/2)*tau) / (sigma*sqrt(tau));
  static_assert(std::is_same<decltype(S*Phi(d1) - exp(-r*tau)*K*Phi(d2)),
    decltype(exp(-r*tau)*K*Phi(-d2) - S*Phi(-d1))>::value, "decltype(call) != decltype(put)");
  if (cp == call)
    return S*Phi(d1) - exp(-r*tau)*K*Phi(d2);
  else
    return exp(-r*tau)*K*Phi(-d2) - S*Phi(-d1);
}

BOOST_AUTO_TEST_SUITE(test_autodiff)

BOOST_AUTO_TEST_CASE(constructors)
{
    constexpr int m = 3;
    constexpr int n = 4;
    // Verify value-initialized instance has all 0 entries.
    const boost::math::differentiation::autodiff::variable<double,m> empty1 = boost::math::differentiation::autodiff::variable<double,m>();
    for (int i=0 ; i<=m ; ++i)
        BOOST_REQUIRE(empty1.derivative(i) == 0.0);
    const auto empty2 = boost::math::differentiation::autodiff::variable<double,m,n>();
    for (int i=0 ; i<=m ; ++i)
        for (int j=0 ; j<=n ; ++j)
            BOOST_REQUIRE(empty2.derivative(i,j) == 0.0);
    // Single variable
    constexpr double cx = 10.0;
    const auto x = boost::math::differentiation::autodiff::variable<double,m>(cx);
    for (int i=0 ; i<=m ; ++i)
        if (i==0)
            BOOST_REQUIRE(x.derivative(i) == cx);
        else if (i==1)
            BOOST_REQUIRE(x.derivative(i) == 1.0);
        else
            BOOST_REQUIRE(x.derivative(i) == 0.0);
    const boost::math::differentiation::autodiff::variable<double,n> xn = x;
    for (int i=0 ; i<=n ; ++i)
        if (i==0)
            BOOST_REQUIRE(xn.derivative(i) == cx);
        else if (i==1)
            BOOST_REQUIRE(xn.derivative(i) == 1.0);
        else
            BOOST_REQUIRE(xn.derivative(i) == 0.0);
    // Second independent variable
    constexpr double cy = 100.0;
    const auto y = boost::math::differentiation::autodiff::variable<double,m,n>(cy);
    for (int i=0 ; i<=m ; ++i)
        for (int j=0 ; j<=n ; ++j)
            if (i==0 && j==0)
                BOOST_REQUIRE(y.derivative(i,j) == cy);
            else if (i==0 && j==1)
                BOOST_REQUIRE(y.derivative(i,j) == 1.0);
            else
                BOOST_REQUIRE(y.derivative(i,j) == 0.0);
}

BOOST_AUTO_TEST_CASE(assignment)
{
    constexpr int m = 3;
    constexpr int n = 4;
    constexpr double cx = 10.0;
    constexpr double cy = 10.0;
    boost::math::differentiation::autodiff::variable<double,m,n> empty; // Uninitialized variable<> may have non-zero values.
    // Single variable
    auto x = boost::math::differentiation::autodiff::variable<double,m>(cx);
    empty = static_cast<decltype(empty)>(x); // Test static_cast of single-variable to double-variable type.
    for (int i=0 ; i<=m ; ++i)
        for (int j=0 ; j<=n ; ++j)
            if (i==0 && j==0)
                BOOST_REQUIRE(empty.derivative(i,j) == cx);
            else if (i==1 && j==0)
                BOOST_REQUIRE(empty.derivative(i,j) == 1.0);
            else
                BOOST_REQUIRE(empty.derivative(i,j) == 0.0);
    auto y = boost::math::differentiation::autodiff::variable<double,m,n>(cy);
    empty = y; // default assignment operator
    for (int i=0 ; i<=m ; ++i)
        for (int j=0 ; j<=n ; ++j)
            if (i==0 && j==0)
                BOOST_REQUIRE(empty.derivative(i,j) == cy);
            else if (i==0 && j==1)
                BOOST_REQUIRE(empty.derivative(i,j) == 1.0);
            else
                BOOST_REQUIRE(empty.derivative(i,j) == 0.0);
    empty = cx; // set a constant
    for (int i=0 ; i<=m ; ++i)
        for (int j=0 ; j<=n ; ++j)
            if (i==0 && j==0)
                BOOST_REQUIRE(empty.derivative(i,j) == cx);
            else
                BOOST_REQUIRE(empty.derivative(i,j) == 0.0);
}

BOOST_AUTO_TEST_CASE(addition_assignment)
{
    constexpr int m = 3;
    constexpr int n = 4;
    constexpr double cx = 10.0;
    auto sum = boost::math::differentiation::autodiff::variable<double,m,n>(); // zero-initialized
    // Single variable
    const auto x = boost::math::differentiation::autodiff::variable<double,m>(cx);
    sum += x;
    for (int i=0 ; i<=m ; ++i)
        for (int j=0 ; j<=n ; ++j)
            if (i==0 && j==0)
                BOOST_REQUIRE(sum.derivative(i,j) == cx);
            else if (i==1 && j==0)
                BOOST_REQUIRE(sum.derivative(i,j) == 1.0);
            else
                BOOST_REQUIRE(sum.derivative(i,j) == 0.0);
    // Arithmetic constant
    constexpr double cy = 11.0;
    sum = 0;
    sum += cy;
    for (int i=0 ; i<=m ; ++i)
        for (int j=0 ; j<=n ; ++j)
            if (i==0 && j==0)
                BOOST_REQUIRE(sum.derivative(i,j) == cy);
            else
                BOOST_REQUIRE(sum.derivative(i,j) == 0.0);
}

BOOST_AUTO_TEST_CASE(subtraction_assignment)
{
    constexpr int m = 3;
    constexpr int n = 4;
    constexpr double cx = 10.0;
    auto sum = boost::math::differentiation::autodiff::variable<double,m,n>(); // zero-initialized
    // Single variable
    const auto x = boost::math::differentiation::autodiff::variable<double,m>(cx);
    sum -= x;
    for (int i=0 ; i<=m ; ++i)
        for (int j=0 ; j<=n ; ++j)
            if (i==0 && j==0)
                BOOST_REQUIRE(sum.derivative(i,j) == -cx);
            else if (i==1 && j==0)
                BOOST_REQUIRE(sum.derivative(i,j) == -1.0);
            else
                BOOST_REQUIRE(sum.derivative(i,j) == 0.0);
    // Arithmetic constant
    constexpr double cy = 11.0;
    sum = 0;
    sum -= cy;
    for (int i=0 ; i<=m ; ++i)
        for (int j=0 ; j<=n ; ++j)
            if (i==0 && j==0)
                BOOST_REQUIRE(sum.derivative(i,j) == -cy);
            else
                BOOST_REQUIRE(sum.derivative(i,j) == 0.0);
}

BOOST_AUTO_TEST_CASE(multiplication_assignment)
{
    constexpr int m = 3;
    constexpr int n = 4;
    constexpr double cx = 10.0;
    auto product = boost::math::differentiation::autodiff::variable<double,m,n>{1}; // unit constant
    // Single variable
    auto x = boost::math::differentiation::autodiff::variable<double,m>(cx);
    product *= x;
    for (int i=0 ; i<=m ; ++i)
        for (int j=0 ; j<=n ; ++j)
            if (i==0 && j==0)
                BOOST_REQUIRE(product.derivative(i,j) == cx);
            else if (i==1 && j==0)
                BOOST_REQUIRE(product.derivative(i,j) == 1.0);
            else
                BOOST_REQUIRE(product.derivative(i,j) == 0.0);
    // Arithmetic constant
    constexpr double cy = 11.0;
    product = 1;
    product *= cy;
    for (int i=0 ; i<=m ; ++i)
        for (int j=0 ; j<=n ; ++j)
            if (i==0 && j==0)
                BOOST_REQUIRE(product.derivative(i,j) == cy);
            else
                BOOST_REQUIRE(product.derivative(i,j) == 0.0);
    // 0 * inf = nan
    x = boost::math::differentiation::autodiff::variable<double,m>(0.0);
    x *= std::numeric_limits<double>::infinity();
    //std::cout << "x = " << x << std::endl;
    for (int i=0 ; i<=m ; ++i)
        if (i==0)
            BOOST_REQUIRE(boost::math::isnan(static_cast<double>(x))); // Correct
            //BOOST_REQUIRE(x.derivative(i) == 0.0); // Wrong. See multiply_assign_by_root_type().
        else if (i==1)
            BOOST_REQUIRE(boost::math::isinf(x.derivative(i)));
        else
            BOOST_REQUIRE(x.derivative(i) == 0.0);
}

BOOST_AUTO_TEST_CASE(division_assignment)
{
    constexpr int m = 3;
    constexpr int n = 4;
    constexpr double cx = 16.0;
    auto quotient = boost::math::differentiation::autodiff::variable<double,m,n>{1}; // unit constant
    // Single variable
    const auto x = boost::math::differentiation::autodiff::variable<double,m>(cx);
    quotient /= x;
    BOOST_REQUIRE(quotient.derivative(0,0) == 1/cx);
    BOOST_REQUIRE(quotient.derivative(1,0) == -1/std::pow(cx,2));
    BOOST_REQUIRE(quotient.derivative(2,0) == 2/std::pow(cx,3));
    BOOST_REQUIRE(quotient.derivative(3,0) == -6/std::pow(cx,4));
    for (int i=0 ; i<=m ; ++i)
        for (int j=1 ; j<=n ; ++j)
            BOOST_REQUIRE(quotient.derivative(i,j) == 0.0);
    // Arithmetic constant
    constexpr double cy = 32.0;
    quotient = 1;
    quotient /= cy;
    for (int i=0 ; i<=m ; ++i)
        for (int j=0 ; j<=n ; ++j)
            if (i==0 && j==0)
                BOOST_REQUIRE(quotient.derivative(i,j) == 1/cy);
            else
                BOOST_REQUIRE(quotient.derivative(i,j) == 0.0);
}

BOOST_AUTO_TEST_CASE(unary_signs)
{
    constexpr int m = 3;
    constexpr int n = 4;
    constexpr double cx = 16.0;
    boost::math::differentiation::autodiff::variable<double,m,n> lhs;
    // Single variable
    const auto x = boost::math::differentiation::autodiff::variable<double,m>(cx);
    lhs = static_cast<decltype(lhs)>(-x);
    for (int i=0 ; i<=m ; ++i)
        for (int j=0 ; j<=n ; ++j)
            if (i==0 && j==0)
                BOOST_REQUIRE(lhs.derivative(i,j) == -cx);
            else if (i==1 && j==0)
                BOOST_REQUIRE(lhs.derivative(i,j) == -1.0);
            else
                BOOST_REQUIRE(lhs.derivative(i,j) == 0.0);
    lhs = static_cast<decltype(lhs)>(+x);
    for (int i=0 ; i<=m ; ++i)
        for (int j=0 ; j<=n ; ++j)
            if (i==0 && j==0)
                BOOST_REQUIRE(lhs.derivative(i,j) == cx);
            else if (i==1 && j==0)
                BOOST_REQUIRE(lhs.derivative(i,j) == 1.0);
            else
                BOOST_REQUIRE(lhs.derivative(i,j) == 0.0);
}

// TODO 3 tests for 3 operator+() definitions.

BOOST_AUTO_TEST_CASE(cast_double)
{
    constexpr double ca = 13.0;
    constexpr int i = 12;
    constexpr int m = 3;
    const auto x = boost::math::differentiation::autodiff::variable<double,m>(ca);
    BOOST_REQUIRE(i < x);
    BOOST_REQUIRE(i*x == i*ca);
}

BOOST_AUTO_TEST_CASE(int_double_casting)
{
    constexpr double ca = 3.0;
    const auto x0 = boost::math::differentiation::autodiff::variable<double,0>(ca);
    BOOST_REQUIRE(static_cast<double>(x0) == ca);
    const auto x1 = boost::math::differentiation::autodiff::variable<double,1>(ca);
    BOOST_REQUIRE(static_cast<double>(x1) == ca);
    const auto x2 = boost::math::differentiation::autodiff::variable<double,2>(ca);
    BOOST_REQUIRE(static_cast<double>(x2) == ca);
}

BOOST_AUTO_TEST_CASE(scalar_addition)
{
    constexpr double ca = 3.0;
    constexpr double cb = 4.0;
    const auto sum0 = boost::math::differentiation::autodiff::variable<double,0>(ca) + boost::math::differentiation::autodiff::variable<double,0>(cb);
    BOOST_REQUIRE(ca+cb == static_cast<double>(sum0));
    const auto sum1 = boost::math::differentiation::autodiff::variable<double,0>(ca) + cb;
    BOOST_REQUIRE(ca+cb == static_cast<double>(sum1));
    const auto sum2 = ca + boost::math::differentiation::autodiff::variable<double,0>(cb);
    BOOST_REQUIRE(ca+cb == static_cast<double>(sum2));
}

BOOST_AUTO_TEST_CASE(power8)
{
    constexpr int n = 8;
    constexpr double ca = 3.0;
    auto x = boost::math::differentiation::autodiff::variable<double,n>(ca);
    // Test operator*=()
    x *= x;
    x *= x;
    x *= x;
    const double power_factorial = boost::math::factorial<double>(n);
    for (int i=0 ; i<=n ; ++i)
        BOOST_REQUIRE(x.derivative(i) == power_factorial/boost::math::factorial<double>(n-i)*std::pow(ca,n-i));
    x = boost::math::differentiation::autodiff::variable<double,n>(ca);
    // Test operator*()
    x = x*x*x*x * x*x*x*x;
    for (int i=0 ; i<=n ; ++i)
        BOOST_REQUIRE(x.derivative(i) == power_factorial/boost::math::factorial<double>(n-i)*std::pow(ca,n-i));
}

BOOST_AUTO_TEST_CASE(dim1_multiplication)
{
    constexpr int m = 2;
    constexpr int n = 3;
    constexpr double cy = 4.0;
    auto y0 = boost::math::differentiation::autodiff::variable<double,m>(cy);
    auto y  = boost::math::differentiation::autodiff::variable<double,n>(cy);
    y *= y0;
    BOOST_REQUIRE(y.derivative(0) == cy*cy);
    BOOST_REQUIRE(y.derivative(1) == 2*cy);
    BOOST_REQUIRE(y.derivative(2) == 2.0);
    BOOST_REQUIRE(y.derivative(3) == 0.0);
    y = y * cy;
    BOOST_REQUIRE(y.derivative(0) == cy*cy*cy);
    BOOST_REQUIRE(y.derivative(1) == 2*cy*cy);
    BOOST_REQUIRE(y.derivative(2) == 2.0*cy);
    BOOST_REQUIRE(y.derivative(3) == 0.0);
}

BOOST_AUTO_TEST_CASE(dim1and2_multiplication)
{
    constexpr int m = 2;
    constexpr int n = 3;
    constexpr double cx = 3.0;
    constexpr double cy = 4.0;
    auto x = boost::math::differentiation::autodiff::variable<double,m>(cx);
    auto y = boost::math::differentiation::autodiff::variable<double,m,n>(cy);
    y *= x;
    BOOST_REQUIRE(y.derivative(0,0) == cx*cy);
    BOOST_REQUIRE(y.derivative(0,1) == cx);
    BOOST_REQUIRE(y.derivative(1,0) == cy);
    BOOST_REQUIRE(y.derivative(1,1) == 1.0);
    for (int i=1 ; i<m ; ++i)
        for (int j=1 ; j<n ; ++j)
            if (i==1 && j==1)
                BOOST_REQUIRE(y.derivative(i,j) == 1.0);
            else
                BOOST_REQUIRE(y.derivative(i,j) == 0.0);
}

BOOST_AUTO_TEST_CASE(dim2_addition)
{
    constexpr int m = 2;
    constexpr int n = 3;
    constexpr double cx = 3.0;
    const auto x = boost::math::differentiation::autodiff::variable<double,m>(cx);
    BOOST_REQUIRE(x.derivative(0) == cx);
    BOOST_REQUIRE(x.derivative(1) == 1.0);
    BOOST_REQUIRE(x.derivative(2) == 0.0);
    constexpr double cy = 4.0;
    const auto y = boost::math::differentiation::autodiff::variable<double,m,n>(cy);
    BOOST_REQUIRE(static_cast<double>(y.derivative(0)) == cy);
    BOOST_REQUIRE(static_cast<double>(y.derivative(1)) == 0.0); // partial of y w.r.t. x.

    BOOST_REQUIRE(y.derivative(0,0) == cy);
    BOOST_REQUIRE(y.derivative(0,1) == 1.0);
    BOOST_REQUIRE(y.derivative(1,0) == 0.0);
    BOOST_REQUIRE(y.derivative(1,1) == 0.0);
    const auto z = x + y;
    BOOST_REQUIRE(z.derivative(0,0) == cx + cy);
    BOOST_REQUIRE(z.derivative(0,1) == 1.0);
    BOOST_REQUIRE(z.derivative(1,0) == 1.0);
    BOOST_REQUIRE(z.derivative(1,1) == 0.0);
    // The following 4 are unnecessarily more expensive than the previous 4.
    BOOST_REQUIRE(z.derivative(0).derivative(0) == cx + cy);
    BOOST_REQUIRE(z.derivative(0).derivative(1) == 1.0);
    BOOST_REQUIRE(z.derivative(1).derivative(0) == 1.0);
    BOOST_REQUIRE(z.derivative(1).derivative(1) == 0.0);
}

BOOST_AUTO_TEST_CASE(dim2_multiplication)
{
    constexpr int m = 3;
    constexpr int n = 4;
    constexpr double cx = 6.0;
    const auto x = boost::math::differentiation::autodiff::variable<double,m>(cx);
    constexpr double cy = 5.0;
    const auto y = boost::math::differentiation::autodiff::variable<double,0,n>(cy);
    const auto z = x*x * y*y*y;
    BOOST_REQUIRE(z.derivative(0,0) == cx*cx * cy*cy*cy); // x^2 * y^3
    BOOST_REQUIRE(z.derivative(0,1) == cx*cx * 3*cy*cy); // x^2 * 3y^2
    BOOST_REQUIRE(z.derivative(0,2) == cx*cx * 6*cy); // x^2 * 6y
    BOOST_REQUIRE(z.derivative(0,3) == cx*cx * 6); // x^2 * 6
    BOOST_REQUIRE(z.derivative(0,4) == 0.0); // x^2 * 0
    BOOST_REQUIRE(z.derivative(1,0) == 2*cx * cy*cy*cy); // 2x * y^3
    BOOST_REQUIRE(z.derivative(1,1) == 2*cx * 3*cy*cy); // 2x * 3y^2
    BOOST_REQUIRE(z.derivative(1,2) == 2*cx * 6*cy); // 2x * 6y
    BOOST_REQUIRE(z.derivative(1,3) == 2*cx * 6); // 2x * 6
    BOOST_REQUIRE(z.derivative(1,4) == 0.0); // 2x * 0
    BOOST_REQUIRE(z.derivative(2,0) == 2 * cy*cy*cy); // 2 * y^3
    BOOST_REQUIRE(z.derivative(2,1) == 2 * 3*cy*cy); // 2 * 3y^2
    BOOST_REQUIRE(z.derivative(2,2) == 2 * 6*cy); // 2 * 6y
    BOOST_REQUIRE(z.derivative(2,3) == 2 * 6); // 2 * 6
    BOOST_REQUIRE(z.derivative(2,4) == 0.0); // 2 * 0
    BOOST_REQUIRE(z.derivative(3,0) == 0.0); // 0 * y^3
    BOOST_REQUIRE(z.derivative(3,1) == 0.0); // 0 * 3y^2
    BOOST_REQUIRE(z.derivative(3,2) == 0.0); // 0 * 6y
    BOOST_REQUIRE(z.derivative(3,3) == 0.0); // 0 * 6
    BOOST_REQUIRE(z.derivative(3,4) == 0.0); // 0 * 0
}

BOOST_AUTO_TEST_CASE(dim2_multiplication_and_subtraction)
{
    constexpr int m = 3;
    constexpr int n = 4;
    constexpr double cx = 6.0;
    const auto x = boost::math::differentiation::autodiff::variable<double,m>(cx);
    constexpr double cy = 5.0;
    const auto y = boost::math::differentiation::autodiff::variable<double,0,n>(cy);
    const auto z = x*x - y*y;
    BOOST_REQUIRE(z.derivative(0,0) == cx*cx - cy*cy);
    BOOST_REQUIRE(z.derivative(0,1) == -2*cy);
    BOOST_REQUIRE(z.derivative(0,2) == -2.0);
    BOOST_REQUIRE(z.derivative(0,3) == 0.0);
    BOOST_REQUIRE(z.derivative(0,4) == 0.0);
    BOOST_REQUIRE(z.derivative(1,0) == 2*cx);
    BOOST_REQUIRE(z.derivative(2,0) == 2.0);
    for (int i=1 ; i<=m ; ++i)
        for (int j=1 ; j<=n ; ++j)
            BOOST_REQUIRE(z.derivative(i,j) == 0.0);
}

BOOST_AUTO_TEST_CASE(inverse)
{
    constexpr int m = 3;
    constexpr double cx = 4.0;
    const auto x = boost::math::differentiation::autodiff::variable<double,m>(cx);
    const auto xinv = x.inverse();
    BOOST_REQUIRE(xinv.derivative(0) == 1/cx);
    BOOST_REQUIRE(xinv.derivative(1) == -1/std::pow(cx,2));
    BOOST_REQUIRE(xinv.derivative(2) == 2/std::pow(cx,3));
    BOOST_REQUIRE(xinv.derivative(3) == -6/std::pow(cx,4));
}

BOOST_AUTO_TEST_CASE(division)
{
    constexpr int m = 3;
    constexpr int n = 4;
    constexpr double cx = 16.0;
    auto x = boost::math::differentiation::autodiff::variable<double,m>(cx);
    constexpr double cy = 4.0;
    auto y = boost::math::differentiation::autodiff::variable<double,1,n>(cy);
    auto z = x*x / (y*y);
    BOOST_REQUIRE(z.derivative(0,0) == cx*cx / (cy*cy)); // x^2 * y^-2
    BOOST_REQUIRE(z.derivative(0,1) == cx*cx * (-2)*std::pow(cy,-3));
    BOOST_REQUIRE(z.derivative(0,2) == cx*cx * (6)*std::pow(cy,-4));
    BOOST_REQUIRE(z.derivative(0,3) == cx*cx * (-24)*std::pow(cy,-5));
    BOOST_REQUIRE(z.derivative(0,4) == cx*cx * (120)*std::pow(cy,-6));
    BOOST_REQUIRE(z.derivative(1,0) == 2*cx / (cy*cy));
    BOOST_REQUIRE(z.derivative(1,1) == 2*cx * (-2)*std::pow(cy,-3));
    BOOST_REQUIRE(z.derivative(1,2) == 2*cx * (6)*std::pow(cy,-4));
    BOOST_REQUIRE(z.derivative(1,3) == 2*cx * (-24)*std::pow(cy,-5));
    BOOST_REQUIRE(z.derivative(1,4) == 2*cx * (120)*std::pow(cy,-6));
    BOOST_REQUIRE(z.derivative(2,0) == 2 / (cy*cy));
    BOOST_REQUIRE(z.derivative(2,1) == 2 * (-2)*std::pow(cy,-3));
    BOOST_REQUIRE(z.derivative(2,2) == 2 * (6)*std::pow(cy,-4));
    BOOST_REQUIRE(z.derivative(2,3) == 2 * (-24)*std::pow(cy,-5));
    BOOST_REQUIRE(z.derivative(2,4) == 2 * (120)*std::pow(cy,-6));
    for (int j=0 ; j<=n ; ++j)
        BOOST_REQUIRE(z.derivative(3,j) == 0.0);

    auto x1 = boost::math::differentiation::autodiff::variable<double,m>(cx);
    auto z1 = x1/cy;
    BOOST_REQUIRE(z1.derivative(0) == cx/cy);
    BOOST_REQUIRE(z1.derivative(1) == 1/cy);
    BOOST_REQUIRE(z1.derivative(2) == 0.0);
    BOOST_REQUIRE(z1.derivative(3) == 0.0);
    auto y2 = boost::math::differentiation::autodiff::variable<double,m,n>(cy);
    auto z2 = cx/y2;
    BOOST_REQUIRE(z2.derivative(0,0) == cx/cy);
    BOOST_REQUIRE(z2.derivative(0,1) == -cx/std::pow(cy,2));
    BOOST_REQUIRE(z2.derivative(0,2) == 2*cx/std::pow(cy,3));
    BOOST_REQUIRE(z2.derivative(0,3) == -6*cx/std::pow(cy,4));
    BOOST_REQUIRE(z2.derivative(0,4) == 24*cx/std::pow(cy,5));
    for (int i=1 ; i<=m ; ++i)
        for (int j=0 ; j<=n ; ++j)
            BOOST_REQUIRE(z2.derivative(i,j) == 0.0);

    const auto z3 = y / x;
    BOOST_REQUIRE(z3.derivative(0,0) == cy / cx);
    BOOST_REQUIRE(z3.derivative(0,1) ==  1 / cx);
    BOOST_REQUIRE(z3.derivative(1,0) == -cy / std::pow(cx,2));
    BOOST_REQUIRE(z3.derivative(1,1) ==  -1 / std::pow(cx,2));
    BOOST_REQUIRE(z3.derivative(2,0) == 2*cy / std::pow(cx,3));
    BOOST_REQUIRE(z3.derivative(2,1) ==    2 / std::pow(cx,3));
    BOOST_REQUIRE(z3.derivative(3,0) == -6*cy / std::pow(cx,4));
    BOOST_REQUIRE(z3.derivative(3,1) ==    -6 / std::pow(cx,4));
    for (int i=0 ; i<=m ; ++i)
        for (int j=2 ; j<=n ; ++j)
            BOOST_REQUIRE(z3.derivative(i,j) == 0.0);
}

BOOST_AUTO_TEST_CASE(equality)
{
    constexpr int m = 3;
    constexpr int n = 4;
    constexpr double cx = 10.0;
    constexpr double cy = 10.0;
    const auto x = boost::math::differentiation::autodiff::variable<double,m>(cx);
    const auto y = boost::math::differentiation::autodiff::variable<double,0,n>(cy);
    BOOST_REQUIRE((x == y));
    BOOST_REQUIRE((x == cy));
    BOOST_REQUIRE((cx == y));
    BOOST_REQUIRE((cy == x));
    BOOST_REQUIRE((y == cx));
}

BOOST_AUTO_TEST_CASE(inequality)
{
    constexpr int m = 3;
    constexpr int n = 4;
    constexpr double cx = 10.0;
    constexpr double cy = 11.0;
    const auto x = boost::math::differentiation::autodiff::variable<double,m>(cx);
    const auto y = boost::math::differentiation::autodiff::variable<double,0,n>(cy);
    BOOST_REQUIRE((x != y));
    BOOST_REQUIRE((x != cy));
    BOOST_REQUIRE((cx != y));
    BOOST_REQUIRE((cy != x));
    BOOST_REQUIRE((y != cx));
}

BOOST_AUTO_TEST_CASE(less_than_or_equal_to)
{
    constexpr int m = 3;
    constexpr int n = 4;
    constexpr double cx = 10.0;
    constexpr double cy = 11.0;
    const auto x = boost::math::differentiation::autodiff::variable<double,m>(cx);
    const auto y = boost::math::differentiation::autodiff::variable<double,0,n>(cy);
    BOOST_REQUIRE((x <= y));
    BOOST_REQUIRE((x <= y-1));
    BOOST_REQUIRE((x < y));
    BOOST_REQUIRE((x <= cy));
    BOOST_REQUIRE((x <= cy-1));
    BOOST_REQUIRE((x < cy));
    BOOST_REQUIRE((cx <= y));
    BOOST_REQUIRE((cx <= y-1));
    BOOST_REQUIRE((cx < y));
}

BOOST_AUTO_TEST_CASE(greater_than_or_equal_to)
{
    constexpr int m = 3;
    constexpr int n = 4;
    constexpr double cx = 11.0;
    constexpr double cy = 10.0;
    const auto x = boost::math::differentiation::autodiff::variable<double,m>(cx);
    const auto y = boost::math::differentiation::autodiff::variable<double,0,n>(cy);
    BOOST_REQUIRE((x >= y));
    BOOST_REQUIRE((x >= y+1));
    BOOST_REQUIRE((x > y));
    BOOST_REQUIRE((x >= cy));
    BOOST_REQUIRE((x >= cy+1));
    BOOST_REQUIRE((x > cy));
    BOOST_REQUIRE((cx >= y));
    BOOST_REQUIRE((cx >= y+1));
    BOOST_REQUIRE((cx > y));
}

BOOST_AUTO_TEST_CASE(abs)
{
    constexpr int m = 3;
    constexpr double cx = 11.0;
    const auto x = boost::math::differentiation::autodiff::variable<double,m>(cx);
    auto abs = boost::math::differentiation::autodiff::abs(x);
    BOOST_REQUIRE(abs.derivative(0) == std::abs(cx));
    BOOST_REQUIRE(abs.derivative(1) == 1.0);
    BOOST_REQUIRE(abs.derivative(2) == 0.0);
    BOOST_REQUIRE(abs.derivative(3) == 0.0);
    abs = boost::math::differentiation::autodiff::abs(-x);
    BOOST_REQUIRE(abs.derivative(0) == std::abs(cx));
    BOOST_REQUIRE(abs.derivative(1) == 1.0); // abs(-x) = abs(x)
    BOOST_REQUIRE(abs.derivative(2) == 0.0);
    BOOST_REQUIRE(abs.derivative(3) == 0.0);
    const auto xneg = boost::math::differentiation::autodiff::variable<double,m>(-cx);
    abs = boost::math::differentiation::autodiff::abs(xneg);
    BOOST_REQUIRE(abs.derivative(0) == std::abs(cx));
    BOOST_REQUIRE(abs.derivative(1) == -1.0);
    BOOST_REQUIRE(abs.derivative(2) == 0.0);
    BOOST_REQUIRE(abs.derivative(3) == 0.0);
    const auto zero = boost::math::differentiation::autodiff::variable<double,m>(0);
    abs = boost::math::differentiation::autodiff::abs(zero);
    for (int i=0 ; i<=m ; ++i)
        BOOST_REQUIRE(abs.derivative(i) == 0.0);
}

BOOST_AUTO_TEST_CASE(ceil_and_floor)
{
    constexpr int m = 3;
    double tests[] { -1.5, 0.0, 1.5 };
    for (unsigned t=0 ; t<sizeof(tests)/sizeof(*tests) ; ++t)
    {
        const auto x = boost::math::differentiation::autodiff::variable<double,m>(tests[t]);
        auto ceil = boost::math::differentiation::autodiff::ceil(x);
        auto floor = boost::math::differentiation::autodiff::floor(x);
        BOOST_REQUIRE(ceil.derivative(0) == std::ceil(tests[t]));
        BOOST_REQUIRE(floor.derivative(0) == std::floor(tests[t]));
        for (int i=1 ; i<=m ; ++i)
        {
            BOOST_REQUIRE(ceil.derivative(i) == 0.0);
            BOOST_REQUIRE(floor.derivative(i) == 0.0);
        }
    }
}

BOOST_AUTO_TEST_CASE(one_over_one_plus_x_squared)
{
    constexpr int m = 4;
    constexpr double cx = 1.0;
    auto f = boost::math::differentiation::autodiff::variable<double,m>(cx);
    //f = ((f *= f) += 1).inverse(); // Microsoft Visual C++ version 14.0: fatal error C1001: An internal error has occurred in the compiler. on call to order_sum() in inverse_apply().
    f = 1 / ((f *= f) += 1);
    BOOST_REQUIRE(f.derivative(0) == 0.5);
    BOOST_REQUIRE(f.derivative(1) == -0.5);
    BOOST_REQUIRE(f.derivative(2) == 0.5);
    BOOST_REQUIRE(f.derivative(3) == 0.0);
    BOOST_REQUIRE(f.derivative(4) == -3.0);
}

BOOST_AUTO_TEST_CASE(exp)
{
    constexpr int m = 4;
    constexpr double cx = 2.0;
    const auto x = boost::math::differentiation::autodiff::variable<double,m>(cx);
    auto y = boost::math::differentiation::autodiff::exp(x);
    for (int i=0 ; i<=m ; ++i)
        BOOST_REQUIRE(y.derivative(i) == std::exp(cx));
}

BOOST_AUTO_TEST_CASE(pow)
{
    constexpr double tolerance = 100e-15; // percent
    constexpr int m = 5;
    constexpr int n = 4;
    constexpr double cx = 2.0;
    constexpr double cy = 3.0;
    const auto x = boost::math::differentiation::autodiff::variable<double,m>(cx);
    const auto y = boost::math::differentiation::autodiff::variable<double,m,n>(cy);
    auto z0 = boost::math::differentiation::autodiff::pow(x,cy);
    BOOST_REQUIRE(z0.derivative(0) == std::pow(cx,cy));
    BOOST_REQUIRE(z0.derivative(1) == cy*std::pow(cx,cy-1));
    BOOST_REQUIRE(z0.derivative(2) == cy*(cy-1)*std::pow(cx,cy-2));
    BOOST_REQUIRE(z0.derivative(3) == cy*(cy-1)*(cy-2)*std::pow(cx,cy-3));
    BOOST_REQUIRE(z0.derivative(4) == 0.0);
    BOOST_REQUIRE(z0.derivative(5) == 0.0);
    auto z1 = boost::math::differentiation::autodiff::pow(cx,y);
    BOOST_REQUIRE_CLOSE(z1.derivative(0,0), std::pow(cx,cy), tolerance);
    for (int j=1 ; j<=n ; ++j)
        BOOST_REQUIRE_CLOSE(z1.derivative(0,j), std::pow(std::log(cx),j)*std::exp(cy*std::log(cx)), tolerance);
    for (int i=1 ; i<=m ; ++i)
        for (int j=0 ; j<=n ; ++j)
            BOOST_REQUIRE(z1.derivative(i,j) == 0.0);
    auto z2 = boost::math::differentiation::autodiff::pow(x,y);
    for (int j=0 ; j<=n ; ++j)
        BOOST_REQUIRE_CLOSE(z2.derivative(0,j), std::pow(cx,cy)*std::pow(std::log(cx),j), tolerance);
    for (int j=0 ; j<=n ; ++j)
        BOOST_REQUIRE_CLOSE(z2.derivative(1,j), std::pow(cx,cy-1)*std::pow(std::log(cx),j-1)*(cy*std::log(cx)+j), tolerance);
    BOOST_REQUIRE_CLOSE(z2.derivative(2,0), std::pow(cx,cy-2)*cy*(cy-1), tolerance);
    BOOST_REQUIRE_CLOSE(z2.derivative(2,1), std::pow(cx,cy-2)*(cy*(cy-1)*std::log(cx)+2*cy-1), tolerance);
    for (int j=2 ; j<=n ; ++j)
        BOOST_REQUIRE_CLOSE(z2.derivative(2,j), std::pow(cx,cy-2)*std::pow(std::log(cx),j-2)*(j*(2*cy-1)*std::log(cx)+(j-1)*j+(cy-1)*cy*std::pow(std::log(cx),2)), tolerance);
    BOOST_REQUIRE_CLOSE(z2.derivative(2,4), std::pow(cx,cy-2)*std::pow(std::log(cx),2)*(4*(2*cy-1)*std::log(cx)+(4-1)*4+(cy-1)*cy*std::pow(std::log(cx),2)), tolerance);
}

BOOST_AUTO_TEST_CASE(sqrt)
{
    constexpr int m = 5;
    constexpr double cx = 4.0;
    auto x = boost::math::differentiation::autodiff::variable<double,m>(cx);
    auto y = boost::math::differentiation::autodiff::sqrt(x);
    BOOST_REQUIRE(y.derivative(0) == std::sqrt(cx));
    BOOST_REQUIRE(y.derivative(1) == 0.5*std::pow(cx,-0.5));
    BOOST_REQUIRE(y.derivative(2) == -0.5*0.5*std::pow(cx,-1.5));
    BOOST_REQUIRE(y.derivative(3) == 0.5*0.5*1.5*std::pow(cx,-2.5));
    BOOST_REQUIRE(y.derivative(4) == -0.5*0.5*1.5*2.5*std::pow(cx,-3.5));
    BOOST_REQUIRE(y.derivative(5) == 0.5*0.5*1.5*2.5*3.5*std::pow(cx,-4.5));
    x = boost::math::differentiation::autodiff::variable<double,m>(0);
    y = boost::math::differentiation::autodiff::sqrt(x);
    //std::cout << "boost::math::differentiation::autodiff::sqrt(0) = " << y << std::endl; // (0,inf,-inf,inf,-inf,inf)
    BOOST_REQUIRE(y.derivative(0) == 0.0);
    for (int i=1; i<=m ; ++i)
        BOOST_REQUIRE(y.derivative(i) == (i&1?1:-1)*std::numeric_limits<double>::infinity());
}

BOOST_AUTO_TEST_CASE(log)
{
    constexpr int m = 5;
    constexpr double cx = 2.0;
    auto x = boost::math::differentiation::autodiff::variable<double,m>(cx);
    auto y = boost::math::differentiation::autodiff::log(x);
    BOOST_REQUIRE(y.derivative(0) == std::log(cx));
    BOOST_REQUIRE(y.derivative(1) == 1/cx);
    BOOST_REQUIRE(y.derivative(2) == -1/std::pow(cx,2));
    BOOST_REQUIRE(y.derivative(3) == 2/std::pow(cx,3));
    BOOST_REQUIRE(y.derivative(4) == -6/std::pow(cx,4));
    BOOST_REQUIRE(y.derivative(5) == 24/std::pow(cx,5));
    x = boost::math::differentiation::autodiff::variable<double,m>(0);
    y = boost::math::differentiation::autodiff::log(x);
    //std::cout << "boost::math::differentiation::autodiff::log(0) = " << y << std::endl; // boost::math::differentiation::autodiff::log(0) = depth(1)(-inf,inf,-inf,inf,-inf,inf)
    for (int i=0; i<=m ; ++i)
        BOOST_REQUIRE(y.derivative(i) == (i&1?1:-1)*std::numeric_limits<double>::infinity());
}

BOOST_AUTO_TEST_CASE(ylogx)
{
    constexpr double tolerance = 100e-15; // percent
    constexpr int m = 5;
    constexpr int n = 4;
    constexpr double cx = 2.0;
    constexpr double cy = 3.0;
    const auto x = boost::math::differentiation::autodiff::variable<double,m>(cx);
    const auto y = boost::math::differentiation::autodiff::variable<double,m,n>(cy);
    auto z = y*boost::math::differentiation::autodiff::log(x);
    BOOST_REQUIRE(z.derivative(0,0) == cy*std::log(cx));
    BOOST_REQUIRE(z.derivative(0,1) == std::log(cx));
    BOOST_REQUIRE(z.derivative(0,2) == 0.0);
    BOOST_REQUIRE(z.derivative(0,3) == 0.0);
    BOOST_REQUIRE(z.derivative(0,4) == 0.0);
    for (size_t i=1 ; i<=m ; ++i)
        BOOST_REQUIRE_CLOSE(z.derivative(i,0), std::pow(-1,i-1)*boost::math::factorial<double>(i-1)*cy/std::pow(cx,i), tolerance);
    for (size_t i=1 ; i<=m ; ++i)
        BOOST_REQUIRE_CLOSE(z.derivative(i,1), std::pow(-1,i-1)*boost::math::factorial<double>(i-1)/std::pow(cx,i), tolerance);
    for (size_t i=1 ; i<=m ; ++i)
        for (size_t j=2 ; j<=n ; ++j)
            BOOST_REQUIRE(z.derivative(i,j) == 0.0);
    auto z1 = boost::math::differentiation::autodiff::exp(z);
    BOOST_REQUIRE_CLOSE(z1.derivative(2,4), std::pow(cx,cy-2)*std::pow(std::log(cx),2)*(4*(2*cy-1)*std::log(cx)+(4-1)*4+(cy-1)*cy*std::pow(std::log(cx),2)), tolerance); // RHS is confirmed by https://www.wolframalpha.com/input/?i=D%5Bx%5Ey,%7Bx,2%7D,%7By,4%7D%5D+%2F.+%7Bx-%3E2.0,+y-%3E3.0%7D
}

BOOST_AUTO_TEST_CASE(frexp)
{
    constexpr int m = 3;
    constexpr double cx = 3.5;
    const auto x = boost::math::differentiation::autodiff::variable<double,m>(cx);
    int exp, testexp;
    auto y = boost::math::differentiation::autodiff::frexp(x,&exp);
    BOOST_REQUIRE(y.derivative(0) == std::frexp(cx,&testexp));
    BOOST_REQUIRE(exp == testexp);
    BOOST_REQUIRE(y.derivative(1) == std::exp2(-exp));
    BOOST_REQUIRE(y.derivative(2) == 0.0);
    BOOST_REQUIRE(y.derivative(3) == 0.0);
}

BOOST_AUTO_TEST_CASE(ldexp)
{
    constexpr int m = 3;
    constexpr double cx = 3.5;
    const auto x = boost::math::differentiation::autodiff::variable<double,m>(cx);
    constexpr int exp = 3;
    auto y = boost::math::differentiation::autodiff::ldexp(x,exp);
    BOOST_REQUIRE(y.derivative(0) == std::ldexp(cx,exp));
    BOOST_REQUIRE(y.derivative(1) == std::exp2(exp));
    BOOST_REQUIRE(y.derivative(2) == 0.0);
    BOOST_REQUIRE(y.derivative(3) == 0.0);
}

BOOST_AUTO_TEST_CASE(cos_and_sin)
{
    constexpr double tolerance = 100e-15; // percent
    constexpr int m = 5;
    constexpr double cx = boost::math::constants::third_pi<double>();
    const auto x = boost::math::differentiation::autodiff::variable<double,m>(cx);
    auto cos5 = boost::math::differentiation::autodiff::cos(x);
    BOOST_REQUIRE(cos5.derivative(0) == std::cos(cx));
    BOOST_REQUIRE_CLOSE(cos5.derivative(1), -std::sin(cx), tolerance);
    BOOST_REQUIRE_CLOSE(cos5.derivative(2), -std::cos(cx), tolerance);
    BOOST_REQUIRE_CLOSE(cos5.derivative(3), std::sin(cx), tolerance);
    BOOST_REQUIRE_CLOSE(cos5.derivative(4), std::cos(cx), tolerance);
    BOOST_REQUIRE_CLOSE(cos5.derivative(5), -std::sin(cx), tolerance);
    auto sin5 = boost::math::differentiation::autodiff::sin(x);
    BOOST_REQUIRE(sin5.derivative(0) == std::sin(cx));
    BOOST_REQUIRE_CLOSE(sin5.derivative(1), std::cos(cx), tolerance);
    BOOST_REQUIRE_CLOSE(sin5.derivative(2), -std::sin(cx), tolerance);
    BOOST_REQUIRE_CLOSE(sin5.derivative(3), -std::cos(cx), tolerance);
    BOOST_REQUIRE_CLOSE(sin5.derivative(4), std::sin(cx), tolerance);
    BOOST_REQUIRE_CLOSE(sin5.derivative(5), std::cos(cx), tolerance);
    // Test Order = 0 for codecov
    auto cos0 = cos(boost::math::differentiation::autodiff::variable<double,0>(cx));
    BOOST_REQUIRE(cos0.derivative(0) == std::cos(cx));
    auto sin0 = sin(boost::math::differentiation::autodiff::variable<double,0>(cx));
    BOOST_REQUIRE(sin0.derivative(0) == std::sin(cx));

}

BOOST_AUTO_TEST_CASE(acos)
{
    constexpr double tolerance = 100e-15; // percent
    constexpr int m = 5;
    constexpr double cx = 0.5;
    auto x = boost::math::differentiation::autodiff::variable<double,m>(cx);
    auto y = boost::math::differentiation::autodiff::acos(x);
    BOOST_REQUIRE(y.derivative(0) == std::acos(cx));
    BOOST_REQUIRE_CLOSE(y.derivative(1), -1/std::sqrt(1-cx*cx), tolerance);
    BOOST_REQUIRE_CLOSE(y.derivative(2), -cx/std::pow(1-cx*cx,1.5), tolerance);
    BOOST_REQUIRE_CLOSE(y.derivative(3), -(2*cx*cx+1)/std::pow(1-cx*cx,2.5), tolerance);
    BOOST_REQUIRE_CLOSE(y.derivative(4), -3*cx*(2*cx*cx+3)/std::pow(1-cx*cx,3.5), tolerance);
    BOOST_REQUIRE_CLOSE(y.derivative(5), -(24*(cx*cx+3)*cx*cx+9)/std::pow(1-cx*cx,4.5), tolerance);
}

BOOST_AUTO_TEST_CASE(asin)
{
    constexpr double tolerance = 100e-15; // percent
    constexpr int m = 5;
    constexpr double cx = 0.5;
    auto x = boost::math::differentiation::autodiff::variable<double,m>(cx);
    auto y = boost::math::differentiation::autodiff::asin(x);
    BOOST_REQUIRE(y.derivative(0) == std::asin(cx));
    BOOST_REQUIRE_CLOSE(y.derivative(1), 1/std::sqrt(1-cx*cx), tolerance);
    BOOST_REQUIRE_CLOSE(y.derivative(2), cx/std::pow(1-cx*cx,1.5), tolerance);
    BOOST_REQUIRE_CLOSE(y.derivative(3), (2*cx*cx+1)/std::pow(1-cx*cx,2.5), tolerance);
    BOOST_REQUIRE_CLOSE(y.derivative(4), 3*cx*(2*cx*cx+3)/std::pow(1-cx*cx,3.5), tolerance);
    BOOST_REQUIRE_CLOSE(y.derivative(5), (24*(cx*cx+3)*cx*cx+9)/std::pow(1-cx*cx,4.5), tolerance);
}

BOOST_AUTO_TEST_CASE(asin_infinity)
{
    constexpr int m = 5;
    auto x = boost::math::differentiation::autodiff::variable<double,m>(1);
    auto y = boost::math::differentiation::autodiff::asin(x);
    //std::cout << "boost::math::differentiation::autodiff::asin(1) = " << y << std::endl; // depth(1)(1.5707963267949,inf,inf,-nan,-nan,-nan)
    BOOST_REQUIRE(y.derivative(0) == boost::math::constants::half_pi<double>());
    BOOST_REQUIRE(y.derivative(1) == std::numeric_limits<double>::infinity());
}

BOOST_AUTO_TEST_CASE(asin_derivative)
{
    constexpr double tolerance = 100e-15; // percent
    constexpr int m = 4;
    constexpr double cx = 0.5;
    auto x = boost::math::differentiation::autodiff::variable<double,m>(cx);
    auto y = 1-x*x;
    BOOST_REQUIRE(y.derivative(0) == 1-cx*cx);
    BOOST_REQUIRE(y.derivative(1) == -2*cx);
    BOOST_REQUIRE(y.derivative(2) == -2);
    BOOST_REQUIRE(y.derivative(3) == 0);
    BOOST_REQUIRE(y.derivative(4) == 0);
    y = boost::math::differentiation::autodiff::sqrt(y);
    BOOST_REQUIRE(y.derivative(0) == std::sqrt(1-cx*cx));
    BOOST_REQUIRE_CLOSE(y.derivative(1), -cx/std::sqrt(1-cx*cx), tolerance);
    BOOST_REQUIRE_CLOSE(y.derivative(2), -1/std::pow(1-cx*cx,1.5), tolerance);
    BOOST_REQUIRE_CLOSE(y.derivative(3), -3*cx/std::pow(1-cx*cx,2.5), tolerance);
    BOOST_REQUIRE_CLOSE(y.derivative(4), -(12*cx*cx+3)/std::pow(1-cx*cx,3.5), tolerance);
    y = y.inverse(); // asin'(x) = 1 / sqrt(1-x*x).
    BOOST_REQUIRE_CLOSE(y.derivative(0), 1/std::sqrt(1-cx*cx), tolerance);
    BOOST_REQUIRE_CLOSE(y.derivative(1), cx/std::pow(1-cx*cx,1.5), tolerance);
    BOOST_REQUIRE_CLOSE(y.derivative(2), (2*cx*cx+1)/std::pow(1-cx*cx,2.5), tolerance);
    BOOST_REQUIRE_CLOSE(y.derivative(3), 3*cx*(2*cx*cx+3)/std::pow(1-cx*cx,3.5), tolerance);
    BOOST_REQUIRE_CLOSE(y.derivative(4), (24*(cx*cx+3)*cx*cx+9)/std::pow(1-cx*cx,4.5), tolerance);
}

BOOST_AUTO_TEST_CASE(tan)
{
    constexpr double tolerance = 200e-15; // percent
    constexpr int m = 5;
    constexpr double cx = boost::math::constants::third_pi<double>();
    const auto x = boost::math::differentiation::autodiff::variable<double,m>(cx);
    auto y = boost::math::differentiation::autodiff::tan(x);
    BOOST_REQUIRE_CLOSE(y.derivative(0), std::sqrt(3), tolerance);
    BOOST_REQUIRE_CLOSE(y.derivative(1), 4.0, tolerance);
    BOOST_REQUIRE_CLOSE(y.derivative(2), 8*std::sqrt(3), tolerance);
    BOOST_REQUIRE_CLOSE(y.derivative(3), 80.0, tolerance);
    BOOST_REQUIRE_CLOSE(y.derivative(4), 352*std::sqrt(3), tolerance);
    BOOST_REQUIRE_CLOSE(y.derivative(5), 5824.0, tolerance);
}

BOOST_AUTO_TEST_CASE(atan)
{
    constexpr int m = 5;
    constexpr double cx = 1.0;
    const auto x = boost::math::differentiation::autodiff::variable<double,m>(cx);
    auto y = boost::math::differentiation::autodiff::atan(x);
    BOOST_REQUIRE(y.derivative(0) == boost::math::constants::pi<double>()/4);
    BOOST_REQUIRE(y.derivative(1) == 0.5);
    BOOST_REQUIRE(y.derivative(2) == -0.5);
    BOOST_REQUIRE(y.derivative(3) == 0.5);
    BOOST_REQUIRE(y.derivative(4) == 0.0);
    BOOST_REQUIRE(y.derivative(5) == -3.0);
}

BOOST_AUTO_TEST_CASE(fmod)
{
    constexpr int m = 3;
    constexpr double cx = 3.25;
    constexpr double cy = 0.5;
    auto x = boost::math::differentiation::autodiff::variable<double,m>(cx);
    auto y = boost::math::differentiation::autodiff::fmod(x,cy);
    BOOST_REQUIRE(y.derivative(0) == 0.25);
    BOOST_REQUIRE(y.derivative(1) == 1.0);
    BOOST_REQUIRE(y.derivative(2) == 0.0);
    BOOST_REQUIRE(y.derivative(3) == 0.0);
}

BOOST_AUTO_TEST_CASE(round_and_trunc)
{
    constexpr int m = 3;
    constexpr double cx = 3.25;
    auto x = boost::math::differentiation::autodiff::variable<double,m>(cx);
    auto y = boost::math::differentiation::autodiff::round(x);
    BOOST_REQUIRE(y.derivative(0) == std::round(cx));
    BOOST_REQUIRE(y.derivative(1) == 0.0);
    BOOST_REQUIRE(y.derivative(2) == 0.0);
    BOOST_REQUIRE(y.derivative(3) == 0.0);
    y = boost::math::differentiation::autodiff::trunc(x);
    BOOST_REQUIRE(y.derivative(0) == std::trunc(cx));
    BOOST_REQUIRE(y.derivative(1) == 0.0);
    BOOST_REQUIRE(y.derivative(2) == 0.0);
    BOOST_REQUIRE(y.derivative(3) == 0.0);
}

BOOST_AUTO_TEST_CASE(iround_and_itrunc)
{
    using namespace boost::math;
    constexpr int m = 3;
    constexpr double cx = 3.25;
    auto x = boost::math::differentiation::autodiff::variable<double,m>(cx);
    int y = iround(x);
    BOOST_REQUIRE(y == iround(cx));
    y = itrunc(x);
    BOOST_REQUIRE(y == itrunc(cx));
}

BOOST_AUTO_TEST_CASE(lround_llround_truncl)
{
    constexpr int m = 3;
    constexpr double cx = 3.25;
    auto x = boost::math::differentiation::autodiff::variable<double,m>(cx);
    long yl = boost::math::differentiation::autodiff::lround(x);
    BOOST_REQUIRE(yl == std::lround(cx));
    long long yll = boost::math::differentiation::autodiff::llround(x);
    BOOST_REQUIRE(yll == std::llround(cx));
    long double yld = boost::math::differentiation::autodiff::truncl(x);
    BOOST_REQUIRE(yld == std::truncl(cx));
}

BOOST_AUTO_TEST_CASE(mixed_partials)
{
    constexpr double tolerance = 100e-12; // percent
    // Derivatives calculated from symbolic differentiation by Mathematica for comparison.
    const double answers[] = {19878.406289804349223,20731.748382749395173,14667.607676239390148,1840.5599364498131187,-9219.3180052370721296,-7272.3006340128117838,-2135.2963700622839242,3095.0810272518467995,4249.0267629086156274,2063.9890610627344166,-885.52841148764960841,-1962.1334204417431580,-1846.8998307870845186,-160.95901276032957552,1091.0394123416339941,452.43955743452299467,666.40139227277049900,-415.64641143336291078,-625.14641790399863613,369.94916697726171101,-24330.896138493893431,-18810.416051756267521,-4890.4061227023590999,8833.0050547689764171,8484.3507396816137478,3097.2041512403988935,-3255.0451367834406121,-4342.7785533321930979,-2407.9872379065234860,861.11739164703000843,2436.7437257633086191,-19.246496107338277838,187.78551488705117144,-1259.4660633352121952,-709.68605239721582613,1423.0005586086045369,484.92081333892339591,763.97468850744531805,-327.41629182280555682,-1122.3377072484945211,23973.060071923469893,8840.5431517787968699,-9082.5710332215493783,-12270.273782892587177,-4320.4340714205998547,3281.3519677072808985,5880.3362630834187672,-1288.4827852197065498,-803.97135376265805266,-2986.3872453316983903,-586.73168598226583063,3929.0731892807393562,1453.7282809838266301,1037.8780716859538297,-1482.7458052774013366,-1877.1347929338288106,-931.71387103692982071,254.65655904203226329,1391.2480647456116638,-431.48205631541379551,16975.340053651795550,19662.603563033417098,15765.851307040200043,3972.1550361959370138,-8681.7485397897205125,-7703.1830424603876567,-3049.7086965695187740,2971.4696859922708762,4370.1964998575500257,2524.6324733574356708,-656.60800002366790717,-2423.4529173252581326,-2074.9876642042632042,-381.22537949881329845,1219.5072457919973510,805.38022398408368773,838.40041900589123805,-390.61251971089838316,-828.20854892982357583,293.89998544549947901,-22965.859858439519778,-20026.691015299296217,-7316.0927450633559965,8632.4661339726146593,8987.0468828704522662,4199.9253995361375411,-2958.4298508960628932,-5665.5638912186240622,-2945.4045522503416159,555.65662724782625247,2936.7964035500791392,651.51916507471100081,444.76294274861551486,-1390.9896717990958013,-1142.8614689467638609,1541.9787231173408435,455.71460632938144702,998.79435039403570373,-204.84855819811212954,-1560.3541154604787861,25278.294506052472235,11873.223371790464699,-8242.1873033688781033,-15939.980564174657519,-5648.8335396980314868,2751.5139261227171185,7349.4320024790771292,194.99725459803711274,-402.81568576826882656,-3518.8719086830633712,-1494.3047934746826191,4640.9275094260800875,1585.7577052032271420,1565.1699924044071379,-1513.2598097335400189,-2974.4378726746800928,-1203.2362926538234416,72.524259498791533840,1871.6252742534199495,-2.4899843373796816664,14462.744235186331026,18367.747409164327117,16565.763244996739614,6054.3152526511029520,-8084.9812719820301461,-7988.3143591282012972,-3989.3193469414926985,2616.7211865346490167,4420.8592709704865621,2973.0335197645479091,-324.14530169827137080,-2843.2420399589692219,-2281.4618061432895177,-642.93532295820559249,1299.2872741769553585,1238.5970833720697622,1021.3340427708481651,-329.05293450692710796,-1046.2543015440520751,134.73430395544806552,-21431.416435076611924,-20856.882814790157847,-9829.2619705919309076,7806.8586470778118280,9319.7000856495681801,5319.8987680257582564,-2387.9548264668417364,-6958.2985251653597607,-3468.5391063919725607,130.41672533427094017,3371.1399302351759874,1569.2326780049081053,750.09121011790652458,-1462.2572096265974522,-1661.5778096302406157,1509.6285286038691333,383.89509025808162595,1248.0510963436380133,17.185695642652602749,-2038.0245980026048531,26118.981320178235148,14943.619434822279033,-6650.6862622761310724,-19519.815295474040679,-6983.1902365008486475,1899.2975028736889830,8715.0036526429634882,2368.1506906818643019,136.89207930934828319,-3954.7327061634171420,-2673.5564402311867864,5078.4839352490435947,1643.4591437212048172,2182.2169795063802937,-1345.8388309636205015,-4309.2853506291084135,-1488.0508699224178177,-228.05849430703437209,2373.3989404257091779,773.84813281039280582,12294.403877378555486,16977.349665718583019,17057.174756225031750,8121.1897585118309359,-7458.4435414062843899,-8134.1311608827380587,-4912.8811586137844196,2030.6531360989337179,4407.4905277094127309,3392.4345688258927524,104.03723558415061987,-3180.8176204844632144,-2460.5239870750694373,-938.22093140691334328,1315.2469055718764567,1735.8623924059921882,1209.7596572231669549,-227.33200545666422971,-1266.1262099919292594,-123.07945723381491568,-19806.907943338346855,-21314.816354405752293,-12317.583844301308050,6349.4186598882814744,9489.8196876965277351,6409.5389484563099944,-1550.2817990131252676,-8109.7111997852175121,-3957.8403302968748777,-404.07965558366678588,3693.6143513011819801,2716.1466583227900648,1094.5910866413989005,-1456.2696455499464209,-2244.3806087356369623,1268.5938915562618711,265.22067303277493466,1496.0915787786394884,354.61373510477227819,-2508.4771100486841292,26517.861408751573247,17922.983877419151441,-4328.2591421276680409,-22704.702459400809491,-8268.6137471737389714,740.40560743926114647,9848.9001828360350810,5213.5983414762103377,801.24629237235082333,-4241.8701339207678459,-4092.2413558685505706,5074.4359092060839438,1607.7653292548209160,2861.1556511165675262,-918.93105463172960902,-5803.2113236460920193,-1767.5418979944773144,-663.06462075200757263,2837.9031946139384145,1976.3196007477977178};
    constexpr int Nw=3;
    constexpr int Nx=2;
    constexpr int Ny=4;
    constexpr int Nz=3;
    const boost::math::differentiation::autodiff::variable<double,Nw> w(11);
    const boost::math::differentiation::autodiff::variable<double,0,Nx> x(12);
    const boost::math::differentiation::autodiff::variable<double,0,0,Ny> y(13);
    const boost::math::differentiation::autodiff::variable<double,0,0,0,Nz> z(14);
    const auto v = mixed_partials_f(w,x,y,z); // auto = boost::math::differentiation::autodiff::variable<double,Nw,Nx,Ny,Nz>
    int ia=0;
    for (int iw=0 ; iw<=Nw ; ++iw)
        for (int ix=0 ; ix<=Nx ; ++ix)
            for (int iy=0 ; iy<=Ny ; ++iy)
                for (int iz=0 ; iz<=Nz ; ++iz)
                    BOOST_REQUIRE_CLOSE(v.derivative(iw,ix,iy,iz), answers[ia++], tolerance);
}

BOOST_AUTO_TEST_CASE(multiprecision)
{
    constexpr double tolerance = 100e-98; // percent
    using cpp_dec_float_100 = boost::multiprecision::cpp_dec_float_100;
    // Calculated from Mathematica symbolic differentiation. See example/multiprecision.nb for script.
    const cpp_dec_float_100 answer("1976.31960074779771777988187529041872090812118921875499076582535951111845769110560421820940516423255314");
    constexpr int Nw=3;
    constexpr int Nx=2;
    constexpr int Ny=4;
    constexpr int Nz=3;
    const boost::math::differentiation::autodiff::variable<cpp_dec_float_100,Nw> w(11);
    const boost::math::differentiation::autodiff::variable<cpp_dec_float_100,0,Nx> x(12);
    const boost::math::differentiation::autodiff::variable<cpp_dec_float_100,0,0,Ny> y(13);
    const boost::math::differentiation::autodiff::variable<cpp_dec_float_100,0,0,0,Nz> z(14);
    const auto v = mixed_partials_f(w,x,y,z); // auto = boost::math::differentiation::autodiff::variable<cpp_dec_float_100,Nw,Nx,Ny,Nz>
    // BOOST_REQUIRE_CLOSE(v.derivative(Nw,Nx,Ny,Nz), answer, tolerance); Doesn't compile on travis-ci trusty.
    using std::fabs;
    const cpp_dec_float_100 relative_error = fabs(v.derivative(Nw,Nx,Ny,Nz)/answer-1);
    BOOST_REQUIRE(100*relative_error.convert_to<double>() < tolerance);
}

BOOST_AUTO_TEST_CASE(black_scholes)
{
  constexpr double tolerance = 100e-14; // percent
  using boost::math::differentiation::autodiff::exp;
  using boost::math::differentiation::autodiff::log;
  using boost::math::differentiation::autodiff::sqrt;
  const double K = 100.0; // Strike price
  const boost::math::differentiation::autodiff::variable<double,3> S(105); // Stock price.
  const boost::math::differentiation::autodiff::variable<double,0,3> sigma(5); // Volatility.
  const boost::math::differentiation::autodiff::variable<double,0,0,1> tau(30.0/365); // Time to expiration in years. (30 days).
  const boost::math::differentiation::autodiff::variable<double,0,0,0,1> r(1.25/100); // Interest rate.
  const auto call_price = black_scholes_option_price(call, K, S, sigma, tau, r);
  const auto put_price  = black_scholes_option_price(put,  K, S, sigma, tau, r);
  // Compare automatically calculated greeks by autodiff with formulas for greeks.
  // https://en.wikipedia.org/wiki/Greeks_(finance)#Formulas_for_European_option_Greeks
  const double d1 = static_cast<double>((log(S/K) + (r+sigma*sigma/2)*tau) / (sigma*sqrt(tau)));
  const double d2 = static_cast<double>((log(S/K) + (r-sigma*sigma/2)*tau) / (sigma*sqrt(tau)));
  const double formula_call_delta = +Phi(+d1);
  const double formula_put_delta  = -Phi(-d1);
  const double formula_vega = static_cast<double>(S*phi(d1)*sqrt(tau));
  const double formula_call_theta = static_cast<double>(-S*phi(d1)*sigma/(2*sqrt(tau))-r*K*exp(-r*tau)*Phi(+d2));
  const double formula_put_theta  = static_cast<double>(-S*phi(d1)*sigma/(2*sqrt(tau))+r*K*exp(-r*tau)*Phi(-d2));
  const double formula_call_rho = static_cast<double>(+K*tau*exp(-r*tau)*Phi(+d2));
  const double formula_put_rho  = static_cast<double>(-K*tau*exp(-r*tau)*Phi(-d2));
  const double formula_gamma = static_cast<double>(phi(d1)/(S*sigma*sqrt(tau)));
  const double formula_vanna = static_cast<double>(-phi(d1)*d2/sigma);
  const double formula_charm = static_cast<double>(phi(d1)*(d2*sigma*sqrt(tau)-2*r*tau)/(2*tau*sigma*sqrt(tau)));
  const double formula_vomma = static_cast<double>(S*phi(d1)*sqrt(tau)*d1*d2/sigma);
  const double formula_veta = static_cast<double>(-S*phi(d1)*sqrt(tau)*(r*d1/(sigma*sqrt(tau))-(1+d1*d2)/(2*tau)));
  const double formula_speed = static_cast<double>(-phi(d1)*(d1/(sigma*sqrt(tau))+1)/(S*S*sigma*sqrt(tau)));
  const double formula_zomma = static_cast<double>(phi(d1)*(d1*d2-1)/(S*sigma*sigma*sqrt(tau)));
  const double formula_color =
    static_cast<double>(-phi(d1)/(2*S*tau*sigma*sqrt(tau))*(1+(2*r*tau-d2*sigma*sqrt(tau))*d1/(sigma*sqrt(tau))));
  const double formula_ultima = -formula_vega*static_cast<double>((d1*d2*(1-d1*d2)+d1*d1+d2*d2)/(sigma*sigma));
  BOOST_REQUIRE_CLOSE( call_price.derivative(1,0,0,0), formula_call_delta, tolerance);
  BOOST_REQUIRE_CLOSE( call_price.derivative(0,1,0,0), formula_vega, tolerance);
  BOOST_REQUIRE_CLOSE(-call_price.derivative(0,0,1,0), formula_call_theta, tolerance); // minus sign from tau = T-time
  BOOST_REQUIRE_CLOSE( call_price.derivative(0,0,0,1), formula_call_rho, tolerance);
  BOOST_REQUIRE_CLOSE(  put_price.derivative(1,0,0,0), formula_put_delta, tolerance);
  BOOST_REQUIRE_CLOSE(  put_price.derivative(0,1,0,0), formula_vega, tolerance);
  BOOST_REQUIRE_CLOSE( -put_price.derivative(0,0,1,0), formula_put_theta, tolerance);
  BOOST_REQUIRE_CLOSE(  put_price.derivative(0,0,0,1), formula_put_rho, tolerance);
  BOOST_REQUIRE_CLOSE( call_price.derivative(2,0,0,0), formula_gamma, tolerance);
  BOOST_REQUIRE_CLOSE(  put_price.derivative(2,0,0,0), formula_gamma, tolerance);
  BOOST_REQUIRE_CLOSE( call_price.derivative(1,1,0,0), formula_vanna, tolerance);
  BOOST_REQUIRE_CLOSE(  put_price.derivative(1,1,0,0), formula_vanna, tolerance);
  BOOST_REQUIRE_CLOSE(-call_price.derivative(1,0,1,0), formula_charm, tolerance);
  BOOST_REQUIRE_CLOSE( -put_price.derivative(1,0,1,0), formula_charm, tolerance);
  BOOST_REQUIRE_CLOSE( call_price.derivative(0,2,0,0), formula_vomma, tolerance);
  BOOST_REQUIRE_CLOSE(  put_price.derivative(0,2,0,0), formula_vomma, tolerance);
  BOOST_REQUIRE_CLOSE( call_price.derivative(0,1,1,0), formula_veta, tolerance);
  BOOST_REQUIRE_CLOSE(  put_price.derivative(0,1,1,0), formula_veta, tolerance);
  BOOST_REQUIRE_CLOSE( call_price.derivative(3,0,0,0), formula_speed, tolerance);
  BOOST_REQUIRE_CLOSE(  put_price.derivative(3,0,0,0), formula_speed, tolerance);
  BOOST_REQUIRE_CLOSE( call_price.derivative(2,1,0,0), formula_zomma, tolerance);
  BOOST_REQUIRE_CLOSE(  put_price.derivative(2,1,0,0), formula_zomma, tolerance);
  BOOST_REQUIRE_CLOSE( call_price.derivative(2,0,1,0), formula_color, tolerance);
  BOOST_REQUIRE_CLOSE(  put_price.derivative(2,0,1,0), formula_color, tolerance);
  BOOST_REQUIRE_CLOSE( call_price.derivative(0,3,0,0), formula_ultima, tolerance);
  BOOST_REQUIRE_CLOSE(  put_price.derivative(0,3,0,0), formula_ultima, tolerance);
}

BOOST_AUTO_TEST_SUITE_END()
