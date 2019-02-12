//           Copyright Matthew Pulver 2018 - 2019.
// Distributed under the Boost Software License, Version 1.0.
//      (See accompanying file LICENSE_1_0.txt or copy at
//           https://www.boost.org/LICENSE_1_0.txt)

#include <boost/fusion/include/algorithm.hpp>
#include <boost/fusion/include/vector.hpp>
#include <boost/math/differentiation/autodiff.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/special_functions/fpclassify.hpp> // isnan
#include <boost/math/special_functions/round.hpp> // iround
#include <boost/math/special_functions/trunc.hpp> // itrunc
#include <boost/multiprecision/cpp_bin_float.hpp>

#define BOOST_TEST_MODULE test_autodiff
#include <boost/test/included/unit_test.hpp>

#include <sstream>

//boost::fusion::vector<float,double,long double,boost::multiprecision::cpp_bin_float_50> bin_float_types;
boost::fusion::vector<float,double,long double> bin_float_types; // cpp_bin_float_50 is fixed in boost 1.70
// cpp_dec_float_50 cannot be used with close_at_tolerance
//boost::fusion::vector<boost::multiprecision::cpp_dec_float_50>
boost::fusion::vector<> multiprecision_float_types;

using namespace boost::math::differentiation;

template<typename W,typename X,typename Y,typename Z>
promote<W,X,Y,Z> mixed_partials_f(const W& w, const X& x, const Y& y, const Z& z)
{
    return exp(w*sin(x*log(y)/z) + sqrt(w*z/(x*y))) + w*w/tan(z);
}

// Equations and function/variable names are from
// https://en.wikipedia.org/wiki/Greeks_(finance)#Formulas_for_European_option_Greeks
//
// Standard normal probability density function
template<typename T>
T phi(const T& x)
{
  return boost::math::constants::one_div_root_two_pi<T>()*exp(-0.5*x*x);
}

// Standard normal cumulative distribution function
template<typename T>
T Phi(const T& x)
{
  return 0.5*erfc(-boost::math::constants::one_div_root_two<T>()*x);
}

enum CP { call, put };

// Assume zero annual dividend yield (q=0).
template<typename Price,typename Sigma,typename Tau,typename Rate>
promote<Price,Sigma,Tau,Rate>
    black_scholes_option_price(CP cp, double K, const Price& S, const Sigma& sigma, const Tau& tau, const Rate& r)
{
  const auto d1 = (log(S/K) + (r+sigma*sigma/2)*tau) / (sigma*sqrt(tau));
  const auto d2 = (log(S/K) + (r-sigma*sigma/2)*tau) / (sigma*sqrt(tau));
  static_assert(std::is_same<decltype(S*Phi(d1) - exp(-r*tau)*K*Phi(d2)),
    decltype(exp(-r*tau)*K*Phi(-d2) - S*Phi(-d1))>::value, "decltype(call) != decltype(put)");
  if (cp == call)
    return S*Phi(d1) - exp(-r*tau)*K*Phi(d2);
  else
    return exp(-r*tau)*K*Phi(-d2) - S*Phi(-d1);
}

template<typename T>
T uncast_return(const T& x)
{
    return x == 0 ? 0 : 1;
}

BOOST_AUTO_TEST_SUITE(test_autodiff)

struct constructors_test
{
  template<typename T>
  void operator()(const T&) const
  {
    constexpr int m = 3;
    constexpr int n = 4;
    // Verify value-initialized instance has all 0 entries.
    const autodiff_fvar<T,m> empty1 = autodiff_fvar<T,m>();
    for (int i=0 ; i<=m ; ++i)
        BOOST_REQUIRE(empty1.derivative(i) == 0.0);
    const auto empty2 = autodiff_fvar<T,m,n>();
    for (int i=0 ; i<=m ; ++i)
        for (int j=0 ; j<=n ; ++j)
            BOOST_REQUIRE(empty2.derivative(i,j) == 0.0);
    // Single variable
    constexpr float cx = 10.0;
    const auto x = make_fvar<T,m>(cx);
    for (int i=0 ; i<=m ; ++i)
        if (i==0)
            BOOST_REQUIRE(x.derivative(i) == cx);
        else if (i==1)
            BOOST_REQUIRE(x.derivative(i) == 1.0);
        else
            BOOST_REQUIRE(x.derivative(i) == 0.0);
    const autodiff_fvar<T,n> xn = x;
    for (int i=0 ; i<=n ; ++i)
        if (i==0)
            BOOST_REQUIRE(xn.derivative(i) == cx);
        else if (i==1)
            BOOST_REQUIRE(xn.derivative(i) == 1.0);
        else
            BOOST_REQUIRE(xn.derivative(i) == 0.0);
    // Second independent variable
    constexpr float cy = 100.0;
    const auto y = make_fvar<T,m,n>(cy);
    for (int i=0 ; i<=m ; ++i)
        for (int j=0 ; j<=n ; ++j)
            if (i==0 && j==0)
                BOOST_REQUIRE(y.derivative(i,j) == cy);
            else if (i==0 && j==1)
                BOOST_REQUIRE(y.derivative(i,j) == 1.0);
            else
                BOOST_REQUIRE(y.derivative(i,j) == 0.0);
  }
};

BOOST_AUTO_TEST_CASE(constructors)
{
    boost::fusion::for_each(bin_float_types, constructors_test());
    boost::fusion::for_each(multiprecision_float_types, constructors_test());
}

struct implicit_constructors_test
{
  template<typename T>
  void operator()(const T&) const
  {
    constexpr int m = 3;
    const autodiff_fvar<T,m> x = 3;
    const autodiff_fvar<T,m> one = uncast_return(x);
    const autodiff_fvar<T,m> two_and_a_half = 2.5;
    BOOST_REQUIRE(static_cast<T>(x) == 3.0);
    BOOST_REQUIRE(static_cast<T>(one) == 1.0);
    BOOST_REQUIRE(static_cast<T>(two_and_a_half) == 2.5);
  }
};

BOOST_AUTO_TEST_CASE(implicit_constructors)
{
    boost::fusion::for_each(bin_float_types, implicit_constructors_test());
    boost::fusion::for_each(multiprecision_float_types, implicit_constructors_test());
}

struct assignment_test
{
  template<typename T>
  void operator()(const T&) const
  {
    constexpr int m = 3;
    constexpr int n = 4;
    constexpr float cx = 10.0;
    constexpr float cy = 10.0;
    autodiff_fvar<T,m,n> empty; // Uninitialized variable<> may have non-zero values.
    // Single variable
    auto x = make_fvar<T,m>(cx);
    empty = static_cast<decltype(empty)>(x); // Test static_cast of single-variable to double-variable type.
    for (int i=0 ; i<=m ; ++i)
        for (int j=0 ; j<=n ; ++j)
            if (i==0 && j==0)
                BOOST_REQUIRE(empty.derivative(i,j) == cx);
            else if (i==1 && j==0)
                BOOST_REQUIRE(empty.derivative(i,j) == 1.0);
            else
                BOOST_REQUIRE(empty.derivative(i,j) == 0.0);
    auto y = make_fvar<T,m,n>(cy);
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
};

BOOST_AUTO_TEST_CASE(assignment)
{
    boost::fusion::for_each(bin_float_types, assignment_test());
    boost::fusion::for_each(multiprecision_float_types, assignment_test());
}

struct ostream_test
{
  template<typename T>
  void operator()(const T&) const
  {
    constexpr int m = 3;
    const T cx = 10;
    const auto x = make_fvar<T,m>(cx);
    std::ostringstream ss;
    ss << "x = " << x;
    BOOST_REQUIRE(ss.str() == "x = depth(1)(10,1,0,0)");
  }
};

BOOST_AUTO_TEST_CASE(ostream)
{
    boost::fusion::for_each(bin_float_types, ostream_test());
    boost::fusion::for_each(multiprecision_float_types, ostream_test());
}

struct addition_assignment_test
{
  template<typename T>
  void operator()(const T&) const
  {
    constexpr int m = 3;
    constexpr int n = 4;
    constexpr float cx = 10.0;
    auto sum = autodiff_fvar<T,m,n>(); // zero-initialized
    // Single variable
    const auto x = make_fvar<T,m>(cx);
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
    constexpr float cy = 11.0;
    sum = 0;
    sum += cy;
    for (int i=0 ; i<=m ; ++i)
        for (int j=0 ; j<=n ; ++j)
            if (i==0 && j==0)
                BOOST_REQUIRE(sum.derivative(i,j) == cy);
            else
                BOOST_REQUIRE(sum.derivative(i,j) == 0.0);
  }
};

BOOST_AUTO_TEST_CASE(addition_assignment)
{
    boost::fusion::for_each(bin_float_types, addition_assignment_test());
    boost::fusion::for_each(multiprecision_float_types, addition_assignment_test());
}

struct subtraction_assignment_test
{
  template<typename T>
  void operator()(const T&) const
  {
    constexpr int m = 3;
    constexpr int n = 4;
    constexpr float cx = 10.0;
    auto sum = autodiff_fvar<T,m,n>(); // zero-initialized
    // Single variable
    const auto x = make_fvar<T,m>(cx);
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
    constexpr float cy = 11.0;
    sum = 0;
    sum -= cy;
    for (int i=0 ; i<=m ; ++i)
        for (int j=0 ; j<=n ; ++j)
            if (i==0 && j==0)
                BOOST_REQUIRE(sum.derivative(i,j) == -cy);
            else
                BOOST_REQUIRE(sum.derivative(i,j) == 0.0);
  }
};

BOOST_AUTO_TEST_CASE(subtraction_assignment)
{
    boost::fusion::for_each(bin_float_types, subtraction_assignment_test());
    boost::fusion::for_each(multiprecision_float_types, subtraction_assignment_test());
}

// Try explicit bracing based on feedback. Doesn't add very much except 26 extra lines.
struct multiplication_assignment_test
{
  template<typename T>
  void operator()(const T&) const
  {
    constexpr int m = 3;
    constexpr int n = 4;
    constexpr float cx = 10.0;
    auto product = autodiff_fvar<T,m,n>(1); // unit constant
    // Single variable
    auto x = make_fvar<T,m>(cx);
    product *= x;
    for (int i=0 ; i<=m ; ++i)
    {
        for (int j=0 ; j<=n ; ++j)
        {
            if (i==0 && j==0)
            {
                BOOST_REQUIRE(product.derivative(i,j) == cx);
            }
            else if (i==1 && j==0)
            {
                BOOST_REQUIRE(product.derivative(i,j) == 1.0);
            }
            else
            {
                BOOST_REQUIRE(product.derivative(i,j) == 0.0);
            }
        }
    }
    // Arithmetic constant
    constexpr float cy = 11.0;
    product = 1;
    product *= cy;
    for (int i=0 ; i<=m ; ++i)
    {
        for (int j=0 ; j<=n ; ++j)
        {
            if (i==0 && j==0)
            {
                BOOST_REQUIRE(product.derivative(i,j) == cy);
            }
            else
            {
                BOOST_REQUIRE(product.derivative(i,j) == 0.0);
            }
        }
    }
    // 0 * inf = nan
    x = make_fvar<T,m>(0.0);
    x *= std::numeric_limits<T>::infinity();
    //std::cout << "x = " << x << std::endl;
    for (int i=0 ; i<=m ; ++i)
    {
        if (i==0)
        {
            BOOST_REQUIRE(boost::math::isnan(static_cast<T>(x))); // Correct
            //BOOST_REQUIRE(x.derivative(i) == 0.0); // Wrong. See multiply_assign_by_root_type().
        }
        else if (i==1)
        {
            BOOST_REQUIRE(boost::math::isinf(x.derivative(i)));
        }
        else
        {
            BOOST_REQUIRE(x.derivative(i) == 0.0);
        }
    }
  }
};

BOOST_AUTO_TEST_CASE(multiplication_assignment)
{
    boost::fusion::for_each(bin_float_types, multiplication_assignment_test());
    boost::fusion::for_each(multiprecision_float_types, multiplication_assignment_test());
}

struct division_assignment_test
{
  template<typename T>
  void operator()(const T&) const
  {
    constexpr int m = 3;
    constexpr int n = 4;
    constexpr float cx = 16.0;
    auto quotient = autodiff_fvar<T,m,n>(1); // unit constant
    // Single variable
    const auto x = make_fvar<T,m>(cx);
    quotient /= x;
    BOOST_REQUIRE(quotient.derivative(0,0) == 1/cx);
    BOOST_REQUIRE(quotient.derivative(1,0) == -1/std::pow(cx,2));
    BOOST_REQUIRE(quotient.derivative(2,0) == 2/std::pow(cx,3));
    BOOST_REQUIRE(quotient.derivative(3,0) == -6/std::pow(cx,4));
    for (int i=0 ; i<=m ; ++i)
        for (int j=1 ; j<=n ; ++j)
            BOOST_REQUIRE(quotient.derivative(i,j) == 0.0);
    // Arithmetic constant
    constexpr float cy = 32.0;
    quotient = 1;
    quotient /= cy;
    for (int i=0 ; i<=m ; ++i)
        for (int j=0 ; j<=n ; ++j)
            if (i==0 && j==0)
                BOOST_REQUIRE(quotient.derivative(i,j) == 1/cy);
            else
                BOOST_REQUIRE(quotient.derivative(i,j) == 0.0);
  }
};

BOOST_AUTO_TEST_CASE(division_assignment)
{
    boost::fusion::for_each(bin_float_types, division_assignment_test());
    boost::fusion::for_each(multiprecision_float_types, division_assignment_test());
}

struct unary_signs_test
{
  template<typename T>
  void operator()(const T&) const
  {
    constexpr int m = 3;
    constexpr int n = 4;
    constexpr float cx = 16.0;
    autodiff_fvar<T,m,n> lhs;
    // Single variable
    const auto x = make_fvar<T,m>(cx);
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
};

BOOST_AUTO_TEST_CASE(unary_signs)
{
    boost::fusion::for_each(bin_float_types, unary_signs_test());
    boost::fusion::for_each(multiprecision_float_types, unary_signs_test());
}

// TODO 3 tests for 3 operator+() definitions.

struct cast_double_test
{
  template<typename T>
  void operator()(const T&) const
  {
    constexpr float ca = 13.0;
    constexpr int i = 12;
    constexpr int m = 3;
    const auto x = make_fvar<T,m>(ca);
    BOOST_REQUIRE(i < x);
    BOOST_REQUIRE(i*x == i*ca);
  }
};

BOOST_AUTO_TEST_CASE(cast_double)
{
    boost::fusion::for_each(bin_float_types, cast_double_test());
    boost::fusion::for_each(multiprecision_float_types, cast_double_test());
}

struct int_double_casting_test
{
  template<typename T>
  void operator()(const T&) const
  {
    constexpr float ca = 3.0;
    const auto x0 = make_fvar<T,0>(ca);
    BOOST_REQUIRE(static_cast<T>(x0) == ca);
    const auto x1 = make_fvar<T,1>(ca);
    BOOST_REQUIRE(static_cast<T>(x1) == ca);
    const auto x2 = make_fvar<T,2>(ca);
    BOOST_REQUIRE(static_cast<T>(x2) == ca);
  }
};

BOOST_AUTO_TEST_CASE(int_double_casting)
{
    boost::fusion::for_each(bin_float_types, int_double_casting_test());
    boost::fusion::for_each(multiprecision_float_types, int_double_casting_test());
}

struct scalar_addition_test
{
  template<typename T>
  void operator()(const T&) const
  {
    constexpr float ca = 3.0;
    constexpr float cb = 4.0;
    const auto sum0 = autodiff_fvar<T,0>(ca) + autodiff_fvar<T,0>(cb);
    BOOST_REQUIRE(ca+cb == static_cast<T>(sum0));
    const auto sum1 = autodiff_fvar<T,0>(ca) + cb;
    BOOST_REQUIRE(ca+cb == static_cast<T>(sum1));
    const auto sum2 = ca + autodiff_fvar<T,0>(cb);
    BOOST_REQUIRE(ca+cb == static_cast<T>(sum2));
  }
};

BOOST_AUTO_TEST_CASE(scalar_addition)
{
    boost::fusion::for_each(bin_float_types, scalar_addition_test());
    boost::fusion::for_each(multiprecision_float_types, scalar_addition_test());
}

struct power8_test
{
  template<typename T>
  void operator()(const T&) const
  {
    constexpr int n = 8;
    constexpr float ca = 3.0;
    auto x = make_fvar<T,n>(ca);
    // Test operator*=()
    x *= x;
    x *= x;
    x *= x;
    const T power_factorial = boost::math::factorial<T>(n);
    for (int i=0 ; i<=n ; ++i)
        BOOST_CHECK(static_cast<T>(x.derivative(i)) == static_cast<T>(power_factorial/boost::math::factorial<T>(n-i)*std::pow(ca,n-i)));
    x = make_fvar<T,n>(ca);
    // Test operator*()
    x = x*x*x*x * x*x*x*x;
    for (int i=0 ; i<=n ; ++i)
        BOOST_REQUIRE(x.derivative(i) == power_factorial/boost::math::factorial<T>(n-i)*std::pow(ca,n-i));
  }
};

BOOST_AUTO_TEST_CASE(power8)
{
    boost::fusion::for_each(bin_float_types, power8_test());
    boost::fusion::for_each(multiprecision_float_types, power8_test());
}

struct dim1_multiplication_test
{
  template<typename T>
  void operator()(const T&) const
  {
    constexpr int m = 2;
    constexpr int n = 3;
    constexpr float cy = 4.0;
    auto y0 = make_fvar<T,m>(cy);
    auto y  = make_fvar<T,n>(cy);
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
};

BOOST_AUTO_TEST_CASE(dim1_multiplication)
{
    boost::fusion::for_each(bin_float_types, dim1_multiplication_test());
    boost::fusion::for_each(multiprecision_float_types, dim1_multiplication_test());
}

struct dim1and2_multiplication_test
{
  template<typename T>
  void operator()(const T&) const
  {
    constexpr int m = 2;
    constexpr int n = 3;
    constexpr float cx = 3.0;
    constexpr float cy = 4.0;
    auto x = make_fvar<T,m>(cx);
    auto y = make_fvar<T,m,n>(cy);
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
};

BOOST_AUTO_TEST_CASE(dim1and2_multiplication)
{
    boost::fusion::for_each(bin_float_types, dim1and2_multiplication_test());
    boost::fusion::for_each(multiprecision_float_types, dim1and2_multiplication_test());
}

struct dim2_addition_test
{
  template<typename T>
  void operator()(const T&) const
  {
    constexpr int m = 2;
    constexpr int n = 3;
    constexpr float cx = 3.0;
    const auto x = make_fvar<T,m>(cx);
    BOOST_REQUIRE(x.derivative(0) == cx);
    BOOST_REQUIRE(x.derivative(1) == 1.0);
    BOOST_REQUIRE(x.derivative(2) == 0.0);
    constexpr float cy = 4.0;
    const auto y = make_fvar<T,m,n>(cy);
    BOOST_REQUIRE(static_cast<T>(y.derivative(0)) == cy);
    BOOST_REQUIRE(static_cast<T>(y.derivative(1)) == 0.0); // partial of y w.r.t. x.

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
};

BOOST_AUTO_TEST_CASE(dim2_addition)
{
    boost::fusion::for_each(bin_float_types, dim2_addition_test());
    boost::fusion::for_each(multiprecision_float_types, dim2_addition_test());
}

struct dim2_multiplication_test
{
  template<typename T>
  void operator()(const T&) const
  {
    constexpr int m = 3;
    constexpr int n = 4;
    constexpr float cx = 6.0;
    const auto x = make_fvar<T,m>(cx);
    constexpr float cy = 5.0;
    const auto y = make_fvar<T,0,n>(cy);
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
};

BOOST_AUTO_TEST_CASE(dim2_multiplication)
{
    boost::fusion::for_each(bin_float_types, dim2_multiplication_test());
    boost::fusion::for_each(multiprecision_float_types, dim2_multiplication_test());
}

struct dim2_multiplication_and_subtraction_test
{
  template<typename T>
  void operator()(const T&) const
  {
    constexpr int m = 3;
    constexpr int n = 4;
    constexpr float cx = 6.0;
    const auto x = make_fvar<T,m>(cx);
    constexpr float cy = 5.0;
    const auto y = make_fvar<T,0,n>(cy);
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
};

BOOST_AUTO_TEST_CASE(dim2_multiplication_and_subtraction)
{
    boost::fusion::for_each(bin_float_types, dim2_multiplication_and_subtraction_test());
    boost::fusion::for_each(multiprecision_float_types, dim2_multiplication_and_subtraction_test());
}

struct inverse_test
{
  template<typename T>
  void operator()(const T&) const
  {
    constexpr int m = 3;
    constexpr float cx = 4.0;
    const auto x = make_fvar<T,m>(cx);
    const auto xinv = x.inverse();
    BOOST_REQUIRE(xinv.derivative(0) == 1/cx);
    BOOST_REQUIRE(xinv.derivative(1) == -1/std::pow(cx,2));
    BOOST_REQUIRE(xinv.derivative(2) == 2/std::pow(cx,3));
    BOOST_REQUIRE(xinv.derivative(3) == -6/std::pow(cx,4));
    const auto zero = make_fvar<T,m>(0);
    const auto inf = zero.inverse();
    for (int i=0 ; i<=m ; ++i)
        BOOST_REQUIRE(inf.derivative(i) == (i&1?-1:1)*std::numeric_limits<T>::infinity());
  }
};

BOOST_AUTO_TEST_CASE(inverse)
{
    boost::fusion::for_each(bin_float_types, inverse_test());
    boost::fusion::for_each(multiprecision_float_types, inverse_test());
}

struct division_test
{
  template<typename T>
  void operator()(const T&) const
  {
    constexpr int m = 3;
    constexpr int n = 4;
    constexpr float cx = 16.0;
    auto x = make_fvar<T,m>(cx);
    constexpr float cy = 4.0;
    auto y = make_fvar<T,1,n>(cy);
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

    auto x1 = make_fvar<T,m>(cx);
    auto z1 = x1/cy;
    BOOST_REQUIRE(z1.derivative(0) == cx/cy);
    BOOST_REQUIRE(z1.derivative(1) == 1/cy);
    BOOST_REQUIRE(z1.derivative(2) == 0.0);
    BOOST_REQUIRE(z1.derivative(3) == 0.0);
    auto y2 = make_fvar<T,m,n>(cy);
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
};

BOOST_AUTO_TEST_CASE(division)
{
    boost::fusion::for_each(bin_float_types, division_test());
    boost::fusion::for_each(multiprecision_float_types, division_test());
}

struct equality_test
{
  template<typename T>
  void operator()(const T&) const
  {
    constexpr int m = 3;
    constexpr int n = 4;
    constexpr float cx = 10.0;
    constexpr float cy = 10.0;
    const auto x = make_fvar<T,m>(cx);
    const auto y = make_fvar<T,0,n>(cy);
    BOOST_REQUIRE((x == y));
    BOOST_REQUIRE((x == cy));
    BOOST_REQUIRE((cx == y));
    BOOST_REQUIRE((cy == x));
    BOOST_REQUIRE((y == cx));
  }
};

BOOST_AUTO_TEST_CASE(equality)
{
    boost::fusion::for_each(bin_float_types, equality_test());
    boost::fusion::for_each(multiprecision_float_types, equality_test());
}

struct inequality_test
{
  template<typename T>
  void operator()(const T&) const
  {
    constexpr int m = 3;
    constexpr int n = 4;
    constexpr float cx = 10.0;
    constexpr float cy = 11.0;
    const auto x = make_fvar<T,m>(cx);
    const auto y = make_fvar<T,0,n>(cy);
    BOOST_REQUIRE((x != y));
    BOOST_REQUIRE((x != cy));
    BOOST_REQUIRE((cx != y));
    BOOST_REQUIRE((cy != x));
    BOOST_REQUIRE((y != cx));
  }
};

BOOST_AUTO_TEST_CASE(inequality)
{
    boost::fusion::for_each(bin_float_types, inequality_test());
    boost::fusion::for_each(multiprecision_float_types, inequality_test());
}

struct less_than_or_equal_to_test
{
  template<typename T>
  void operator()(const T&) const
  {
    constexpr int m = 3;
    constexpr int n = 4;
    constexpr float cx = 10.0;
    constexpr float cy = 11.0;
    const auto x = make_fvar<T,m>(cx);
    const auto y = make_fvar<T,0,n>(cy);
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
};

BOOST_AUTO_TEST_CASE(less_than_or_equal_to)
{
    boost::fusion::for_each(bin_float_types, less_than_or_equal_to_test());
    boost::fusion::for_each(multiprecision_float_types, less_than_or_equal_to_test());
}

struct greater_than_or_equal_to_test
{
  template<typename T>
  void operator()(const T&) const
  {
    constexpr int m = 3;
    constexpr int n = 4;
    constexpr float cx = 11.0;
    constexpr float cy = 10.0;
    const auto x = make_fvar<T,m>(cx);
    const auto y = make_fvar<T,0,n>(cy);
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
};

BOOST_AUTO_TEST_CASE(greater_than_or_equal_to)
{
    boost::fusion::for_each(bin_float_types, greater_than_or_equal_to_test());
    boost::fusion::for_each(multiprecision_float_types, greater_than_or_equal_to_test());
}

struct abs_test_test
{
  template<typename T>
  void operator()(const T&) const
  {
    constexpr int m = 3;
    constexpr float cx = 11.0;
    const auto x = make_fvar<T,m>(cx);
    auto a = abs(x);
    BOOST_REQUIRE(a.derivative(0) == std::abs(cx));
    BOOST_REQUIRE(a.derivative(1) == 1.0);
    BOOST_REQUIRE(a.derivative(2) == 0.0);
    BOOST_REQUIRE(a.derivative(3) == 0.0);
    a = abs(-x);
    BOOST_REQUIRE(a.derivative(0) == std::abs(cx));
    BOOST_REQUIRE(a.derivative(1) == 1.0); // abs(-x) = abs(x)
    BOOST_REQUIRE(a.derivative(2) == 0.0);
    BOOST_REQUIRE(a.derivative(3) == 0.0);
    const auto xneg = make_fvar<T,m>(-cx);
    a = abs(xneg);
    BOOST_REQUIRE(a.derivative(0) == std::abs(cx));
    BOOST_REQUIRE(a.derivative(1) == -1.0);
    BOOST_REQUIRE(a.derivative(2) == 0.0);
    BOOST_REQUIRE(a.derivative(3) == 0.0);
    const auto zero = make_fvar<T,m>(0);
    a = abs(zero);
    for (int i=0 ; i<=m ; ++i)
        BOOST_REQUIRE(a.derivative(i) == 0.0);
  }
};

BOOST_AUTO_TEST_CASE(abs_test)
{
    boost::fusion::for_each(bin_float_types, abs_test_test());
    boost::fusion::for_each(multiprecision_float_types, abs_test_test());
}

struct ceil_and_floor_test
{
  template<typename T>
  void operator()(const T&) const
  {
    constexpr int m = 3;
    float tests[] { -1.5, 0.0, 1.5 };
    for (unsigned t=0 ; t<sizeof(tests)/sizeof(*tests) ; ++t)
    {
        const auto x = make_fvar<T,m>(tests[t]);
        auto c = ceil(x);
        auto f = floor(x);
        BOOST_REQUIRE(c.derivative(0) == std::ceil(tests[t]));
        BOOST_REQUIRE(f.derivative(0) == std::floor(tests[t]));
        for (int i=1 ; i<=m ; ++i)
        {
            BOOST_REQUIRE(c.derivative(i) == 0.0);
            BOOST_REQUIRE(f.derivative(i) == 0.0);
        }
    }
  }
};

BOOST_AUTO_TEST_CASE(ceil_and_floor)
{
    boost::fusion::for_each(bin_float_types, ceil_and_floor_test());
    boost::fusion::for_each(multiprecision_float_types, ceil_and_floor_test());
}

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
    const T eps = 100*std::numeric_limits<T>::epsilon(); // percent
    constexpr int m = 5;
    const T cx = boost::math::constants::third_pi<T>();
    const auto x = make_fvar<T,m>(cx);
    auto cos5 = cos(x);
    BOOST_REQUIRE(cos5.derivative(0) == cos(cx));
    BOOST_REQUIRE_CLOSE(cos5.derivative(1), -sin(cx), eps);
    BOOST_REQUIRE_CLOSE(cos5.derivative(2), -cos(cx), eps);
    BOOST_REQUIRE_CLOSE(cos5.derivative(3), sin(cx), eps);
    BOOST_REQUIRE_CLOSE(cos5.derivative(4), cos(cx), eps);
    BOOST_REQUIRE_CLOSE(cos5.derivative(5), -sin(cx), eps);
    auto sin5 = sin(x);
    BOOST_REQUIRE(sin5.derivative(0) == sin(cx));
    BOOST_REQUIRE_CLOSE(sin5.derivative(1), cos(cx), eps);
    BOOST_REQUIRE_CLOSE(sin5.derivative(2), -sin(cx), eps);
    BOOST_REQUIRE_CLOSE(sin5.derivative(3), -cos(cx), eps);
    BOOST_REQUIRE_CLOSE(sin5.derivative(4), sin(cx), eps);
    BOOST_REQUIRE_CLOSE(sin5.derivative(5), cos(cx), eps);
    // Test Order = 0 for codecov
    auto cos0 = cos(make_fvar<T,0>(cx));
    BOOST_REQUIRE(cos0.derivative(0) == cos(cx));
    auto sin0 = sin(make_fvar<T,0>(cx));
    BOOST_REQUIRE(sin0.derivative(0) == sin(cx));
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
    BOOST_REQUIRE(y.derivative(0) == acos(cx));
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
    BOOST_REQUIRE(y.derivative(0) == asin(cx));
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
    BOOST_REQUIRE(y.derivative(0) == asinh(cx));
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
    constexpr int m = 5;
    const T cx = 1;
    auto x = make_fvar<T,m>(cx);
    auto s = sinh(x);
    auto c = cosh(x);
    BOOST_REQUIRE(s.derivative(0) == sinh(static_cast<T>(x)));
    BOOST_REQUIRE(c.derivative(0) == cosh(static_cast<T>(x)));
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

struct lambert_w0_test_test
{
  template<typename T>
  void operator()(const T&) const
  {
    const T eps = 1000*std::numeric_limits<T>::epsilon(); // percent
    constexpr int m = 10;
    const T cx = 3;
    // Mathematica: N[Table[D[ProductLog[x], {x, n}], {n, 0, 10}] /. x -> 3, 52]
    const char* const answers[m+1] {
        "1.049908894964039959988697070552897904589466943706341",
        "0.1707244807388472968312949774415522047470762509741737",
        "-0.04336545501146252734105411312976167858858970875797718",
        "0.02321456264324789334313200360870492961288748451791104",
        "-0.01909049778427783072663170526188353869136655225133878",
        "0.02122935002563637629500975949987796094687564718834156",
        "-0.02979093848448877259041971538394953658978044986784643",
        "0.05051290266216717699803334605370337985567016837482099",
        "-0.1004503154972645060971099914384090562800544486549660",
        "0.2292464437392250211967939182075930820454464472006425",
        "-0.5905839053125614593682763387470654123192290838719517"};
    auto x = make_fvar<T,m>(cx);
    auto y = lambert_w0(x);
    for (int i=0 ; i<=m ; ++i)
    {
        const T answer = boost::lexical_cast<T>(answers[i]);
        BOOST_REQUIRE_CLOSE(y.derivative(i), answer, eps);
    }
    //const T cx0 = -1 / boost::math::constants::e<T>();
    //auto edge = lambert_w0(make_fvar<T,m>(cx0));
    //std::cout << "edge = " << edge << std::endl;
    //edge = depth(1)(-1,inf,-inf,inf,-inf,inf,-inf,inf,-inf,inf,-inf)
    //edge = depth(1)(-1,inf,-inf,inf,-inf,inf,-inf,inf,-inf,inf,-inf)
    //edge = depth(1)(-1,3.68935e+19,-9.23687e+57,4.62519e+96,-2.89497e+135,2.02945e+174,-1.52431e+213,1.19943e+252,-9.75959e+290,8.14489e+329,-6.93329e+368)
  }
};

BOOST_AUTO_TEST_CASE(lambert_w0_test)
{
    boost::fusion::for_each(bin_float_types, lambert_w0_test_test());
    boost::fusion::for_each(multiprecision_float_types, lambert_w0_test_test());
}

struct lround_llround_truncl_test
{
  template<typename T>
  void operator()(const T&) const
  {
    using std::lround;
    using std::llround;
    //using std::truncl; // truncl not supported by boost::multiprecision types.
    constexpr int m = 3;
    const T cx = 3.25;
    auto x = make_fvar<T,m>(cx);
    long yl = lround(x);
    BOOST_REQUIRE(yl == lround(cx));
    long long yll = llround(x);
    BOOST_REQUIRE(yll == llround(cx));
    //long double yld = truncl(x);
    //BOOST_REQUIRE(yld == truncl(cx));
  }
};

BOOST_AUTO_TEST_CASE(lround_llround_truncl)
{
    boost::fusion::for_each(bin_float_types, lround_llround_truncl_test());
    boost::fusion::for_each(multiprecision_float_types, lround_llround_truncl_test());
}

struct mixed_partials_test
{
  template<typename T>
  void operator()(const T&) const
  {
    const T eps = 20000e2*std::numeric_limits<T>::epsilon(); // percent
    // Derivatives calculated from symbolic differentiation by Mathematica for comparison.
    const char* const answers[] = {"19878.40628980434922342465374997798674242532797789489","20731.74838274939517275508122761443159515217855975002","14667.60767623939014840117674691707821648144188283774","1840.559936449813118734351750381849294157477519107602","-9219.318005237072129605008516120710807803827373819700","-7272.300634012811783845589472196110804386170683300081","-2135.296370062283924160196772166043360841114107521292","3095.081027251846799545897828297310835169325417217168","4249.026762908615627428402369471953790564918480025345","2063.989061062734416582172072883742097425754355167541","-885.5284114876496084068555333811894392182458751895290","-1962.133420441743158021558423645064067562765178375508","-1846.899830787084518564013512948598850243350915531775","-160.9590127603295755195950112199107484483554942817846","1091.039412341633994110997652976585409621806446647794","452.4395574345229946707651998323417632800605985181691","666.4013922727704990031159406121675703174518834914461","-415.6464114333629107803309520898363153301435468382605","-625.1464179039986361267627631122900331946746137220517","369.9491669772617110087494756677334192842413470837587","-24330.89613849389343130420303653062335840497802221681","-18810.41605175626752065686192937776868736029049989926","-4890.406122702359099863022925593448420259414896197252","8833.005054768976417065486877649473665597894570245307","8484.350739681613747819854384228795938450532463850094","3097.204151240398893507362023543393154680147349049848","-3255.045136783440612110181337652522080890693968833148","-4342.778553332193097878812792875447018366988006584840","-2407.987237906523486012534085031032446996713414362131","861.1173916470300084261504495377425043024739914571554","2436.743725763308619092960749816106318692933687303014","-19.24649610733827783846392798978023489104363382129689","187.7855148870511714395275130898958731897480766620821","-1259.466063335212195169531010871023748854744563232277","-709.6860523972158261343923419671629587637051060458295","1423.000558608604536932163648918899935569543711292466","484.9208133389233959103861107714757012185008046446372","763.9746885074453180462508029718247316712990115789154","-327.4162918228055568224139277603073169658358026440432","-1122.337707248494521123614369562896901904418640152220","23973.06007192346989337502250398494874845408708506720","8840.543151778796869949670401421984604862699128880003","-9082.571033221549378277312292526023838132689941236879","-12270.27378289258717737657881957466807305650429436397","-4320.434071420599854743576892819691675331049612545664","3281.351967707280898543984556670710235259118405463698","5880.336263083418767219493592767818708317492833223933","-1288.482785219706549809211085113790275109642879331959","-803.9713537626580526627976840414468844364935388365037","-2986.387245331698390346145949708414455858834967096376","-586.7316859822658306283656047992829723003491823675739","3929.073189280739356198769778905960586080418779863615","1453.728280983826630077825553258703050898056317382483","1037.878071685953829685046234106860743366780050925514","-1482.745805277401336553926171580259185140208053329753","-1877.134792933828810602377451370316364621357891989679","-931.7138710369298207131581126980851620513905805624544","254.6565590420322632851077818917210811815919344882311","1391.248064745611663849820246430123214796614030838600","-431.4820563154137955051720207563800896297257103310465","16975.34005365179555009050533000516107937041784876054","19662.60356303341709846238790020024593550984564081068","15765.85130704020004301064240357947656083104783442825","3972.155036195937013764185795634749937308876197976202","-8681.748539789720512499473840242996096730194203989543","-7703.183042460387656743498394861780784700076575106134","-3049.708696569518774040135942468704911634779352213044","2971.469685992270876159892302788930292108129670398058","4370.196499857550025657084783894747734031876677385611","2524.632473357435670756946837415389227139966527203701","-656.6080000236679071742450437463693211275208125750923","-2423.452917325258132591368397957959217829861665178601","-2074.987664204263204162199830716851483704870169031179","-381.2253794988132984501358802316138392247470857452486","1219.507245791997351017860252538035146744682380716428","805.3802239840836877339667281819652171888443003165988","838.4004190058912380470543219448821914235443115661655","-390.6125197108983831575656956558201636111305409512701","-828.2085489298235758253219930356006757081473789845849","293.8999854454994790079171865082094494146506490533363","-22965.85985843951977785883587223006628792405076928067","-20026.69101529929621743747554537576887048069629325374","-7316.092745063355996548975300169565482331369744607021","8632.466133972614659252310985982644793465043032940318","8987.046882870452266200748127338744248816756004290490","4199.925399536137541108783465785304128965582292174062","-2958.429850896062893179851696175634522187021390095560","-5665.563891218624062243686482808197054863235184904433","-2945.404552250341615883104643651287431663294281737652","555.6566272478262524735403145861484390537770707372992","2936.796403550079139218970638242013974322758744804216","651.5191650747110008135060635556227666232180743487328","444.7629427486155148584918602702161457622049333694568","-1390.989671799095801316658971275073184600067187023729","-1142.861468946763860859271224968631944511098747155437","1541.978723117340843491920690654997335632919116206279","455.7146063293814470171599782651235242129856311098151","998.7943503940357037260061331795191352937661538946216","-204.8485581981121295383497187536442450324011940647949","-1560.354115460478786113711476250386112014306509906244","25278.29450605247223516529112562423587288781657290275","11873.22337179046469888005044109378787446671408425048","-8242.187303368878103323785658604027555126374435611949","-15939.98056417465751946455567789306872745912255628512","-5648.833539698031486810309720694416837861242341227280","2751.513926122717118525029734574022921057261239749143","7349.432002479077129245930487320138527887196396579062","194.9972545980371127390142753318206783334452047502143","-402.8156857682688265622049800462325595907987257153782","-3518.871908683063371167722463713374376552181380727802","-1494.304793474682619087166400375396721307777439607909","4640.927509426080087451995953783429589632369803588940","1585.757705203227141964561144798400703219894640413562","1565.169992404407137888592924342582799362959736185298","-1513.259809733540018859089666188672238777297615451800","-2974.437872674680092826212901753475972242208819679978","-1203.236292653823441598437153564865951527142648802876","72.52425949879153384040698301599842998884036742649047","1871.625274253419949517250818647194858608124560073483","-2.489984337379681666361341362948045621969765070197429","14462.74423518633102580192225823524237502860825596609","18367.74740916432711689913219912502810575714860430297","16565.76324499673961400925630526921000337443450249297","6054.315252651102952034254100792777051580892954459740","-8084.981271982030146065497115893934803061545998433631","-7988.314359128201297240919364015959817416101519999194","-3989.319346941492698525859335371231602272119870228687","2616.721186534649016680934493970036169897788778926434","4420.859270970486562095630193355634655337290952862363","2973.033519764547909146474824627687039969488363657908","-324.1453016982713707989332262410969595194473127209825","-2843.242039958969221918101261762794653424879358390111","-2281.461806143289517702658392470195144560150025832652","-642.9353229582055924928927665183236308235598082837497","1299.287274176955358490409470855361289523321919337117","1238.597083372069762230817383681570828675426312803376","1021.334042770848165110529668635291528449691525937968","-329.0529345069271079573348500899329811170455711610811","-1046.254301544052075124857362060924818517694048905299","134.7343039554480655186788228552325941588620079791654","-21431.41643507661192392650726158493697457993678274754","-20856.88281479015784660571401663659059349708627445067","-9829.261970591930907585958999196966814861251125275804","7806.858647077811827981774785577363365546600234846335","9319.700085649568180114405924685286453652118439999060","5319.898768025758256383579171601100187435481641933401","-2387.954826466841736373447020403170264502066930376059","-6958.298525165359760665355886221309296550746152109847","-3468.539106391972560670887295398968213297736424267559","130.4167253342709401698825285623058661085645012029873","3371.139930235175987370940343096776588915600470241960","1569.232678004908105313880673484968847566948896728142","750.0912101179065245750415609380442359608197763310413","-1462.257209626597452197736652121394535208578921869658","-1661.577809630240615684355192771059515041884351493459","1509.628528603869133250456671040505284128185908768108","383.8950902580816259502239917715884779698864996879279","1248.051096343638013308778159911906703363730187986273","17.18569564265260274901760034571610990094333217519021","-2038.024598002604853054532645991188063394308018947374","26118.98132017823514803387529120810044029492871875474","14943.61943482227903328457116850255971625430735856355","-6650.686262276131072415580833374348889422387492668440","-19519.81529547404067945704333355155941895199228108631","-6983.190236500848647457042860591724089812405118922223","1899.297502873688983038424995203515277346497811783168","8715.003652642963488202943622358986745434720576722170","2368.150690681864301926962120618658083737878227231428","136.8920793093482831910443246272238406481527839521448","-3954.732706163417141961077488373290331419627965482785","-2673.556440231186786375595871506657802723673830409989","5078.483935249043594670125721926702845818403229980691","1643.459143721204817182772630730123271413273760820347","2182.216979506380293664703833586468523416961563720645","-1345.838830963620501537777318021157952722412472356094","-4309.285350629108413525304135326225818270616857298235","-1488.050869922417817689426519211523527088509094291312","-228.0584943070343720919835603886532454450555855354340","2373.398940425709177876367020236623713151456855728138","773.8481328103928058186643458500631723389600248582833","12294.40387737855548614823173849184004455244840062464","16977.34966571858301862913845572077593071467784570724","17057.17475622503175013658695220988017704387344177727","8121.189758511830935868344768490586007624092305459885","-7458.443541406284389918808653948439156033975014107187","-8134.131160882738058651976911725365291142418949378248","-4912.881158613784419581465435995807691111897279859302","2030.653136098933717888434825960516061206391833398177","4407.490527709412730881592594976776779312299897714205","3392.434568825892752350943548729559313328141534290860","104.0372355841506198680609232049783930050635078746762","-3180.817620484463214391157460812371170723810181051096","-2460.523987075069437321629265332968914260047631079537","-938.2209314069133432825590545267820890922150850657831","1315.246905571876456706320919211807375254975062430487","1735.862392405992188189147617586418269768276241147998","1209.759657223166954850207025399731503326968841680649","-227.3320054566642297128407910803774238020746116287390","-1266.126209991929259396966729664100401813091860201682","-123.0794572338149156803989321165094334755661021559442","-19806.90794333834685506732819834090525250045748665845","-21314.81635440575229337844631555492486744407550254908","-12317.58384430130805020250005527399703840208659666608","6349.418659888281474363154227419204673663621492760982","9489.819687696527735093973063679592839666155440941289","6409.538948456309994399374417972222747225748405617373","-1550.281799013125267606263057621300789555474258987989","-8109.711199785217512061886243157800006692908759687186","-3957.840330296874877742767473517819198882831790006004","-404.0796555836667858753163727999380679499192203780272","3693.614351301181980145006883746936633676934626580499","2716.146658322790064799415509615557123789406209068981","1094.591086641398900496318896947912437274250932576747","-1456.269645549946420883827817869876763706452982413420","-2244.380608735636962338392373719455877272151458411079","1268.593891556261871090883000459505759446497182073132","265.2206730327749346649809229271069944357537135668622","1496.091578778639488439197917198148587432113387871024","354.6137351047722781932932090799444060236757625488818","-2508.477110048684129181005769771219369377836598443263","26517.86140875157324686379805134248778305979287686214","17922.98387741915144079932445041215068937644694653527","-4328.259142127668040873054918170572859673703425721293","-22704.70245940080949074466622805971940616027152354999","-8268.613747173738971390434576274225941735552759965376","740.4056074392611464740778308961471299437619012164253","9848.900182836035080973766381422758538530595451048714","5213.598341476210337710365441072904970861063876340963","801.2462923723508233330997243930793458484750729415321","-4241.870133920767845856621968904769727964770527614244","-4092.241355868550570635569815488217469506874233892269","5074.435909206083943809967780457349942315503368249477","1607.765329254820915989772546102530187884674235100928","2861.155651116567526208762405651011317435252198548496","-918.9310546317296090214320737728927500362088478158839","-5803.211323646092019259074499814222806376618363553826","-1767.541897994477314401145980308432268207111761980100","-663.0646207520075726320417301262932382663072876188661","2837.903194613938414496183429129769829434890424213252","1976.319600747797717779881875290418720908121189218755"};
    constexpr int Nw=3;
    constexpr int Nx=2;
    constexpr int Ny=4;
    constexpr int Nz=3;
    const auto w = make_fvar<T,Nw>(11);
    const auto x = make_fvar<T,0,Nx>(12);
    const auto y = make_fvar<T,0,0,Ny>(13);
    const auto z = make_fvar<T,0,0,0,Nz>(14);
    const auto v = mixed_partials_f(w,x,y,z); // auto = autodiff_fvar<double,Nw,Nx,Ny,Nz>
    int ia=0;
    for (int iw=0 ; iw<=Nw ; ++iw)
        for (int ix=0 ; ix<=Nx ; ++ix)
            for (int iy=0 ; iy<=Ny ; ++iy)
                for (int iz=0 ; iz<=Nz ; ++iz)
                {
                    const T answer = boost::lexical_cast<T>(answers[ia++]);
                    BOOST_REQUIRE_CLOSE(v.derivative(iw,ix,iy,iz), answer, eps);
                }
  }
};

BOOST_AUTO_TEST_CASE(mixed_partials)
{
    boost::fusion::for_each(bin_float_types, mixed_partials_test());
}

struct multiprecision_test
{
  template<typename T>
  void operator()(const T&) const
  {
    const T eps = 600*std::numeric_limits<T>::epsilon(); // percent
    constexpr int Nw=3;
    constexpr int Nx=2;
    constexpr int Ny=4;
    constexpr int Nz=3;
    const auto w = make_fvar<T,Nw>(11);
    const auto x = make_fvar<T,0,Nx>(12);
    const auto y = make_fvar<T,0,0,Ny>(13);
    const auto z = make_fvar<T,0,0,0,Nz>(14);
    const auto v = mixed_partials_f(w,x,y,z); // auto = autodiff_fvar<T,Nw,Nx,Ny,Nz>
    // Calculated from Mathematica symbolic differentiation.
    const T answer = boost::lexical_cast<T>("1976.3196007477977177798818752904187209081211892187"
        "5499076582535951111845769110560421820940516423255314");
    // BOOST_REQUIRE_CLOSE(v.derivative(Nw,Nx,Ny,Nz), answer, eps); // Doesn't work for cpp_dec_float
    using std::fabs;
    const double relative_error = static_cast<double>(fabs(v.derivative(Nw,Nx,Ny,Nz)/answer-1));
    BOOST_REQUIRE(100*relative_error < eps);
  }
};

BOOST_AUTO_TEST_CASE(multiprecision)
{
    //multiprecision_test()(boost::multiprecision::cpp_bin_float_50());
    boost::fusion::for_each(bin_float_types, mixed_partials_test());
}

struct black_scholes_test
{
  template<typename T>
  void operator()(const T&) const
  {
    //const T eps = 2725*std::numeric_limits<T>::epsilon(); // percent
    const T eps = 2600e2*std::numeric_limits<T>::epsilon(); // percent - requied by OSX
    const double K = 100.0; // Strike price
    const auto S     = make_fvar<T,3>(105); // Stock price.
    const auto sigma = make_fvar<T,0,3>(5); // Volatility.
    const auto tau   = make_fvar<T,0,0,1>(T(30.0)/365); // Time to expiration in years. (30 days).
    const auto r     = make_fvar<T,0,0,0,1>(T(1.25)/100); // Interest rate.
    const auto call_price = black_scholes_option_price(call, K, S, sigma, tau, r);
    const auto put_price  = black_scholes_option_price(put,  K, S, sigma, tau, r);
    // Compare automatically calculated greeks by autodiff with formulas for greeks.
    // https://en.wikipedia.org/wiki/Greeks_(finance)#Formulas_for_European_option_Greeks
    const T d1 = static_cast<T>((log(S/K) + (r+sigma*sigma/2)*tau) / (sigma*sqrt(tau)));
    const T d2 = static_cast<T>((log(S/K) + (r-sigma*sigma/2)*tau) / (sigma*sqrt(tau)));
    const T Phi_pd1 = Phi(d1);
    // intermediate cpp_dec_float calculation can't go to template function as it can't be implicitly cast back to T.
    const T Phi_nd1 = Phi(static_cast<T>(-d1));
    const T Phi_pd2 = Phi(d2);
    const T Phi_nd2 = Phi(static_cast<T>(-d2));
    //const T formula_call_delta = +Phi(+d1);
    const T formula_call_delta = +Phi_pd1;
    //const T formula_put_delta  = -Phi(-d1);
    const T formula_put_delta  = -Phi_nd1;
    const T formula_vega = static_cast<T>(S*phi(d1)*sqrt(tau));
    //const T formula_call_theta = static_cast<T>(-S*phi(d1)*sigma/(2*sqrt(tau))-r*K*exp(-r*tau)*Phi(+d2));
    const T formula_call_theta = static_cast<T>(-S*phi(d1)*sigma/(2*sqrt(tau))-r*K*exp(-r*tau)*Phi_pd2);
    //const T formula_put_theta  = static_cast<T>(-S*phi(d1)*sigma/(2*sqrt(tau))+r*K*exp(-r*tau)*Phi(-d2));
    const T formula_put_theta  = static_cast<T>(-S*phi(d1)*sigma/(2*sqrt(tau))+r*K*exp(-r*tau)*Phi_nd2);
    //const T formula_call_rho = static_cast<T>(+K*tau*exp(-r*tau)*Phi(+d2));
    const T formula_call_rho = static_cast<T>(+K*tau*exp(-r*tau)*Phi_pd2);
    //const T formula_put_rho  = static_cast<T>(-K*tau*exp(-r*tau)*Phi(-d2));
    const T formula_put_rho  = static_cast<T>(-K*tau*exp(-r*tau)*Phi_nd2);
    const T formula_gamma = static_cast<T>(phi(d1)/(S*sigma*sqrt(tau)));
    const T formula_vanna = static_cast<T>(-phi(d1)*d2/sigma);
    const T formula_charm = static_cast<T>(phi(d1)*(d2*sigma*sqrt(tau)-2*r*tau)/(2*tau*sigma*sqrt(tau)));
    const T formula_vomma = static_cast<T>(S*phi(d1)*sqrt(tau)*d1*d2/sigma);
    const T formula_veta = static_cast<T>(-S*phi(d1)*sqrt(tau)*(r*d1/(sigma*sqrt(tau))-(1+d1*d2)/(2*tau)));
    const T formula_speed = static_cast<T>(-phi(d1)*(d1/(sigma*sqrt(tau))+1)/(S*S*sigma*sqrt(tau)));
    const T formula_zomma = static_cast<T>(phi(d1)*(d1*d2-1)/(S*sigma*sigma*sqrt(tau)));
    const T formula_color =
      static_cast<T>(-phi(d1)/(2*S*tau*sigma*sqrt(tau))*(1+(2*r*tau-d2*sigma*sqrt(tau))*d1/(sigma*sqrt(tau))));
    const T formula_ultima = -formula_vega*static_cast<T>((d1*d2*(1-d1*d2)+d1*d1+d2*d2)/(sigma*sigma));
    BOOST_REQUIRE_CLOSE( call_price.derivative(1,0,0,0), formula_call_delta, eps);
    BOOST_REQUIRE_CLOSE( call_price.derivative(0,1,0,0), formula_vega, eps);
    BOOST_REQUIRE_CLOSE(-call_price.derivative(0,0,1,0), formula_call_theta, eps); // minus sign from tau = T-time
    BOOST_REQUIRE_CLOSE( call_price.derivative(0,0,0,1), formula_call_rho, eps);
    BOOST_REQUIRE_CLOSE(  put_price.derivative(1,0,0,0), formula_put_delta, eps);
    BOOST_REQUIRE_CLOSE(  put_price.derivative(0,1,0,0), formula_vega, eps);
    BOOST_REQUIRE_CLOSE( -put_price.derivative(0,0,1,0), formula_put_theta, eps);
    BOOST_REQUIRE_CLOSE(  put_price.derivative(0,0,0,1), formula_put_rho, eps);
    BOOST_REQUIRE_CLOSE( call_price.derivative(2,0,0,0), formula_gamma, eps);
    BOOST_REQUIRE_CLOSE(  put_price.derivative(2,0,0,0), formula_gamma, eps);
    BOOST_REQUIRE_CLOSE( call_price.derivative(1,1,0,0), formula_vanna, eps);
    BOOST_REQUIRE_CLOSE(  put_price.derivative(1,1,0,0), formula_vanna, eps);
    BOOST_REQUIRE_CLOSE(-call_price.derivative(1,0,1,0), formula_charm, eps);
    BOOST_REQUIRE_CLOSE( -put_price.derivative(1,0,1,0), formula_charm, eps);
    BOOST_REQUIRE_CLOSE( call_price.derivative(0,2,0,0), formula_vomma, eps);
    BOOST_REQUIRE_CLOSE(  put_price.derivative(0,2,0,0), formula_vomma, eps);
    BOOST_REQUIRE_CLOSE( call_price.derivative(0,1,1,0), formula_veta, eps);
    BOOST_REQUIRE_CLOSE(  put_price.derivative(0,1,1,0), formula_veta, eps);
    BOOST_REQUIRE_CLOSE( call_price.derivative(3,0,0,0), formula_speed, eps);
    BOOST_REQUIRE_CLOSE(  put_price.derivative(3,0,0,0), formula_speed, eps);
    BOOST_REQUIRE_CLOSE( call_price.derivative(2,1,0,0), formula_zomma, eps);
    BOOST_REQUIRE_CLOSE(  put_price.derivative(2,1,0,0), formula_zomma, eps);
    BOOST_REQUIRE_CLOSE( call_price.derivative(2,0,1,0), formula_color, eps);
    BOOST_REQUIRE_CLOSE(  put_price.derivative(2,0,1,0), formula_color, eps);
    BOOST_REQUIRE_CLOSE( call_price.derivative(0,3,0,0), formula_ultima, eps);
    BOOST_REQUIRE_CLOSE(  put_price.derivative(0,3,0,0), formula_ultima, eps);
  }
};

BOOST_AUTO_TEST_CASE(black_scholes)
{
    boost::fusion::for_each(bin_float_types, black_scholes_test());
}

/*
// Compilation tests for boost special functions.
struct boost_special_functions_test
{
  template<typename T>
  void operator()(const T&) const
  {
    using namespace boost;
    constexpr int m = 3;
    BOOST_REQUIRE(math::acosh(make_fvar<T,m>(1.5)) == math::acosh(static_cast<T>(1.5)));
    // Policy parameter prevents proper ADL for autodiff_fvar objects. E.g. iround(v,pol) instead of iround(v).
    // In cyl_bessel_j_imp() call is made to iround(v, pol) with v of type autodiff_fvar. It it were just iround(v)
    // then autodiff's iround would properly be called via ADL.
    //BOOST_REQUIRE(math::airy_ai(make_fvar<T,m>(1)) == math::airy_ai(static_cast<T>(1)));
    //BOOST_REQUIRE(math::airy_bi(make_fvar<T,m>(1)) == math::airy_bi(static_cast<T>(1)));
    //BOOST_REQUIRE(math::airy_ai_prime(make_fvar<T,m>(1)) == math::airy_ai_prime(static_cast<T>(1)));
    //BOOST_REQUIRE(math::airy_bi_prime(make_fvar<T,m>(1)) == math::airy_bi_prime(static_cast<T>(1)));
    BOOST_REQUIRE(math::asinh(make_fvar<T,m>(0.5)) == math::asinh(static_cast<T>(0.5)));
    BOOST_REQUIRE(math::atanh(make_fvar<T,m>(0.5)) == math::atanh(static_cast<T>(0.5)));
    // Policy parameter prevents ADL.
    //BOOST_REQUIRE(math::cyl_bessel_j(0,make_fvar<T,m>(0.5)) == math::cyl_bessel_j(0,static_cast<T>(0.5)));
    //BOOST_REQUIRE(math::cyl_neumann(0,make_fvar<T,m>(0.5)) == math::cyl_neumann(0,static_cast<T>(0.5)));
    //BOOST_REQUIRE(math::cyl_bessel_j_zero(make_fvar<T,m>(0.5),0) == math::cyl_bessel_j_zero(static_cast<T>(0.5),0));
    //BOOST_REQUIRE(math::cyl_neumann_zero(make_fvar<T,m>(0.5),0) == math::cyl_neumann_zero(static_cast<T>(0.5),0));
    // Required sinh() (added) but then has policy parameter ADL issue.
    //BOOST_REQUIRE(math::cyl_bessel_i(0,make_fvar<T,m>(0.5)) == math::cyl_bessel_i(0,static_cast<T>(0.5)));
    BOOST_REQUIRE(math::cyl_bessel_k(0,make_fvar<T,m>(0.5)) == math::cyl_bessel_k(0,static_cast<T>(0.5)));
    // Policy parameter prevents ADL.
    //BOOST_REQUIRE(math::sph_bessel(0,make_fvar<T,m>(0.5)) == math::sph_bessel(0,static_cast<T>(0.5)));
    // Required fmod() but then has policy parameter ADL issue.
    //BOOST_REQUIRE(math::sph_neumann(0,make_fvar<T,m>(0.5)) == math::sph_neumann(0,static_cast<T>(0.5)));
    // Policy parameter prevents ADL.
    //BOOST_REQUIRE(math::cyl_bessel_j_prime(0,make_fvar<T,m>(0.5)) == math::cyl_bessel_j_prime(0,static_cast<T>(0.5)));
    //BOOST_REQUIRE(math::cyl_neumann_prime(0,make_fvar<T,m>(0.5)) == math::cyl_neumann_prime(0,static_cast<T>(0.5)));
    //BOOST_REQUIRE(math::cyl_bessel_i_prime(0,make_fvar<T,m>(0.5)) == math::cyl_bessel_i_prime(0,static_cast<T>(0.5)));
    BOOST_REQUIRE(math::cyl_bessel_k_prime(0,make_fvar<T,m>(0.5)) == math::cyl_bessel_k_prime(0,static_cast<T>(0.5)));
    // Policy parameter prevents ADL.
    //BOOST_REQUIRE(math::sph_bessel_prime(0,make_fvar<T,m>(0.5)) == math::sph_bessel_prime(0,static_cast<T>(0.5)));
    //BOOST_REQUIRE(math::sph_neumann_prime(0,make_fvar<T,m>(0.5)) == math::sph_neumann_prime(0,static_cast<T>(0.5)));
    // Per docs: "the functions can only be instantiated on types float, double and long double."
    //BOOST_REQUIRE(math::cyl_hankel_1(0,make_fvar<T,m>(0.5)).real() == math::cyl_hankel_1(0,static_cast<T>(0.5)).real());
    //BOOST_REQUIRE(math::cyl_hankel_2(0,make_fvar<T,m>(0.5)).real() == math::cyl_hankel_2(0,static_cast<T>(0.5)).real());
    //BOOST_REQUIRE(math::sph_hankel_1(0,make_fvar<T,m>(0.5)).real() == math::sph_hankel_1(0,static_cast<T>(0.5)).real());
    //BOOST_REQUIRE(math::sph_hankel_2(0,make_fvar<T,m>(0.5)).real() == math::sph_hankel_2(0,static_cast<T>(0.5)).real());
    // Compiles, but compares 0.74868571757768376251 == 0.74868571757768354047 which is false.
    // BOOST_REQUIRE(static_cast<T>(math::beta(make_fvar<T,m>(1.1),make_fvar<T,m>(1.2))) == static_cast<T>(math::beta(static_cast<T>(1.1),static_cast<T>(1.2))));
    // Skipped other beta functions.
    std::cout.precision(20);
    // Compiles, but compares 0.7937005259840996807 == 0.79370052598409979172 which is false.
    //BOOST_REQUIRE(math::cbrt(make_fvar<T,m>(0.5)) == math::cbrt(static_cast<T>(0.5)));
    //Skipping other Basic Functions
    BOOST_REQUIRE(math::chebyshev_next(make_fvar<T,m>(0.5),make_fvar<T,m>(0.5),make_fvar<T,m>(0.5)) ==
        math::chebyshev_next(static_cast<T>(0.5),static_cast<T>(0.5),static_cast<T>(0.5)));
    // Requires acosh() (added)
    BOOST_REQUIRE(math::chebyshev_t(0,make_fvar<T,m>(0.5)) == math::chebyshev_t(0,static_cast<T>(0.5)));
    BOOST_REQUIRE(math::chebyshev_u(0,make_fvar<T,m>(0.5)) == math::chebyshev_u(0,static_cast<T>(0.5)));
    BOOST_REQUIRE(math::chebyshev_t_prime(0,make_fvar<T,m>(0.5)) == math::chebyshev_t_prime(0,static_cast<T>(0.5)));
    {
        std::array<double,4> c{14.2, -13.7, 82.3, 96};
        // /usr/include/boost/math/special_functions/chebyshev.hpp:164:40: error: cannot convert ‘boost::math::differentiation::autodiff_v1::detail::fvar<double, 3>’ to ‘double’ in return
        //BOOST_REQUIRE(math::chebyshev_clenshaw_recurrence(c.data(),c.size(),make_fvar<T,m>(0.5)) == math::chebyshev_clenshaw_recurrence(c.data(),c.size(),static_cast<T>(0.5)));
    }
    // TODO Continue alphabetically through <boost/math/special_functions.hpp>
  }
};

BOOST_AUTO_TEST_CASE(boost_special_functions)
{
    boost_special_functions_test()(static_cast<double>(0));
    //boost::fusion::for_each(bin_float_types, boost_special_functions_test());
    //boost::fusion::for_each(multiprecision_float_types, boost_special_functions_test());
}
*/

BOOST_AUTO_TEST_SUITE_END()
