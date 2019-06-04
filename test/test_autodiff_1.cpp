//           Copyright Matthew Pulver 2018 - 2019.
// Distributed under the Boost Software License, Version 1.0.
//      (See accompanying file LICENSE_1_0.txt or copy at
//           https://www.boost.org/LICENSE_1_0.txt)

#include "test_autodiff.hpp"

BOOST_AUTO_TEST_SUITE(test_autodiff_1)

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

BOOST_AUTO_TEST_SUITE_END()
