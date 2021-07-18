//  (C) Copyright Matt Borland 2021.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cmath>
#include <limits>
#include <type_traits>
#include <boost/math/ccmath/sqrt.hpp>
#include <boost/core/lightweight_test.hpp>

template <typename Real>
void test_float_sqrt()
{
    using std::abs;
    
    constexpr Real tol = 2*std::numeric_limits<Real>::epsilon();
    
    constexpr Real test_val = boost::math::ccmath::sqrt(Real(2));
    constexpr Real dummy = 1;
    static_assert(test_val > dummy, "Not constexpr");

    Real known_val = std::sqrt(Real(2));
    BOOST_TEST(abs(test_val - known_val) < tol);

    // 1000 eps
    constexpr Real test_1000 = boost::math::ccmath::sqrt(1000*std::numeric_limits<Real>::epsilon());
    Real known_1000 = std::sqrt(1000*std::numeric_limits<Real>::epsilon());
    BOOST_TEST(abs(test_1000 - known_1000) < tol);

    // inf
    constexpr Real test_inf = boost::math::ccmath::sqrt(std::numeric_limits<Real>::infinity());
    Real known_inf = std::sqrt(std::numeric_limits<Real>::infinity());
    BOOST_TEST_EQ(test_inf, known_inf);

    // NAN
    constexpr Real test_nan = boost::math::ccmath::sqrt(std::numeric_limits<Real>::quiet_NaN());
    Real known_nan = std::sqrt(std::numeric_limits<Real>::quiet_NaN());
    BOOST_TEST(std::isnan(test_nan) && std::isnan(known_nan));

    // 100'000'000
    constexpr Real test_100m = boost::math::ccmath::sqrt(100000000);
    BOOST_TEST_EQ(test_100m, Real(10000));

    // MAX / 2
    // Only tests float since double and long double will exceed maximum template depth
    if constexpr (std::is_same_v<float, Real>)
    {
        constexpr Real test_max = boost::math::ccmath::sqrt(std::numeric_limits<Real>::max() / 2);
        Real known_max = std::sqrt(std::numeric_limits<Real>::max() / 2);
        BOOST_TEST(abs(test_max - known_max) < tol);
    }
}

template <typename Z>
void test_int_sqrt()
{
    using std::abs;

    constexpr double tol = 2*std::numeric_limits<double>::epsilon();

    constexpr double test_val = boost::math::ccmath::sqrt(Z(2));
    constexpr double dummy = 1;
    static_assert(test_val > dummy, "Not constexpr");

    double known_val = std::sqrt(2.0);

    BOOST_TEST(abs(test_val - known_val) < tol);
}

int main()
{
    test_float_sqrt<float>();
    test_float_sqrt<double>();
    
    #ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
    test_float_sqrt<long double>();
    #endif

    test_int_sqrt<int>();
    test_int_sqrt<unsigned>();
    
    return boost::report_errors();
}
