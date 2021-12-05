//  (C) Copyright Matt Borland 2021.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
//  Each test is pulled from the larger ccmath_"function"_test.cpp
//  for minimal testing to make sure the module works as intended

#include <cmath>
#include <cfloat>
#include <limits>

import boost.math.ccmath;

template <typename T>
inline constexpr T base_helper(const T val)
{
    int i = 0;
    const T ans = boost::math::ccmath::frexp(val, &i);

    return ans;
}

template <typename T>
inline constexpr int exp_helper(const T val)
{
    int i = 0;
    boost::math::ccmath::frexp(val, &i);

    return i;
}

template <typename T>
inline constexpr T floating_point_value(const T val)
{
    T i = 0;
    const T ans = boost::math::ccmath::modf(val, &i);

    return ans;
}

template <typename T>
inline constexpr T integral_value(const T val)
{
    T i = 0;
    boost::math::ccmath::modf(val, &i);

    return i;
}

template <typename T>
void float_type_tests()
{
    static_assert(boost::math::ccmath::abs(T(-3)) == 3);
    static_assert(boost::math::ccmath::ceil(T(2.4)) == T(3));
    static_assert(boost::math::ccmath::copysign(T(1), T(-2)) == T(-1));
    static_assert(boost::math::ccmath::floor(T(2.9)) == T(2));
    static_assert(boost::math::ccmath::fmod(T(7.0/3), T(2.0) == T(1.0/3)));
    static_assert(boost::math::ccmath::fpclassify(T(0)) == FP_ZERO);
    static_assert(boost::math::ccmath::hypot(T(1), T(2)) == boost::math::ccmath::sqrt(T(5)));
    static_assert(boost::math::ccmath::isfinite(T(0)), "Wrong response to 0");
    static_assert(boost::math::ccmath::isinf(std::numeric_limits<T>::infinity()));
    static_assert(boost::math::ccmath::remainder(T(3.0/2), T(1.0) == T(3.0/2)));
    static_assert(boost::math::ccmath::round(T(2.3)) == T(2));
    static_assert(boost::math::ccmath::scalbln(T(1), 2l) == T(4));
    static_assert(boost::math::ccmath::scalbn(T(1.2), 10) == T(1228.8));
    static_assert(boost::math::ccmath::trunc(T(2.4)) == T(2));

    if constexpr (std::numeric_limits<T>::has_quiet_NaN)
    {
        static_assert(boost::math::ccmath::isnan(std::numeric_limits<T>::quiet_NaN()), "Quiet NAN failed");
        static_assert(!boost::math::ccmath::isnormal(std::numeric_limits<T>::quiet_NaN()), "Wrong response to quiet NAN");
    }

    // N[125/32, 30]
    // 3.90625000000000000000000000000
    // 0.976562500000000000000000000000 * 2^2
    constexpr T test_base = base_helper(T(125.0/32));
    static_assert(test_base == T(0.9765625));
    constexpr int test_exp = exp_helper(T(125.0/32));
    static_assert(test_exp == 2);

    // 123.45 = 1.92891 * 2^6
    constexpr int test_exp_ilogb = boost::math::ccmath::ilogb(T(123.45));
    static_assert(test_exp_ilogb == 6);

    // 1 * 2^2 = 4
    static_assert(boost::math::ccmath::ldexp(T(1), 2) == T(4));

    // 123.45 = 1.92891 * 2^6
    constexpr T test_exp_logb = boost::math::ccmath::logb(T(123.45));
    static_assert(test_exp_logb == T(6));

    // The returned value is exact, the current rounding mode is ignored
    // The return value and *iptr each have the same type and sign as x
    static_assert(integral_value(T(123.45)) == 123);
    static_assert(integral_value(T(-234.56)) == -234);
    static_assert(floating_point_value(T(1.0/2)) == T(1.0/2));
    static_assert(floating_point_value(T(-1.0/3)) == T(-1.0/3));

    // sqrt
    constexpr T tol = 2 * std::numeric_limits<T>::epsilon();
    constexpr T test_val = boost::math::ccmath::sqrt(T(2));
    constexpr T sqrt2 = T(1.4142135623730950488016887l);
    constexpr T abs_test_error = (test_val - sqrt2) > 0 ? (test_val - sqrt2) : (sqrt2 - test_val);
    static_assert(abs_test_error < tol, "Out of tolerance");
}

template <typename Z>
void integer_type_tests()
{
    constexpr auto test_val = boost::math::ccmath::div(Z(1'000'000), Z(3));
    static_assert(test_val.quot == Z(333'333));
    static_assert(test_val.rem == Z(1));
}

int main()
{
    float_type_tests<float>();
    float_type_tests<double>();
    float_type_tests<long double>();

    integer_type_tests<int>();
    integer_type_tests<long>();
    integer_type_tests<long long>();
    
    return 0;
}
