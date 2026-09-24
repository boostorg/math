//  (C) Copyright Matt Borland 2024.
//  (C) Copyrigh Fancidev 2024.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/math/special_functions/pow1p.hpp>
#include <boost/math/concepts/real_concept.hpp>
#include <exception>
#include <random>

#if __has_include(<stdfloat>) && !defined(BOOST_MATH_HAS_GPU_SUPPORT)
#  include <stdfloat>
#endif

#include "math_unit_test.hpp"

#ifndef BOOST_MATH_STANDALONE
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>

// Expression-template types must instantiate.  Arguments are kept moderate because
// pow1p calls pow(1+x, y), and multiprecision pow loses accuracy for huge y.
template <typename T>
void test_multiprecision()
{
    using big = boost::multiprecision::cpp_bin_float_100;
    CHECK_ULP_CLOSE(static_cast<T>(pow(big(1.5), big(2.25))), boost::math::pow1p(T(0.5), T(2.25)), 10);
    CHECK_ULP_CLOSE(static_cast<T>(pow(big(0.75), big(-7.5))), boost::math::pow1p(T(-0.25), T(-7.5)), 10);
    CHECK_ULP_CLOSE(static_cast<T>(pow(1 + big(T(0.001)), big(1000.5))), boost::math::pow1p(T(0.001), T(1000.5)), 10);
}
#endif

template <typename T>
void test()
{
    using std::pow;

    // First we hit all the special cases
    // pow(x, +/-0)
    CHECK_EQUAL(boost::math::pow1p(T(1), T(0)), T(1));
    
    // pow(0, y)
    #ifndef BOOST_MATH_NO_EXCEPTIONS
    CHECK_THROW(boost::math::pow1p(T(-1), T(-1)), std::domain_error);
    #endif
    CHECK_EQUAL(boost::math::pow1p(T(-1), T(1)), T(0));

    // pow(-1, inf)
    CHECK_EQUAL(boost::math::pow1p(T(-2), boost::math::numeric_limits<T>::infinity()), T(1));
    
    // pow(1, y)
    CHECK_EQUAL(boost::math::pow1p(T(0), T(2)), T(1));

    // pow(x, +/-inf)
    BOOST_MATH_IF_CONSTEXPR (boost::math::numeric_limits<T>::has_infinity)
    {
        CHECK_EQUAL(boost::math::pow1p(T(5), -boost::math::numeric_limits<T>::infinity()), T(0));
        CHECK_EQUAL(boost::math::pow1p(T(5), boost::math::numeric_limits<T>::infinity()), boost::math::numeric_limits<T>::infinity());
        // |1+x| < 1 reverses the limits:
        CHECK_EQUAL(boost::math::pow1p(T(-0.5), boost::math::numeric_limits<T>::infinity()), T(0));
        CHECK_EQUAL(boost::math::pow1p(T(-0.5), -boost::math::numeric_limits<T>::infinity()), boost::math::numeric_limits<T>::infinity());
        CHECK_EQUAL(boost::math::pow1p(T(-1.5), boost::math::numeric_limits<T>::infinity()), T(0));
        CHECK_EQUAL(boost::math::pow1p(T(-3), boost::math::numeric_limits<T>::infinity()), boost::math::numeric_limits<T>::infinity());
        CHECK_EQUAL(boost::math::pow1p(T(-3), -boost::math::numeric_limits<T>::infinity()), T(0));
        BOOST_MATH_IF_CONSTEXPR (boost::math::numeric_limits<T>::has_quiet_NaN)
        {
            CHECK_EQUAL(boost::math::isnan(boost::math::pow1p(boost::math::numeric_limits<T>::quiet_NaN(), boost::math::numeric_limits<T>::infinity())), true);
        }

        // Tiny x, huge finite y: the exponent y*log1p(x) is far out of range, and the
        // result must overflow or underflow cleanly rather than become inf*0 = NaN.
        using std::ldexp;
        const T tiny = ldexp(T(1), -boost::math::tools::digits<T>() - 7);
        const T huge = boost::math::tools::max_value<T>() / 4;
        CHECK_EQUAL(boost::math::pow1p(tiny, huge), boost::math::numeric_limits<T>::infinity());
        CHECK_EQUAL(boost::math::pow1p(tiny, -huge), T(0));
        CHECK_EQUAL(boost::math::pow1p(-tiny, huge), T(0));
        CHECK_EQUAL(boost::math::pow1p(-tiny, -huge), boost::math::numeric_limits<T>::infinity());

    // pow(+/-inf, y)
        CHECK_EQUAL(boost::math::pow1p(boost::math::numeric_limits<T>::infinity(), T(2)), boost::math::numeric_limits<T>::infinity());
        CHECK_EQUAL(boost::math::pow1p(-boost::math::numeric_limits<T>::infinity(), T(2)), boost::math::numeric_limits<T>::infinity());
    }

    // NANs for x and y
    BOOST_MATH_IF_CONSTEXPR (boost::math::numeric_limits<T>::has_quiet_NaN)
    {
        CHECK_EQUAL(boost::math::isnan(boost::math::pow1p(boost::math::numeric_limits<T>::quiet_NaN(), T(1))), true);
        CHECK_EQUAL(boost::math::isnan(boost::math::pow1p(T(1), boost::math::numeric_limits<T>::quiet_NaN())), true);
    }

    // pow(x, +/-1)
    CHECK_ULP_CLOSE(boost::math::pow1p(T(2), T(1)), pow(T(3), T(1)), 10);
    CHECK_ULP_CLOSE(boost::math::pow1p(T(2), T(-1)), pow(T(3), T(-1)), 10);
    
    // (1+x) < 0
    CHECK_ULP_CLOSE(boost::math::pow1p(T(-3), T(2)), pow(T(-2), T(2)), 10);

    // Tiny x, huge y: 1+x is inexact and the double-T branch is taken.
    // Reference values are exp(y*log1p(x)) computed at 100 digits.
    {
        using std::ldexp;
        const T x = ldexp(T(3), -40);
        const T y = ldexp(T(5), 42);
        CHECK_ULP_CLOSE(T(1.1420073897222058133288923989033585274659e+26L), boost::math::pow1p(x, y), 10);
        CHECK_ULP_CLOSE(T(8.7565107619797603254946543126347414974562e-27L), boost::math::pow1p(-x, y), 10);
        CHECK_ULP_CLOSE(T(8.7565107634132803515388737801922127069804e-27L), boost::math::pow1p(x, -y), 10);
        CHECK_ULP_CLOSE(T(1.1420073899091627540050136756953893069185e+26L), boost::math::pow1p(-x, -y), 10);
    }
    if (boost::math::tools::log_max_value<T>() > 480)
    {
        using std::ldexp;
        const T x = ldexp(T(3), -70);
        const T y = ldexp(T(5), 75);
        CHECK_ULP_CLOSE(T(2.8930191842539452447704692545247205720025e+208L), boost::math::pow1p(x, y), 10);
        CHECK_ULP_CLOSE(T(3.4565965045886172534873245069347178738856e-209L), boost::math::pow1p(-x, y), 10);
        CHECK_ULP_CLOSE(T(3.4565965045886172577034301265755950639994e-209L), boost::math::pow1p(x, -y), 10);
        CHECK_ULP_CLOSE(T(2.8930191842539452482991641497113860091743e+208L), boost::math::pow1p(-x, -y), 10);
    }
    // As above, but with full 53-bit mantissas so that the double-T products are not trivially exact.
    // The inputs are double literals, hence exact for T at least as wide as double.
    if ((boost::math::tools::digits<T>() >= 53) && (boost::math::tools::log_max_value<T>() > 662))
    {
        CHECK_ULP_CLOSE(T(7.6653726844726161792796892810864958094044e+261L), boost::math::pow1p(T(5.100075373549856e-19), T(1.1823578638298747e+21)), 10);
        CHECK_ULP_CLOSE(T(1.9141657450800789253871830037475352857803e-261L), boost::math::pow1p(T(-2.1413450492090156e-22), T(2.8034969310889015e+24)), 10);
        CHECK_ULP_CLOSE(T(1.5108017556765851783235726673009126604833e+287L), boost::math::pow1p(T(6.3377799876470245e-17), T(1.0433536087590337e+19)), 10);
    }

    // x < 0
    std::mt19937_64 rng;
    std::uniform_real_distribution<double> dist (-1, 0);
    std::uniform_real_distribution<double> dist_y (0, 10);
    constexpr int N = 1024;
    for (int i = 0; i < N; ++i)
    {
        const auto x = static_cast<T>(dist(rng));
        const auto y = static_cast<T>(dist_y(rng));

        CHECK_ULP_CLOSE(boost::math::pow1p(x, y), pow(x + 1, y), 100);
    }

    // 0 < x < 1
    std::uniform_real_distribution<double> dist_x_1(0, 1);
    for (int i = 0; i < N; ++i)
    {
        const auto x = static_cast<T>(dist_x_1(rng));
        const auto y = static_cast<T>(dist_y(rng));

        CHECK_ULP_CLOSE(boost::math::pow1p(x, y), pow(x + 1, y), 100);
    }

    // Else
    std::uniform_real_distribution<double> dist_other_x(1, 1000);
    for (int i = 0; i < N; ++i)
    {
        const auto x = static_cast<T>(dist_other_x(rng));
        const auto y = static_cast<T>(dist_y(rng));

        CHECK_ULP_CLOSE(boost::math::pow1p(x, y), pow(x + 1, y), 100);
    }
}

int main()
{
    #ifdef __STDCPP_FLOAT32_T__
    test<std::float32_t>();
    #else
    test<float>();
    #endif

    #ifdef __STDCPP_FLOAT64_T__
    test<std::float64_t>();
    #else
    test<double>();
    #endif

    #ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
    test<long double>();
    #endif

    #ifndef BOOST_MATH_NO_REAL_CONCEPT_TESTS
    test<boost::math::concepts::real_concept>();
    #endif

    #ifndef BOOST_MATH_STANDALONE
    test_multiprecision<boost::multiprecision::cpp_bin_float_50>();
    test_multiprecision<boost::multiprecision::cpp_dec_float_50>();
    #endif

    return boost::math::test::report_errors();
}
