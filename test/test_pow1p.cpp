//  (C) Copyright Matt Borland 2024.
//  (C) Copyrigh Fancidev 2024.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/math/special_functions/pow1p.hpp>
#include <boost/math/concepts/real_concept.hpp>
#include <exception>
#include <random>
#include "math_unit_test.hpp"

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
    test<float>();
    test<double>();

    #ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
    test<long double>();
    #endif

    #ifndef BOOST_MATH_NO_REAL_CONCEPT_TESTS
    test<boost::math::concepts::real_concept>();
    #endif

    return boost::math::test::report_errors();
}
