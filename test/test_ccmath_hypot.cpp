//  (C) Copyright Matt Borland 2021.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cmath>
#include <cfloat>
#include <cstdint>
#include <limits>
#include <type_traits>
#include <boost/math/ccmath/hypot.hpp>
#include <boost/math/ccmath/isnan.hpp>
#include <boost/math/ccmath/isinf.hpp>
#include <boost/math/ccmath/sqrt.hpp>

#ifdef BOOST_HAS_FLOAT128
#include <boost/multiprecision/float128.hpp>
#endif

#if !defined(BOOST_MATH_NO_CONSTEXPR_DETECTION) && !defined(BOOST_MATH_USING_BUILTIN_CONSTANT_P)
template <typename T>
constexpr void test()
{
    // Error Handling
    if constexpr (std::numeric_limits<T>::has_quiet_NaN)
    {
        static_assert(boost::math::ccmath::isnan(boost::math::ccmath::hypot(std::numeric_limits<T>::quiet_NaN(), T(1))), "If x is NaN, NaN is returned");
        static_assert(boost::math::ccmath::isnan(boost::math::ccmath::hypot(T(1), std::numeric_limits<T>::quiet_NaN())), "If y is NaN, NaN is returned");
    }

    static_assert(boost::math::ccmath::isinf(boost::math::ccmath::hypot(std::numeric_limits<T>::infinity(), T(1))));
    static_assert(boost::math::ccmath::isinf(boost::math::ccmath::hypot(-std::numeric_limits<T>::infinity(), T(1))));
    static_assert(boost::math::ccmath::isinf(boost::math::ccmath::hypot(T(1), std::numeric_limits<T>::infinity())));
    static_assert(boost::math::ccmath::isinf(boost::math::ccmath::hypot(T(1), -std::numeric_limits<T>::infinity())));

    // Correct promoted types
    if constexpr (!std::is_same_v<T, float>)
    {
        constexpr auto test_type = boost::math::ccmath::hypot(T(1), 1.0f);
        static_assert(std::is_same_v<T, std::remove_cv_t<decltype(test_type)>>);
    }
    else
    {
        constexpr auto test_type = boost::math::ccmath::hypot(1.0f, 1);
        static_assert(std::is_same_v<double, std::remove_cv_t<decltype(test_type)>>);
    }

    // Functionality
    static_assert(boost::math::ccmath::hypot(T(1), T(1)) == boost::math::ccmath::sqrt(T(2)));
    static_assert(boost::math::ccmath::hypot(T(-1), T(1)) == boost::math::ccmath::sqrt(T(2)));
    static_assert(boost::math::ccmath::hypot(T(-1), T(-1)) == boost::math::ccmath::sqrt(T(2)));
    static_assert(boost::math::ccmath::hypot(T(1), T(-1)) == boost::math::ccmath::sqrt(T(2)));
    static_assert(boost::math::ccmath::hypot(T(1), T(2)) == boost::math::ccmath::sqrt(T(5)));
    static_assert(boost::math::ccmath::hypot(T(2), T(2)) == boost::math::ccmath::sqrt(T(8)));
}

template <typename T>
constexpr void test_3_arg()
{
    constexpr T tol = 5 * std::numeric_limits<T>::epsilon();
    
    // Error Handling
    if constexpr (std::numeric_limits<T>::has_quiet_NaN)
    {
        static_assert(boost::math::ccmath::isnan(boost::math::ccmath::hypot(std::numeric_limits<T>::quiet_NaN(), T(1), T(1))), "If x is NaN, NaN is returned");
        static_assert(boost::math::ccmath::isnan(boost::math::ccmath::hypot(T(1), std::numeric_limits<T>::quiet_NaN(), T(1))), "If y is NaN, NaN is returned");
        static_assert(boost::math::ccmath::isnan(boost::math::ccmath::hypot(T(1), T(1), std::numeric_limits<T>::quiet_NaN())), "If z is NaN, NaN is returned");
    }

    static_assert(boost::math::ccmath::isinf(boost::math::ccmath::hypot(std::numeric_limits<T>::infinity(), T(1), T(1))));
    static_assert(boost::math::ccmath::isinf(boost::math::ccmath::hypot(-std::numeric_limits<T>::infinity(), T(1), T(1))));
    static_assert(boost::math::ccmath::isinf(boost::math::ccmath::hypot(T(1), std::numeric_limits<T>::infinity(), T(1))));
    static_assert(boost::math::ccmath::isinf(boost::math::ccmath::hypot(T(1), -std::numeric_limits<T>::infinity(), T(1))));
    static_assert(boost::math::ccmath::isinf(boost::math::ccmath::hypot(T(1), T(1), std::numeric_limits<T>::infinity())));
    static_assert(boost::math::ccmath::isinf(boost::math::ccmath::hypot(T(1), T(1), -std::numeric_limits<T>::infinity())));
    
    // Correct promoted types
    if constexpr (!std::is_same_v<T, float>)
    {
        constexpr auto test_type = boost::math::ccmath::hypot(T(1), T(1), 1.0f);
        static_assert(std::is_same_v<T, std::remove_cv_t<decltype(test_type)>>);
    }
    else
    {
        constexpr auto test_type = boost::math::ccmath::hypot(1.0f, 1, 1);
        static_assert(std::is_same_v<double, std::remove_cv_t<decltype(test_type)>>);
    }
    
    // Functionality
    static_assert(boost::math::ccmath::hypot(T(1), T(1), T(1)) == boost::math::ccmath::sqrt(T(3)));
    static_assert(boost::math::ccmath::hypot(T(-1), T(1), T(1)) == boost::math::ccmath::sqrt(T(3)));
    static_assert(boost::math::ccmath::hypot(T(2), T(2), T(1)) == T(3));
    static_assert(boost::math::ccmath::hypot(T(2), T(-2), T(1)) == T(3));
    static_assert(boost::math::ccmath::abs(boost::math::ccmath::hypot(T(1), T(2), T(3)) - boost::math::ccmath::sqrt(T(14))) < tol);

    // Zeros: no 0/0
    static_assert(boost::math::ccmath::hypot(T(0), T(0), T(0)) == T(0));
    static_assert(boost::math::ccmath::hypot(T(-0.0), T(-0.0), T(-0.0)) == T(0));
    static_assert(boost::math::ccmath::hypot(T(0), T(0), T(-3)) == T(3));

    // inf wins over NaN, as for the two argument version
    if constexpr (std::numeric_limits<T>::has_quiet_NaN)
    {
        static_assert(boost::math::ccmath::isinf(boost::math::ccmath::hypot(std::numeric_limits<T>::infinity(), std::numeric_limits<T>::quiet_NaN(), T(1))));
    }

    // No spurious overflow when the squares overflow
    constexpr T big = (std::numeric_limits<T>::max)() / 4;
    static_assert(boost::math::ccmath::abs(boost::math::ccmath::hypot(big, big, big) / (big * boost::math::ccmath::sqrt(T(3))) - 1) < tol);
}

template <typename T>
bool test_3_arg_runtime()
{
    // Not constant evaluated, so this goes to std::hypot
    volatile T x = 2;
    volatile T y = 2;
    volatile T z = 1;
    return boost::math::ccmath::hypot(T(x), T(y), T(z)) == T(3);
}

int main()
{
    test<float>();
    test<double>();

    #ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
    test<long double>();
    #endif
    
    #ifdef BOOST_HAS_FLOAT128
    test<boost::multiprecision::float128>();
    #endif

    test_3_arg<float>();
    test_3_arg<double>();

    #ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
    test_3_arg<long double>();
    #endif

    if (!test_3_arg_runtime<float>() || !test_3_arg_runtime<double>())
    {
        return 1;
    }

    return 0;
}
#else
int main()
{
    return 0;
}
#endif
