//  (C) Copyright Matt Borland 2021.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cmath>
#include <cfloat>
#include <cstdint>
#include <limits>
#include <type_traits>
#include <boost/math/ccmath/fmod.hpp>
#include <boost/math/ccmath/isnan.hpp>
#include <boost/math/ccmath/isinf.hpp>

#ifdef BOOST_HAS_FLOAT128
#include <boost/multiprecision/float128.hpp>
#endif

template <typename T>
constexpr void test()
{
    // Error Handling
    if constexpr (std::numeric_limits<T>::has_quiet_NaN)
    {
        static_assert(boost::math::ccmath::isnan(boost::math::ccmath::fmod(std::numeric_limits<T>::quiet_NaN(), T(1))), "If x is NaN, NaN is returned");
        static_assert(boost::math::ccmath::isnan(boost::math::ccmath::fmod(T(1), std::numeric_limits<T>::quiet_NaN())), "If y is NaN, NaN is returned");
    }

    static_assert(boost::math::ccmath::fmod(T(0), T(1.0)) == 0);
    static_assert(boost::math::ccmath::fmod(T(-0), T(1.0)) == -0);
    static_assert(boost::math::ccmath::isnan(boost::math::ccmath::fmod(std::numeric_limits<T>::infinity(), T(1.0))));
    static_assert(boost::math::ccmath::isnan(boost::math::ccmath::fmod(-std::numeric_limits<T>::infinity(), T(1.0))));

    static_assert(boost::math::ccmath::isnan(boost::math::ccmath::fmod(T(1), T(0))));
    static_assert(boost::math::ccmath::isnan(boost::math::ccmath::fmod(T(1), T(-0))));
    static_assert(boost::math::ccmath::fmod(T(1), std::numeric_limits<T>::infinity()) == T(1));
    static_assert(boost::math::ccmath::fmod(T(1), -std::numeric_limits<T>::infinity()) == T(1));

    // Functionality
    static_assert(boost::math::ccmath::fmod(T(3.0/2), T(1.0) == T(3.0/2)));
    static_assert(boost::math::ccmath::fmod(T(7.0/3), T(2.0) == T(1.0/3)));
    static_assert(boost::math::ccmath::fmod(T(-8.0/3), T(2.0) == T(-2.0/3)));

    // Correct promoted types
    if constexpr (!std::is_same_v<T, float>)
    {
        constexpr auto test_type = boost::math::ccmath::fmod(T(1), 1.0f);
        static_assert(std::is_same_v<T, std::remove_cv_t<decltype(test_type)>>);
    }
    else
    {
        constexpr auto test_type = boost::math::ccmath::fmod(1.0f, 1);
        static_assert(std::is_same_v<double, std::remove_cv_t<decltype(test_type)>>);
    }
}

#if !defined(BOOST_MATH_NO_CONSTEXPR_DETECTION) && !defined(BOOST_MATH_USING_BUILTIN_CONSTANT_P)
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

    return 0;
}
#else
int main()
{
    return 0;
}
#endif
