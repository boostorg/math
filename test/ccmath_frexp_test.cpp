//  (C) Copyright Matt Borland 2021.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cmath>
#include <cfloat>
#include <cstdint>
#include <limits>
#include <type_traits>
#include <boost/math/ccmath/frexp.hpp>
#include <boost/math/ccmath/isnan.hpp>
#include <boost/math/ccmath/isinf.hpp>

#ifdef BOOST_HAS_FLOAT128
#include <boost/multiprecision/float128.hpp>
#endif

template <typename T>
void test()
{
    static int i;
    if constexpr (std::numeric_limits<T>::has_quiet_NaN)
    {
        static_assert(boost::math::ccmath::isnan(boost::math::ccmath::frexp(std::numeric_limits<T>::quiet_NaN(), &i)));
    }

    static_assert(boost::math::ccmath::frexp(0, &i) == 0);
    static_assert(boost::math::ccmath::isinf(boost::math::ccmath::frexp(std::numeric_limits<T>::infinity(), &i)));
    static_assert(boost::math::ccmath::frexp(T(1024), &i) == T(0.5));
}

#ifndef BOOST_MATH_NO_CONSTEXPR_DETECTION
int main()
{
    test<float>();
    test<double>();

    #ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
    test<long double>();
    #endif
    
    #if defined(BOOST_HAS_FLOAT128) && !defined(BOOST_MATH_USING_BUILTIN_CONSTANT_P)
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
