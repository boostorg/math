//  (C) Copyright Matt Borland 2021.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <limits>
#include <boost/math/ccmath/abs.hpp>
#include <boost/math/tools/config.hpp>

#ifdef BOOST_HAS_FLOAT128
#include <boost/multiprecision/float128.hpp>
#endif

template <typename T>
void test()
{
    static_assert(boost::math::ccmath::abs(T(3)) == 3);
    //static_assert(boost::math::ccmath::abs(T(-3)) == 3);
    //static_assert(boost::math::ccmath::abs(T(-0)) == 0);
    //static_assert(boost::math::ccmath::abs(-std::numeric_limits<T>::infinity()));
    //static_assert(boost::math::ccmath::abs(-std::numeric_limits<T>::quiet_NaN()));
}

// Only test on platforms that provide BOOST_MATH_IS_CONSTANT_EVALUATED
#ifdef BOOST_MATH_IS_CONSTANT_EVALUATED
int main()
{
    test<float>();
    /*
    test<double>();
    
    #ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
    test<long double>();
    #endif

    #if defined(BOOST_HAS_FLOAT128) && !defined(BOOST_MATH_USING_BUILTIN_CONSTANT_P)
    test<boost::multiprecision::float128>();
    #endif

    test<int>();
    test<long>();
    test<std::int32_t>();
    test<std::int64_t>();
    */
    return 0;
}
#else
int main()
{
    return 0;
}
#endif
