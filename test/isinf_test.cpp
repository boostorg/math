//  (C) Copyright Matt Borland 2021.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cmath>
#include <cfloat>
#include <limits>
#include <boost/math/ccmath/isinf.hpp>
#include <boost/core/lightweight_test.hpp>

template <typename T>
void test()
{
    constexpr bool test_val = boost::math::ccmath::isinf(T(0));
    static_assert(!test_val, "Not constexpr");

    static_assert(!boost::math::ccmath::isinf(std::numeric_limits<T>::quiet_NaN()));
    static_assert(boost::math::ccmath::isinf(std::numeric_limits<T>::infinity()));
    BOOST_TEST(boost::math::ccmath::isinf(std::exp(T(10'000'000))));
}

// Only test on platforms that provide BOOST_MATH_IS_CONSTANT_EVALUATED
#ifdef BOOST_MATH_IS_CONSTANT_EVALUATED
int main()
{
    test<float>();
    test<double>();

    #ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
    test<long double>();
    #endif
    
    return boost::report_errors();
}
#else
int main()
{
    return 0;
}
#endif
