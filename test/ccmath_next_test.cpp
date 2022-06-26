//  (C) Copyright John Maddock 2008 - 2022.
//  (C) Copyright Matt Borland 2022.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <iomanip>
#include <limits>
#include <boost/math/tools/precision.hpp>
#include <boost/math/ccmath/next.hpp>
#include <boost/math/ccmath/fpclassify.hpp>

#if !defined(BOOST_MATH_NO_CONSTEXPR_DETECTION) && !defined(BOOST_MATH_USING_BUILTIN_CONSTANT_P)
template <typename T>
void test_next()
{
    // NaN handling
    static_assert(boost::math::ccmath::isnan(boost::math::ccmath::nextafter(std::numeric_limits<T>::quiet_NaN(), T(0))));
    static_assert(boost::math::ccmath::isnan(boost::math::ccmath::nextafter(T(0), std::numeric_limits<T>::quiet_NaN())));

    // Handling of 0
    static_assert(boost::math::ccmath::nextafter(T(-0.0), T(0.0)) == T(0.0));
    static_assert(boost::math::ccmath::nextafter(T(0.0), T(-0.0)) == T(-0.0));

    // val = 1
    //static_assert(boost::math::ccmath::detail::float_distance(boost::math::ccmath::nextafter(1, boost::math::tools::max_value<T>()), T(1)) == -1);
}

int main(void)
{
    test_next<float>();
    test_next<double>();

    #ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
    test_next<long double>();
    #endif

    return 0;
}
#else
int main(void)
{
    return 0;
}
#endif
