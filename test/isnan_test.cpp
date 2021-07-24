//  (C) Copyright Matt Borland 2021.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cmath>
#include <cfloat>
#include <limits>
#include <boost/math/ccmath/isnan.hpp>
#include <boost/core/lightweight_test.hpp>

template <typename T>
void test()
{
    constexpr bool test_val = boost::math::ccmath::isnan(T(0));
    static_assert(!test_val, "Not constexpr");

    static_assert(boost::math::ccmath::isnan(std::numeric_limits<T>::quiet_NaN()), "Quiet NAN failed");
    static_assert(boost::math::ccmath::isnan(std::numeric_limits<T>::signaling_NaN()), "Signaling NAN failed");
    static_assert(!boost::math::ccmath::isnan(std::numeric_limits<T>::infinity()), "Infininty failed");
    static_assert(!boost::math::ccmath::isnan(T(0)), "Real 0 failed");
}

int main()
{
    test<float>();
    test<double>();

    #ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
    test<long double>();
    #endif

    return boost::report_errors();
}