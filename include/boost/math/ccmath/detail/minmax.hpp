//  (C) Copyright Matt Borland 2022.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_CCMATH_MINMAX_HPP
#define BOOST_MATH_CCMATH_MINMAX_HPP

#include <utility>

namespace boost::math::ccmath::detail {

template <typename T>
constexpr T min(T a, T b)
{
    return a > b ? b : a;
}

template <typename T>
constexpr T max(T a, T b)
{
    return a > b ? a : b;
}

}

#endif // BOOST_MATH_CCMATH_MINMAX_HPP
