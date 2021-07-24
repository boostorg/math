//  (C) Copyright Matt Borland 2021.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_CCMATH_ISINF
#define BOOST_MATH_CCMATH_ISINF

#include <limits>

namespace boost::math::ccmath {

template <typename T>
inline constexpr bool isinf(T x)
{
    return x == std::numeric_limits<T>::infinity() || -x == std::numeric_limits<T>::infinity();
}

}

#endif // BOOST_MATH_CCMATH_ISINF