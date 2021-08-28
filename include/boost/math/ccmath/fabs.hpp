//  (C) Copyright Matt Borland 2021.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
//  Constepxr implementation of fabs (see c.math.abs secion 26.8.2 of the ISO standard)

#ifndef BOOST_MATH_CCMATH_FABS
#define BOOST_MATH_CCMATH_FABS

#include <boost/math/ccmath/abs.hpp>

namespace boost::math::ccmath {

template <typename T>
inline constexpr auto fabs(T x)
{
    return boost::math::ccmath::abs(x);
}

}

#endif // BOOST_MATH_CCMATH_FABS
