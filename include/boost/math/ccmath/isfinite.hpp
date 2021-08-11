//  (C) Copyright Matt Borland 2021.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_CCMATH_ISFINITE
#define BOOST_MATH_CCMATH_ISFINITE

#include <boost/math/ccmath/isinf.hpp>
#include <boost/math/ccmath/isnan.hpp>

namespace boost::math::ccmath {

template <typename T>
inline constexpr bool isfinite(T x)
{
    return !boost::math::ccmath::isinf(x) && !boost::math::ccmath::isnan(x);
}

}

#endif // BOOST_MATH_CCMATH_ISFINITE
