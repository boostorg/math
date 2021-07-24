//  (C) Copyright Matt Borland 2021.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_CCMATH_ISINF
#define BOOST_MATH_CCMATH_ISINF

#include <cmath>
#include <limits>
#include <boost/math/tools/is_constant_evaluated.hpp>

namespace boost::math::ccmath {

template <typename T>
inline constexpr bool isinf(T x)
{
    if(BOOST_MATH_IS_CONSTANT_EVALUATED(x))
    {
        return x == std::numeric_limits<T>::infinity() || -x == std::numeric_limits<T>::infinity();
    }
    else
    {
        using std::isinf;
        return isinf(x);
    }
}

}

#endif // BOOST_MATH_CCMATH_ISINF
