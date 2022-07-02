//  (C) Copyright Matt Borland 2022.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_CCMATH_SIGNBIT_HPP
#define BOOST_MATH_CCMATH_SIGNBIT_HPP

#include <cmath>
#include <limits>
#include <type_traits>
#include <boost/math/tools/is_constant_evaluated.hpp>
#include <boost/math/ccmath/isnan.hpp>

namespace boost::math::ccmath {

namespace detail {

// Typical implementations of signbit involve type punning via union and manipulating
// overflow (see libc++ or musl). Neither of these are allowed in constexpr contexts
// (technically type punning via union in general is UB in c++ but well defined in C) 
// therefore NANs and 0s are treated as positive

template <typename T>
constexpr bool signbit_impl(T arg)
{
    if (boost::math::ccmath::isnan(arg))
    {
        return false;
    }
    
    return arg < static_cast<T>(0);
}

}

// Return value: true if arg is negative, false if arg is 0, NAN, or positive
template <typename Real, std::enable_if_t<!std::is_integral_v<Real>, bool> = true>
constexpr bool signbit(Real arg)
{
    if (BOOST_MATH_IS_CONSTANT_EVALUATED(arg))
    {
        return boost::math::ccmath::detail::signbit_impl(arg);
    }
    else
    {
        using std::signbit;
        return signbit(arg);
    }
}

template <typename Z, std::enable_if_t<std::is_integral_v<Z>, bool> = true>
constexpr bool signbit(Z arg)
{
    return boost::math::ccmath::signbit(static_cast<double>(arg));
}

} // Namespaces

#endif // BOOST_MATH_CCMATH_SIGNBIT_HPP
