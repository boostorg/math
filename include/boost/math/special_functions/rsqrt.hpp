//  (C) Copyright Nick Thompson 2020.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_SPECIAL_FUNCTIONS_RSQRT_HPP
#define BOOST_MATH_SPECIAL_FUNCTIONS_RSQRT_HPP
#include <cmath>
#include <type_traits>
#include <limits>
#include <boost/math/ccmath/sqrt.hpp>

namespace boost { namespace math {

template <typename Real, typename std::enable_if<std::is_same<Real, float>::value ||
                                                 std::is_same<Real, double>::value ||
                                                 std::is_same<Real, long double>::value, bool>::type = true>
inline constexpr Real rsqrt(Real const & x)
{
    return 1 / boost::math::ccmath::sqrt(x);
}

template<typename Real, typename std::enable_if<!std::is_same<Real, float>::value &&
                                                !std::is_same<Real, double>::value &&
                                                !std::is_same<Real, long double>::value, bool>::type = true>
inline Real rsqrt(Real const & x)
{
    using std::sqrt;

    // if it's so tiny it rounds to 0 as long double,
    // no performance gains are possible:
    if (x < std::numeric_limits<long double>::denorm_min() || x > (std::numeric_limits<long double>::max)()) {
        return 1/sqrt(x);
    }
    Real x0 = 1/sqrt(static_cast<long double>(x));
    // Divide by 512 for leeway:
    Real s = sqrt(std::numeric_limits<Real>::epsilon())*x0/512;
    Real x1 = x0 + x0*(1-x*x0*x0)/2;
    while(abs(x1 - x0) > s) {
        x0 = x1;
        x1 = x0 + x0*(1-x*x0*x0)/2;
    }
    // Final iteration get ~2ULPs:
    return  x1 + x1*(1-x*x1*x1)/2;
}


}}
#endif
