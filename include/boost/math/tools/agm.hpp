//  (C) Copyright Nick Thompson 2020.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_TOOLS_AGM_HPP
#define BOOST_MATH_TOOLS_AGM_HPP
#include <boost/math/special_functions/next.hpp>
#include <cmath>

namespace boost::math::tools {

namespace detail {
    #pragma GCC diagnostic ignored "-Wstrict-aliasing"
    int32_t fast_float_distance(float x, float y) {
        static_assert(sizeof(float) == sizeof(int32_t), "float is incorrect size.");
        int32_t xi = *reinterpret_cast<int32_t*>(&x);
        int32_t yi = *reinterpret_cast<int32_t*>(&y);
        return yi - xi;
    }

    int64_t fast_float_distance(double x, double y) {
        static_assert(sizeof(double) == sizeof(int64_t), "double is incorrect size.");
        int64_t xi = *reinterpret_cast<int64_t*>(&x);
        int64_t yi = *reinterpret_cast<int64_t*>(&y);
        return yi - xi;
    }
    #pragma GCC diagnostic pop
}
template<typename Real>
Real agm(Real a, Real g)
{
    if (a < g)
    {
        throw std::domain_error("a >= g is required");
        return std::numeric_limits<Real>::quiet_NaN();
    }
    // Use: M(rx, ry) = rM(x,y)
    if (a <= 0 || g <= 0) {
        if (a < 0 || g < 0) {
             throw std::domain_error("a > 0, g > 0 is required");
             return std::numeric_limits<Real>::quiet_NaN();
        }
        return Real(0);
    }

    // a > g:
    if constexpr (std::is_same_v<float, Real> || std::is_same_v<double, Real>)
    {
        while (detail::fast_float_distance(g, a) > 2000)
        {
            Real anp1 = (a + g)/2;
            g = sqrt(a*g);
            a = anp1;
        }
    }
    else
    {
        while (boost::math::float_distance(g, a) > 2000)
        {
            Real anp1 = (a + g)/2;
            g = sqrt(a*g);
            a = anp1;
        }
    }

    // Final cleanup iteration recovers down to ~2ULPs:
    return (a + g)/2;
}



}
#endif
