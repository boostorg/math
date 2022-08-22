//  (C) Copyright Matt Borland 2022.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_SF_FAST_FLOAT_DISTANCE
#define BOOST_MATH_SF_FAST_FLOAT_DISTANCE

#include <boost/math/special_functions/next.hpp>
#include <boost/math/tools/throw_exception.hpp>
#include <stdexcept>

#if defined(BOOST_MATH_USE_FLOAT128) && !defined(BOOST_MATH_STANDALONE)
#include <boost/multiprecision/float128.hpp>
#include <boost/multiprecision/detail/standalone_config.hpp>
#define BOOST_MATH_USE_FAST_FLOAT128
#endif

namespace boost { namespace math { 

// https://randomascii.wordpress.com/2012/01/23/stupid-float-tricks-2/
// https://blog.regehr.org/archives/959
inline std::int32_t fast_float_distance(float a, float b)
{
    return boost::math::float_distance(a, b);
}

inline std::int64_t fast_float_distance(double a, double b)
{
    return boost::math::float_distance(a, b);
}

#ifdef BOOST_MATH_USE_FAST_FLOAT128
boost::multiprecision::int128_type fast_float_distance(boost::multiprecision::float128_type a, boost::multiprecision::float128_type b)
{
    using std::abs;
    using std::isfinite;

    constexpr boost::multiprecision::float128_type tol = BOOST_MP_QUAD_MIN;

    // 0, very small, and large magnitude distances all need special handling
    if (abs(a) == 0 || abs(b) == 0)
    {
        return 0;
    }
    else if (abs(a) < tol || abs(b) < tol)
    {
        BOOST_MATH_THROW_EXCEPTION(std::domain_error("special handling is required for tiny distances. Please use boost::math::float_distance for a slower but safe solution"));
    }

    if (!(isfinite)(a))
    {  
        BOOST_MATH_THROW_EXCEPTION(std::domain_error("Both arguments to fast_float_distnace must be finite"));
    }
    else if (!(isfinite)(b))
    {
        BOOST_MATH_THROW_EXCEPTION(std::domain_error("Both arguments to fast_float_distnace must be finite"));
    }

    static_assert(sizeof(boost::multiprecision::int128_type) == sizeof(boost::multiprecision::float128_type));

    boost::multiprecision::int128_type ai;
    boost::multiprecision::int128_type bi;
    std::memcpy(&ai, &a, sizeof(boost::multiprecision::float128_type));
    std::memcpy(&bi, &b, sizeof(boost::multiprecision::float128_type));

    boost::multiprecision::int128_type result = bi - ai;

    if (ai < 0 || bi < 0)
    {
        result = -result;
    }

    return result;
}

#endif

}} // Namespaces

#endif
