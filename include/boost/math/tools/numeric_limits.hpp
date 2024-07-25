//  Copyright (c) 2024 Matt Borland
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
//  Regular use of std::numeric_limits functions can not be used on 
//  GPU platforms like CUDA since they are missing the __device__ marker
//  and libcu++ does not provide something analogous.
//  Rather than using giant if else blocks make our own version of numeric limits

#include <boost/math/tools/config.hpp>
#include <type_traits>
#include <limits>
#include <climits>

namespace boost {
namespace math {

template <typename T>
struct numeric_limits : public std::numeric_limits<T> {};

#ifdef BOOST_MATH_HAS_GPU_SUPPORT

template <>
struct numeric_limits<float>
{
    BOOST_MATH_STATIC constexpr bool is_specialized = std::numeric_limits<float>::is_specialized;
    BOOST_MATH_STATIC constexpr bool is_signed = std::numeric_limits<float>::is_signed;
    BOOST_MATH_STATIC constexpr bool is_integer = std::numeric_limits<float>::is_integer;
    BOOST_MATH_STATIC constexpr bool is_exact = std::numeric_limits<float>::is_exact;
    BOOST_MATH_STATIC constexpr bool has_infinity = std::numeric_limits<float>::has_infinity;
    BOOST_MATH_STATIC constexpr bool has_quiet_NaN = std::numeric_limits<float>::has_quiet_NaN;
    BOOST_MATH_STATIC constexpr bool has_signaling_NaN = std::numeric_limits<float>::has_signaling_NaN;

    BOOST_MATH_STATIC constexpr std::float_round_style round_style = std::numeric_limits<float>::round_style;
    BOOST_MATH_STATIC constexpr bool is_iec559 = std::numeric_limits<float>::is_iec559;
    BOOST_MATH_STATIC constexpr bool is_bounded = std::numeric_limits<float>::is_bounded;
    BOOST_MATH_STATIC constexpr bool is_modulo = std::numeric_limits<float>::is_modulo;
    BOOST_MATH_STATIC constexpr int digits = std::numeric_limits<float>::digits;
    BOOST_MATH_STATIC constexpr int digits10 = std::numeric_limits<float>::digits10;
    BOOST_MATH_STATIC constexpr int max_digits10 = std::numeric_limits<float>::max_digits10;
    BOOST_MATH_STATIC constexpr int radix = std::numeric_limits<float>::radix;
    BOOST_MATH_STATIC constexpr int min_exponent = std::numeric_limits<float>::min_exponent;
    BOOST_MATH_STATIC constexpr int min_exponent10 = std::numeric_limits<float>::min_exponent10;
    BOOST_MATH_STATIC constexpr int max_exponent = std::numeric_limits<float>::max_exponent;
    BOOST_MATH_STATIC constexpr int max_exponent10 = std::numeric_limits<float>::max_exponent10;
    BOOST_MATH_STATIC constexpr bool traps = std::numeric_limits<float>::traps;
    BOOST_MATH_STATIC constexpr bool tinyness_before = std::numeric_limits<float>::tinyness_before;

    // Member Functions
    BOOST_MATH_GPU_ENABLED BOOST_MATH_STATIC constexpr float (min)         () { return FLT_MIN; }
    BOOST_MATH_GPU_ENABLED BOOST_MATH_STATIC constexpr float (max)         () { return FLT_MAX; }
    BOOST_MATH_GPU_ENABLED BOOST_MATH_STATIC constexpr float lowest        () { return -FLT_MAX; }
    BOOST_MATH_GPU_ENABLED BOOST_MATH_STATIC constexpr float epsilon       () { return FLT_EPSILON; }
    BOOST_MATH_GPU_ENABLED BOOST_MATH_STATIC constexpr float round_error   () { return 0.5F; }
    BOOST_MATH_GPU_ENABLED BOOST_MATH_STATIC constexpr float infinity      () { return static_cast<float>(INFINITY); }
    BOOST_MATH_GPU_ENABLED BOOST_MATH_STATIC constexpr float quiet_NaN     () { return static_cast<float>(NAN); }
    BOOST_MATH_GPU_ENABLED BOOST_MATH_STATIC constexpr float signaling_NaN () 
    { 
        #ifdef FLT_SNAN
        return FLT_SNAN;
        #else
        return static_cast<float>(NAN);
        #endif
    }
    BOOST_MATH_GPU_ENABLED BOOST_MATH_STATIC constexpr float denorm_min    () { return FLT_TRUE_MIN; }
};

template <>
struct numeric_limits<double>
{
    BOOST_MATH_STATIC constexpr bool is_specialized = std::numeric_limits<double>::is_specialized;
    BOOST_MATH_STATIC constexpr bool is_signed = std::numeric_limits<double>::is_signed;
    BOOST_MATH_STATIC constexpr bool is_integer = std::numeric_limits<double>::is_integer;
    BOOST_MATH_STATIC constexpr bool is_exact = std::numeric_limits<double>::is_exact;
    BOOST_MATH_STATIC constexpr bool has_infinity = std::numeric_limits<double>::has_infinity;
    BOOST_MATH_STATIC constexpr bool has_quiet_NaN = std::numeric_limits<double>::has_quiet_NaN;
    BOOST_MATH_STATIC constexpr bool has_signaling_NaN = std::numeric_limits<double>::has_signaling_NaN;

    BOOST_MATH_STATIC constexpr std::float_round_style round_style = std::numeric_limits<double>::round_style;
    BOOST_MATH_STATIC constexpr bool is_iec559 = std::numeric_limits<double>::is_iec559;
    BOOST_MATH_STATIC constexpr bool is_bounded = std::numeric_limits<double>::is_bounded;
    BOOST_MATH_STATIC constexpr bool is_modulo = std::numeric_limits<double>::is_modulo;
    BOOST_MATH_STATIC constexpr int digits = std::numeric_limits<double>::digits;
    BOOST_MATH_STATIC constexpr int digits10 = std::numeric_limits<double>::digits10;
    BOOST_MATH_STATIC constexpr int max_digits10 = std::numeric_limits<double>::max_digits10;
    BOOST_MATH_STATIC constexpr int radix = std::numeric_limits<double>::radix;
    BOOST_MATH_STATIC constexpr int min_exponent = std::numeric_limits<double>::min_exponent;
    BOOST_MATH_STATIC constexpr int min_exponent10 = std::numeric_limits<double>::min_exponent10;
    BOOST_MATH_STATIC constexpr int max_exponent = std::numeric_limits<double>::max_exponent;
    BOOST_MATH_STATIC constexpr int max_exponent10 = std::numeric_limits<double>::max_exponent10;
    BOOST_MATH_STATIC constexpr bool traps = std::numeric_limits<double>::traps;
    BOOST_MATH_STATIC constexpr bool tinyness_before = std::numeric_limits<double>::tinyness_before;

    // Member Functions
    BOOST_MATH_GPU_ENABLED BOOST_MATH_STATIC constexpr double (min)         () { return DBL_MIN; }
    BOOST_MATH_GPU_ENABLED BOOST_MATH_STATIC constexpr double (max)         () { return DBL_MAX; }
    BOOST_MATH_GPU_ENABLED BOOST_MATH_STATIC constexpr double lowest        () { return -DBL_MAX; }
    BOOST_MATH_GPU_ENABLED BOOST_MATH_STATIC constexpr double epsilon       () { return DBL_EPSILON; }
    BOOST_MATH_GPU_ENABLED BOOST_MATH_STATIC constexpr double round_error   () { return 0.5; }
    BOOST_MATH_GPU_ENABLED BOOST_MATH_STATIC constexpr double infinity      () { return static_cast<double>(INFINITY); }
    BOOST_MATH_GPU_ENABLED BOOST_MATH_STATIC constexpr double quiet_NaN     () { return static_cast<double>(NAN); }
    BOOST_MATH_GPU_ENABLED BOOST_MATH_STATIC constexpr double signaling_NaN () 
    { 
        #ifdef DBL_SNAN
        return DBL_SNAN;
        #else
        return static_cast<double>(NAN);
        #endif
    }
    BOOST_MATH_GPU_ENABLED BOOST_MATH_STATIC constexpr double denorm_min    () { return DBL_TRUE_MIN; }
};

#endif

} // namespace math
} // namespace boost
