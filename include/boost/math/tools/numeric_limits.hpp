//  Copyright (c) 2024 Matt Borland
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
//  Regular use of std::numeric_limits functions can not be used on 
//  GPU platforms like CUDA since they are missing the __device__ marker
//  and libcu++ does not provide something analogous.
//  Rather than using giant if else blocks make our own version of numeric limits

#ifndef BOOST_MATH_TOOLS_NUMERIC_LIMITS_HPP
#define BOOST_MATH_TOOLS_NUMERIC_LIMITS_HPP

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

template <>
struct numeric_limits<short>
{
    BOOST_MATH_STATIC constexpr bool is_specialized = std::numeric_limits<short>::is_specialized;
    BOOST_MATH_STATIC constexpr bool is_signed = std::numeric_limits<short>::is_signed;
    BOOST_MATH_STATIC constexpr bool is_integer = std::numeric_limits<short>::is_integer;
    BOOST_MATH_STATIC constexpr bool is_exact = std::numeric_limits<short>::is_exact;
    BOOST_MATH_STATIC constexpr bool has_infinity = std::numeric_limits<short>::has_infinity;
    BOOST_MATH_STATIC constexpr bool has_quiet_NaN = std::numeric_limits<short>::has_quiet_NaN;
    BOOST_MATH_STATIC constexpr bool has_signaling_NaN = std::numeric_limits<short>::has_signaling_NaN;

    BOOST_MATH_STATIC constexpr std::float_round_style round_style = std::numeric_limits<short>::round_style;
    BOOST_MATH_STATIC constexpr bool is_iec559 = std::numeric_limits<short>::is_iec559;
    BOOST_MATH_STATIC constexpr bool is_bounded = std::numeric_limits<short>::is_bounded;
    BOOST_MATH_STATIC constexpr bool is_modulo = std::numeric_limits<short>::is_modulo;
    BOOST_MATH_STATIC constexpr int digits = std::numeric_limits<short>::digits;
    BOOST_MATH_STATIC constexpr int digits10 = std::numeric_limits<short>::digits10;
    BOOST_MATH_STATIC constexpr int max_digits10 = std::numeric_limits<short>::max_digits10;
    BOOST_MATH_STATIC constexpr int radix = std::numeric_limits<short>::radix;
    BOOST_MATH_STATIC constexpr int min_exponent = std::numeric_limits<short>::min_exponent;
    BOOST_MATH_STATIC constexpr int min_exponent10 = std::numeric_limits<short>::min_exponent10;
    BOOST_MATH_STATIC constexpr int max_exponent = std::numeric_limits<short>::max_exponent;
    BOOST_MATH_STATIC constexpr int max_exponent10 = std::numeric_limits<short>::max_exponent10;
    BOOST_MATH_STATIC constexpr bool traps = std::numeric_limits<short>::traps;
    BOOST_MATH_STATIC constexpr bool tinyness_before = std::numeric_limits<short>::tinyness_before;

    // Member Functions
    BOOST_MATH_GPU_ENABLED BOOST_MATH_STATIC constexpr short (min)         () { return SHRT_MIN; }
    BOOST_MATH_GPU_ENABLED BOOST_MATH_STATIC constexpr short (max)         () { return SHRT_MAX; }
    BOOST_MATH_GPU_ENABLED BOOST_MATH_STATIC constexpr short lowest        () { return SHRT_MIN; }
    BOOST_MATH_GPU_ENABLED BOOST_MATH_STATIC constexpr short epsilon       () { return 0; }
    BOOST_MATH_GPU_ENABLED BOOST_MATH_STATIC constexpr short round_error   () { return 0; }
    BOOST_MATH_GPU_ENABLED BOOST_MATH_STATIC constexpr short infinity      () { return 0; }
    BOOST_MATH_GPU_ENABLED BOOST_MATH_STATIC constexpr short quiet_NaN     () { return 0; }
    BOOST_MATH_GPU_ENABLED BOOST_MATH_STATIC constexpr short signaling_NaN () { return 0; }
    BOOST_MATH_GPU_ENABLED BOOST_MATH_STATIC constexpr short denorm_min    () { return 0; }
};

template <>
struct numeric_limits<unsigned short>
{
    BOOST_MATH_STATIC constexpr bool is_specialized = std::numeric_limits<unsigned short>::is_specialized;
    BOOST_MATH_STATIC constexpr bool is_signed = std::numeric_limits<unsigned short>::is_signed;
    BOOST_MATH_STATIC constexpr bool is_integer = std::numeric_limits<unsigned short>::is_integer;
    BOOST_MATH_STATIC constexpr bool is_exact = std::numeric_limits<unsigned short>::is_exact;
    BOOST_MATH_STATIC constexpr bool has_infinity = std::numeric_limits<unsigned short>::has_infinity;
    BOOST_MATH_STATIC constexpr bool has_quiet_NaN = std::numeric_limits<unsigned short>::has_quiet_NaN;
    BOOST_MATH_STATIC constexpr bool has_signaling_NaN = std::numeric_limits<unsigned short>::has_signaling_NaN;

    BOOST_MATH_STATIC constexpr std::float_round_style round_style = std::numeric_limits<unsigned short>::round_style;
    BOOST_MATH_STATIC constexpr bool is_iec559 = std::numeric_limits<unsigned short>::is_iec559;
    BOOST_MATH_STATIC constexpr bool is_bounded = std::numeric_limits<unsigned short>::is_bounded;
    BOOST_MATH_STATIC constexpr bool is_modulo = std::numeric_limits<unsigned short>::is_modulo;
    BOOST_MATH_STATIC constexpr int digits = std::numeric_limits<unsigned short>::digits;
    BOOST_MATH_STATIC constexpr int digits10 = std::numeric_limits<unsigned short>::digits10;
    BOOST_MATH_STATIC constexpr int max_digits10 = std::numeric_limits<unsigned short>::max_digits10;
    BOOST_MATH_STATIC constexpr int radix = std::numeric_limits<unsigned short>::radix;
    BOOST_MATH_STATIC constexpr int min_exponent = std::numeric_limits<unsigned short>::min_exponent;
    BOOST_MATH_STATIC constexpr int min_exponent10 = std::numeric_limits<unsigned short>::min_exponent10;
    BOOST_MATH_STATIC constexpr int max_exponent = std::numeric_limits<unsigned short>::max_exponent;
    BOOST_MATH_STATIC constexpr int max_exponent10 = std::numeric_limits<unsigned short>::max_exponent10;
    BOOST_MATH_STATIC constexpr bool traps = std::numeric_limits<unsigned short>::traps;
    BOOST_MATH_STATIC constexpr bool tinyness_before = std::numeric_limits<unsigned short>::tinyness_before;

    // Member Functions
    BOOST_MATH_GPU_ENABLED BOOST_MATH_STATIC constexpr unsigned short (min)         () { return 0; }
    BOOST_MATH_GPU_ENABLED BOOST_MATH_STATIC constexpr unsigned short (max)         () { return USHRT_MAX; }
    BOOST_MATH_GPU_ENABLED BOOST_MATH_STATIC constexpr unsigned short lowest        () { return 0; }
    BOOST_MATH_GPU_ENABLED BOOST_MATH_STATIC constexpr unsigned short epsilon       () { return 0; }
    BOOST_MATH_GPU_ENABLED BOOST_MATH_STATIC constexpr unsigned short round_error   () { return 0; }
    BOOST_MATH_GPU_ENABLED BOOST_MATH_STATIC constexpr unsigned short infinity      () { return 0; }
    BOOST_MATH_GPU_ENABLED BOOST_MATH_STATIC constexpr unsigned short quiet_NaN     () { return 0; }
    BOOST_MATH_GPU_ENABLED BOOST_MATH_STATIC constexpr unsigned short signaling_NaN () { return 0; }
    BOOST_MATH_GPU_ENABLED BOOST_MATH_STATIC constexpr unsigned short denorm_min    () { return 0; }
};

template <>
struct numeric_limits<int>
{
    BOOST_MATH_STATIC constexpr bool is_specialized = std::numeric_limits<int>::is_specialized;
    BOOST_MATH_STATIC constexpr bool is_signed = std::numeric_limits<int>::is_signed;
    BOOST_MATH_STATIC constexpr bool is_integer = std::numeric_limits<int>::is_integer;
    BOOST_MATH_STATIC constexpr bool is_exact = std::numeric_limits<int>::is_exact;
    BOOST_MATH_STATIC constexpr bool has_infinity = std::numeric_limits<int>::has_infinity;
    BOOST_MATH_STATIC constexpr bool has_quiet_NaN = std::numeric_limits<int>::has_quiet_NaN;
    BOOST_MATH_STATIC constexpr bool has_signaling_NaN = std::numeric_limits<int>::has_signaling_NaN;

    BOOST_MATH_STATIC constexpr std::float_round_style round_style = std::numeric_limits<int>::round_style;
    BOOST_MATH_STATIC constexpr bool is_iec559 = std::numeric_limits<int>::is_iec559;
    BOOST_MATH_STATIC constexpr bool is_bounded = std::numeric_limits<int>::is_bounded;
    BOOST_MATH_STATIC constexpr bool is_modulo = std::numeric_limits<int>::is_modulo;
    BOOST_MATH_STATIC constexpr int digits = std::numeric_limits<int>::digits;
    BOOST_MATH_STATIC constexpr int digits10 = std::numeric_limits<int>::digits10;
    BOOST_MATH_STATIC constexpr int max_digits10 = std::numeric_limits<int>::max_digits10;
    BOOST_MATH_STATIC constexpr int radix = std::numeric_limits<int>::radix;
    BOOST_MATH_STATIC constexpr int min_exponent = std::numeric_limits<int>::min_exponent;
    BOOST_MATH_STATIC constexpr int min_exponent10 = std::numeric_limits<int>::min_exponent10;
    BOOST_MATH_STATIC constexpr int max_exponent = std::numeric_limits<int>::max_exponent;
    BOOST_MATH_STATIC constexpr int max_exponent10 = std::numeric_limits<int>::max_exponent10;
    BOOST_MATH_STATIC constexpr bool traps = std::numeric_limits<int>::traps;
    BOOST_MATH_STATIC constexpr bool tinyness_before = std::numeric_limits<int>::tinyness_before;

    // Member Functions
    BOOST_MATH_GPU_ENABLED BOOST_MATH_STATIC constexpr int (min)         () { return INT_MIN; }
    BOOST_MATH_GPU_ENABLED BOOST_MATH_STATIC constexpr int (max)         () { return INT_MAX; }
    BOOST_MATH_GPU_ENABLED BOOST_MATH_STATIC constexpr int lowest        () { return INT_MIN; }
    BOOST_MATH_GPU_ENABLED BOOST_MATH_STATIC constexpr int epsilon       () { return 0; }
    BOOST_MATH_GPU_ENABLED BOOST_MATH_STATIC constexpr int round_error   () { return 0; }
    BOOST_MATH_GPU_ENABLED BOOST_MATH_STATIC constexpr int infinity      () { return 0; }
    BOOST_MATH_GPU_ENABLED BOOST_MATH_STATIC constexpr int quiet_NaN     () { return 0; }
    BOOST_MATH_GPU_ENABLED BOOST_MATH_STATIC constexpr int signaling_NaN () { return 0; }
    BOOST_MATH_GPU_ENABLED BOOST_MATH_STATIC constexpr int denorm_min    () { return 0; }
};

template <>
struct numeric_limits<unsigned int>
{
    BOOST_MATH_STATIC constexpr bool is_specialized = std::numeric_limits<unsigned int>::is_specialized;
    BOOST_MATH_STATIC constexpr bool is_signed = std::numeric_limits<unsigned int>::is_signed;
    BOOST_MATH_STATIC constexpr bool is_integer = std::numeric_limits<unsigned int>::is_integer;
    BOOST_MATH_STATIC constexpr bool is_exact = std::numeric_limits<unsigned int>::is_exact;
    BOOST_MATH_STATIC constexpr bool has_infinity = std::numeric_limits<unsigned int>::has_infinity;
    BOOST_MATH_STATIC constexpr bool has_quiet_NaN = std::numeric_limits<unsigned int>::has_quiet_NaN;
    BOOST_MATH_STATIC constexpr bool has_signaling_NaN = std::numeric_limits<unsigned int>::has_signaling_NaN;

    BOOST_MATH_STATIC constexpr std::float_round_style round_style = std::numeric_limits<unsigned int>::round_style;
    BOOST_MATH_STATIC constexpr bool is_iec559 = std::numeric_limits<unsigned int>::is_iec559;
    BOOST_MATH_STATIC constexpr bool is_bounded = std::numeric_limits<unsigned int>::is_bounded;
    BOOST_MATH_STATIC constexpr bool is_modulo = std::numeric_limits<unsigned int>::is_modulo;
    BOOST_MATH_STATIC constexpr int digits = std::numeric_limits<unsigned int>::digits;
    BOOST_MATH_STATIC constexpr int digits10 = std::numeric_limits<unsigned int>::digits10;
    BOOST_MATH_STATIC constexpr int max_digits10 = std::numeric_limits<unsigned int>::max_digits10;
    BOOST_MATH_STATIC constexpr int radix = std::numeric_limits<unsigned int>::radix;
    BOOST_MATH_STATIC constexpr int min_exponent = std::numeric_limits<unsigned int>::min_exponent;
    BOOST_MATH_STATIC constexpr int min_exponent10 = std::numeric_limits<unsigned int>::min_exponent10;
    BOOST_MATH_STATIC constexpr int max_exponent = std::numeric_limits<unsigned int>::max_exponent;
    BOOST_MATH_STATIC constexpr int max_exponent10 = std::numeric_limits<unsigned int>::max_exponent10;
    BOOST_MATH_STATIC constexpr bool traps = std::numeric_limits<unsigned int>::traps;
    BOOST_MATH_STATIC constexpr bool tinyness_before = std::numeric_limits<unsigned int>::tinyness_before;

    // Member Functions
    BOOST_MATH_GPU_ENABLED BOOST_MATH_STATIC constexpr unsigned int (min)         () { return 0; }
    BOOST_MATH_GPU_ENABLED BOOST_MATH_STATIC constexpr unsigned int (max)         () { return UINT_MAX; }
    BOOST_MATH_GPU_ENABLED BOOST_MATH_STATIC constexpr unsigned int lowest        () { return 0; }
    BOOST_MATH_GPU_ENABLED BOOST_MATH_STATIC constexpr unsigned int epsilon       () { return 0; }
    BOOST_MATH_GPU_ENABLED BOOST_MATH_STATIC constexpr unsigned int round_error   () { return 0; }
    BOOST_MATH_GPU_ENABLED BOOST_MATH_STATIC constexpr unsigned int infinity      () { return 0; }
    BOOST_MATH_GPU_ENABLED BOOST_MATH_STATIC constexpr unsigned int quiet_NaN     () { return 0; }
    BOOST_MATH_GPU_ENABLED BOOST_MATH_STATIC constexpr unsigned int signaling_NaN () { return 0; }
    BOOST_MATH_GPU_ENABLED BOOST_MATH_STATIC constexpr unsigned int denorm_min    () { return 0; }
};

template <>
struct numeric_limits<long>
{
    BOOST_MATH_STATIC constexpr bool is_specialized = std::numeric_limits<long>::is_specialized;
    BOOST_MATH_STATIC constexpr bool is_signed = std::numeric_limits<long>::is_signed;
    BOOST_MATH_STATIC constexpr bool is_integer = std::numeric_limits<long>::is_integer;
    BOOST_MATH_STATIC constexpr bool is_exact = std::numeric_limits<long>::is_exact;
    BOOST_MATH_STATIC constexpr bool has_infinity = std::numeric_limits<long>::has_infinity;
    BOOST_MATH_STATIC constexpr bool has_quiet_NaN = std::numeric_limits<long>::has_quiet_NaN;
    BOOST_MATH_STATIC constexpr bool has_signaling_NaN = std::numeric_limits<long>::has_signaling_NaN;

    BOOST_MATH_STATIC constexpr std::float_round_style round_style = std::numeric_limits<long>::round_style;
    BOOST_MATH_STATIC constexpr bool is_iec559 = std::numeric_limits<long>::is_iec559;
    BOOST_MATH_STATIC constexpr bool is_bounded = std::numeric_limits<long>::is_bounded;
    BOOST_MATH_STATIC constexpr bool is_modulo = std::numeric_limits<long>::is_modulo;
    BOOST_MATH_STATIC constexpr int digits = std::numeric_limits<long>::digits;
    BOOST_MATH_STATIC constexpr int digits10 = std::numeric_limits<long>::digits10;
    BOOST_MATH_STATIC constexpr int max_digits10 = std::numeric_limits<long>::max_digits10;
    BOOST_MATH_STATIC constexpr int radix = std::numeric_limits<long>::radix;
    BOOST_MATH_STATIC constexpr int min_exponent = std::numeric_limits<long>::min_exponent;
    BOOST_MATH_STATIC constexpr int min_exponent10 = std::numeric_limits<long>::min_exponent10;
    BOOST_MATH_STATIC constexpr int max_exponent = std::numeric_limits<long>::max_exponent;
    BOOST_MATH_STATIC constexpr int max_exponent10 = std::numeric_limits<long>::max_exponent10;
    BOOST_MATH_STATIC constexpr bool traps = std::numeric_limits<long>::traps;
    BOOST_MATH_STATIC constexpr bool tinyness_before = std::numeric_limits<long>::tinyness_before;

    // Member Functions
    BOOST_MATH_GPU_ENABLED BOOST_MATH_STATIC constexpr long (min)         () { return LONG_MIN; }
    BOOST_MATH_GPU_ENABLED BOOST_MATH_STATIC constexpr long (max)         () { return LONG_MAX; }
    BOOST_MATH_GPU_ENABLED BOOST_MATH_STATIC constexpr long lowest        () { return LONG_MIN; }
    BOOST_MATH_GPU_ENABLED BOOST_MATH_STATIC constexpr long epsilon       () { return 0; }
    BOOST_MATH_GPU_ENABLED BOOST_MATH_STATIC constexpr long round_error   () { return 0; }
    BOOST_MATH_GPU_ENABLED BOOST_MATH_STATIC constexpr long infinity      () { return 0; }
    BOOST_MATH_GPU_ENABLED BOOST_MATH_STATIC constexpr long quiet_NaN     () { return 0; }
    BOOST_MATH_GPU_ENABLED BOOST_MATH_STATIC constexpr long signaling_NaN () { return 0; }
    BOOST_MATH_GPU_ENABLED BOOST_MATH_STATIC constexpr long denorm_min    () { return 0; }
};

template <>
struct numeric_limits<unsigned long>
{
    BOOST_MATH_STATIC constexpr bool is_specialized = std::numeric_limits<unsigned long>::is_specialized;
    BOOST_MATH_STATIC constexpr bool is_signed = std::numeric_limits<unsigned long>::is_signed;
    BOOST_MATH_STATIC constexpr bool is_integer = std::numeric_limits<unsigned long>::is_integer;
    BOOST_MATH_STATIC constexpr bool is_exact = std::numeric_limits<unsigned long>::is_exact;
    BOOST_MATH_STATIC constexpr bool has_infinity = std::numeric_limits<unsigned long>::has_infinity;
    BOOST_MATH_STATIC constexpr bool has_quiet_NaN = std::numeric_limits<unsigned long>::has_quiet_NaN;
    BOOST_MATH_STATIC constexpr bool has_signaling_NaN = std::numeric_limits<unsigned long>::has_signaling_NaN;

    BOOST_MATH_STATIC constexpr std::float_round_style round_style = std::numeric_limits<unsigned long>::round_style;
    BOOST_MATH_STATIC constexpr bool is_iec559 = std::numeric_limits<unsigned long>::is_iec559;
    BOOST_MATH_STATIC constexpr bool is_bounded = std::numeric_limits<unsigned long>::is_bounded;
    BOOST_MATH_STATIC constexpr bool is_modulo = std::numeric_limits<unsigned long>::is_modulo;
    BOOST_MATH_STATIC constexpr int digits = std::numeric_limits<unsigned long>::digits;
    BOOST_MATH_STATIC constexpr int digits10 = std::numeric_limits<unsigned long>::digits10;
    BOOST_MATH_STATIC constexpr int max_digits10 = std::numeric_limits<unsigned long>::max_digits10;
    BOOST_MATH_STATIC constexpr int radix = std::numeric_limits<unsigned long>::radix;
    BOOST_MATH_STATIC constexpr int min_exponent = std::numeric_limits<unsigned long>::min_exponent;
    BOOST_MATH_STATIC constexpr int min_exponent10 = std::numeric_limits<unsigned long>::min_exponent10;
    BOOST_MATH_STATIC constexpr int max_exponent = std::numeric_limits<unsigned long>::max_exponent;
    BOOST_MATH_STATIC constexpr int max_exponent10 = std::numeric_limits<unsigned long>::max_exponent10;
    BOOST_MATH_STATIC constexpr bool traps = std::numeric_limits<unsigned long>::traps;
    BOOST_MATH_STATIC constexpr bool tinyness_before = std::numeric_limits<unsigned long>::tinyness_before;

    // Member Functions
    BOOST_MATH_GPU_ENABLED BOOST_MATH_STATIC constexpr unsigned long (min)         () { return 0; }
    BOOST_MATH_GPU_ENABLED BOOST_MATH_STATIC constexpr unsigned long (max)         () { return ULONG_MAX; }
    BOOST_MATH_GPU_ENABLED BOOST_MATH_STATIC constexpr unsigned long lowest        () { return 0; }
    BOOST_MATH_GPU_ENABLED BOOST_MATH_STATIC constexpr unsigned long epsilon       () { return 0; }
    BOOST_MATH_GPU_ENABLED BOOST_MATH_STATIC constexpr unsigned long round_error   () { return 0; }
    BOOST_MATH_GPU_ENABLED BOOST_MATH_STATIC constexpr unsigned long infinity      () { return 0; }
    BOOST_MATH_GPU_ENABLED BOOST_MATH_STATIC constexpr unsigned long quiet_NaN     () { return 0; }
    BOOST_MATH_GPU_ENABLED BOOST_MATH_STATIC constexpr unsigned long signaling_NaN () { return 0; }
    BOOST_MATH_GPU_ENABLED BOOST_MATH_STATIC constexpr unsigned long denorm_min    () { return 0; }
};

template <>
struct numeric_limits<long long>
{
    BOOST_MATH_STATIC constexpr bool is_specialized = std::numeric_limits<long long>::is_specialized;
    BOOST_MATH_STATIC constexpr bool is_signed = std::numeric_limits<long long>::is_signed;
    BOOST_MATH_STATIC constexpr bool is_integer = std::numeric_limits<long long>::is_integer;
    BOOST_MATH_STATIC constexpr bool is_exact = std::numeric_limits<long long>::is_exact;
    BOOST_MATH_STATIC constexpr bool has_infinity = std::numeric_limits<long long>::has_infinity;
    BOOST_MATH_STATIC constexpr bool has_quiet_NaN = std::numeric_limits<long long>::has_quiet_NaN;
    BOOST_MATH_STATIC constexpr bool has_signaling_NaN = std::numeric_limits<long long>::has_signaling_NaN;

    BOOST_MATH_STATIC constexpr std::float_round_style round_style = std::numeric_limits<long long>::round_style;
    BOOST_MATH_STATIC constexpr bool is_iec559 = std::numeric_limits<long long>::is_iec559;
    BOOST_MATH_STATIC constexpr bool is_bounded = std::numeric_limits<long long>::is_bounded;
    BOOST_MATH_STATIC constexpr bool is_modulo = std::numeric_limits<long long>::is_modulo;
    BOOST_MATH_STATIC constexpr int digits = std::numeric_limits<long long>::digits;
    BOOST_MATH_STATIC constexpr int digits10 = std::numeric_limits<long long>::digits10;
    BOOST_MATH_STATIC constexpr int max_digits10 = std::numeric_limits<long long>::max_digits10;
    BOOST_MATH_STATIC constexpr int radix = std::numeric_limits<long long>::radix;
    BOOST_MATH_STATIC constexpr int min_exponent = std::numeric_limits<long long>::min_exponent;
    BOOST_MATH_STATIC constexpr int min_exponent10 = std::numeric_limits<long long>::min_exponent10;
    BOOST_MATH_STATIC constexpr int max_exponent = std::numeric_limits<long long>::max_exponent;
    BOOST_MATH_STATIC constexpr int max_exponent10 = std::numeric_limits<long long>::max_exponent10;
    BOOST_MATH_STATIC constexpr bool traps = std::numeric_limits<long long>::traps;
    BOOST_MATH_STATIC constexpr bool tinyness_before = std::numeric_limits<long long>::tinyness_before;

    // Member Functions
    BOOST_MATH_GPU_ENABLED BOOST_MATH_STATIC constexpr long long (min)         () { return LLONG_MIN; }
    BOOST_MATH_GPU_ENABLED BOOST_MATH_STATIC constexpr long long (max)         () { return LLONG_MAX; }
    BOOST_MATH_GPU_ENABLED BOOST_MATH_STATIC constexpr long long lowest        () { return LLONG_MIN; }
    BOOST_MATH_GPU_ENABLED BOOST_MATH_STATIC constexpr long long epsilon       () { return 0; }
    BOOST_MATH_GPU_ENABLED BOOST_MATH_STATIC constexpr long long round_error   () { return 0; }
    BOOST_MATH_GPU_ENABLED BOOST_MATH_STATIC constexpr long long infinity      () { return 0; }
    BOOST_MATH_GPU_ENABLED BOOST_MATH_STATIC constexpr long long quiet_NaN     () { return 0; }
    BOOST_MATH_GPU_ENABLED BOOST_MATH_STATIC constexpr long long signaling_NaN () { return 0; }
    BOOST_MATH_GPU_ENABLED BOOST_MATH_STATIC constexpr long long denorm_min    () { return 0; }
};

template <>
struct numeric_limits<unsigned long long>
{
    BOOST_MATH_STATIC constexpr bool is_specialized = std::numeric_limits<unsigned long long>::is_specialized;
    BOOST_MATH_STATIC constexpr bool is_signed = std::numeric_limits<unsigned long long>::is_signed;
    BOOST_MATH_STATIC constexpr bool is_integer = std::numeric_limits<unsigned long long>::is_integer;
    BOOST_MATH_STATIC constexpr bool is_exact = std::numeric_limits<unsigned long long>::is_exact;
    BOOST_MATH_STATIC constexpr bool has_infinity = std::numeric_limits<unsigned long long>::has_infinity;
    BOOST_MATH_STATIC constexpr bool has_quiet_NaN = std::numeric_limits<unsigned long long>::has_quiet_NaN;
    BOOST_MATH_STATIC constexpr bool has_signaling_NaN = std::numeric_limits<unsigned long long>::has_signaling_NaN;

    BOOST_MATH_STATIC constexpr std::float_round_style round_style = std::numeric_limits<unsigned long long>::round_style;
    BOOST_MATH_STATIC constexpr bool is_iec559 = std::numeric_limits<unsigned long long>::is_iec559;
    BOOST_MATH_STATIC constexpr bool is_bounded = std::numeric_limits<unsigned long long>::is_bounded;
    BOOST_MATH_STATIC constexpr bool is_modulo = std::numeric_limits<unsigned long long>::is_modulo;
    BOOST_MATH_STATIC constexpr int digits = std::numeric_limits<unsigned long long>::digits;
    BOOST_MATH_STATIC constexpr int digits10 = std::numeric_limits<unsigned long long>::digits10;
    BOOST_MATH_STATIC constexpr int max_digits10 = std::numeric_limits<unsigned long long>::max_digits10;
    BOOST_MATH_STATIC constexpr int radix = std::numeric_limits<unsigned long long>::radix;
    BOOST_MATH_STATIC constexpr int min_exponent = std::numeric_limits<unsigned long long>::min_exponent;
    BOOST_MATH_STATIC constexpr int min_exponent10 = std::numeric_limits<unsigned long long>::min_exponent10;
    BOOST_MATH_STATIC constexpr int max_exponent = std::numeric_limits<unsigned long long>::max_exponent;
    BOOST_MATH_STATIC constexpr int max_exponent10 = std::numeric_limits<unsigned long long>::max_exponent10;
    BOOST_MATH_STATIC constexpr bool traps = std::numeric_limits<unsigned long long>::traps;
    BOOST_MATH_STATIC constexpr bool tinyness_before = std::numeric_limits<unsigned long long>::tinyness_before;

    // Member Functions
    BOOST_MATH_GPU_ENABLED BOOST_MATH_STATIC constexpr unsigned long long (min)         () { return 0; }
    BOOST_MATH_GPU_ENABLED BOOST_MATH_STATIC constexpr unsigned long long (max)         () { return ULLONG_MAX; }
    BOOST_MATH_GPU_ENABLED BOOST_MATH_STATIC constexpr unsigned long long lowest        () { return 0; }
    BOOST_MATH_GPU_ENABLED BOOST_MATH_STATIC constexpr unsigned long long epsilon       () { return 0; }
    BOOST_MATH_GPU_ENABLED BOOST_MATH_STATIC constexpr unsigned long long round_error   () { return 0; }
    BOOST_MATH_GPU_ENABLED BOOST_MATH_STATIC constexpr unsigned long long infinity      () { return 0; }
    BOOST_MATH_GPU_ENABLED BOOST_MATH_STATIC constexpr unsigned long long quiet_NaN     () { return 0; }
    BOOST_MATH_GPU_ENABLED BOOST_MATH_STATIC constexpr unsigned long long signaling_NaN () { return 0; }
    BOOST_MATH_GPU_ENABLED BOOST_MATH_STATIC constexpr unsigned long long denorm_min    () { return 0; }
};

#endif // BOOST_MATH_HAS_GPU_SUPPORT

} // namespace math
} // namespace boost

#endif
