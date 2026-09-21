//  (C) Copyright Matt Borland 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
//  Portable bit manipulation helpers usable from host and device code.
//  Everything here works on 64-bit unsigned values and never relies on
//  function-local static tables, so it is safe under CUDA and SYCL.

#ifndef BOOST_MATH_TOOLS_BIT_HPP
#define BOOST_MATH_TOOLS_BIT_HPP

#include <boost/math/tools/config.hpp>
#include <boost/math/tools/cstdint.hpp>

#ifndef BOOST_MATH_HAS_NVRTC
#ifndef BOOST_MATH_BUILD_MODULE
#if defined(__has_include)
#  if __has_include(<version>)
#    include <version>
#  endif
#endif
#if defined(__cpp_lib_bitops) && (__cpp_lib_bitops >= 201907L)
#  include <bit>
#endif
#if defined(_MSC_VER) && !defined(__clang__)
#  include <intrin.h>
#endif
#endif // BOOST_MATH_BUILD_MODULE
#endif // BOOST_MATH_HAS_NVRTC

namespace boost { namespace math { namespace tools {

// Number of trailing zero bits in x. Returns 64 when x is zero.
BOOST_MATH_GPU_ENABLED inline int countr_zero(boost::math::uint64_t x) noexcept
{
#if defined(__CUDA_ARCH__)
    return x == 0 ? 64 : __ffsll(static_cast<long long>(x)) - 1;
#elif defined(__cpp_lib_bitops) && (__cpp_lib_bitops >= 201907L) && !defined(BOOST_MATH_HAS_GPU_SUPPORT)
    return std::countr_zero(static_cast<unsigned long long>(x));
#elif defined(__GNUC__) || defined(__clang__)
    return x == 0 ? 64 : __builtin_ctzll(static_cast<unsigned long long>(x));
#elif defined(_MSC_VER) && (defined(_M_X64) || defined(_M_ARM64))
    unsigned long index {};
    return _BitScanForward64(&index, static_cast<unsigned __int64>(x)) ? static_cast<int>(index) : 64;
#else
    if (x == 0)
    {
        return 64;
    }
    int n {0};
    while ((x & 1u) == 0)
    {
        x >>= 1;
        ++n;
    }
    return n;
#endif
}

// Number of leading zero bits in x. Returns 64 when x is zero.
BOOST_MATH_GPU_ENABLED inline int countl_zero(boost::math::uint64_t x) noexcept
{
#if defined(__CUDA_ARCH__)
    return __clzll(static_cast<long long>(x));
#elif defined(__cpp_lib_bitops) && (__cpp_lib_bitops >= 201907L) && !defined(BOOST_MATH_HAS_GPU_SUPPORT)
    return std::countl_zero(static_cast<unsigned long long>(x));
#elif defined(__GNUC__) || defined(__clang__)
    return x == 0 ? 64 : __builtin_clzll(static_cast<unsigned long long>(x));
#elif defined(_MSC_VER) && (defined(_M_X64) || defined(_M_ARM64))
    unsigned long index {};
    return _BitScanReverse64(&index, static_cast<unsigned __int64>(x)) ? 63 - static_cast<int>(index) : 64;
#else
    if (x == 0)
    {
        return 64;
    }
    int n {0};
    while ((x >> 63) == 0)
    {
        x <<= 1;
        ++n;
    }
    return n;
#endif
}

// Number of set bits in x.
BOOST_MATH_GPU_ENABLED inline int popcount(boost::math::uint64_t x) noexcept
{
#if defined(__CUDA_ARCH__)
    return __popcll(static_cast<unsigned long long>(x));
#elif defined(__cpp_lib_bitops) && (__cpp_lib_bitops >= 201907L) && !defined(BOOST_MATH_HAS_GPU_SUPPORT)
    return std::popcount(static_cast<unsigned long long>(x));
#elif defined(__GNUC__) || defined(__clang__)
    return __builtin_popcountll(static_cast<unsigned long long>(x));
#else
    // Parallel bit count (Hacker's Delight, figure 5-2)
    x = x - ((x >> 1) & 0x5555555555555555ULL);
    x = (x & 0x3333333333333333ULL) + ((x >> 2) & 0x3333333333333333ULL);
    x = (x + (x >> 4)) & 0x0f0f0f0f0f0f0f0fULL;
    return static_cast<int>((x * 0x0101010101010101ULL) >> 56);
#endif
}

// Largest power of two not exceeding x. Returns 0 for x == 0.
BOOST_MATH_GPU_ENABLED inline boost::math::uint64_t floor_pow2(boost::math::uint64_t x) noexcept
{
    return x == 0 ? 0 : boost::math::uint64_t(1) << (63 - countl_zero(x));
}

// Exact integer square root: the largest r with r * r <= n.
BOOST_MATH_GPU_ENABLED inline boost::math::uint64_t isqrt(boost::math::uint64_t n) noexcept
{
    if (n < 2)
    {
        return n;
    }
    // Newton iteration from a power of two above the root converges monotonically down.
    const int bits {64 - countl_zero(n)};
    boost::math::uint64_t x {boost::math::uint64_t(1) << ((bits + 1) / 2)};
    while (true)
    {
        const boost::math::uint64_t y {(x + n / x) / 2};
        if (y >= x)
        {
            return x;
        }
        x = y;
    }
}

// Reads eight bytes as a little-endian 64-bit value regardless of host byte order.
BOOST_MATH_GPU_ENABLED inline boost::math::uint64_t load_le64(const unsigned char* p) noexcept
{
    return  static_cast<boost::math::uint64_t>(p[0])
         | (static_cast<boost::math::uint64_t>(p[1]) << 8)
         | (static_cast<boost::math::uint64_t>(p[2]) << 16)
         | (static_cast<boost::math::uint64_t>(p[3]) << 24)
         | (static_cast<boost::math::uint64_t>(p[4]) << 32)
         | (static_cast<boost::math::uint64_t>(p[5]) << 40)
         | (static_cast<boost::math::uint64_t>(p[6]) << 48)
         | (static_cast<boost::math::uint64_t>(p[7]) << 56);
}

}}} // namespace boost::math::tools

#endif // BOOST_MATH_TOOLS_BIT_HPP
