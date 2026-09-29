//  (C) Copyright Matt Borland 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
//  Best-effort detection of the data cache sizes that drive the segment geometry.
//  Falls back to 32 KiB L1 and 1 MiB L2 per core where nothing can be queried.

#ifndef BOOST_MATH_SF_DETAIL_PRIME_SIEVE_CPU_CACHE_HPP
#define BOOST_MATH_SF_DETAIL_PRIME_SIEVE_CPU_CACHE_HPP

#include <boost/math/tools/config.hpp>

#ifndef BOOST_MATH_BUILD_MODULE
#include <cstddef>
#include <cstdint>
#endif

#if defined(__APPLE__)
#  include <sys/types.h>
#  include <sys/sysctl.h>
#elif defined(__linux__)
#  include <unistd.h>
#endif

namespace boost::math::detail::prime_sieve {

struct cache_info
{
    std::size_t l1d_bytes {32768};
    std::size_t l2_bytes {1048576};   // per core (L2 size divided by the cores sharing it)
};

inline cache_info detect_cache_info() noexcept
{
    cache_info info {};

#if defined(__APPLE__)
    std::int64_t value {0};
    std::size_t size {sizeof(value)};
    if (sysctlbyname("hw.perflevel0.l1dcachesize", &value, &size, nullptr, 0) == 0 && value > 0)
    {
        info.l1d_bytes = static_cast<std::size_t>(value);
    }
    else
    {
        value = 0;
        size = sizeof(value);
        if (sysctlbyname("hw.l1dcachesize", &value, &size, nullptr, 0) == 0 && value > 0)
        {
            info.l1d_bytes = static_cast<std::size_t>(value);
        }
    }
    value = 0;
    size = sizeof(value);
    std::int64_t l2 {0};
    if (sysctlbyname("hw.perflevel0.l2cachesize", &value, &size, nullptr, 0) == 0 && value > 0)
    {
        l2 = value;
        std::int64_t sharing {0};
        size = sizeof(sharing);
        if (sysctlbyname("hw.perflevel0.cpusperl2", &sharing, &size, nullptr, 0) == 0 && sharing > 1)
        {
            l2 /= sharing;
        }
    }
    else
    {
        value = 0;
        size = sizeof(value);
        if (sysctlbyname("hw.l2cachesize", &value, &size, nullptr, 0) == 0 && value > 0)
        {
            l2 = value;
        }
    }
    if (l2 > 0)
    {
        info.l2_bytes = static_cast<std::size_t>(l2);
    }
#elif defined(__linux__) && defined(_SC_LEVEL1_DCACHE_SIZE)
    const long l1 {sysconf(_SC_LEVEL1_DCACHE_SIZE)};
    if (l1 > 0)
    {
        info.l1d_bytes = static_cast<std::size_t>(l1);
    }
#  if defined(_SC_LEVEL2_CACHE_SIZE)
    const long l2 {sysconf(_SC_LEVEL2_CACHE_SIZE)};
    if (l2 > 0)
    {
        info.l2_bytes = static_cast<std::size_t>(l2);
    }
#  endif
#endif

    // Guard against nonsense from virtual machines
    if (info.l1d_bytes < 4096 || info.l1d_bytes > (1u << 20))
    {
        info.l1d_bytes = 32768;
    }
    if (info.l2_bytes < info.l1d_bytes || info.l2_bytes > (std::size_t(1) << 30))
    {
        info.l2_bytes = 1048576;
    }
    return info;
}

inline const cache_info& cached_cache_info() noexcept
{
    static const cache_info info {detect_cache_info()};
    return info;
}

} // namespace boost::math::detail::prime_sieve

#endif // BOOST_MATH_SF_DETAIL_PRIME_SIEVE_CPU_CACHE_HPP
