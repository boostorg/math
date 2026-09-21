//  (C) Copyright Matt Borland 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
//  User-facing tuning options and the derived segment geometry.

#ifndef BOOST_MATH_SF_DETAIL_PRIME_SIEVE_OPTIONS_HPP
#define BOOST_MATH_SF_DETAIL_PRIME_SIEVE_OPTIONS_HPP

#include <boost/math/tools/config.hpp>
#include <boost/math/tools/bit.hpp>
#include <boost/math/special_functions/detail/prime_sieve/cpu_cache.hpp>

#ifndef BOOST_MATH_BUILD_MODULE
#include <cstddef>
#include <cstdint>
#include <algorithm>
#endif

namespace boost::math {

// Strategy for ranges that fit in 64 bits: sieve to sqrt(upper) or test survivors of a
// shallow sieve individually (better for short intervals at large magnitudes).
BOOST_MATH_EXPORT enum class prime_range_strategy
{
    automatic,
    full_sieve,
    test_survivors
};

BOOST_MATH_EXPORT struct prime_sieve_options
{
    std::size_t l1d_bytes {0};          // 0 = detect
    std::size_t l2_bytes {0};           // 0 = detect (per core)
    std::size_t sieve_bytes {0};        // 0 = derive from sqrt(upper) and the caches; clamped to [16 KiB, 8 MiB]
    unsigned max_threads {0};           // 0 = std::thread::hardware_concurrency()
    std::size_t chunk_primes {1u << 19};  // target primes per parallel chunk when storing output
    prime_range_strategy range_strategy {prime_range_strategy::automatic};
    bool probable_prime_only {false};   // beyond 2^64: skip the deterministic pseudosquares test
};

namespace detail::prime_sieve {

inline constexpr std::size_t min_sieve_bytes {16u * 1024u};
inline constexpr std::size_t max_sieve_bytes {8u * 1024u * 1024u};   // 23-bit multiple index
inline constexpr double factor_small {0.2};
inline constexpr double factor_medium {3.0};
inline constexpr double factor_sieve_size {2.0};

struct sieve_geometry
{
    std::uint64_t start {7};      // >= 7
    std::uint64_t stop {7};       // inclusive
    std::uint64_t sqrt_stop {2};
    std::size_t sieve_bytes {min_sieve_bytes};   // multiple of 8, power of two when big primes exist
    std::size_t l1_bytes {32768};                // erat_small chunk
    std::uint64_t max_small {0};  // primes <= max_small go to erat_small
    std::uint64_t max_medium {0}; // primes <= max_medium go to erat_medium, larger to erat_big
    bool has_big {false};
};

inline std::size_t round_down_multiple(std::size_t value, std::size_t unit) noexcept
{
    return unit == 0 ? value : (value / unit) * unit;
}

inline sieve_geometry make_geometry(std::uint64_t start, std::uint64_t stop, const prime_sieve_options& options) noexcept
{
    sieve_geometry g {};
    g.start = start < 7 ? 7 : start;
    g.stop = stop;
    g.sqrt_stop = boost::math::tools::isqrt(stop);

    const cache_info& cache {cached_cache_info()};
    const std::size_t l1 {options.l1d_bytes != 0 ? options.l1d_bytes : cache.l1d_bytes};
    const std::size_t l2 {options.l2_bytes != 0 ? options.l2_bytes : cache.l2_bytes};

    std::size_t sieve_bytes {};
    if (options.sieve_bytes != 0)
    {
        sieve_bytes = options.sieve_bytes;
    }
    else
    {
        // Scale with sqrt(stop) so that medium primes have several multiples per segment
        const double wanted {factor_sieve_size * static_cast<double>(g.sqrt_stop)};
        sieve_bytes = wanted > static_cast<double>(max_sieve_bytes) ? max_sieve_bytes : static_cast<std::size_t>(wanted);
        if (sieve_bytes > l1)
        {
            sieve_bytes = round_down_multiple(sieve_bytes, l1);
        }
        // Keep the segment inside the per-core L2 and within 16 L1 sizes
        std::size_t cap {static_cast<std::size_t>(boost::math::tools::floor_pow2(l2))};
        cap = (std::min)(cap, 16 * l1);
        cap = (std::max)(cap, l1);
        sieve_bytes = (std::max)(sieve_bytes, l1);
        sieve_bytes = (std::min)(sieve_bytes, cap);
    }
    sieve_bytes = (std::max)(sieve_bytes, min_sieve_bytes);
    sieve_bytes = (std::min)(sieve_bytes, max_sieve_bytes);

    // A range that fits in one segment does not need more memory than that, but the segment
    // must stay large enough that no prime falls into the big class if it can be avoided
    if (stop - g.start < 30u * static_cast<std::uint64_t>(sieve_bytes))
    {
        std::size_t shrunk {static_cast<std::size_t>((stop - g.start) / 30 + 2)};
        const std::size_t no_big {static_cast<std::size_t>(static_cast<double>(g.sqrt_stop) / factor_medium) + 1};
        shrunk = (std::max)(shrunk, no_big);
        sieve_bytes = (std::min)(sieve_bytes, shrunk);
        sieve_bytes = (std::max)(sieve_bytes, std::size_t(8));
    }
    sieve_bytes = (sieve_bytes + 7) / 8 * 8;

    // Big primes need a power-of-two segment for cheap bucket indexing. The bucket sieve
    // is not cache bound, and the largest segment measured fastest, so use the maximum
    // unless the caller fixed the size.
    if (static_cast<double>(g.sqrt_stop) > factor_medium * static_cast<double>(sieve_bytes))
    {
        sieve_bytes = options.sieve_bytes != 0 ? static_cast<std::size_t>(boost::math::tools::floor_pow2(sieve_bytes)) : max_sieve_bytes;
        sieve_bytes = (std::max)(sieve_bytes, min_sieve_bytes);
    }

    g.sieve_bytes = sieve_bytes;
    g.l1_bytes = (std::min)(l1, sieve_bytes);
    g.max_small = static_cast<std::uint64_t>(factor_small * static_cast<double>(g.l1_bytes));
    g.max_medium = static_cast<std::uint64_t>(factor_medium * static_cast<double>(sieve_bytes));
    g.max_small = (std::min)(g.max_small, g.sqrt_stop);
    g.max_medium = (std::min)(g.max_medium, g.sqrt_stop);
    g.has_big = g.sqrt_stop > g.max_medium;
    return g;
}

} // namespace detail::prime_sieve
} // namespace boost::math

#endif // BOOST_MATH_SF_DETAIL_PRIME_SIEVE_OPTIONS_HPP
