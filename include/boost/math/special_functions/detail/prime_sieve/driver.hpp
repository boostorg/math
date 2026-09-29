//  (C) Copyright Matt Borland 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
//  Sieving prime generation and the single-threaded 64-bit driver.

#ifndef BOOST_MATH_SF_DETAIL_PRIME_SIEVE_DRIVER_HPP
#define BOOST_MATH_SF_DETAIL_PRIME_SIEVE_DRIVER_HPP

#include <boost/math/tools/config.hpp>
#include <boost/math/tools/bit.hpp>
#include <boost/math/special_functions/detail/prime_sieve/engine.hpp>
#include <boost/math/special_functions/detail/prime_sieve/sinks.hpp>

#ifndef BOOST_MATH_BUILD_MODULE
#include <cstdint>
#include <cstddef>
#include <cmath>
#include <vector>
#include <algorithm>
#endif

namespace boost::math::detail::prime_sieve {

// Upper bound on the number of primes <= x (Dusart 2010 for x >= 60184).
inline std::uint64_t prime_count_upper_bound(std::uint64_t x) noexcept
{
    if (x < 2)
    {
        return 0;
    }
    if (x < 60184)
    {
        return static_cast<std::uint64_t>(1.25 * static_cast<double>(x) / std::log(static_cast<double>(x))) + 1;
    }
    const double d {static_cast<double>(x)};
    return static_cast<std::uint64_t>(d / (std::log(d) - 1.1)) + 1;
}

// Upper bound on the number of primes in [lower, upper].
inline std::uint64_t prime_count_upper_bound(std::uint64_t lower, std::uint64_t upper) noexcept
{
    if (upper < lower)
    {
        return 0;
    }
    if (lower < 60184)
    {
        return prime_count_upper_bound(upper);
    }
    // pi(upper) - pi(lower) <= (upper - lower) / (ln lower - 1.1) since density decreases
    const double width {static_cast<double>(upper - lower) + 1.0};
    return static_cast<std::uint64_t>(width / (std::log(static_cast<double>(lower)) - 1.1)) + 1;
}

// Odd-only byte sieve for the primes in [167, n] with n <= 2^32 - 1 (used for n up to 65535
// to seed the segmented generator, but correct for any n that fits memory).
inline std::vector<std::uint32_t> simple_primes_from_167(std::uint32_t n)
{
    std::vector<std::uint32_t> out;
    if (n < 167)
    {
        return out;
    }
    std::vector<std::uint8_t> composite((n / 2) + 1, 0);
    for (std::uint32_t i {3}; static_cast<std::uint64_t>(i) * i <= n; i += 2)
    {
        if (!composite[i / 2])
        {
            for (std::uint64_t j {static_cast<std::uint64_t>(i) * i}; j <= n; j += 2u * i)
            {
                composite[static_cast<std::size_t>(j / 2)] = 1;
            }
        }
    }
    out.reserve(static_cast<std::size_t>(prime_count_upper_bound(n)));
    for (std::uint32_t i {167}; i <= n; i += 2)
    {
        if (!composite[i / 2])
        {
            out.push_back(i);
        }
    }
    return out;
}

// Primes in [167, n] for n < 2^32, produced by the segmented engine itself.
inline std::vector<std::uint32_t> sieving_primes_upto(std::uint64_t n)
{
    std::vector<std::uint32_t> out;
    if (n < 167)
    {
        return out;
    }
    if (n < 100000)
    {
        return simple_primes_from_167(static_cast<std::uint32_t>(n));
    }
    const std::uint32_t root {static_cast<std::uint32_t>(boost::math::tools::isqrt(n))};
    const std::vector<std::uint32_t> seed {simple_primes_from_167(root)};

    prime_sieve_options options {};
    const sieve_geometry g {make_geometry(167, n, options)};
    extract_sink_u32 sink {out, static_cast<std::size_t>(prime_count_upper_bound(n)) + extract_sink<append_u32>::slack};
    segment_sieve engine {g, seed};
    engine.run(sink);
    return out;
}

// Sieves [start, stop] (start >= 7) sequentially into sink.
template <class Sink>
void run_u64(std::uint64_t start, std::uint64_t stop, const prime_sieve_options& options, Sink& sink)
{
    if (stop < start)
    {
        sink.flush();
        return;
    }
    const sieve_geometry g {make_geometry(start, stop, options)};
    const std::vector<std::uint32_t> primes {sieving_primes_upto(g.sqrt_stop)};
    segment_sieve engine {g, primes};
    engine.run(sink);
}

// Counts the primes in [start, stop] (start >= 7) sequentially.
inline std::uint64_t count_u64(std::uint64_t start, std::uint64_t stop, const prime_sieve_options& options)
{
    count_sink sink {};
    run_u64(start, stop, options, sink);
    return sink.count;
}

} // namespace boost::math::detail::prime_sieve

#endif // BOOST_MATH_SF_DETAIL_PRIME_SIEVE_DRIVER_HPP
