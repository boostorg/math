//  (C) Copyright Matt Borland 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
//  Ranges beyond 2^64 and short intervals at large magnitudes: sieve a window of odd
//  numbers by the primes up to a modest depth, then test each survivor individually
//  (deterministic Miller-Rabin below 2^64, the pseudosquares test while the table
//  covers n / depth, Baillie-PSW beyond). Sorenson, "The pseudosquares prime sieve" (2006).

#ifndef BOOST_MATH_SF_DETAIL_PRIME_SIEVE_BIG_RANGE_HPP
#define BOOST_MATH_SF_DETAIL_PRIME_SIEVE_BIG_RANGE_HPP

#include <boost/math/tools/config.hpp>

#ifndef BOOST_MATH_HAS_NVRTC

#include <boost/math/tools/bit.hpp>
#include <boost/math/special_functions/detail/prime_sieve/options.hpp>
#include <boost/math/special_functions/detail/prime_sieve/integer_traits.hpp>
#include <boost/math/special_functions/detail/prime_sieve/primality.hpp>
#include <boost/math/special_functions/detail/prime_sieve/driver.hpp>

#ifndef BOOST_MATH_BUILD_MODULE
#include <cstdint>
#include <cstddef>
#include <cmath>
#include <vector>
#include <algorithm>
#include <limits>
#include <type_traits>
#endif

namespace boost::math::detail::prime_sieve {

inline constexpr std::uint64_t window_min_depth {100000};
inline constexpr std::uint64_t window_max_depth {100000000};
inline constexpr std::size_t window_max_bits {std::size_t(1) << 24};   // odd candidates per window (2 MiB)

// The short-interval heuristic for 64-bit ranges: when (width * ln(stop))^2 < stop the
// per-candidate tests are cheaper than generating pi(sqrt(stop)) sieving primes.
inline bool prefer_test_path(std::uint64_t start, std::uint64_t stop, const prime_sieve_options& options) noexcept
{
    if (options.range_strategy == prime_range_strategy::full_sieve)
    {
        return false;
    }
    if (options.range_strategy == prime_range_strategy::test_survivors)
    {
        return true;
    }
    const double width {static_cast<double>(stop - start) + 1.0};
    const double t {width * std::log(static_cast<double>(stop))};
    return t * t < static_cast<double>(stop);
}

// Sieving depth s balancing pi(s) modular reductions against the tests on the survivors,
// whose density after sieving to s is about 0.56 / ln(s) (Mertens' theorem).
// test_cost is the cost of one primality test relative to one modular reduction.
inline std::uint64_t choose_window_depth(double window_width, double test_cost, std::uint64_t sqrt_upper) noexcept
{
    double s {window_min_depth};
    for (int i {0}; i < 8; ++i)
    {
        const double ln_s {std::log(s)};
        s = 0.5615 * window_width * test_cost / ln_s;
        s = (std::max)(s, static_cast<double>(window_min_depth));
        s = (std::min)(s, static_cast<double>(window_max_depth));
    }
    std::uint64_t depth {static_cast<std::uint64_t>(s)};
    depth = (std::min)(depth, sqrt_upper);
    depth = (std::max)(depth, std::uint64_t(7));
    return depth;
}

// All odd primes up to depth.
inline std::vector<std::uint32_t> window_sieving_primes(std::uint64_t depth)
{
    std::vector<std::uint32_t> primes;
    for (std::size_t i {1}; i < 71 && small_primes_to_353[i] <= presieve_max_prime && small_primes_to_353[i] <= depth; ++i)
    {
        primes.push_back(small_primes_to_353[i]);
    }
    const std::vector<std::uint32_t> rest {sieving_primes_upto(depth)};
    primes.insert(primes.end(), rest.begin(), rest.end());
    return primes;
}

// Residue of the (odd) window base modulo small primes, specialized per base representation.
struct base_mod_u64
{
    std::uint64_t base;

    std::uint32_t operator()(std::uint32_t p) const noexcept
    {
        return static_cast<std::uint32_t>(base % p);
    }
};

// base = high * 2^64 + low, reduced with 64-bit arithmetic only.
struct base_mod_u128
{
    std::uint64_t high;
    std::uint64_t low;

    std::uint32_t operator()(std::uint32_t p) const noexcept
    {
        const std::uint64_t two64_mod_p {((std::numeric_limits<std::uint64_t>::max)() % p + 1) % p};
        const std::uint64_t r {((high % p) * two64_mod_p + low % p) % p};
        return static_cast<std::uint32_t>(r);
    }
};

template <class Integer>
struct base_mod_generic
{
    const Integer& base;

    std::uint32_t operator()(std::uint32_t p) const
    {
        return static_cast<std::uint32_t>(base % Integer(p));
    }
};

// Marks composites among the odd numbers base + 2 i, i < bits, using the given primes.
// A prime equal to a candidate is left standing.
template <class BaseMod>
void sieve_window(std::vector<std::uint8_t>& composite, std::size_t bits, const BaseMod& base_mod,
                  const std::vector<std::uint32_t>& primes, bool base_is_small)
{
    composite.assign(bits, 0);
    for (const std::uint32_t p : primes)
    {
        const std::uint32_t r {base_mod(p)};
        // first multiple >= base, made odd
        std::uint64_t offset {r == 0 ? 0 : static_cast<std::uint64_t>(p - r)};
        if (offset & 1u)
        {
            offset += p;
        }
        if (base_is_small && offset == 0)
        {
            // base itself could be the prime p; skip to the next odd multiple in that case
            // (the caller only sets base_is_small when the base is at most the depth)
            offset = 2u * p;
        }
        for (std::uint64_t i {offset / 2}; i < bits; i += p)
        {
            composite[static_cast<std::size_t>(i)] = 1;
        }
    }
}

// Sieves [start, stop] (start >= 7, both < 2^64) by the window method and passes batches
// of primes to consume(const std::uint64_t*, std::size_t).
template <class Consumer>
void test_range_u64(std::uint64_t start, std::uint64_t stop, Consumer& consume)
{
    if (stop < start)
    {
        return;
    }
    const std::uint64_t root {boost::math::tools::isqrt(stop)};
    const double width {static_cast<double>(stop - start) + 1.0};
    const std::uint64_t depth {choose_window_depth((std::min)(width, 2.0 * static_cast<double>(window_max_bits)), 100.0, root)};
    const std::vector<std::uint32_t> primes {window_sieving_primes(depth)};
    std::vector<std::uint8_t> composite;
    std::vector<std::uint64_t> batch;
    batch.reserve(4096);

    std::uint64_t base {start | 1u};
    while (base <= stop)
    {
        const std::uint64_t remaining {(stop - base) / 2 + 1};
        const std::size_t bits {static_cast<std::size_t>((std::min)(remaining, static_cast<std::uint64_t>(window_max_bits)))};
        sieve_window(composite, bits, base_mod_u64 {base}, primes, base <= depth);
        for (std::size_t i {0}; i < bits; ++i)
        {
            if (!composite[i])
            {
                const std::uint64_t n {base + 2u * static_cast<std::uint64_t>(i)};
                if (is_prime_u64(n))
                {
                    batch.push_back(n);
                    if (batch.size() == batch.capacity())
                    {
                        consume(batch.data(), batch.size());
                        batch.clear();
                    }
                }
            }
        }
        if (bits == window_max_bits && base <= stop - 2u * static_cast<std::uint64_t>(window_max_bits))
        {
            base += 2u * static_cast<std::uint64_t>(window_max_bits);
        }
        else
        {
            break;
        }
    }
    if (!batch.empty())
    {
        consume(batch.data(), batch.size());
    }
}

// Classifies a survivor above 2^64.
template <class Integer>
bool survivor_is_prime(const Integer& n, std::size_t pss_index, std::uint64_t depth, bool probable_only)
{
    if (!probable_only && pss_index < 49)
    {
        return pseudosquares_prime_test(n, pss_index, depth);
    }
    return is_probable_prime_bpsw(n);
}

// Primes in [lower, upper] (inclusive, lower >= 7, upper >= 2^64) passed one at a time to
// emit(const Integer&).
template <class Integer, class Emit>
void big_range_impl(const Integer& lower, const Integer& upper, const prime_sieve_options& options, Emit&& emit)
{
    if (upper < lower)
    {
        return;
    }
    const Integer two64 {Integer(std::numeric_limits<std::uint64_t>::max()) + Integer(1)};
    const bool fits_128 {upper < two64 * two64};

    const Integer total_width {upper - lower};
    const double width {total_width > Integer(2u * window_max_bits) ? 2.0 * static_cast<double>(window_max_bits) : static_cast<double>(total_width) + 1.0};
    const double test_cost {options.probable_prime_only ? 2000.0 : 40000.0};
    const std::uint64_t depth {choose_window_depth(width, test_cost, window_max_depth)};
    const std::vector<std::uint32_t> primes {window_sieving_primes(depth)};
    std::vector<std::uint8_t> composite;

    Integer base {lower % 2 == 0 ? lower + Integer(1) : lower};
    while (base <= upper)
    {
        const Integer remaining {(upper - base) / 2 + 1};
        const std::size_t bits {remaining > Integer(window_max_bits) ? window_max_bits : static_cast<std::size_t>(static_cast<std::uint64_t>(remaining))};
        const Integer window_high {base + Integer(2u) * Integer(bits - 1)};

        if (fits_128)
        {
            const std::uint64_t high {static_cast<std::uint64_t>((base / two64) % two64)};
            const std::uint64_t low {static_cast<std::uint64_t>(base % two64)};
            sieve_window(composite, bits, base_mod_u128 {high, low}, primes, false);
        }
        else
        {
            sieve_window(composite, bits, base_mod_generic<Integer> {base}, primes, false);
        }

        const Integer n_over_s {window_high / Integer(depth)};
        const std::size_t pss_index {pseudosquare_index_for(n_over_s)};
        for (std::size_t i {0}; i < bits; ++i)
        {
            if (!composite[i])
            {
                const Integer n {base + Integer(2u) * Integer(i)};
                if (fits_u64(n))
                {
                    if (is_prime_u64(to_u64(n)))
                    {
                        emit(n);
                    }
                }
                else if (survivor_is_prime(n, pss_index, depth, options.probable_prime_only))
                {
                    emit(n);
                }
            }
        }
        if (bits < window_max_bits)
        {
            break;
        }
        base += Integer(2u) * Integer(window_max_bits);
    }
}

// Primes in [lower, upper] (inclusive) for bounds that do not fit 64 bits. Builtin types
// can never reach this point, so they get an empty instantiation.
template <class Integer, class OutputIterator>
OutputIterator big_range(const Integer& lower, const Integer& upper, OutputIterator out, const prime_sieve_options& options)
{
    if constexpr (std::is_integral<Integer>::value)
    {
        (void)lower;
        (void)upper;
        (void)options;
        return out;
    }
    else
    {
        big_range_impl(lower, upper, options, [&](const Integer& p) { *out = p; ++out; });
        return out;
    }
}

template <class Integer>
std::uint64_t big_count(const Integer& lower, const Integer& upper, const prime_sieve_options& options)
{
    if constexpr (std::is_integral<Integer>::value)
    {
        (void)lower;
        (void)upper;
        (void)options;
        return 0;
    }
    else
    {
        std::uint64_t count {0};
        big_range_impl(lower, upper, options, [&](const Integer&) { ++count; });
        return count;
    }
}

} // namespace boost::math::detail::prime_sieve

#endif // BOOST_MATH_HAS_NVRTC
#endif // BOOST_MATH_SF_DETAIL_PRIME_SIEVE_BIG_RANGE_HPP
