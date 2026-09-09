//  (C) Copyright Matt Borland 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
//  Pre-sieving: multiples of the primes 7..163 are removed by ANDing periodic byte patterns
//  into each segment instead of crossing them off one by one. The 16 patterns (products of
//  two or three primes, 6 to 10 KB each) are generated on first use.

#ifndef BOOST_MATH_SF_DETAIL_PRIME_SIEVE_PRESIEVE_HPP
#define BOOST_MATH_SF_DETAIL_PRIME_SIEVE_PRESIEVE_HPP

#include <boost/math/tools/config.hpp>
#include <boost/math/special_functions/detail/prime_sieve/layout.hpp>

#ifndef BOOST_MATH_BUILD_MODULE
#include <cstdint>
#include <cstddef>
#include <array>
#include <vector>
#include <algorithm>
#endif

namespace boost::math::detail::prime_sieve {

inline constexpr unsigned presieve_max_prime {163};
inline constexpr std::size_t presieve_table_count {16};

// The 35 primes handled by pre-sieving, 7..163.
inline constexpr std::uint8_t presieve_primes[35] =
{
    7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97,
    101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163
};

// Groups whose products give tables of 6 to 10 KB each.
inline constexpr std::uint8_t presieve_groups[presieve_table_count][3] =
{
    {7, 23, 37}, {11, 19, 31}, {13, 17, 29}, {41, 163, 0}, {43, 157, 0}, {47, 151, 0}, {53, 149, 0}, {59, 139, 0},
    {61, 137, 0}, {67, 131, 0}, {71, 127, 0}, {73, 113, 0}, {79, 109, 0}, {83, 107, 0}, {89, 103, 0}, {97, 101, 0}
};

// Crosses the multiples of p off a pattern of period bytes covering the integers
// [0, 30 * period + 1]. Usable at compile time.
constexpr void presieve_mark(std::uint8_t* pattern, std::size_t period, std::uint64_t p) noexcept
{
    constexpr std::uint8_t residues[8] = {1, 7, 11, 13, 17, 19, 23, 29};
    const std::uint64_t limit {30 * static_cast<std::uint64_t>(period) + 1};
    for (std::uint64_t base {0}; ; base += 30)
    {
        for (unsigned k {0}; k < 8; ++k)
        {
            const std::uint64_t n {p * (base + residues[k])};
            if (n > limit)
            {
                return;
            }
            pattern[static_cast<std::size_t>((n - 7) / 30)] = static_cast<std::uint8_t>(pattern[(n - 7) / 30] & unset_bit[bit_of_residue(static_cast<unsigned>(n % 30))]);
        }
    }
}

template <std::size_t Period>
struct presieve_pattern
{
    std::uint8_t bytes[Period];
};

// The pattern for one group of primes; the third prime may be 1 (absent).
template <unsigned A, unsigned B, unsigned C>
constexpr presieve_pattern<A * B * C> make_presieve_pattern() noexcept
{
    presieve_pattern<A * B * C> t {};
    for (std::size_t i {0}; i < A * B * C; ++i)
    {
        t.bytes[i] = 0xff;
    }
    presieve_mark(t.bytes, A * B * C, A);
    presieve_mark(t.bytes, A * B * C, B);
    if (C != 1)
    {
        presieve_mark(t.bytes, A * B * C, C);
    }
    return t;
}

// View of the 16 patterns used by presieve_segment.
struct presieve_tables
{
    const std::uint8_t* table[presieve_table_count];
    std::size_t period[presieve_table_count];
};

// The patterns are constant expressions on compilers whose constexpr step limits allow it.
// MSVC (100000 steps by default) and the nvcc front end fall back to generation at first use.
#if !defined(BOOST_MATH_PRIME_SIEVE_RUNTIME_PRESIEVE) && !defined(_MSC_VER) && !defined(__CUDACC__)
#  define BOOST_MATH_PRIME_SIEVE_CONSTEXPR_PRESIEVE
#endif

#ifdef BOOST_MATH_PRIME_SIEVE_CONSTEXPR_PRESIEVE

namespace presieve_detail {

inline constexpr auto pattern_0 = make_presieve_pattern<7, 23, 37>();
inline constexpr auto pattern_1 = make_presieve_pattern<11, 19, 31>();
inline constexpr auto pattern_2 = make_presieve_pattern<13, 17, 29>();
inline constexpr auto pattern_3 = make_presieve_pattern<41, 163, 1>();
inline constexpr auto pattern_4 = make_presieve_pattern<43, 157, 1>();
inline constexpr auto pattern_5 = make_presieve_pattern<47, 151, 1>();
inline constexpr auto pattern_6 = make_presieve_pattern<53, 149, 1>();
inline constexpr auto pattern_7 = make_presieve_pattern<59, 139, 1>();
inline constexpr auto pattern_8 = make_presieve_pattern<61, 137, 1>();
inline constexpr auto pattern_9 = make_presieve_pattern<67, 131, 1>();
inline constexpr auto pattern_10 = make_presieve_pattern<71, 127, 1>();
inline constexpr auto pattern_11 = make_presieve_pattern<73, 113, 1>();
inline constexpr auto pattern_12 = make_presieve_pattern<79, 109, 1>();
inline constexpr auto pattern_13 = make_presieve_pattern<83, 107, 1>();
inline constexpr auto pattern_14 = make_presieve_pattern<89, 103, 1>();
inline constexpr auto pattern_15 = make_presieve_pattern<97, 101, 1>();

inline constexpr presieve_tables tables
{
    {
        pattern_0.bytes, pattern_1.bytes, pattern_2.bytes, pattern_3.bytes, pattern_4.bytes, pattern_5.bytes, pattern_6.bytes, pattern_7.bytes,
        pattern_8.bytes, pattern_9.bytes, pattern_10.bytes, pattern_11.bytes, pattern_12.bytes, pattern_13.bytes, pattern_14.bytes, pattern_15.bytes
    },
    {
        sizeof(pattern_0.bytes), sizeof(pattern_1.bytes), sizeof(pattern_2.bytes), sizeof(pattern_3.bytes), sizeof(pattern_4.bytes), sizeof(pattern_5.bytes), sizeof(pattern_6.bytes), sizeof(pattern_7.bytes),
        sizeof(pattern_8.bytes), sizeof(pattern_9.bytes), sizeof(pattern_10.bytes), sizeof(pattern_11.bytes), sizeof(pattern_12.bytes), sizeof(pattern_13.bytes), sizeof(pattern_14.bytes), sizeof(pattern_15.bytes)
    }
};

} // namespace presieve_detail

inline const presieve_tables& get_presieve_tables() noexcept
{
    return presieve_detail::tables;
}

#else

namespace presieve_detail {

struct runtime_tables
{
    std::vector<std::uint8_t> storage[presieve_table_count];
    presieve_tables view {};

    runtime_tables()
    {
        for (std::size_t k {0}; k < presieve_table_count; ++k)
        {
            std::size_t period {1};
            for (const std::uint8_t p : presieve_groups[k])
            {
                if (p != 0)
                {
                    period *= p;
                }
            }
            storage[k].assign(period, 0xff);
            for (const std::uint8_t p : presieve_groups[k])
            {
                if (p != 0)
                {
                    presieve_mark(storage[k].data(), period, p);
                }
            }
            view.table[k] = storage[k].data();
            view.period[k] = period;
        }
    }
};

} // namespace presieve_detail

inline const presieve_tables& get_presieve_tables()
{
    static const presieve_detail::runtime_tables tables {};
    return tables.view;
}

#endif

// Applies four tables to sieve[0, bytes) starting at segment_low; the first pass assigns,
// later passes AND. Runs are cut at table wrap points so the inner loops auto-vectorize.
template <bool Assign>
inline void presieve_pass(std::uint8_t* sieve, std::size_t bytes, std::uint64_t segment_low,
                          const presieve_tables& t, std::size_t first_table)
{
    const std::uint8_t* tab[4] {};
    std::size_t pos[4] {};
    std::size_t per[4] {};
    for (std::size_t k {0}; k < 4; ++k)
    {
        tab[k] = t.table[first_table + k];
        per[k] = t.period[first_table + k];
        pos[k] = static_cast<std::size_t>((segment_low / 30) % per[k]);
    }

    std::size_t i {0};
    while (i < bytes)
    {
        std::size_t run {bytes - i};
        for (std::size_t k {0}; k < 4; ++k)
        {
            run = (std::min)(run, per[k] - pos[k]);
        }
        const std::uint8_t* t0 {tab[0] + pos[0]};
        const std::uint8_t* t1 {tab[1] + pos[1]};
        const std::uint8_t* t2 {tab[2] + pos[2]};
        const std::uint8_t* t3 {tab[3] + pos[3]};
        std::uint8_t* s {sieve + i};
        if (Assign)
        {
            for (std::size_t j {0}; j < run; ++j)
            {
                s[j] = static_cast<std::uint8_t>(t0[j] & t1[j] & t2[j] & t3[j]);
            }
        }
        else
        {
            for (std::size_t j {0}; j < run; ++j)
            {
                s[j] = static_cast<std::uint8_t>(s[j] & t0[j] & t1[j] & t2[j] & t3[j]);
            }
        }
        i += run;
        for (std::size_t k {0}; k < 4; ++k)
        {
            pos[k] += run;
            if (pos[k] == per[k])
            {
                pos[k] = 0;
            }
        }
    }
}

// Initializes sieve[0, bytes) for the segment at segment_low with multiples of 7..163 removed,
// then restores the bits of those primes themselves when the segment contains them.
inline void presieve_segment(std::uint8_t* sieve, std::size_t bytes, std::uint64_t segment_low)
{
    const presieve_tables& t {get_presieve_tables()};
    presieve_pass<true>(sieve, bytes, segment_low, t, 0);
    presieve_pass<false>(sieve, bytes, segment_low, t, 4);
    presieve_pass<false>(sieve, bytes, segment_low, t, 8);
    presieve_pass<false>(sieve, bytes, segment_low, t, 12);

    if (segment_low <= presieve_max_prime)
    {
        for (const std::uint8_t p : presieve_primes)
        {
            if (p >= segment_low + 7)
            {
                const std::size_t byte {static_cast<std::size_t>((p - segment_low - 7) / 30)};
                if (byte < bytes)
                {
                    sieve[byte] = static_cast<std::uint8_t>(sieve[byte] | (1u << bit_of_residue(p % 30)));
                }
            }
        }
    }
}

} // namespace boost::math::detail::prime_sieve

#endif // BOOST_MATH_SF_DETAIL_PRIME_SIEVE_PRESIEVE_HPP
