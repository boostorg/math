//  (C) Copyright Matt Borland 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
//  Bit layout of the sieve array and the wheel tables shared by every backend.
//
//  One byte represents 30 consecutive integers: bit k of byte i is the integer
//  segment_low + 30 * i + wheel_offsets[k] with wheel_offsets = {7, 11, 13, 17, 19, 23, 29, 31}.
//  segment_low is always a multiple of 30 and the first byte of a range starting at s is the
//  one holding s, i.e. segment_low = 30 * ((s - 7) / 30). Multiples of 2, 3 and 5 have no bit.
//
//  A sieving prime p = 30 * p30 + r advances through its multiples p * q where q runs over the
//  residues coprime to the wheel modulus. The wheel tables give, per (residue class of p,
//  residue class of q), the bit to clear, the distance to the next q, and the byte correction.
//
//  The tables are host objects; the CUDA backend copies what it needs into a kernel argument.

#ifndef BOOST_MATH_SF_DETAIL_PRIME_SIEVE_LAYOUT_HPP
#define BOOST_MATH_SF_DETAIL_PRIME_SIEVE_LAYOUT_HPP

#include <boost/math/tools/config.hpp>
#include <boost/math/tools/bit.hpp>
#include <boost/math/tools/cstdint.hpp>

#ifndef BOOST_MATH_BUILD_MODULE
#include <cstdint>
#include <cstddef>
#include <array>
#endif

// Prime extraction with AVX-512 byte compression (Ice Lake and later, or Zen 4), selected at
// compile time; define BOOST_MATH_PRIME_SIEVE_NO_SIMD to force the portable loop.
#if defined(__AVX512F__) && defined(__AVX512VBMI__) && defined(__AVX512VBMI2__) && !defined(__CUDACC__) && !defined(BOOST_MATH_PRIME_SIEVE_NO_SIMD)
#  define BOOST_MATH_PRIME_SIEVE_AVX512_EXTRACT
#  ifndef BOOST_MATH_BUILD_MODULE
#    include <immintrin.h>
#  endif
#endif

namespace boost::math::detail::prime_sieve {

inline constexpr std::uint8_t wheel_offsets[8] = {7, 11, 13, 17, 19, 23, 29, 31};
inline constexpr std::uint8_t unset_bit[8] = {0xfe, 0xfd, 0xfb, 0xf7, 0xef, 0xdf, 0xbf, 0x7f};

// Bit index of a residue modulo 30 that is coprime to 30, or 0xff otherwise.
// The same mapping orders the residue classes of the sieving primes.
BOOST_MATH_GPU_ENABLED constexpr std::uint8_t bit_of_residue(unsigned r) noexcept
{
    return r == 7 ? 0 : r == 11 ? 1 : r == 13 ? 2 : r == 17 ? 3 : r == 19 ? 4 : r == 23 ? 5 : r == 29 ? 6 : r == 1 ? 7 : 0xff;
}

// Residue modulo 30 of the sieving primes in class c.
BOOST_MATH_GPU_ENABLED constexpr unsigned class_residue(unsigned c) noexcept
{
    return c == 7 ? 1u : wheel_offsets[c];
}

namespace layout_detail {

constexpr std::array<std::uint8_t, 65> make_bit_values() noexcept
{
    std::array<std::uint8_t, 65> a {};
    for (std::size_t b {0}; b < 64; ++b)
    {
        a[b] = static_cast<std::uint8_t>(30 * (b / 8) + wheel_offsets[b % 8]);
    }
    a[64] = 0;
    return a;
}

// keep_low[d]: bits of the first byte whose offset is >= d (d = start - segment_low in [7, 36])
constexpr std::array<std::uint8_t, 37> make_keep_low() noexcept
{
    std::array<std::uint8_t, 37> a {};
    for (unsigned d {0}; d < 37; ++d)
    {
        std::uint8_t m {0};
        for (unsigned k {0}; k < 8; ++k)
        {
            if (wheel_offsets[k] >= d)
            {
                m = static_cast<std::uint8_t>(m | (1u << k));
            }
        }
        a[d] = m;
    }
    return a;
}

// keep_high[e]: bits of the last byte whose offset is <= e (e = stop - byte_base in [0, 36])
constexpr std::array<std::uint8_t, 37> make_keep_high() noexcept
{
    std::array<std::uint8_t, 37> a {};
    for (unsigned e {0}; e < 37; ++e)
    {
        std::uint8_t m {0};
        for (unsigned k {0}; k < 8; ++k)
        {
            if (wheel_offsets[k] <= e)
            {
                m = static_cast<std::uint8_t>(m | (1u << k));
            }
        }
        a[e] = m;
    }
    return a;
}

constexpr unsigned gcd_u(unsigned a, unsigned b) noexcept
{
    while (b != 0)
    {
        const unsigned t {a % b};
        a = b;
        b = t;
    }
    return a;
}

} // namespace layout_detail

// Offset within a 64-bit word (240 integers) of each bit; entry 64 absorbs countr_zero(0).
inline constexpr std::array<std::uint8_t, 65> bit_values = layout_detail::make_bit_values();
inline constexpr std::array<std::uint8_t, 37> keep_low = layout_detail::make_keep_low();
inline constexpr std::array<std::uint8_t, 37> keep_high = layout_detail::make_keep_high();

// Wheel tables for a modulus with Size residues coprime to it (30/8 or 210/48).
template <unsigned Modulo, unsigned Size>
struct wheel_tables
{
    struct init_entry
    {
        std::uint8_t next_multiple_factor;  // distance from q to the next coprime multiplier
        std::uint8_t wheel_index;           // phase of that multiplier
    };

    struct element
    {
        std::uint8_t unset_bit;             // mask clearing the bit of p * q
        std::uint8_t next_multiple_factor;  // multiplier increment to the next coprime q
        std::uint8_t correct;               // byte correction beyond factor * (p / 30)
        std::uint16_t next;                 // following state
    };

    static constexpr unsigned modulo {Modulo};
    static constexpr unsigned size {Size};
    static constexpr unsigned states {8 * Size};

    std::array<std::uint8_t, Size> residues {};
    std::array<init_entry, Modulo> init {};
    std::array<element, 8 * Size> wheel {};
};

template <unsigned Modulo, unsigned Size>
constexpr wheel_tables<Modulo, Size> make_wheel_tables() noexcept
{
    wheel_tables<Modulo, Size> t {};

    unsigned n {0};
    for (unsigned r {0}; r < Modulo; ++r)
    {
        if (layout_detail::gcd_u(r, Modulo) == 1)
        {
            t.residues[n++] = static_cast<std::uint8_t>(r);
        }
    }

    for (unsigned x {0}; x < Modulo; ++x)
    {
        unsigned y {x};
        unsigned distance {0};
        while (layout_detail::gcd_u(y % Modulo, Modulo) != 1)
        {
            ++y;
            ++distance;
        }
        unsigned phase {0};
        while (t.residues[phase] != y % Modulo)
        {
            ++phase;
        }
        t.init[x].next_multiple_factor = static_cast<std::uint8_t>(distance);
        t.init[x].wheel_index = static_cast<std::uint8_t>(phase);
    }

    for (unsigned c {0}; c < 8; ++c)
    {
        const unsigned rp {class_residue(c)};
        for (unsigned i {0}; i < Size; ++i)
        {
            const unsigned rq {t.residues[i]};
            const unsigned product {(rp * rq) % 30};
            const unsigned bit {bit_of_residue(product)};
            const unsigned i_next {(i + 1) % Size};
            const unsigned factor {(t.residues[i_next] + Modulo - t.residues[i]) % Modulo};
            const unsigned correct {(wheel_offsets[bit] - 6u + factor * rp) / 30u};

            auto& e = t.wheel[c * Size + i];
            e.unset_bit = unset_bit[bit];
            e.next_multiple_factor = static_cast<std::uint8_t>(factor);
            e.correct = static_cast<std::uint8_t>(correct);
            e.next = static_cast<std::uint16_t>(c * Size + i_next);
        }
    }
    return t;
}

using wheel30_t = wheel_tables<30, 8>;
using wheel210_t = wheel_tables<210, 48>;

inline constexpr wheel30_t wheel30 = make_wheel_tables<30, 8>();
inline constexpr wheel210_t wheel210 = make_wheel_tables<210, 48>();

// Packed sieving prime: 23 bits of byte index, 9 bits of wheel state, and p / 30.
struct sieving_prime
{
    static constexpr std::uint32_t max_multiple_index {(1u << 23) - 1};

    std::uint32_t indexes;
    std::uint32_t prime30;

    constexpr std::size_t multiple_index() const noexcept
    {
        return indexes & max_multiple_index;
    }

    constexpr unsigned wheel_index() const noexcept
    {
        return indexes >> 23;
    }

    constexpr void set(std::size_t multiple_index, unsigned wheel_index) noexcept
    {
        indexes = static_cast<std::uint32_t>(multiple_index) | (static_cast<std::uint32_t>(wheel_index) << 23);
    }
};

// Computes the byte index (relative to segment_low) and wheel state of the first multiple of
// prime that is at least max(prime * prime, segment_low + 7) and coprime to the wheel modulus.
// Returns false when no such multiple is <= stop; the outputs are then meaningless.
// Host only: device code uses its own copy of the tables (see cuda.hpp).
template <class Tables>
inline bool first_multiple(const Tables& t, std::uint64_t prime, std::uint64_t segment_low,
                                                  std::uint64_t stop, std::uint64_t& multiple_index, unsigned& wheel_index) noexcept
{
    // Offsets 7..31 shifted by 6 divide cleanly into bytes 0..
    // Written without early returns: for short windows at large magnitudes about half of the
    // primes have no multiple in range and a data-dependent branch here would mispredict.
    const std::uint64_t low6 {segment_low + 6};
    std::uint64_t q {low6 / prime + 1};
    q = q < prime ? prime : q;
    // prime * q lies in (low6, low6 + prime] or equals prime * prime (< 2^64), so a wrapped
    // product is below low6: one comparison replaces a second 64-bit division
    const std::uint64_t multiple {prime * q};
    const bool in_range {static_cast<bool>(static_cast<int>(multiple >= low6) & static_cast<int>(multiple <= stop))};
    const auto e = t.init[q % Tables::modulo];
    const std::uint64_t advance {prime * e.next_multiple_factor};
    const bool fits {static_cast<bool>(static_cast<int>(in_range) & static_cast<int>(advance <= stop - multiple))};
    multiple_index = (multiple + advance - low6) / 30;
    wheel_index = bit_of_residue(static_cast<unsigned>(prime % 30)) * Tables::size + e.wheel_index;
    return fits;
}

#ifdef BOOST_MATH_PRIME_SIEVE_AVX512_EXTRACT

// Writes the integers whose bits are set in word (base value low) to out and returns how many.
// The 64 byte offsets are compressed by the bit mask so the set bits' offsets become the first
// count bytes; each group of eight is then spread into 64-bit lanes and added to the base.
// Writes whole groups of eight, so out needs 7 slots of slack beyond the count.
inline std::size_t extract_word(std::uint64_t word, std::uint64_t low, std::uint64_t* out) noexcept
{
    const std::size_t count {static_cast<std::size_t>(boost::math::tools::popcount(word))};
    if (count == 0)
    {
        return 0;
    }
    const __m512i offsets {_mm512_loadu_si512(bit_values.data())};
    const __m512i compressed {_mm512_maskz_compress_epi8(static_cast<__mmask64>(word), offsets)};
    const __m512i base {_mm512_set1_epi64(static_cast<long long>(low))};
    __m512i index {_mm512_set_epi64(7, 6, 5, 4, 3, 2, 1, 0)};
    const __m512i eight {_mm512_set1_epi64(8)};
    for (std::size_t j {0}; j < count; j += 8)
    {
        // byte 0 of each 64-bit lane receives compressed[8 * group + lane]
        const __m512i values {_mm512_maskz_permutexvar_epi8(static_cast<__mmask64>(0x0101010101010101ULL), index, compressed)};
        _mm512_storeu_si512(out + j, _mm512_add_epi64(values, base));
        index = _mm512_add_epi64(index, eight);
    }
    return count;
}

// Same for 32-bit outputs: sixteen values per group, so 15 slots of slack are needed.
inline std::size_t extract_word_u32(std::uint64_t word, std::uint64_t low, std::uint32_t* out) noexcept
{
    const std::size_t count {static_cast<std::size_t>(boost::math::tools::popcount(word))};
    if (count == 0)
    {
        return 0;
    }
    const __m512i offsets {_mm512_loadu_si512(bit_values.data())};
    const __m512i compressed {_mm512_maskz_compress_epi8(static_cast<__mmask64>(word), offsets)};
    const __m512i base {_mm512_set1_epi32(static_cast<int>(static_cast<std::uint32_t>(low)))};
    __m512i index {_mm512_set_epi32(15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0)};
    const __m512i sixteen {_mm512_set1_epi32(16)};
    for (std::size_t j {0}; j < count; j += 16)
    {
        const __m512i values {_mm512_maskz_permutexvar_epi8(static_cast<__mmask64>(0x1111111111111111ULL), index, compressed)};
        _mm512_storeu_si512(out + j, _mm512_add_epi32(values, base));
        index = _mm512_add_epi32(index, sixteen);
    }
    return count;
}

#else

// Writes the integers whose bits are set in word (base value low) to out and returns how many.
// Writes four values per iteration unconditionally, so out needs 3 slots of slack beyond
// the count (64 bits at most, so 67 slots always suffice).
inline std::size_t extract_word(std::uint64_t word, std::uint64_t low, std::uint64_t* out) noexcept
{
    const std::size_t count {static_cast<std::size_t>(boost::math::tools::popcount(word))};
    std::size_t j {0};
    while (j < count)
    {
        out[j] = low + bit_values[boost::math::tools::countr_zero(word)];
        word &= word - 1;
        out[j + 1] = low + bit_values[boost::math::tools::countr_zero(word)];
        word &= word - 1;
        out[j + 2] = low + bit_values[boost::math::tools::countr_zero(word)];
        word &= word - 1;
        out[j + 3] = low + bit_values[boost::math::tools::countr_zero(word)];
        word &= word - 1;
        j += 4;
    }
    return count;
}

// Same as extract_word for values known to fit 32 bits (sieving primes).
inline std::size_t extract_word_u32(std::uint64_t word, std::uint64_t low, std::uint32_t* out) noexcept
{
    const std::size_t count {static_cast<std::size_t>(boost::math::tools::popcount(word))};
    const std::uint32_t base {static_cast<std::uint32_t>(low)};
    std::size_t j {0};
    while (j < count)
    {
        out[j] = base + bit_values[boost::math::tools::countr_zero(word)];
        word &= word - 1;
        out[j + 1] = base + bit_values[boost::math::tools::countr_zero(word)];
        word &= word - 1;
        out[j + 2] = base + bit_values[boost::math::tools::countr_zero(word)];
        word &= word - 1;
        out[j + 3] = base + bit_values[boost::math::tools::countr_zero(word)];
        word &= word - 1;
        j += 4;
    }
    return count;
}

#endif // BOOST_MATH_PRIME_SIEVE_AVX512_EXTRACT

} // namespace boost::math::detail::prime_sieve

#endif // BOOST_MATH_SF_DETAIL_PRIME_SIEVE_LAYOUT_HPP
