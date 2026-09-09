//  (C) Copyright Matt Borland 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
//  Sinks consume finished segments: count set bits, or extract the primes into a
//  64-bit buffer that a consumer converts to the user's integer type.

#ifndef BOOST_MATH_SF_DETAIL_PRIME_SIEVE_SINKS_HPP
#define BOOST_MATH_SF_DETAIL_PRIME_SIEVE_SINKS_HPP

#include <boost/math/tools/config.hpp>
#include <boost/math/tools/bit.hpp>
#include <boost/math/special_functions/detail/prime_sieve/layout.hpp>

#ifndef BOOST_MATH_BUILD_MODULE
#include <cstdint>
#include <cstddef>
#include <vector>
#endif

namespace boost::math::detail::prime_sieve {

// Sink contract: segment(words, n_words, low) receives the final bits of one segment whose
// first byte represents low + {7, ..., 31}; bits outside the requested range are already clear.
struct count_sink
{
    std::uint64_t count {0};

    void segment(const std::uint64_t* words, std::size_t n_words, std::uint64_t) noexcept
    {
        count += count_bits(words, n_words);
    }

    // Without a hardware popcount instruction the compiler expands popcount into a library
    // call or a long bit trick per word; a Harley-Seal carry-save adder over 8 words needs
    // about a third of the operations (Hacker's Delight, 2nd edition, chapter 5).
    static std::uint64_t count_bits(const std::uint64_t* words, std::size_t n) noexcept
    {
#if defined(__x86_64__) && !defined(__POPCNT__)
        std::uint64_t total {0};
        std::uint64_t ones {0};
        std::uint64_t twos {0};
        std::uint64_t fours {0};
        std::size_t i {0};
        for (; i + 8 <= n; i += 8)
        {
            std::uint64_t twos_a {};
            std::uint64_t twos_b {};
            std::uint64_t fours_a {};
            std::uint64_t fours_b {};
            csa(ones, twos_a, ones, words[i], words[i + 1]);
            csa(ones, twos_b, ones, words[i + 2], words[i + 3]);
            csa(twos, fours_a, twos, twos_a, twos_b);
            csa(ones, twos_a, ones, words[i + 4], words[i + 5]);
            csa(ones, twos_b, ones, words[i + 6], words[i + 7]);
            csa(twos, fours_b, twos, twos_a, twos_b);
            std::uint64_t eights {};
            csa(fours, eights, fours, fours_a, fours_b);
            total += static_cast<std::uint64_t>(swar_popcount(eights));
        }
        total = 8 * total + 4 * swar_popcount(fours) + 2 * swar_popcount(twos) + swar_popcount(ones);
        for (; i < n; ++i)
        {
            total += static_cast<std::uint64_t>(swar_popcount(words[i]));
        }
        return total;
#else
        std::uint64_t c {0};
        for (std::size_t i {0}; i < n; ++i)
        {
            c += static_cast<std::uint64_t>(boost::math::tools::popcount(words[i]));
        }
        return c;
#endif
    }

    static void csa(std::uint64_t& sum, std::uint64_t& carry, std::uint64_t a, std::uint64_t b, std::uint64_t c) noexcept
    {
        const std::uint64_t u {a ^ b};
        sum = u ^ c;
        carry = (a & b) | (u & c);
    }

    static std::uint64_t swar_popcount(std::uint64_t x) noexcept
    {
        x = x - ((x >> 1) & 0x5555555555555555ULL);
        x = (x & 0x3333333333333333ULL) + ((x >> 2) & 0x3333333333333333ULL);
        x = (x + (x >> 4)) & 0x0f0f0f0f0f0f0f0fULL;
        return (x * 0x0101010101010101ULL) >> 56;
    }

    void flush() noexcept
    {
    }
};

// Extracts primes into a buffer of batch entries and hands full batches to Consumer,
// which is called as consume(const std::uint64_t* primes, std::size_t n).
template <class Consumer>
class extract_sink
{
public:
    static constexpr std::size_t slack {68};

    explicit extract_sink(Consumer& consume, std::size_t batch = 65536) : consume_(consume), batch_(batch)
    {
        buffer_.resize(batch + slack);
    }

    void segment(const std::uint64_t* words, std::size_t n_words, std::uint64_t low)
    {
        const unsigned char* bytes {reinterpret_cast<const unsigned char*>(words)};
        for (std::size_t w {0}; w < n_words; ++w)
        {
            if (size_ >= batch_)
            {
                flush();
            }
            const std::uint64_t word {boost::math::tools::load_le64(bytes + 8 * w)};
            size_ += extract_word(word, low + 240u * w, buffer_.data() + size_);
        }
    }

    void flush()
    {
        if (size_ != 0)
        {
            consume_(buffer_.data(), size_);
            size_ = 0;
        }
    }

private:
    Consumer& consume_;
    std::vector<std::uint64_t> buffer_;
    std::size_t size_ {0};
    std::size_t batch_;
};

// Consumer that converts to Integer and writes through an output iterator.
template <class Integer, class OutputIterator>
struct output_converter
{
    OutputIterator out;

    void operator()(const std::uint64_t* primes, std::size_t n)
    {
        for (std::size_t i {0}; i < n; ++i)
        {
            *out = static_cast<Integer>(primes[i]);
            ++out;
        }
    }
};

// Writes 32-bit primes straight into a pre-sized vector (sieving prime generation), avoiding
// the 64-bit batch buffer and the conversion pass.
class extract_sink_u32
{
public:
    // capacity must be an upper bound on the number of primes plus 68 slots of slack
    extract_sink_u32(std::vector<std::uint32_t>& out, std::size_t capacity) : out_(out)
    {
        out_.resize(capacity);
    }

    void segment(const std::uint64_t* words, std::size_t n_words, std::uint64_t low) noexcept
    {
        const unsigned char* bytes {reinterpret_cast<const unsigned char*>(words)};
        std::uint32_t* dst {out_.data()};
        for (std::size_t w {0}; w < n_words; ++w)
        {
            const std::uint64_t word {boost::math::tools::load_le64(bytes + 8 * w)};
            size_ += extract_word_u32(word, low + 240u * w, dst + size_);
        }
    }

    void flush()
    {
        out_.resize(size_);
    }

private:
    std::vector<std::uint32_t>& out_;
    std::size_t size_ {0};
};

// Consumer that appends 32-bit values to a vector (sieving primes).
struct append_u32
{
    std::vector<std::uint32_t>& out;

    void operator()(const std::uint64_t* primes, std::size_t n)
    {
        const std::size_t old_size {out.size()};
        out.resize(old_size + n);
        std::uint32_t* dst {out.data() + old_size};
        for (std::size_t i {0}; i < n; ++i)
        {
            dst[i] = static_cast<std::uint32_t>(primes[i]);
        }
    }
};

// Consumer that appends 64-bit values to a vector.
struct append_u64
{
    std::vector<std::uint64_t>& out;

    void operator()(const std::uint64_t* primes, std::size_t n)
    {
        out.insert(out.end(), primes, primes + n);
    }
};

} // namespace boost::math::detail::prime_sieve

#endif // BOOST_MATH_SF_DETAIL_PRIME_SIEVE_SINKS_HPP
