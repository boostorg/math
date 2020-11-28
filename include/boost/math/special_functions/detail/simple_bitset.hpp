// Copyright 2020 John Maddock and Matt Borland
//
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_SPECIAL_FUNCTIONS_DETAIL_SIMPLE_BITSET_HPP
#define BOOST_MATH_SPECIAL_FUNCTIONS_DETAIL_SIMPLE_BITSET_HPP

#include <array>
#include <memory>
#include <type_traits>
#include <cmath>
#include <cstring>
#include <cstdint>
#include <climits>
#include <cstddef>

namespace boost::math::detail::prime_sieve
{
template <typename I = std::uint64_t>
class simple_bitset
{
private:   
    static constexpr std::size_t ln2(std::size_t n) noexcept
    {
        return n <= 1 ? 0 : 1 + ln2(n >> 1);
    }
    
    // https://www.chessprogramming.org/BitScan#De_Bruijn_Multiplication
    static constexpr std::array<std::uint_fast8_t, 64> index64
    {
         0, 47,  1, 56, 48, 27,  2, 60,
        57, 49, 41, 37, 28, 16,  3, 61,
        54, 58, 35, 52, 50, 42, 21, 44,
        38, 32, 29, 23, 17, 11,  4, 62,
        46, 55, 26, 59, 40, 36, 15, 53,
        34, 51, 20, 43, 31, 22, 10, 45,
        25, 39, 14, 33, 19, 30,  9, 24,
        13, 18,  8, 12,  7,  6,  5, 63
    };

    static constexpr std::uint64_t debruijn64 {0x03f79d71b4cb0a89};
    
    static constexpr auto num_bits {sizeof(I) * CHAR_BIT};
    static constexpr I mask {num_bits - 1};
    static constexpr std::size_t shift {ln2(num_bits)};
    
    std::unique_ptr<I[]> bits;
    std::size_t m_size;
    
public:
    simple_bitset() noexcept : bits {nullptr}, m_size {} {};
    
    explicit simple_bitset(std::size_t n) noexcept : bits {nullptr}, m_size {n}
    {
        resize(n);
        reset();
    }

    // Initialize with specific pattern:
    simple_bitset(std::size_t n, const I* pattern, std::size_t len) noexcept
    {
        const std::size_t block_count = n / num_bits + (n % num_bits ? 1 : 0);
        if (block_count <= len)
        {
            std::memcpy(bits.get(), pattern, block_count * sizeof(I));
        }
        else
        {
            I* p = bits.get();
            std::memcpy(p, pattern, len * sizeof(I));
            I* base = p;
            p += len;
            
            while (len <= block_count - len)
            {
                std::memcpy(p, base, len * sizeof(I));
                p += len;
                len *= 2;
            }
            
            if (block_count > len)
            {
                std::memcpy(p, base, (block_count - len) * sizeof(I));
            }
        }
    }

    ~simple_bitset() = default;
    
    inline I* limbs() noexcept
    {
        return bits.get();
    }

    inline I test(std::size_t n) const noexcept
    {
        return bits[n >> shift] & (I(1u) << (n & mask));
    }

    inline void clear(std::size_t n) noexcept
    {
        bits[n >> shift] &= ~I(I(1u) << (n & mask));
    }

    inline void set(std::size_t n) noexcept
    {
        bits[n >> shift] |= (I(1u) << (n & mask));
    }

    inline std::size_t size() const noexcept { return m_size; }

    void resize(std::size_t n) noexcept
    {
        if(n != m_size)
        {
            m_size = n;
            bits.reset(new I[n / num_bits + (n % num_bits ? 1 : 0)]);
        }
    }

    void reset() noexcept 
    { 
        std::memset(bits.get(), 0xff, m_size / CHAR_BIT + (m_size % CHAR_BIT ? 1 : 0)); 
    }
        
    void clear_all() noexcept
    { 
        std::memset(bits.get(), 0, m_size / CHAR_BIT + (m_size % CHAR_BIT ? 1 : 0)); 
    }

    inline I operator[](const std::size_t n) const noexcept
    {
        return test(n);
    }

    // https://bisqwit.iki.fi/source/misc/bitcounting/
    // WP3 - Uses hardcoded constants if type is U64
    I count() const noexcept
    {
        static constexpr I m1 {std::is_same_v<I, std::uint64_t> ? 0x5555555555555555 : (~static_cast<I>(0)) / 3};
        static constexpr I m2 {std::is_same_v<I, std::uint64_t> ? 0x3333333333333333 : (~static_cast<I>(0)) / 5};
        static constexpr I m4 {std::is_same_v<I, std::uint64_t> ? 0x0f0f0f0f0f0f0f0f : (~static_cast<I>(0)) / 17};
        static constexpr I h1 {std::is_same_v<I, std::uint64_t> ? 0x0101010101010101 : (~static_cast<I>(0)) / 255};

        I counter {};

        for(std::size_t i {}; i < m_size; i += num_bits)
        {
            I x = bits[std::floor(i / num_bits)];
            x -= (x >> 1) & m1;
            x = (x & m2) + ((x >> 2) & m2);
            x = (x + (x >> 4)) & m4;

            counter += (x * h1) >> 56;
        }

        return counter;
    }

    template<std::enable_if_t<std::is_same_v<I, std::uint64_t>, bool> = true>
    std::size_t bit_scan_forward(std::size_t pos) const noexcept
    {
        pos = std::ceil(pos / 64.0);
        
        while(bits[pos] == 0 && pos < m_size)
        {
            ++pos;
        }
        
        const std::uint64_t temp = bits[pos];
        if(temp == 0)
        {
            return m_size;
        }
        else
        {
            return index64[((temp ^ (temp-1)) * debruijn64) >> 58] + pos * 64;
        }
    }

    template<std::enable_if_t<std::is_same_v<I, std::uint64_t>, bool> = true>
    std::size_t bit_scan_reverse(std::size_t pos) const noexcept
    {
        pos = std::floor(pos / 64.0) - 1;
        
        while(bits[pos] == 0 && pos >= 0)
        {
            --pos;
        }

        std::uint64_t temp = bits[pos];
        if(temp == 0)
        {
            return 0;
        }
        else
        {
            temp |= temp >> 1;
            temp |= temp >> 1; 
            temp |= temp >> 2;
            temp |= temp >> 4;
            temp |= temp >> 8;
            temp |= temp >> 16;
            temp |= temp >> 32;
            
            return index64[(temp * debruijn64) >> 58] + pos * 64;
        }
    }
};
}

#endif // BOOST_MATH_SPECIAL_FUNCTIONS_DETAIL_SIMPLE_BITSET_HPP
