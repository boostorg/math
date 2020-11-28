// Copyright 2020 John Maddock and Matt Borland
//
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_SPECIAL_FUNCTIONS_DETAIL_SIMPLE_BITSET_HPP
#define BOOST_MATH_SPECIAL_FUNCTIONS_DETAIL_SIMPLE_BITSET_HPP

#include <memory>
#include <cstring>
#include <cstdint>
#include <climits>
#include <cstddef>
#include <type_traits>

namespace boost::math::detail::prime_sieve
{
template <typename I = std::uint64_t>
class simple_bitset
{
private:
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
        const std::size_t block_count = n / (sizeof(I) * CHAR_BIT) + (n % (sizeof(I) * CHAR_BIT) ? 1 : 0);
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
    
    static constexpr std::size_t ln2(std::size_t n) noexcept
    {
        return n <= 1 ? 0 : 1 + ln2(n >> 1);
    }
    
    I* limbs() noexcept
    {
        return bits.get();
    }

    I test(std::size_t n) const noexcept
    {
        constexpr I mask = (sizeof(I) * CHAR_BIT) - 1;
        constexpr std::size_t shift = ln2(sizeof(I) * CHAR_BIT);
        return bits[n >> shift] & (I(1u) << (n & mask));
    }

    void clear(std::size_t n) noexcept
    {
        constexpr I mask = (sizeof(I) * CHAR_BIT) - 1;
        constexpr std::size_t shift = ln2(sizeof(I) * CHAR_BIT);
        bits[n >> shift] &= ~I(I(1u) << (n & mask));
    }

    void set(std::size_t n) noexcept
    {
        constexpr I mask = (sizeof(I) * CHAR_BIT) - 1;
        constexpr std::size_t shift = ln2(sizeof(I) * CHAR_BIT);
        bits[n >> shift] |= (I(1u) << (n & mask));
    }

    std::size_t size() const noexcept { return m_size; }

    void resize(std::size_t n) noexcept
    {
        if(n != m_size)
        {
            m_size = n;
            bits.reset(new I[n / (sizeof(I) * CHAR_BIT) + (n % (sizeof(I) * CHAR_BIT) ? 1 : 0)]);
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

        for(std::size_t i {}; i < m_size; ++i)
        {
            I x = bits[i];
            x -= (x >> 1) & m1;
            x = (x & m2) + ((x >> 2) & m2);
            x = (x + (x >> 4)) & m4;

            counter += (x * h1) >> 56;
        }

        return counter;
    }
};
}

#endif // BOOST_MATH_SPECIAL_FUNCTIONS_DETAIL_SIMPLE_BITSET_HPP
