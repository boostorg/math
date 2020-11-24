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

namespace boost::math::detail::prime_sieve
{
template <typename I = std::uint64_t>
class simple_bitset
{
private:
    std::unique_ptr<I[]> bits;
    std::size_t m_size;
    
public:
    simple_bitset() = delete;
    
    explicit simple_bitset(std::size_t n) noexcept : bits(new I[n / (sizeof(I) * CHAR_BIT) + (n % (sizeof(I) * CHAR_BIT) ? 1 : 0)]), m_size(n)
    {
        std::memset(bits.get(), 0xff, n / CHAR_BIT + (n % CHAR_BIT ? 1 : 0));
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
};
}

#endif // BOOST_MATH_SPECIAL_FUNCTIONS_DETAIL_SIMPLE_BITSET_HPP
