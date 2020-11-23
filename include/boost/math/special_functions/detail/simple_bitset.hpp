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
    explicit simple_bitset(std::size_t n) : bits(new I[n / (sizeof(I) * CHAR_BIT) + (n % (sizeof(I) * CHAR_BIT) ? 1 : 0)]), m_size(n)
    {
        std::memset(bits.get(), 0xff, n / CHAR_BIT + (n % CHAR_BIT ? 1 : 0));
    }
    ~simple_bitset() = default;
    
    static constexpr std::size_t ln2(std::size_t n) noexcept
    {
        return n <= 1 ? 0 : 1 + ln2(n >> 1);
    }
    
    I test(const std::size_t n) const noexcept
    {
        constexpr I masks[] = { static_cast<I>(1uLL), static_cast<I>(2uLL), static_cast<I>(4uLL), static_cast<I>(8uLL), static_cast<I>(16uLL), 
            static_cast<I>(32uLL), static_cast<I>(64uLL), static_cast<I>(128uLL), static_cast<I>(256uLL),
            static_cast<I>(1uLL << 9),  static_cast<I>(1uLL << 10), static_cast<I>(1uLL << 11), static_cast<I>(1uLL << 12), static_cast<I>(1uLL << 13), 
            static_cast<I>(1uLL << 14), static_cast<I>(1uLL << 15), static_cast<I>(1uLL << 16), static_cast<I>(1uLL << 17), static_cast<I>(1uLL << 18), 
            static_cast<I>(1uLL << 19), static_cast<I>(1uLL << 20), static_cast<I>(1uLL << 21), static_cast<I>(1uLL << 22), static_cast<I>(1uLL << 23), 
            static_cast<I>(1uLL << 24), static_cast<I>(1uLL << 25), static_cast<I>(1uLL << 26), static_cast<I>(1uLL << 27), static_cast<I>(1uLL << 28), 
            static_cast<I>(1uLL << 29), static_cast<I>(1uLL << 30), static_cast<I>(1uLL << 31), static_cast<I>(1uLL << 32), static_cast<I>(1uLL << 33), 
            static_cast<I>(1uLL << 34), static_cast<I>(1uLL << 35), static_cast<I>(1uLL << 36), static_cast<I>(1uLL << 37), static_cast<I>(1uLL << 38), 
            static_cast<I>(1uLL << 39), static_cast<I>(1uLL << 40), static_cast<I>(1uLL << 41), static_cast<I>(1uLL << 42), static_cast<I>(1uLL << 43), 
            static_cast<I>(1uLL << 44), static_cast<I>(1uLL << 45), static_cast<I>(1uLL << 46), static_cast<I>(1uLL << 47), static_cast<I>(1uLL << 48),
            static_cast<I>(1uLL << 49), static_cast<I>(1uLL << 50), static_cast<I>(1uLL << 51), static_cast<I>(1uLL << 52), static_cast<I>(1uLL << 53), 
            static_cast<I>(1uLL << 54), static_cast<I>(1uLL << 55), static_cast<I>(1uLL << 56), static_cast<I>(1uLL << 57), static_cast<I>(1uLL << 58), 
            static_cast<I>(1uLL << 59), static_cast<I>(1uLL << 60), static_cast<I>(1uLL << 61), static_cast<I>(1uLL << 62), static_cast<I>(1uLL << 63)
        };
        constexpr I mask = (sizeof(I) * CHAR_BIT) - 1;
        constexpr std::size_t shift = ln2(sizeof(I) * CHAR_BIT);

        return bits[n >> shift] & masks[n & mask];
    }
    
    void clear(const std::size_t n) noexcept
    {
        constexpr I masks[] = { static_cast<I>(~1uLL), static_cast<I>(~2uLL), static_cast<I>(~4uLL), static_cast<I>(~8uLL), static_cast<I>(~16uLL), 
            static_cast<I>(~32uLL), static_cast<I>(~64uLL), static_cast<I>(~128uLL), static_cast<I>(~256uLL),
            static_cast<I>(~(1uLL << 9)),  static_cast<I>(~(1uLL << 10)), static_cast<I>(~(1uLL << 11)), static_cast<I>(~(1uLL << 12)), 
            static_cast<I>(~(1uLL << 13)), static_cast<I>(~(1uLL << 14)), static_cast<I>(~(1uLL << 15)), static_cast<I>(~(1uLL << 16)),
            static_cast<I>(~(1uLL << 17)), static_cast<I>(~(1uLL << 18)), static_cast<I>(~(1uLL << 19)), static_cast<I>(~(1uLL << 20)), 
            static_cast<I>(~(1uLL << 21)), static_cast<I>(~(1uLL << 22)), static_cast<I>(~(1uLL << 23)), static_cast<I>(~(1uLL << 24)),
            static_cast<I>(~(1uLL << 25)), static_cast<I>(~(1uLL << 26)), static_cast<I>(~(1uLL << 27)), static_cast<I>(~(1uLL << 28)), 
            static_cast<I>(~(1uLL << 29)), static_cast<I>(~(1uLL << 30)), static_cast<I>(~(1uLL << 31)), static_cast<I>(~(1uLL << 32)),
            static_cast<I>(~(1uLL << 33)), static_cast<I>(~(1uLL << 34)), static_cast<I>(~(1uLL << 35)), static_cast<I>(~(1uLL << 36)), 
            static_cast<I>(~(1uLL << 37)), static_cast<I>(~(1uLL << 38)), static_cast<I>(~(1uLL << 39)), static_cast<I>(~(1uLL << 40)),
            static_cast<I>(~(1uLL << 41)), static_cast<I>(~(1uLL << 42)), static_cast<I>(~(1uLL << 43)), static_cast<I>(~(1uLL << 44)), 
            static_cast<I>(~(1uLL << 45)), static_cast<I>(~(1uLL << 46)), static_cast<I>(~(1uLL << 47)), static_cast<I>(~(1uLL << 48)),
            static_cast<I>(~(1uLL << 49)), static_cast<I>(~(1uLL << 50)), static_cast<I>(~(1uLL << 51)), static_cast<I>(~(1uLL << 52)), 
            static_cast<I>(~(1uLL << 53)), static_cast<I>(~(1uLL << 54)), static_cast<I>(~(1uLL << 55)), static_cast<I>(~(1uLL << 56)),
            static_cast<I>(~(1uLL << 57)), static_cast<I>(~(1uLL << 58)), static_cast<I>(~(1uLL << 59)), static_cast<I>(~(1uLL << 60)), 
            static_cast<I>(~(1uLL << 61)), static_cast<I>(~(1uLL << 62)), static_cast<I>(~(1uLL << 63))
        };
        constexpr I mask = (sizeof(I) * CHAR_BIT) - 1;
        constexpr std::size_t shift = ln2(sizeof(I) * CHAR_BIT);

        bits[n >> shift] &= masks[n & mask];
    }
    
    std::size_t size() const noexcept { return m_size; }
    
    void reset() 
    { 
        std::memset(bits.get(), 0xff, m_size / CHAR_BIT + (m_size % CHAR_BIT ? 1 : 0)); 
    }

    inline I operator[](const std::size_t n) const noexcept
    {
        return test(n);
    }
};
}

#endif // BOOST_MATH_SPECIAL_FUNCTIONS_DETAIL_SIMPLE_BITSET_HPP
