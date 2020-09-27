// Copyright 2020 Matt Borland
//
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_SPECIAL_FUNCTIONS_DETAIL_LINEAR_PRIME_SIEVE_HPP
#define BOOST_MATH_SPECIAL_FUNCTIONS_DETAIL_LINEAR_PRIME_SIEVE_HPP

#include <memory>

namespace boost::math::detail::prime_sieve 
{
template<typename Integer, typename OutputIterator>
decltype(auto) linear_sieve(const Integer upper_bound, OutputIterator resultant_primes)
{    
    const std::size_t masks_size {static_cast<std::size_t>(upper_bound / 2 + 1)};
    std::unique_ptr<bool[]> masks {new bool[masks_size]};
    memset(masks.get(), true, sizeof(*masks.get()) * (masks_size));

    *resultant_primes++ = 2;

    for(std::size_t index {1}; index < masks_size; ++index)
    {
        if(masks[index])
        {
            *resultant_primes++ = static_cast<Integer>(2 * index + 1);
            for(std::size_t clear {index * 3 + 1}; clear < masks_size; clear += index * 2 + 1)
            {
                masks[clear] = false;
            }
        }
    }

    return resultant_primes;
}

// 4'096 is where benchmarked performance of linear_sieve begins to diverge
template<typename Integer>
static const Integer linear_sieve_limit = Integer(4'096); // Constexpr does not work with boost::multiprecision types
}

#endif // BOOST_MATH_SPECIAL_FUNCTIONS_DETAIL_LINEAR_PRIME_SIEVE_HPP
