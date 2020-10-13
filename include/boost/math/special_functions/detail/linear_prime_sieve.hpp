// Copyright 2020 Matt Borland
//
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_SPECIAL_FUNCTIONS_DETAIL_LINEAR_PRIME_SIEVE_HPP
#define BOOST_MATH_SPECIAL_FUNCTIONS_DETAIL_LINEAR_PRIME_SIEVE_HPP

#include <boost/dynamic_bitset.hpp>
#include <memory>
#include <algorithm>
#include <array>
#include <cmath>

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

// Stepanov Sieve - From Mathematics to Generic Programming Chap 3
template<typename RandomAccessIterator, typename Integer>
void mark_sieve(RandomAccessIterator first, RandomAccessIterator last, Integer factor)
{
    *first = false;
    while(last - first > factor)
    {
        first = first + factor;
        *first = false;
    }
}

template<typename Bitset, typename Integer>
inline void mark_sieve(Bitset& bits, const Integer factor)
{
    for(Integer i {factor * factor}; i < bits.size(); i += factor)
    {
        bits[static_cast<std::size_t>(i)] = 0;
    }
}

template<typename RandomAccessIterator, typename Integer>
void sift(RandomAccessIterator first, Integer n)
{
    const auto last {std::next(first, static_cast<Integer>(n))};
    std::fill(first, last, true);
    Integer i {0};
    Integer index_square {3};
    Integer factor {3};
    
    for(; index_square < n; index_square += factor + factor - 2)
    {
        if(first[i])
        {
            mark_sieve(first + index_square, last, factor);
        }

        ++i;
        factor += 2;
    }
}

// TODO(mborland): Pass in a more efficient data structure (likely dynamic_bitset) to sift and post-process
template<typename Integer, typename OutputIterator>
inline decltype(auto) stepanov_sieve(Integer upper_bound, OutputIterator resultant_primes)
{
    if(upper_bound == 2)
    {
        return resultant_primes;
    }
    
    sift(resultant_primes, upper_bound);
    return resultant_primes;
}

// TODO(mborland): Pass in execution policy. mark_sieve can readily be converted to std::for_each, but dynamic_bitset would need replaced with something
//                 that has iterators
template<typename Integer, typename OutputIterator>
decltype(auto) wheel_sieve_of_eratosthenes(const Integer upper_bound, OutputIterator resultant_primes)
{
    if(upper_bound == 2)
    {
        *resultant_primes++ = static_cast<Integer>(2);
        return resultant_primes;
    }

    const Integer sqrt_upper_bound {static_cast<Integer>(std::floor(std::sqrt(static_cast<double>(upper_bound)))) + 1};
    boost::dynamic_bitset<> trial(static_cast<std::size_t>(upper_bound));
    trial.set();
    std::array<Integer, 3> primes {2, 3, 5}; // Wheel basis
    std::array<Integer, 8> wheel {7, 11, 13, 17, 19, 23, 29, 31}; // MOD 30 wheel
    const Integer wheel_mod {30};

    for(std::size_t i {}; i < primes.size(); ++i)
    {
        mark_sieve(trial, primes[i]);
        *resultant_primes++ = primes[i];
    }

    // Last value in the wheel is the starting point for the next step
    for(std::size_t i {}; i < wheel.size(); ++i)
    {
        mark_sieve(trial, wheel[i]);
        *resultant_primes++ = wheel[i];
    }

    Integer i {wheel_mod};
    for(; (i + wheel.front()) < sqrt_upper_bound; i += wheel_mod)
    {
        for(std::size_t j {}; j < wheel.size(); ++j)
        {
            Integer spoke {i + wheel[j]};
            if(trial[static_cast<std::size_t>(spoke)])
            {
                mark_sieve(trial, spoke);
                *resultant_primes++ = std::move(spoke);
            }
        }
    }

    for(; (i + wheel.front()) < upper_bound; i += wheel_mod)
    {
        for(std::size_t j {}; j < wheel.size(); ++j)
        {
            Integer spoke {i + wheel[j]};
            if(trial[static_cast<std::size_t>(spoke)])
            {
                *resultant_primes++ = std::move(spoke); 
            }
        }
    }

    return resultant_primes;
}
}

#endif // BOOST_MATH_SPECIAL_FUNCTIONS_DETAIL_LINEAR_PRIME_SIEVE_HPP
