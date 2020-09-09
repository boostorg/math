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
// https://mathworld.wolfram.com/SieveofEratosthenes.html
// https://www.cs.utexas.edu/users/misra/scannedPdf.dir/linearSieve.pdf
template <typename Integer, typename ForwardIterator>
void linear_sieve(const Integer upper_bound, const ForwardIterator first, const ForwardIterator last) noexcept
{
    const std::size_t least_divisors_size{static_cast<std::size_t>(upper_bound + 1)};
    std::unique_ptr<Integer[]> least_divisors{new Integer[least_divisors_size]{0}};
    auto current {first};

    for (std::size_t i{2}; i < upper_bound; ++i)
    {
        if (least_divisors[i] == 0)
        {
            least_divisors[i] = i;
            *current++ = i;
        }

        for (std::size_t j{}; (first + j) < last && i * *(first + j) <= upper_bound && *(first + j) <= least_divisors[i] && j < least_divisors_size; ++j)
        {
            least_divisors[i * static_cast<std::size_t>(*(first + j))] = *(first + j);
        }
    }
}

// 4'096 is where benchmarked performance of linear_sieve begins to diverge
template<typename Integer>
static const Integer linear_sieve_limit = Integer(4'096); // Constexpr does not work with boost::multiprecision types
}

#endif // BOOST_MATH_SPECIAL_FUNCTIONS_DETAIL_LINEAR_PRIME_SIEVE_HPP
