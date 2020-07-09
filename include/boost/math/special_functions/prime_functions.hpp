// Copyright 2020 Matt Borland
//
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_SPECIAL_FUNCTIONS_PRIME_FUNCTIONS_HPP
#define BOOST_MATH_SPECIAL_FUNCTIONS_PRIME_FUNCTIONS_HPP

#include <boost/math/special_functions/prime.hpp>
#include <boost/assert.hpp>
#include <deque>
#include <vector>
#include <iterator>
#include <iostream>

namespace boost { namespace math
{

// https://mathworld.wolfram.com/SieveofEratosthenes.html
// https://www.cs.utexas.edu/users/misra/scannedPdf.dir/linearSieve.pdf
template<class Z, class OutputIterator>
auto prime_sieve(Z lower_bound, Z upper_bound, OutputIterator output) -> decltype(output)
{
    static_assert(std::is_integral<Z>::value, "No primes for floating point types");
    std::vector<Z> least_divisors;
    std::deque<Z> primes;

    try
    {
        least_divisors.reserve(upper_bound + 1);
        for (size_t i{}; i < upper_bound + 1; ++i)
        {
            least_divisors.emplace_back(0);
        }
    }

    catch (const std::exception &e)
    {
        // If exception is thrown it is most likely std::bad_alloc
        std::cerr << e.what() << '\n';
        throw;
    }


    for (Z i{2}; i <= upper_bound; ++i)
    {
        if (least_divisors[i] == 0)
        {
            least_divisors[i] = i;
            primes.emplace_back(i);
        }

        for (size_t j{}; j < least_divisors.size(); ++j)
        {
            if (j >= primes.size())
            {
                break;
            }

            if (primes[j] > least_divisors[i])
            {
                break;
            }

            if (i * primes[j] > upper_bound)
            {
                break;
            }

            least_divisors[i * primes[j]] = primes[j];
        }
    }

    auto it{primes.begin()};
    while (*it < lower_bound && it != primes.end())
    {
        primes.pop_front();
        ++it;
    }

    return std::move(primes.begin(), primes.end(), output);
}

template<class Z, class OutputIterator>
auto prime_range(Z lower_bound, Z upper_bound, OutputIterator output) -> decltype(output)
{
    if (upper_bound <= 104729)
    {
        Z i{2};
        unsigned counter {};
        std::deque<Z> primes;
        while (i <= upper_bound)
        {
            if (i >= lower_bound)
            {
                primes.emplace_back(i);
            }
            ++counter;
            i = static_cast<Z>(boost::math::prime(counter));
        }

        return std::move(primes.begin(), primes.end(), output);
    } else
    {
        return prime_sieve(lower_bound, upper_bound, output);
    }
}

template<class Z, class OutputIterator>
inline auto prime_range(Z upper_bound, OutputIterator output) -> decltype(output)
{
    return prime_range(static_cast<Z>(2), upper_bound, output);
}
}}

#endif //BOOST_MATH_SPECIAL_FUNCTIONS_PRIME_FUNCTIONS_HPP
