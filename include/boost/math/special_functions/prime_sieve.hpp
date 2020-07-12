// Copyright 2020 Matt Borland
//
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_SPECIAL_FUNCTIONS_PRIME_SIEVE_HPP
#define BOOST_MATH_SPECIAL_FUNCTIONS_PRIME_SIEVE_HPP

#include <boost/math/special_functions/prime.hpp>
#include <boost/assert.hpp>
#include <deque>
#include <vector>
#include <iterator>
#include <cmath>
#include <thread>

namespace boost { namespace math { namespace detail
{
// https://mathworld.wolfram.com/SieveofEratosthenes.html
// https://www.cs.utexas.edu/users/misra/scannedPdf.dir/linearSieve.pdf
template<class Z, class Container>
void linear_sieve(Z upper_bound, Container &c)
{
    Z least_divisors_size{upper_bound + 1};
    Z *least_divisors{new Z[least_divisors_size]{0}};

    for (Z i{2}; i <= upper_bound; ++i)
    {
        if (least_divisors[i] == 0)
        {
            least_divisors[i] = i;
            c.emplace_back(i);
        }

        for (size_t j{}; j < least_divisors_size; ++j)
        {
            if (j >= c.size())
            {
                break;
            }

            else if (c[j] > least_divisors[i])
            {
                break;
            }

            else if (i * c[j] > upper_bound)
            {
                break;
            }

            else
            {
                least_divisors[i * c[j]] = c[j];
            }
        }
    }

    delete[] least_divisors;
}

template<class Z, class Container>
void prime_table(Z upper_bound, Container &c)
{
    Z i{2};
    unsigned counter{};

    while (i <= upper_bound && counter < 9999) // 10k elements are in the lookup table
    {
        c.emplace_back(i);
        ++counter;
        i = static_cast<Z>(boost::math::prime(counter));
    }
}

template<class Z, class Container>
void mask_sieve(Z lower_bound, Z upper_bound, Container &c)
{
    Z limit{static_cast<Z>(std::floor(std::sqrt(static_cast<double>(upper_bound)))) + 1};
    std::vector<Z> primes;
    primes.reserve(limit / std::log(limit));

    boost::math::detail::linear_sieve(limit, primes);

    const Z n{upper_bound - lower_bound + 1};
    bool *mask{new bool[n + 1]{false}};

    for (size_t i{}; i < primes.size(); ++i)
    {
        Z lower_limit = std::floor(lower_bound / primes[i]) * primes[i];

        if (lower_limit < lower_bound)
        {
            lower_limit += primes[i];
        }

        if (lower_limit == primes[i])
        {
            lower_limit += primes[i];
        }

        for (Z j{lower_limit}; j <= upper_bound; j += primes[i])
        {
            mask[j - lower_bound] = true;
        }
    }

    // Numbers which are not masked in range, are prime
    for (Z i{lower_bound}; i <= upper_bound; i++)
    {
        if (!mask[i - lower_bound])
        {
            if (i >= lower_bound)
            {
                c.emplace_back(i);
            }
        }
    }

    delete[] mask;
}
} // End namespace detail

template<typename Z, class OutputIterator>
auto prime_sieve(Z upper_bound, OutputIterator output) -> decltype(output)
{
    static_assert(std::is_integral<Z>::value, "No primes for floating point types");
    BOOST_ASSERT_MSG(upper_bound + 1 < std::numeric_limits<Z>::max(), "Type Overflow");

    std::vector<Z> primes;
    primes.reserve(upper_bound / std::log(upper_bound));

    if (upper_bound <= 104729)
    {
        boost::math::detail::prime_table(upper_bound, primes);
    }

    else
    {
        std::vector<Z> small_primes;
        small_primes.reserve(1000);

        // Spilt into two vectors and merge after joined to avoid data races
        std::thread t1([upper_bound, &small_primes]{boost::math::detail::prime_table(static_cast<Z>(104729), small_primes);});
        std::thread t2([upper_bound, &primes]{boost::math::detail::mask_sieve(static_cast<Z>(104729), upper_bound, primes);});

        t1.join();
        t2.join();
        primes.insert(primes.begin(), small_primes.begin(), small_primes.end());
    }

    return std::move(primes.begin(), primes.end(), output);
}

template<class Z, class OutputIterator>
auto prime_range(Z lower_bound, Z upper_bound, OutputIterator output) -> decltype(output)
{
    std::vector<Z> primes;
    primes.reserve(upper_bound / std::log(static_cast<double>(upper_bound)));

    boost::math::prime_sieve(upper_bound, std::back_inserter(primes));

    auto it{primes.begin()};
    while(*it < lower_bound && it != primes.end())
    {
        ++it;
    }

    return std::move(it, primes.end(), output);
}

template<class Z, class OutputIterator>
inline auto prime_range(Z upper_bound, OutputIterator output) -> decltype(output)
{
    return prime_range(static_cast<Z>(2), upper_bound, output);
}
}}

#endif //BOOST_MATH_SPECIAL_FUNCTIONS_PRIME_SIEVE_HPP
