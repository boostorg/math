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

#if __has_include(<execution>)
#include <thread>
#include <execution>
#endif

namespace boost { namespace math { namespace detail
{
// https://mathworld.wolfram.com/SieveofEratosthenes.html
// https://www.cs.utexas.edu/users/misra/scannedPdf.dir/linearSieve.pdf
template<class Z, class Container>
void linear_sieve(Z upper_bound, Container &c)
{
    size_t least_divisors_size{static_cast<size_t>(upper_bound + 1)};
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

#if __has_include(<execution>)
template<class ExecutionPolicy, typename Z, class OutputIterator>
auto prime_sieve(ExecutionPolicy&& policy, Z upper_bound, OutputIterator output) -> decltype(output)
{
    static_assert(std::is_integral<Z>::value, "No primes for floating point types");
    BOOST_ASSERT_MSG(upper_bound + 1 < std::numeric_limits<Z>::max(), "Type Overflow");

    std::vector<Z> primes;
    primes.reserve(upper_bound / std::log(static_cast<double>(upper_bound)));

    // Range for when the linear sieve is no longer faster than threading
    if (upper_bound < 8192)
    {
        boost::math::detail::linear_sieve(upper_bound, primes);
    }


    else if (std::is_same_v<decltype(policy), std::execution::sequenced_policy>)
    {
        boost::math::detail::mask_sieve(static_cast<Z>(2), upper_bound, primes);
    }

    else
    {
        unsigned processor_count {std::thread::hardware_concurrency()};

        // May return 0 when unable to detect
        if(processor_count == 0)
        {
            processor_count = 2;
        }

        std::vector<Z> small_primes;
        small_primes.reserve(1000);

        // Threshold for when 2 thread performance begins to be non-linear, or when the system can only support two threads
        if(upper_bound <= 16777216 || processor_count == 2)
        {
            // Split into two vectors and merge after joined to avoid data races
            std::thread t1([&small_primes] {
                boost::math::detail::linear_sieve(static_cast<Z>(8192), small_primes);
            });
            std::thread t2([upper_bound, &primes] {
                boost::math::detail::mask_sieve(static_cast<Z>(8192), upper_bound, primes);
            });

            t1.join();
            t2.join();
            primes.insert(primes.begin(), small_primes.begin(), small_primes.end());
        }

        //If sufficiently large upper bound spawn as many threads as the system has processors for
        else
        {
            std::vector<std::thread> thread_manager;
            std::vector<std::vector<Z>> prime_vectors(processor_count - 1);
            const Z range_per_thread = upper_bound / (processor_count - 1);
            Z current_lower_bound {8192};
            Z current_upper_bound {current_lower_bound + range_per_thread};
            Z primes_in_range {static_cast<Z>(current_upper_bound / std::log(static_cast<double>(current_upper_bound)) -
                               current_lower_bound / std::log(static_cast<double>(current_lower_bound)))};

            std::thread t1([upper_bound, &small_primes] {
                boost::math::detail::linear_sieve(static_cast<Z>(8192), small_primes);
            });
            thread_manager.push_back(std::move(t1));

            for(size_t i{1}; i < processor_count - 1; ++i)
            {
                std::vector<Z> temp;
                temp.reserve(primes_in_range);
                prime_vectors.emplace_back(temp);
                std::thread t([current_lower_bound, current_upper_bound, &prime_vectors, i] {
                    boost::math::detail::mask_sieve(current_lower_bound, current_upper_bound, prime_vectors[i]);
                });

                thread_manager.push_back(std::move(t));

                current_lower_bound = current_upper_bound;
                current_upper_bound += range_per_thread;
                primes_in_range = static_cast<Z>(current_upper_bound / std::log(static_cast<double>(current_upper_bound)) -
                                                 current_lower_bound / std::log(static_cast<double>(current_lower_bound)));
            }

            std::vector<Z> temp;
            temp.reserve(primes_in_range);
            prime_vectors.emplace_back(temp);
            std::thread t([current_lower_bound, upper_bound, &prime_vectors] {
                boost::math::detail::mask_sieve(current_lower_bound, upper_bound, prime_vectors.back());
            });

            thread_manager.push_back(std::move(t));

            for(auto &thread : thread_manager)
            {
                thread.join();
            }

            primes.insert(primes.begin(), small_primes.begin(), small_primes.end());
            for(auto &v : prime_vectors)
            {
                primes.insert(primes.begin(), v.begin(), v.end());
            }
        }
    }

    return std::move(primes.begin(), primes.end(), output);
}

template<class ExecutionPolicy, class Z, class OutputIterator>
auto prime_range(ExecutionPolicy&& policy, Z lower_bound, Z upper_bound, OutputIterator output) -> decltype(output)
{
    std::vector<Z> primes;
    primes.reserve(upper_bound / std::log(static_cast<double>(upper_bound)));

    boost::math::prime_sieve(policy, upper_bound, std::back_inserter(primes));

    auto it{primes.begin()};
    while(*it < lower_bound && it != primes.end())
    {
        ++it;
    }

    return std::move(it, primes.end(), output);
}
#endif //__has_include(<execution>)

template<typename Z, class OutputIterator>
auto prime_sieve(Z upper_bound, OutputIterator output) -> decltype(output)
{
    static_assert(std::is_integral<Z>::value, "No primes for floating point types");
    BOOST_ASSERT_MSG(upper_bound + 1 < std::numeric_limits<Z>::max(), "Type Overflow");

    std::vector<Z> primes;
    primes.reserve(upper_bound / std::log(upper_bound));

    if (upper_bound <= 8192)
    {
        boost::math::detail::linear_sieve(upper_bound, primes);
    }

    else
    {
        boost::math::detail::mask_sieve(static_cast<Z>(2), upper_bound, primes);
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
}}

#endif //BOOST_MATH_SPECIAL_FUNCTIONS_PRIME_SIEVE_HPP
