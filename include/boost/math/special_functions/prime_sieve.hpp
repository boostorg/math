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
#include <memory>
#include <future>

#ifdef _MSVC_LANG
#if _MSVC_LANG >= 201703 // _MSVC_LANG == __cplusplus: https://devblogs.microsoft.com/cppblog/msvc-now-correctly-reports-__cplusplus/
#include <execution>
#endif
#else
#if __cplusplus >= 201703
#include <execution>
#endif
#endif

namespace boost { namespace math { namespace detail
{
// https://mathworld.wolfram.com/SieveofEratosthenes.html
// https://www.cs.utexas.edu/users/misra/scannedPdf.dir/linearSieve.pdf
template<class Z, class Container>
void linear_sieve(Z upper_bound, Container &c)
{
    size_t least_divisors_size{static_cast<size_t>(upper_bound + 1)};
    std::unique_ptr<Z[]> least_divisors{new Z[least_divisors_size]{0}};

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

//https://core.ac.uk/download/pdf/62440589.pdf
template<typename Z, class PrimeContainer, class Container>
void mask_sieve(Z lower_bound, Z upper_bound, const PrimeContainer& primes, Container &c)
{
    Z limit {static_cast<Z>(std::floor(std::sqrt(static_cast<double>(upper_bound)))) + 1};

    size_t primes_size {};
    auto it{primes.begin()};
    while(it != primes.end() && *it < limit)
    {
        ++primes_size;
        ++it;
    }

    const size_t n {static_cast<size_t>(upper_bound - lower_bound + 1)};
    std::unique_ptr<bool[]> is_prime {new bool[n]};
    memset(is_prime.get(), true, sizeof(*is_prime.get()) * (n));

    for(size_t i{}; i < primes_size; ++i)
    {
        Z current_prime{primes[i]};
        for(Z j {std::max(current_prime * current_prime, (lower_bound + current_prime - 1) / current_prime * current_prime)};
            j <= upper_bound; j += current_prime)
        {
            is_prime[j - lower_bound] = false;
        }
    }

    if(lower_bound == 1)
    {
        is_prime[0] = false;
    }

    for(Z i{lower_bound}; i <= upper_bound; ++i)
    {
        if(is_prime[i - lower_bound])
        {
            c.emplace_back(i);
        }
    }
}

template<typename Z, class Container>
void mask_sieve(Z lower_bound, Z upper_bound, Container &c)
{
    Z limit{static_cast<Z>(std::floor(std::sqrt(static_cast<double>(upper_bound)))) + 1};
    std::vector<Z> primes {};
    primes.reserve(limit / std::log(limit));

    boost::math::detail::linear_sieve(limit, primes);

    boost::math::detail::mask_sieve(lower_bound, upper_bound, primes, c);
}

template <typename Z, class PrimesContainer, class Container>
void segmented_sieve(Z lower_bound, Z upper_bound, const PrimesContainer &primes, Container &c)
{
    const Z L1_SIZE {32648};
    const Z interval {L1_SIZE * 4};
    Z current_lower_bound{lower_bound};
    Z current_upper_bound{current_lower_bound + interval};

    if(current_upper_bound > upper_bound)
    {
        current_upper_bound = upper_bound;
    }

    size_t ranges {static_cast<size_t>((upper_bound - lower_bound) / interval)};

    std::vector<std::vector<Z>> prime_vectors(ranges + 1);
    std::vector<std::future<void>> future_manager(ranges);

    Z primes_in_range {static_cast<Z>(current_upper_bound / std::log(static_cast<double>(current_upper_bound)) -
                                      current_lower_bound / std::log(static_cast<double>(current_lower_bound)))};

    for(size_t i {}; i < ranges; ++i)
    {
        prime_vectors[i].reserve(primes_in_range);

        future_manager.emplace_back(std::async([current_lower_bound, current_upper_bound, &primes, &prime_vectors, i]{
            boost::math::detail::mask_sieve(current_lower_bound, current_upper_bound, primes, prime_vectors[i]);
        }));

        current_lower_bound = current_upper_bound + 1;
        current_upper_bound += interval;
    }

    prime_vectors[ranges].reserve(primes_in_range);
    future_manager.emplace_back(std::async([current_lower_bound, upper_bound, &primes, &prime_vectors]{
        boost::math::detail::mask_sieve(current_lower_bound, upper_bound, primes, prime_vectors.back());
    }));

    for(auto &&future : future_manager)
    {
        if(future.valid())
        {
            future.get();
        }
    }

    for(auto &v : prime_vectors)
    {
        c.insert(c.end(), v.begin(), v.end());
    }
}

template <typename Z, class Container>
void segmented_sieve(Z lower_bound, Z upper_bound, Container &c)
{
    Z limit{static_cast<Z>(std::floor(std::sqrt(static_cast<double>(upper_bound)))) + 1};
    std::vector<Z> primes {};
    primes.reserve(limit / std::log(limit));

    // Prepare for max value so you do not have to calculate this again
    if(limit < 8192)
    {
        boost::math::detail::linear_sieve(limit, primes);
    }

    else
    {
        boost::math::detail::mask_sieve(static_cast<Z>(2), limit, primes);
    }

    boost::math::detail::segmented_sieve(lower_bound, upper_bound, primes, c);
}
} // End namespace detail

#if __cplusplus >= 201703 || _MSVC_LANG >= 201703
template<class ExecutionPolicy, typename Z, class OutputIterator>
auto prime_sieve(ExecutionPolicy&& policy, Z upper_bound, OutputIterator output) -> decltype(output)
{
    static_assert(std::is_integral<Z>::value, "No primes for floating point types");
    BOOST_ASSERT_MSG(upper_bound + 1 < std::numeric_limits<Z>::max(), "Type Overflow");

    --upper_bound; // Not inclusive, but several methods in boost::math::detail need to be

    std::vector<Z> primes {};
    primes.reserve(upper_bound / std::log(static_cast<double>(upper_bound)));

    // Range for when the linear sieve is no longer faster than threading
    if (upper_bound < 8192)
    {
        boost::math::detail::linear_sieve(upper_bound, primes);
    }


    else if (std::is_same_v<decltype(policy), std::execution::sequenced_policy>)
    {
        boost::math::detail::segmented_sieve(static_cast<Z>(2), upper_bound, primes);
    }

    else
    {
        unsigned processor_count {std::thread::hardware_concurrency()};

        // May return 0 when unable to detect
        if(processor_count == 0)
        {
            processor_count = 2;
        }

        std::vector<Z> small_primes {};
        small_primes.reserve(1000);

        // Threshold for when 2 thread performance begins to be non-linear, or when the system can only support two threads
        if(upper_bound < 16777216 || processor_count == 2)
        {
            // Split into two vectors and merge after joined to avoid data races
            std::thread t1([&small_primes] {
                boost::math::detail::linear_sieve(static_cast<Z>(8192), small_primes);
            });
            std::thread t2([upper_bound, &primes] {
                boost::math::detail::segmented_sieve(static_cast<Z>(8192), upper_bound, primes);
            });

            t1.join();
            t2.join();
            primes.insert(primes.begin(), small_primes.begin(), small_primes.end());
        }

        //If sufficiently large upper bound spawn as many threads as the system has processors for
        else
        {
            //Pre-generate all of the primes so that each thread does not have to
            Z limit{static_cast<Z>(std::floor(std::sqrt(static_cast<double>(upper_bound)))) + 1};
            std::vector<Z> pre_generated_primes {};
            pre_generated_primes.reserve(limit / std::log(limit));

            if(limit < 8192)
            {
                boost::math::detail::linear_sieve(limit, pre_generated_primes);
            }

            else
            {
                std::thread t1([&small_primes] {
                    boost::math::detail::linear_sieve(static_cast<Z>(8192), small_primes);
                });
                std::thread t2([limit, &pre_generated_primes] {
                    boost::math::detail::segmented_sieve(static_cast<Z>(8192), limit, pre_generated_primes);
                });

                t1.join();
                t2.join();
                pre_generated_primes.insert(pre_generated_primes.begin(), small_primes.begin(), small_primes.end());
            }

            std::vector<std::thread> thread_manager {};
            std::vector<std::vector<Z>> prime_vectors(processor_count);
            const Z range_per_thread = upper_bound / (processor_count);
            Z current_lower_bound {8192};
            Z current_upper_bound {current_lower_bound + range_per_thread};
            Z primes_in_range {static_cast<Z>(current_upper_bound / std::log(static_cast<double>(current_upper_bound)) -
                               current_lower_bound / std::log(static_cast<double>(current_lower_bound)))};

            for(size_t i{}; i < processor_count - 1; ++i)
            {
                prime_vectors[i].reserve(primes_in_range);

                std::thread t([current_lower_bound, current_upper_bound, &prime_vectors, i, &pre_generated_primes] {
                    boost::math::detail::segmented_sieve(current_lower_bound, current_upper_bound, pre_generated_primes,
                                                         prime_vectors[i]);
                });

                thread_manager.push_back(std::move(t));

                current_lower_bound = current_upper_bound;
                current_upper_bound += range_per_thread;
                primes_in_range = static_cast<Z>(current_upper_bound / std::log(static_cast<double>(current_upper_bound)) -
                                                 current_lower_bound / std::log(static_cast<double>(current_lower_bound)));
            }

            prime_vectors.back().reserve(primes_in_range);

            std::thread t([current_lower_bound, upper_bound, &prime_vectors, &pre_generated_primes] {
                boost::math::detail::segmented_sieve(current_lower_bound, upper_bound, pre_generated_primes, prime_vectors.back());
            });
            thread_manager.push_back(std::move(t));

            for(auto &thread : thread_manager)
            {
                thread.join();
            }

            primes.insert(primes.begin(), small_primes.begin(), small_primes.end());
            for(auto &v : prime_vectors)
            {
                primes.insert(primes.end(), v.begin(), v.end());
            }
        }
    }

    return std::move(primes.begin(), primes.end(), output);
}

template<class ExecutionPolicy, class Z, class OutputIterator>
auto prime_range(ExecutionPolicy&& policy, Z lower_bound, Z upper_bound, OutputIterator output) -> decltype(output)
{
    --upper_bound; // Not inclusive, but several methods in boost::math::detail need to be
    std::vector<Z> primes {};
    primes.reserve(upper_bound / std::log(static_cast<double>(upper_bound)));

    boost::math::prime_sieve(policy, upper_bound, std::back_inserter(primes));

    auto it{primes.begin()};
    while(*it < lower_bound && it != primes.end())
    {
        ++it;
    }

    return std::move(it, primes.end(), output);
}
#endif //__cplusplus >= 201703

template<typename Z, class OutputIterator>
auto prime_sieve(Z upper_bound, OutputIterator output) -> decltype(output)
{
    static_assert(std::is_integral<Z>::value, "No primes for floating point types");
    BOOST_ASSERT_MSG(upper_bound + 1 < std::numeric_limits<Z>::max(), "Type Overflow");

    --upper_bound; // Not inclusive, but several methods in boost::math::detail need to be
    std::vector<Z> primes{};
    primes.reserve(upper_bound / std::log(upper_bound));

    if (upper_bound <= 8192)
    {
        boost::math::detail::linear_sieve(upper_bound, primes);
    }

    else
    {
        boost::math::detail::segmented_sieve(static_cast<Z>(2), upper_bound, primes);
    }

    return std::move(primes.begin(), primes.end(), output);
}

template<class Z, class OutputIterator>
auto prime_range(Z lower_bound, Z upper_bound, OutputIterator output) -> decltype(output)
{
    --upper_bound; // Not inclusive, but several methods in boost::math::detail need to be
    std::vector<Z> primes{};
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
