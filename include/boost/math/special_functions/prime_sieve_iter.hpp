// Copyright 2020 Matt Borland
//
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_SPECIAL_FUNCTIONS_PRIME_SEIVE_ITER_HPP
#define BOOST_MATH_SPECIAL_FUNCTIONS_PRIME_SEIVE_ITER_HPP

#include <boost/math/special_functions/detail/linear_prime_sieve.hpp>
#include <boost/math/special_functions/detail/interval_prime_sieve.hpp>
#include <boost/math/special_functions/prime_approximation.hpp>
#include <boost/math/special_functions/prime_sieve.hpp>
#include <execution>
#include <cstdint>
#include <vector>
#include <future>
#include <thread>

namespace boost::math::detail::prime_sieve
{
inline std::size_t L1D_SIZE {32'768};

template<typename Integer, typename OutputIterator>
decltype(auto) sequential_segmented_sieve(const Integer lower_bound, const Integer upper_bound, OutputIterator resultant_primes)
{
    const Integer interval {static_cast<Integer>(L1D_SIZE * 8)};
    Integer current_lower_bound {lower_bound};
    Integer current_upper_bound {current_lower_bound + interval};

    if (current_upper_bound > upper_bound)
    {
        current_upper_bound = upper_bound;
    }

    std::size_t ranges {static_cast<std::size_t>((upper_bound - lower_bound) / interval)};

    IntervalSieve sieve(current_lower_bound, current_upper_bound, resultant_primes);

    for(std::size_t i {}; i < ranges; ++i)
    {
        current_lower_bound = current_upper_bound;
        current_upper_bound += interval;
        if(current_upper_bound > upper_bound)
        {
            current_upper_bound = upper_bound;
        }
        resultant_primes = sieve.NewRange(current_lower_bound, current_upper_bound);
    }
    
    return resultant_primes;
}

template<typename Integer, typename OutputIterator>
decltype(auto) segmented_sieve(const Integer lower_bound, const Integer upper_bound, OutputIterator resultant_primes)
{
    const auto num_threads {std::thread::hardware_concurrency() > 0 ? std::thread::hardware_concurrency() : 2u};
    const Integer thread_range {(upper_bound - lower_bound) / static_cast<Integer>(num_threads)};

    std::vector<std::vector<Integer>> prime_vectors(num_threads);
    std::vector<std::future<void>> future_manager;
    future_manager.reserve(num_threads);

    Integer current_lower_bound {lower_bound};
    Integer current_upper_bound {current_lower_bound + thread_range};

    for(std::size_t i {}; i < num_threads - 1; ++i)
    {
        prime_vectors[i].resize(static_cast<std::size_t>(prime_approximation(current_lower_bound, current_upper_bound)));

        future_manager.emplace_back(std::async(std::launch::async, [current_lower_bound, current_upper_bound, &prime_vectors, i]{
            sequential_segmented_sieve(current_lower_bound, current_upper_bound, prime_vectors[i].begin());
        }));

        current_lower_bound = current_upper_bound;
        current_upper_bound += thread_range;
    }

    prime_vectors.back().resize(static_cast<std::size_t>(prime_approximation(current_lower_bound, upper_bound)));
    future_manager.emplace_back(std::async(std::launch::async, [current_lower_bound, upper_bound, &prime_vectors]{
            sequential_segmented_sieve(current_lower_bound, upper_bound, prime_vectors.back().begin());
    }));

    std::size_t i {};
    for(auto&& future : future_manager)
    {
        future.get(); // Blocks to maintain proper sorting
        
        for(auto& val : prime_vectors[i])
        {
            *resultant_primes++ = std::move(val);
        }

        ++i;
    }
    
    return resultant_primes;
}
}

namespace boost::math
{
template<typename ExecutionPolicy, typename Integer, typename OutputIterator>
decltype(auto) prime_sieve_iter(ExecutionPolicy&& policy, const Integer upper_bound, OutputIterator resultant_primes)
{
    if (upper_bound == 2)
    {
        return resultant_primes;
    }

    else if (upper_bound <= detail::prime_sieve::linear_sieve_limit<Integer>)
    {
        detail::prime_sieve::linear_sieve(upper_bound, resultant_primes);
        return resultant_primes;
    }
    
    resultant_primes = detail::prime_sieve::linear_sieve(detail::prime_sieve::linear_sieve_limit<Integer>, resultant_primes);

    if constexpr (std::is_same_v<std::remove_reference_t<decltype(policy)>, decltype(std::execution::seq)> 
                  #if __cpp_lib_execution > 201900
                  || std::is_same_v<std::remove_reference_t<decltype(policy)>, decltype(std::execution::unseq)>
                  #endif
                  )
    {
        resultant_primes = detail::prime_sieve::sequential_segmented_sieve(detail::prime_sieve::linear_sieve_limit<Integer>, upper_bound, resultant_primes);
    }

    else
    {
        resultant_primes = detail::prime_sieve::segmented_sieve(detail::prime_sieve::linear_sieve_limit<Integer>, upper_bound, resultant_primes);
    }
    
    return resultant_primes;
}

template<typename Integer, typename OutputIterator>
inline decltype(auto) prime_sieve_iter(const Integer upper_bound, OutputIterator resultant_primes)
{
    return prime_sieve_iter(std::execution::seq, upper_bound, resultant_primes);
}

template<typename Integer>
inline void set_l1d_size(const Integer size)
{
    detail::prime_sieve::L1D_SIZE = static_cast<std::size_t>(size);
}
}

#endif // BOOST_MATH_SPECIAL_FUNCTIONS_PRIME_SEIVE_ITER_HPP
