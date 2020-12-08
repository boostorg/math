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
#include <boost/math/special_functions/detail/small_primes.hpp>
#include <boost/math/special_functions/prime_approximation.hpp>
#include <execution>
#include <vector>
#include <future>
#include <thread>
#include <limits>
#include <type_traits>
#include <iterator>
#include <cstdint>

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

template<typename ExecutionPolicy, typename Integer, typename OutputIterator>
decltype(auto) prime_sieve_iter_impl(ExecutionPolicy&& policy, const Integer upper_bound, OutputIterator resultant_primes)
{
    if (upper_bound == 2)
    {
        return resultant_primes;
    }

    else if (upper_bound <= small_prime_limit<Integer>)
    {
        small_primes(upper_bound, resultant_primes);
        return resultant_primes;
    }
    
    resultant_primes = small_primes(small_prime_limit<Integer>, resultant_primes);

    if constexpr (std::is_same_v<std::remove_reference_t<decltype(policy)>, decltype(std::execution::seq)> 
                  #if __cpp_lib_execution > 201900
                  || std::is_same_v<std::remove_reference_t<decltype(policy)>, decltype(std::execution::unseq)>
                  #endif
                  )
    {
        resultant_primes = sequential_segmented_sieve(small_prime_limit<Integer>, upper_bound, resultant_primes);
    }

    else
    {
        resultant_primes = segmented_sieve(small_prime_limit<Integer>, upper_bound, resultant_primes);
    }
    
    return resultant_primes;
}

template<typename Integer, typename OutputIterator>
inline decltype(auto) prime_sieve_iter_impl(const Integer upper_bound, OutputIterator resultant_primes)
{
    return prime_sieve_iter_impl(std::execution::seq, upper_bound, resultant_primes);
}

template<typename ExecutionPolicy, typename Integer, typename OutputIterator>
decltype(auto) prime_range_iter_impl(ExecutionPolicy&& policy, const Integer lower_bound, const Integer upper_bound, OutputIterator resultant_primes)
{
    if(lower_bound <= small_prime_limit<Integer> || upper_bound <= small_prime_limit<Integer>)
    {
        std::vector<Integer> primes;
        primes.resize(static_cast<std::size_t>(prime_approximation(upper_bound)));
        small_primes(upper_bound, primes.begin());

        auto it {primes.begin()};
        while(*it < lower_bound)
        {
            ++it;
        }

        while(it != primes.end() && *it <= upper_bound)
        {
            *resultant_primes++ = std::move(*it++);
        }
    }

    if(lower_bound <= small_prime_limit<Integer> && upper_bound > small_prime_limit<Integer>)
    {
        if constexpr (std::is_same_v<std::remove_reference_t<decltype(policy)>, decltype(std::execution::seq)> 
                      #if __cpp_lib_execution > 201900
                      || std::is_same_v<std::remove_reference_t<decltype(policy)>, decltype(std::execution::unseq)>
                      #endif
                     )
        {
            resultant_primes = detail::prime_sieve::sequential_segmented_sieve(small_prime_limit<Integer>, upper_bound, resultant_primes);
        }

        else
        {
            resultant_primes = detail::prime_sieve::segmented_sieve(small_prime_limit<Integer>, upper_bound, resultant_primes);
        }
    }

    else
    {
        if constexpr (std::is_same_v<std::remove_reference_t<decltype(policy)>, decltype(std::execution::seq)> 
                      #if __cpp_lib_execution > 201900
                      || std::is_same_v<std::remove_reference_t<decltype(policy)>, decltype(std::execution::unseq)>
                      #endif
                     )
        {
            resultant_primes = detail::prime_sieve::sequential_segmented_sieve(lower_bound, upper_bound, resultant_primes);
        }

        else
        {
            resultant_primes = detail::prime_sieve::segmented_sieve(lower_bound, upper_bound, resultant_primes);
        }
    }

    return resultant_primes;
}

template<typename Integer, typename OutputIterator>
inline decltype(auto) prime_range_iter_impl(const Integer lower_bound, const Integer upper_bound, OutputIterator resultant_primes)
{
    return prime_range_iter_impl(std::execution::seq, lower_bound, upper_bound, resultant_primes);
}

// SFINAE for dual interface
template<typename T>
class is_container
{
private:
    using yes = char; 
    struct no { char x[2]; };

    template<typename U>  
    static yes test( decltype(&U::size) );
    
    template<typename U> 
    static no test(...);    

public:
    enum { type = sizeof(test<T>(0)) == sizeof(char) };
};

// Add C++14/17esqe interface to is_container
template<typename T>
static constexpr auto is_container_t = is_container<T>::type;
}

namespace boost::math
{
template<typename ExecutionPolicy, typename Integer, typename T>
inline decltype(auto) prime_sieve_iter(ExecutionPolicy&& exec, Integer upper_bound, T output)
{
    if constexpr (detail::prime_sieve::is_container_t<std::remove_pointer_t<T>>)
    {
        if constexpr (std::is_integral_v<Integer>)
        {
            return detail::prime_sieve::prime_sieve_iter_impl(exec, upper_bound, std::begin(*output));
        }
        else
        {
            if (upper_bound < static_cast<Integer>(std::numeric_limits<std::size_t>::max()))
            {
                return detail::prime_sieve::prime_sieve_iter_impl(exec, static_cast<std::size_t>(upper_bound), std::begin(*output));
            }
            else
            {
                return detail::prime_sieve::prime_sieve_iter_impl(exec, upper_bound, std::begin(*output));
            }
        }
    }
    else
    {
        if constexpr (std::is_integral_v<Integer>)
        {
            return detail::prime_sieve::prime_sieve_iter_impl(exec, upper_bound, output);
        }
        else
        {
            if (upper_bound < static_cast<Integer>(std::numeric_limits<std::size_t>::max()))
            {
                return detail::prime_sieve::prime_sieve_iter_impl(exec, static_cast<std::size_t>(upper_bound), output);
            }
            else
            {
                return detail::prime_sieve::prime_sieve_iter_impl(exec, upper_bound, output);
            }
        }
    }
}

template<typename Integer, typename T>
inline decltype(auto) prime_sieve_iter(Integer upper_bound, T output)
{
    return prime_sieve_iter(std::execution::seq, upper_bound, output);
}

template<typename ExecutionPolicy, typename Integer, typename T>
inline decltype(auto) prime_range_iter(ExecutionPolicy&& exec, Integer upper_bound, Integer lower_bound, T output)
{
    if constexpr (detail::prime_sieve::is_container_t<std::remove_pointer_t<T>>)
    {
        if constexpr (std::is_integral_v<Integer>)
        {
            return detail::prime_sieve::prime_range_iter_impl(exec, lower_bound, upper_bound, std::begin(*output));
        }
        else
        {
            if (upper_bound < static_cast<Integer>(std::numeric_limits<std::size_t>::max()))
            {
                return detail::prime_sieve::prime_range_iter_impl(exec, static_cast<std::size_t>(lower_bound), 
                                                                  static_cast<std::size_t>(upper_bound), std::begin(*output));
            }
            else
            {
                return detail::prime_sieve::prime_range_iter_impl(exec, lower_bound, upper_bound, std::begin(*output));
            }
        }
    }
    else
    {
        if constexpr (std::is_integral_v<Integer>)
        {
            return detail::prime_sieve::prime_range_iter_impl(exec, lower_bound, upper_bound, output);
        }
        else
        {
            if (upper_bound < static_cast<Integer>(std::numeric_limits<std::size_t>::max()))
            {
                return detail::prime_sieve::prime_range_iter_impl(exec, static_cast<std::size_t>(lower_bound), static_cast<std::size_t>(upper_bound), output);
            }
            else
            {
                return detail::prime_sieve::prime_range_iter_impl(exec, lower_bound, upper_bound, output);
            }
        }
    }
}

template<typename Integer, typename T>
inline decltype(auto) prime_range_iter(Integer lower_bound, Integer upper_bound, T output)
{
    return prime_range_iter(std::execution::seq, lower_bound, upper_bound, output);
}

template<typename Integer>
inline void set_l1d_size(const Integer size)
{
    detail::prime_sieve::L1D_SIZE = static_cast<std::size_t>(size);
}
}

#endif // BOOST_MATH_SPECIAL_FUNCTIONS_PRIME_SEIVE_ITER_HPP
