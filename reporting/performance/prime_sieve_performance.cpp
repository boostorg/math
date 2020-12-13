// Copyright 2020 Matt Borland
//
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/math/special_functions/prime_sieve.hpp>
#include <boost/math/special_functions/detail/linear_prime_sieve.hpp>
#include <boost/math/special_functions/detail/small_primes.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/gmp.hpp>
#include <benchmark/benchmark.h>
#include <primesieve.hpp>
#include <vector>

// Individual Algos
template<typename Integer>
void linear_sieve(benchmark::State& state)
{
    Integer upper {static_cast<Integer>(state.range(0))};
    std::vector<Integer> primes(boost::math::prime_approximation(upper));
    
    for(auto _ : state)
    {
        benchmark::DoNotOptimize(boost::math::detail::prime_sieve::linear_sieve(upper, primes.begin()));
    }
    state.SetComplexityN(state.range(0));
}

template<typename Integer>
void small_primes(benchmark::State& state)
{
    Integer upper {static_cast<Integer>(state.range(0))};
    std::vector<Integer> primes(boost::math::prime_approximation(upper));
    
    for(auto _ : state)
    {
        benchmark::DoNotOptimize(boost::math::detail::prime_sieve::small_primes(upper, primes.begin()));
    }
    state.SetComplexityN(state.range(0));
}

template<typename Integer, typename OutputIterator>
inline OutputIterator interval_sieve_helper(Integer lower_bound, Integer upper_bound, OutputIterator out)
{
    boost::math::detail::prime_sieve::IntervalSieve sieve(lower_bound, upper_bound, out);
    return out;
}

template<typename Integer>
void interval_sieve(benchmark::State& state)
{
    // In practice the lower bound is never going to be less than the limit for small primes
    Integer lower {boost::math::detail::prime_sieve::small_prime_limit<Integer>()};
    Integer upper {static_cast<Integer>(state.range(0))};
    std::vector<Integer> primes(boost::math::prime_approximation(upper));

    for(auto _ : state)
    {
        benchmark::DoNotOptimize(interval_sieve_helper(lower, upper, primes.begin()));
    }
    state.SetComplexityN(state.range(0));
}

// Composite Algos
template<typename Integer>
void prime_sieve(benchmark::State& state)
{
    Integer upper {static_cast<Integer>(state.range(0))};
    std::vector<Integer> primes(boost::math::prime_approximation(upper));

    for(auto _ : state)
    {
        benchmark::DoNotOptimize(boost::math::prime_sieve(std::execution::par, upper, primes.begin()));
    }
    state.SetComplexityN(state.range(0));
}

// Reference Algo
template <typename Integer>
inline auto kimwalish_primes_helper(Integer upper, std::vector<Integer> primes) -> std::vector<Integer>
{
    primesieve::generate_primes(upper, &primes);
    return primes;
}

template <typename Integer>
void kimwalish_primes(benchmark::State& state)
{
    Integer upper {static_cast<Integer>(state.range(0))};
    std::vector<Integer> primes;
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(kimwalish_primes_helper(upper, primes));
    }
    state.SetComplexityN(state.range(0));
}

// Invidiual Implementations
// Linear
//BENCHMARK_TEMPLATE(linear_sieve, int32_t)->RangeMultiplier(2)->Range(1 << 1, 1 << 16)->Complexity(benchmark::oN);
//BENCHMARK_TEMPLATE(linear_sieve, int64_t)->RangeMultiplier(2)->Range(1 << 1, 1 << 16)->Complexity(benchmark::oN);
//BENCHMARK_TEMPLATE(linear_sieve, uint32_t)->RangeMultiplier(2)->Range(1 << 1, 1 << 16)->Complexity(benchmark::oN);
//BENCHMARK_TEMPLATE(small_primes, int64_t)->RangeMultiplier(2)->Range(1 << 3, 1 << 10)->Complexity(benchmark::oN)->UseRealTime();
//BENCHMARK_TEMPLATE(linear_sieve, int64_t)->RangeMultiplier(2)->Range(1 << 3, 1 << 10)->Complexity(benchmark::oN)->UseRealTime();
BENCHMARK_TEMPLATE(interval_sieve, int64_t)->RangeMultiplier(2)->Range(1 << 14, 2 << 26)->Complexity();

// Complete Implemenations
//BENCHMARK_TEMPLATE(prime_sieve, int32_t)->RangeMultiplier(2)->Range(1 << 1, 1 << 30)->Complexity(benchmark::oN)->UseRealTime();
//BENCHMARK_TEMPLATE(prime_sieve, int64_t)->RangeMultiplier(2)->Range(1 << 1, 1 << 30)->Complexity(benchmark::oN)->UseRealTime();
//BENCHMARK_TEMPLATE(prime_sieve_wrapper, int64_t)->RangeMultiplier(2)->Range(1 << 1, 1 << 30)->Complexity(benchmark::oN)->UseRealTime();
BENCHMARK_TEMPLATE(prime_sieve, int64_t)->RangeMultiplier(2)->Range(1 << 1, 1 << 30)->Complexity(benchmark::oN)->UseRealTime();
BENCHMARK_TEMPLATE(kimwalish_primes, int64_t)->RangeMultiplier(2)->Range(1 << 1, 1 << 30)->Complexity(benchmark::oN)->UseRealTime(); // Benchmark
//BENCHMARK_TEMPLATE(prime_sieve, uint32_t)->RangeMultiplier(2)->Range(1 << 1, 1 << 30)->Complexity(benchmark::oN)->UseRealTime();
//BENCHMARK_TEMPLATE(prime_sieve, boost::multiprecision::cpp_int)->RangeMultiplier(2)->Range(1 << 1, 1 << 30)->Complexity(benchmark::oN)->UseRealTime();
//BENCHMARK_TEMPLATE(prime_sieve, boost::multiprecision::mpz_int)->RangeMultiplier(2)->Range(1 << 1, 1 << 30)->Complexity(benchmark::oN)->UseRealTime();

BENCHMARK_MAIN();
