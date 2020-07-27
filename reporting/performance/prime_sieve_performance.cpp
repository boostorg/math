// Copyright 2020 Matt Borland
//
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "../../include/boost/math/special_functions/prime_sieve.hpp"
//#include <boost/math/special_functions/prime_sieve.hpp>
#include <benchmark/benchmark.h>
#include <primesieve.hpp>
#include <vector>

// Individual Algos
// Linear
template<class Z>
inline auto linear_sieve_helper(Z upper_bound, std::vector<Z> c) -> std::vector<Z>
{
    boost::math::detail::linear_sieve(upper_bound, c);
    return c;
}

template<class Z>
void linear_sieve(benchmark::State& state)
{
    Z upper = static_cast<Z>(state.range(0));
    for(auto _ : state)
    {
        std::vector<Z> primes;
        benchmark::DoNotOptimize(linear_sieve_helper(upper, primes));
    }
    state.SetComplexityN(state.range(0));
}

template<class Z>
inline auto sub_linear_sieve_helper(Z upper_bound, std::vector<Z> c) -> std::vector<Z>
{
    boost::math::detail::sub_linear_wheel_sieve(upper_bound, c);
    return c;
}

template<class Z>
void sub_linear_sieve(benchmark::State& state)
{
    Z upper = static_cast<Z>(state.range(0));
    for(auto _ : state)
    {
        std::vector<Z> primes;
        benchmark::DoNotOptimize(sub_linear_sieve_helper(upper, primes));
    }
    state.SetComplexityN(state.range(0));
}

// Segmented
template<class Z>
inline auto mask_sieve_helper(Z lower_bound, Z upper_bound, std::vector<Z> c) -> std::vector<Z>
{
    boost::math::detail::mask_sieve(lower_bound, upper_bound, c);
    return c;
}

template<class Z>
void mask_sieve(benchmark::State& state)
{
    Z lower = static_cast<Z>(2);
    Z upper = static_cast<Z>(state.range(0));
    for(auto _ : state)
    {
        std::vector<Z> primes;
        benchmark::DoNotOptimize(mask_sieve_helper(lower, upper, primes));
    }
    state.SetComplexityN(state.range(0));
}

template<class Z>
inline auto segmented_wheel_sieve_helper(Z lower_bound, Z upper_bound, std::vector<Z> c) -> std::vector<Z>
{
    boost::math::detail::linear_segmented_wheel_sieve(lower_bound, upper_bound, c);
    return c;
}

template<class Z>
void segmented_wheel_sieve(benchmark::State& state)
{
    Z lower = static_cast<Z>(2);
    Z upper = static_cast<Z>(state.range(0));
    for(auto _ : state)
    {
        std::vector<Z> primes;
        benchmark::DoNotOptimize(segmented_wheel_sieve_helper(lower, upper, primes));
    }
    state.SetComplexityN(state.range(0));
}

// Complete Implementations
template <class Z>
void prime_sieve(benchmark::State& state)
{
    Z upper = static_cast<Z>(state.range(0));
    for(auto _ : state)
    {
        std::vector<Z> primes;
        benchmark::DoNotOptimize(boost::math::prime_sieve(std::execution::par, upper, std::back_inserter(primes)));
    }
    state.SetComplexityN(state.range(0));
}

template <class Z>
void prime_sieve_partial_range(benchmark::State& state)
{
    Z upper = static_cast<Z>(state.range(0));
    Z lower = static_cast<Z>(state.range(0)) > 2 ? static_cast<Z>(state.range(0)) : 2;
    for(auto _ : state)
    {
        std::vector<Z> primes;
        benchmark::DoNotOptimize(boost::math::prime_range(std::execution::par, lower, upper, std::back_inserter(primes)));
    }
    state.SetComplexityN(state.range(0));
}

template <class Z>
void kimwalish_primes(benchmark::State& state)
{

    Z upper = static_cast<Z>(state.range(0));
    for (auto _ : state)
    {
        std::vector<Z> primes;
        primesieve::generate_primes(upper, &primes);
        benchmark::DoNotOptimize(primes.back());
    }
    state.SetComplexityN(state.range(0));
}

// Individual Algos

// Linear
BENCHMARK_TEMPLATE(linear_sieve, int64_t)->RangeMultiplier(2)->Range(1 << 1, 1 << 16)->Complexity(benchmark::oN);
BENCHMARK_TEMPLATE(sub_linear_sieve, int64_t)->RangeMultiplier(2)->Range(1 << 1, 1 << 16)->Complexity();


// Segmented
BENCHMARK_TEMPLATE(mask_sieve, int64_t)->RangeMultiplier(2)->Range(1 << 2, 2 << 22)->Complexity(benchmark::oNLogN);
BENCHMARK_TEMPLATE(segmented_wheel_sieve, int64_t)->RangeMultiplier(2)->Range(1 << 2, 2 << 22)->Complexity();

/*
// Complete Implementations
BENCHMARK_TEMPLATE(prime_sieve, int32_t)->RangeMultiplier(2)->Range(1 << 1, 1 << 22)->Complexity()->UseRealTime();
BENCHMARK_TEMPLATE(prime_sieve, int64_t)->RangeMultiplier(2)->Range(1 << 1, 1 << 30)->Complexity(benchmark::oN)->UseRealTime();
BENCHMARK_TEMPLATE(kimwalish_primes, int64_t)->RangeMultiplier(2)->Range(1 << 1, 1 << 30)->Complexity(benchmark::oN)->UseRealTime();
BENCHMARK_TEMPLATE(prime_sieve, uint32_t)->RangeMultiplier(2)->Range(1 << 1, 1 << 22)->Complexity()->UseRealTime();
BENCHMARK_TEMPLATE(prime_sieve_partial_range, int32_t)->RangeMultiplier(2)->Range(1 << 1, 1 << 22)->Complexity()->UseRealTime();
BENCHMARK_TEMPLATE(prime_sieve_partial_range, int64_t)->RangeMultiplier(2)->Range(1 << 1, 1 << 22)->Complexity()->UseRealTime();
BENCHMARK_TEMPLATE(prime_sieve_partial_range, uint32_t)->RangeMultiplier(2)->Range(1 << 1, 1 << 22)->Complexity()->UseRealTime();
*/
BENCHMARK_MAIN();
