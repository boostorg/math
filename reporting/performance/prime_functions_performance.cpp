// Copyright 2020 Matt Borland
//
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "../../include/boost/math/special_functions/prime_functions.hpp"

#include <benchmark/benchmark.h>

template <class Z>
void prime_sieve(benchmark::State& state)
{
    Z upper = static_cast<Z>(state.range(0));
    for(auto _ : state)
    {
        std::vector<Z> primes;
        benchmark::DoNotOptimize(boost::math::prime_sieve(static_cast<Z>(2), upper, std::back_inserter(primes)));
    }
}

template <class Z>
void prime_sieve_partial_range(benchmark::State& state)
{
    Z upper = static_cast<Z>(state.range(0));
    Z lower = static_cast<Z>(state.range(0)) > 2 ? static_cast<Z>(state.range(0)) : 2;
    for(auto _ : state)
    {
        std::vector<Z> primes;
        benchmark::DoNotOptimize(boost::math::prime_sieve(lower, upper, std::back_inserter(primes)));
    }
}

BENCHMARK_TEMPLATE(prime_sieve, int32_t)->RangeMultiplier(2)->Range(1 << 1, 1 << 22)->Complexity();
BENCHMARK_TEMPLATE(prime_sieve, int64_t)->RangeMultiplier(2)->Range(1 << 1, 1 << 22)->Complexity();
BENCHMARK_TEMPLATE(prime_sieve, uint32_t)->RangeMultiplier(2)->Range(1 << 1, 1 << 22)->Complexity();
BENCHMARK_TEMPLATE(prime_sieve_partial_range, int32_t)->RangeMultiplier(2)->Range(1 << 1, 1 << 22)->Complexity();
BENCHMARK_TEMPLATE(prime_sieve_partial_range, int64_t)->RangeMultiplier(2)->Range(1 << 1, 1 << 22)->Complexity();
BENCHMARK_TEMPLATE(prime_sieve_partial_range, uint32_t)->RangeMultiplier(2)->Range(1 << 1, 1 << 22)->Complexity();

BENCHMARK_MAIN();
