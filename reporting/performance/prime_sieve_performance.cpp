//  (C) Copyright Matt Borland 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
//  Benchmarks of boost::math::prime_sieve / prime_count against primesieve
//  (https://github.com/kimwalisch/primesieve). Built by hand, for example:
//
//  clang++ -std=c++17 -O3 -march=native -DNDEBUG -DBOOST_MATH_STANDALONE -fexperimental-library \
//      -I<boost-root>/libs/math/include -I<primesieve>/include -I/opt/homebrew/include \
//      prime_sieve_performance.cpp <primesieve>/build/libprimesieve.a \
//      -L/opt/homebrew/lib -lbenchmark -lpthread -o prime_sieve_performance
//
//  Define BOOST_MATH_BENCH_NO_PRIMESIEVE to build without the reference library.

#include <boost/math/special_functions/prime_sieve.hpp>
#include <benchmark/benchmark.h>
#ifndef BOOST_MATH_BENCH_NO_PRIMESIEVE
#include <primesieve.hpp>
#endif
#include <cstdint>
#include <vector>
#include <thread>

namespace bm = boost::math;

#ifdef BOOST_MATH_PRIME_SIEVE_HAS_STD_EXECUTION
// Policies are not default constructible in every standard library, so select them by tag
template <bool Parallel>
constexpr decltype(auto) policy()
{
    if constexpr (Parallel)
    {
        return (std::execution::par);
    }
    else
    {
        return (std::execution::seq);
    }
}

template <bool Parallel>
void count_ours(benchmark::State& state)
{
    const std::uint64_t n {static_cast<std::uint64_t>(state.range(0))};
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(bm::prime_count(policy<Parallel>(), n));
    }
    state.SetItemsProcessed(static_cast<std::int64_t>(state.iterations()) * static_cast<std::int64_t>(n));
}

template <bool Parallel>
void generate_ours(benchmark::State& state)
{
    const std::uint64_t n {static_cast<std::uint64_t>(state.range(0))};
    std::vector<std::uint64_t> primes;
    bm::prime_reserve(n, primes);
    for (auto _ : state)
    {
        primes.clear();
        bm::prime_sieve(policy<Parallel>(), n, primes);
        benchmark::DoNotOptimize(primes.data());
    }
    state.SetItemsProcessed(static_cast<std::int64_t>(state.iterations()) * static_cast<std::int64_t>(n));
}

template <bool Parallel>
void window_ours(benchmark::State& state)
{
    const std::uint64_t lo {static_cast<std::uint64_t>(state.range(0))};
    const std::uint64_t width {static_cast<std::uint64_t>(state.range(1))};
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(bm::prime_count(policy<Parallel>(), lo, lo + width));
    }
}
#endif

#ifndef BOOST_MATH_BENCH_NO_PRIMESIEVE
void count_primesieve(benchmark::State& state)
{
    const std::uint64_t n {static_cast<std::uint64_t>(state.range(0))};
    // range(1) == 0 means all hardware threads
    const int threads {state.range(1) == 0 ? static_cast<int>(std::thread::hardware_concurrency()) : static_cast<int>(state.range(1))};
    primesieve::set_num_threads(threads);
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(primesieve::count_primes(0, n));
    }
    state.SetItemsProcessed(static_cast<std::int64_t>(state.iterations()) * static_cast<std::int64_t>(n));
}

void generate_primesieve(benchmark::State& state)
{
    const std::uint64_t n {static_cast<std::uint64_t>(state.range(0))};
    std::vector<std::uint64_t> primes;
    for (auto _ : state)
    {
        primes.clear();
        primesieve::generate_primes(n, &primes);
        benchmark::DoNotOptimize(primes.data());
    }
    state.SetItemsProcessed(static_cast<std::int64_t>(state.iterations()) * static_cast<std::int64_t>(n));
}

void window_primesieve(benchmark::State& state)
{
    const std::uint64_t lo {static_cast<std::uint64_t>(state.range(0))};
    const std::uint64_t width {static_cast<std::uint64_t>(state.range(1))};
    primesieve::set_num_threads(1);
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(primesieve::count_primes(lo, lo + width - 1));
    }
}
#endif

#ifdef BOOST_MATH_PRIME_SIEVE_HAS_STD_EXECUTION
BENCHMARK_TEMPLATE(count_ours, false)->RangeMultiplier(10)->Range(100000000, 100000000000)->Unit(benchmark::kMillisecond);
BENCHMARK_TEMPLATE(count_ours, true)->RangeMultiplier(10)->Range(100000000, 100000000000)->Unit(benchmark::kMillisecond)->UseRealTime();
BENCHMARK_TEMPLATE(generate_ours, false)->RangeMultiplier(10)->Range(100000000, 10000000000)->Unit(benchmark::kMillisecond);
BENCHMARK_TEMPLATE(generate_ours, true)->RangeMultiplier(10)->Range(100000000, 10000000000)->Unit(benchmark::kMillisecond)->UseRealTime();
BENCHMARK_TEMPLATE(window_ours, false)->Args({1000000000000, 1000000000})->Args({1000000000000000, 1000000000})->Args({1000000000000000000, 1000000000})->Unit(benchmark::kMillisecond);
#endif

#ifndef BOOST_MATH_BENCH_NO_PRIMESIEVE
BENCHMARK(count_primesieve)->ArgsProduct({{100000000, 1000000000, 10000000000, 100000000000}, {1}})->Unit(benchmark::kMillisecond);
BENCHMARK(count_primesieve)->ArgsProduct({{100000000, 1000000000, 10000000000, 100000000000}, {0}})->Unit(benchmark::kMillisecond)->UseRealTime();
BENCHMARK(generate_primesieve)->RangeMultiplier(10)->Range(100000000, 10000000000)->Unit(benchmark::kMillisecond);
BENCHMARK(window_primesieve)->Args({1000000000000, 1000000000})->Args({1000000000000000, 1000000000})->Args({1000000000000000000, 1000000000})->Unit(benchmark::kMillisecond);
#endif

BENCHMARK_MAIN();
