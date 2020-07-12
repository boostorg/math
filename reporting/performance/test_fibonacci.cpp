//  Copyright Madhur Chauhan 2020.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <benchmark/benchmark.h>
#include <boost/math/special_functions/fibonacci.hpp>
#include <boost/multiprecision/gmp.hpp>
#include <utility>

using T = boost::multiprecision::mpz_int;

auto fib_rec(unsigned long long n) -> std::pair<T, T> {
    if (n == 0) return {0, 1};
    auto p = fib_rec(n >> 1);
    T c = p.first * (2 * p.second - p.first);
    T d = p.first * p.first + p.second * p.second;
    return (n & 1) ? std::make_pair(d, c + d) : std::make_pair(c, d);
}

static void recursive_slow(benchmark::State &state) {
    for (auto _ : state)
        benchmark::DoNotOptimize(fib_rec(state.range(0)).first);
}
constexpr int bm_start = 1 << 3, bm_end = 1 << 22;
BENCHMARK(recursive_slow)->Range(bm_start, bm_end);

static void iterative_fast(benchmark::State &state) {
    for (auto _ : state)
        benchmark::DoNotOptimize(boost::math::fibonacci<T>(state.range(0)));
}
BENCHMARK(iterative_fast)->Range(bm_start, bm_end);

BENCHMARK_MAIN();

/*
Expected output:

Run on (4 X 2400 MHz CPU s)
CPU Caches:
  L1 Data 32K (x4)
  L1 Instruction 32K (x4)
  L2 Unified 256K (x4)
  L3 Unified 8192K (x4)
Load Average: 0.53, 0.65, 0.80
-----------------------------------------------------------------
Benchmark                       Time             CPU   Iterations
-----------------------------------------------------------------
recursive_slow/8             3525 ns         3521 ns       153697
recursive_slow/64            6014 ns         6003 ns       118048
recursive_slow/512           8802 ns         8784 ns        78288
recursive_slow/4096         13661 ns        13660 ns        52241
recursive_slow/32768        70636 ns        70632 ns         9678
recursive_slow/262144     1217721 ns      1217522 ns          570
recursive_slow/2097152   20469874 ns     20468924 ns           35
recursive_slow/4194304   43150902 ns     43142089 ns           17
iterative_fast/8             1419 ns         1419 ns       495748
iterative_fast/64            2389 ns         2389 ns       285226
iterative_fast/512           4168 ns         4167 ns       169008
iterative_fast/4096          8288 ns         8288 ns        92749
iterative_fast/32768        64556 ns        64555 ns        11057
iterative_fast/262144     1200102 ns      1200053 ns          576
iterative_fast/2097152   19781647 ns     19781013 ns           35
iterative_fast/4194304   42263450 ns     42195817 ns           17
*/