//  (C) Copyright Matt Borland 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
//  The CUDA prime sieve against the CPU engine: fixed counts, boundaries, random ranges,
//  element-wise output comparison, two back-to-back passes and a device memory leak check.

#include <boost/math/special_functions/prime_sieve.hpp>
#include <cuda_runtime.h>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <random>
#include <chrono>
#include <limits>
#include <execution>

#ifndef BOOST_MATH_HAS_CUDA_PRIME_SIEVE
#error "BOOST_MATH_HAS_CUDA_PRIME_SIEVE must be defined when compiling this test with nvcc and BOOST_MATH_ENABLE_CUDA"
#endif

namespace bm = boost::math;

static int failures {0};

static double now()
{
    return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}

static void check_count(std::uint64_t lo, std::uint64_t hi, std::uint64_t expected, const char* what)
{
    const double t0 {now()};
    const std::uint64_t got {bm::prime_count(bm::execution::cuda, lo, hi)};
    const double t1 {now()};
    if (got != expected)
    {
        std::printf("FAIL %s: count [%llu, %llu) = %llu, expected %llu\n", what, (unsigned long long)lo, (unsigned long long)hi, (unsigned long long)got, (unsigned long long)expected);
        ++failures;
    }
    else
    {
        std::printf("ok   %s: count [%llu, %llu) = %llu in %.3fs\n", what, (unsigned long long)lo, (unsigned long long)hi, (unsigned long long)got, t1 - t0);
    }
}

static void check_range(std::uint64_t lo, std::uint64_t hi)
{
    bm::prime_sieve_options cpu_options {};
    cpu_options.range_strategy = bm::prime_range_strategy::full_sieve;
    std::vector<std::uint64_t> cpu;
    bm::prime_range(lo, hi, cpu, cpu_options);
    std::vector<std::uint64_t> gpu;
    bm::prime_sieve_options gpu_options {};
    gpu_options.range_strategy = bm::prime_range_strategy::full_sieve;
    bm::prime_range(bm::execution::cuda, lo, hi, gpu, gpu_options);
    const std::uint64_t gpu_count {bm::prime_count(bm::execution::cuda, lo, hi, gpu_options)};
    if (cpu != gpu || gpu_count != cpu.size())
    {
        std::printf("FAIL range [%llu, %llu): cpu %zu, gpu %zu, gpu count %llu\n", (unsigned long long)lo, (unsigned long long)hi, cpu.size(), gpu.size(), (unsigned long long)gpu_count);
        for (std::size_t i {0}; i < cpu.size() && i < gpu.size(); ++i)
        {
            if (cpu[i] != gpu[i])
            {
                std::printf("     first difference at index %zu: cpu %llu gpu %llu\n", i, (unsigned long long)cpu[i], (unsigned long long)gpu[i]);
                break;
            }
        }
        ++failures;
    }
}

static void run_pass()
{
    check_count(0, 1000000, 78498, "pi(1e6)");
    check_count(0, 100000000, 5761455, "pi(1e8)");
    check_count(0, 1000000000, 50847534, "pi(1e9)");
    check_count(0, 10000000000ull, 455052511, "pi(1e10)");
    check_count(0, 100000000000ull, 4118054813ull, "pi(1e11)");
    check_count(1000000000000ull, 1000000000000ull + 100000000ull, bm::prime_count(1000000000000ull, 1000000000000ull + 100000000ull), "1e12 window");
    check_count(1000000000000000ull, 1000000000000000ull + 100000000ull, 2893937, "1e15 window");
    check_count(1000000000000000000ull, 1000000000000000000ull + 100000000ull, bm::prime_count(1000000000000000000ull, 1000000000000000000ull + 100000000ull), "1e18 window");

    const std::uint64_t top {(std::numeric_limits<std::uint64_t>::max)()};
    check_count(top - 1000, top, bm::prime_count(top - 1000, top), "top of range");
    check_count(top - 300000000ull, top, bm::prime_count(static_cast<std::uint64_t>(top - 300000000ull), top), "3e8 below 2^64");
    std::vector<std::uint64_t> last;
    bm::prime_range(bm::execution::cuda, static_cast<std::uint64_t>(18446744073709551557ull), top, last);
    if (last.size() != 1 || last[0] != 18446744073709551557ull)
    {
        std::printf("FAIL largest 64-bit prime\n");
        ++failures;
    }

    // boundaries around bytes, segments and chunks
    const std::uint64_t segment_span {32768ull * 30ull};
    const std::uint64_t chunk {std::uint64_t(1) << 32};
    const std::uint64_t lows[] = {0, 1, 2, 7, 8, 29, 30, 31, 36, 37, 59, 60, segment_span - 1, segment_span, segment_span + 1,
                                  chunk - 100000, chunk - 1, chunk, chunk + 1, 3 * chunk - 50000};
    const std::uint64_t widths[] = {1, 2, 30, 1000, 200000};
    for (const std::uint64_t lo : lows)
    {
        for (const std::uint64_t width : widths)
        {
            check_range(lo, lo + width);
        }
    }
    // several interior chunk boundaries (candidates congruent to 1 mod 30 straddle them)
    check_count(chunk - 1000, 3 * chunk + 1000, bm::prime_count(std::execution::par, chunk - 1000, 3 * chunk + 1000), "three chunk boundaries");
    check_range(2, 8);
    check_range(7, 8);
    check_range(0, 2);
    check_range(100, 100);

    // random ranges
    std::mt19937_64 rng {2026};
    for (int i {0}; i < 40; ++i)
    {
        const double e {std::uniform_real_distribution<double>(3.0, 19.0)(rng)};
        std::uint64_t lo {static_cast<std::uint64_t>(std::pow(10.0, e))};
        const std::uint64_t width {rng() % 50000000ull + 1};
        if (lo > top - width)
        {
            lo = top - width;
        }
        check_range(lo, lo + width);
    }
}

int main()
{
    std::size_t free_before {0};
    std::size_t total {0};
    cudaMemGetInfo(&free_before, &total);

    const double t0 {now()};
    run_pass();
    const double t1 {now()};
    std::printf("first pass: %.1fs, %d failures\n", t1 - t0, failures);
    const int first_failures {failures};
    run_pass();
    const double t2 {now()};
    std::printf("second pass: %.1fs, %d failures\n", t2 - t1, failures - first_failures);

    cudaDeviceSynchronize();
    std::size_t free_after {0};
    cudaMemGetInfo(&free_after, &total);
    if (free_after + (16u << 20) < free_before)
    {
        std::printf("FAIL: device memory dropped from %zu to %zu bytes free\n", free_before, free_after);
        ++failures;
    }
    std::printf("%s\n", failures == 0 ? "ALL OK" : "FAILED");
    return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
