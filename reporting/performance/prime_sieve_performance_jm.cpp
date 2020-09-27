// Copyright 2020 Matt Borland
//
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/math/special_functions/prime_sieve.hpp>
#include <boost/math/special_functions/prime_sieve_jm.hpp>
#include <boost/math/special_functions/interval_sieve.hpp>
#include <boost/math/special_functions/detail/linear_prime_sieve.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/dynamic_bitset.hpp>
//#include <boost/multiprecision/gmp.hpp>
#include <benchmark/benchmark.h>
#include <primesieve.hpp>
#include <vector>

template<class Integer>
void linear_sieve(benchmark::State& state)
{
    Integer upper = static_cast<Integer>(state.range(0));
    std::vector<Integer> primes;
    boost::math::prime_reserve(upper, primes);
    for(auto _ : state)
    {
        primes.clear();
        boost::math::detail::linear_sieve(upper, primes);
    }
    state.SetComplexityN(state.range(0));
}

template<class Integer>
void linear_sieve_jm_helper(Integer upper, std::vector<Integer>& primes)
{
   jm::detail::simple_bitset<std::uint64_t> masks((upper + 3) / 2);
   jm::detail::linear_sieve_classical<Integer>(masks, std::back_inserter(primes));
}

template<class Integer>
void linear_sieve_jm(benchmark::State& state)
{
    Integer upper = static_cast<Integer>(state.range(0));
    std::vector<Integer> primes;
    boost::math::prime_reserve(upper, primes);
    for(auto _ : state)
    {
        primes.clear();
        linear_sieve_jm_helper(upper, primes);
    }
    state.SetComplexityN(state.range(0));
}

// Complete Implementations
template <class Integer>
void prime_sieve_seq(benchmark::State& state)
{
    Integer upper = static_cast<Integer>(state.range(0));
    std::vector<Integer> primes;
    boost::math::prime_reserve(upper, primes);
    for(auto _ : state)
    {
        primes.clear();
        boost::math::prime_sieve(upper, primes);
    }
    state.SetComplexityN(state.range(0));
}
template <class Integer>
void prime_sieve_seq_jm(benchmark::State& state)
{
    Integer upper = static_cast<Integer>(state.range(0));
    std::vector<Integer> primes;
    boost::math::prime_reserve(upper, primes);
    for(auto _ : state)
    {
        primes.clear();
        jm::prime_sieve(upper, primes);
    }
    state.SetComplexityN(state.range(0));
}
template <class Integer>
void prime_sieve_seq_jm_oi(benchmark::State& state)
{
    Integer upper = static_cast<Integer>(state.range(0));
    std::vector<Integer> primes;
    boost::math::prime_reserve(upper, primes);
    for(auto _ : state)
    {
        primes.clear();
        jm::prime_sieve(upper, std::back_inserter(primes));
    }
    state.SetComplexityN(state.range(0));
}

template <class Integer>
void kimwalish_primes(benchmark::State& state)
{
    Integer upper = static_cast<Integer>(state.range(0));
    std::vector<Integer> primes;
    for (auto _ : state)
    {
       primes.clear();
       primesieve::generate_primes(upper, &primes);
    }
    state.SetComplexityN(state.range(0));
}

template <class Integer, class I2>
inline Integer kimwalish_prime_factorizer_helper(Integer upper, I2 value)
{
   std::vector<Integer> primes;
   primesieve::generate_primes(upper, &primes);
   for (unsigned i = 0; i < primes.size(); ++i)
      while (value % primes[i] == 0)
         value /= primes[i];
   return value;
}

template <class Integer>
void kimwalish_prime_factorizer(benchmark::State& state)
{
    Integer upper = static_cast<Integer>(state.range(0));
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(kimwalish_prime_factorizer_helper(upper, std::numeric_limits<std::uint32_t>::max()));
    }
    state.SetComplexityN(state.range(0));
}

template<class Integer>
void prime_sieve_seq_oi(benchmark::State& state)
{
   Integer upper = static_cast<Integer>(state.range(0));
   std::vector<Integer> primes;
   boost::math::prime_reserve(upper, primes);
   for (auto _ : state)
   {
      benchmark::DoNotOptimize(boost::math::detail::prime_sieve::linear_sieve(upper, std::back_inserter(primes)));
   }
   state.SetComplexityN(state.range(0));
}

constexpr uint64_t low_range = 4;
constexpr uint64_t high_range = uint64_t(1) << 32;

// Invidiual Implementations
// Linear
//BENCHMARK_TEMPLATE(linear_sieve, int32_t)->RangeMultiplier(2)->Range(1 << 1, 1 << 16)->Complexity(benchmark::oN);
//BENCHMARK_TEMPLATE(linear_sieve, int64_t)->RangeMultiplier(2)->Range(1 << 1, 1 << 16)->Complexity(benchmark::oN);

//BENCHMARK_TEMPLATE(linear_sieve, uint32_t)->RangeMultiplier(2)->Range(low_range, high_range)->Complexity(benchmark::oN);
//BENCHMARK_TEMPLATE(linear_sieve_jm, uint32_t)->RangeMultiplier(2)->Range(low_range, high_range)->Complexity(benchmark::oN);

// Segmented
//BENCHMARK_TEMPLATE(mask_sieve, int32_t)->RangeMultiplier(2)->Range(1 << 2, 2 << 22)->Complexity(benchmark::oNLogN);
//BENCHMARK_TEMPLATE(mask_sieve, int64_t)->RangeMultiplier(2)->Range(1 << 14, 2 << 26)->Complexity(benchmark::oNLogN);
//BENCHMARK_TEMPLATE(interval_sieve, int64_t)->RangeMultiplier(2)->Range(1 << 14, 2 << 26)->Complexity();
//BENCHMARK_TEMPLATE(mask_sieve, uint32_t)->RangeMultiplier(2)->Range(1 << 2, 2 << 22)->Complexity(benchmark::oNLogN);

// Complete Implemenations
//BENCHMARK_TEMPLATE(prime_sieve, int32_t)->RangeMultiplier(2)->Range(1 << 1, 1 << 30)->Complexity(benchmark::oN)->UseRealTime();
//BENCHMARK_TEMPLATE(prime_sieve, int64_t)->RangeMultiplier(2)->Range(1 << 1, 1 << 30)->Complexity(benchmark::oN)->UseRealTime();
BENCHMARK_TEMPLATE(kimwalish_primes, int64_t)->RangeMultiplier(2)->Range(low_range, high_range)->Complexity(benchmark::oN);
BENCHMARK_TEMPLATE(prime_sieve_seq, uint64_t)->RangeMultiplier(2)->Range(low_range, high_range)->Complexity(benchmark::oN);
BENCHMARK_TEMPLATE(prime_sieve_seq_oi, uint64_t)->RangeMultiplier(2)->Range(low_range, high_range)->Complexity(benchmark::oN);
BENCHMARK_TEMPLATE(prime_sieve_seq_jm, uint64_t)->RangeMultiplier(2)->Range(low_range, high_range)->Complexity(benchmark::oN);
BENCHMARK_TEMPLATE(prime_sieve_seq_jm_oi, uint64_t)->RangeMultiplier(2)->Range(low_range, high_range)->Complexity(benchmark::oN);
//BENCHMARK_TEMPLATE(prime_sieve, boost::multiprecision::cpp_int)->RangeMultiplier(2)->Range(1 << 1, 1 << 30)->Complexity(benchmark::oN)->UseRealTime();
//BENCHMARK_TEMPLATE(prime_sieve, boost::multiprecision::mpz_int)->RangeMultiplier(2)->Range(1 << 1, 1 << 30)->Complexity(benchmark::oN)->UseRealTime();

BENCHMARK_MAIN();
