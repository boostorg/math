//  (C) Copyright Nick Thompson 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Build: c++ -std=c++17 -O2 -I../../include fresnel_performance.cpp -lbenchmark -o fresnel_performance
//
// One benchmark per evaluation path of fresnel_c, and fresnel_cis on the series.

#include <random>
#include <vector>
#include <benchmark/benchmark.h>
#include <boost/math/special_functions/fresnel.hpp>

namespace
{

std::vector<double> uniform(double a, double b)
{
   std::mt19937_64 gen(12345);
   std::uniform_real_distribution<double> dis(a, b);
   std::vector<double> v(1024);
   for (auto& x : v)
   {
      x = dis(gen);
   }
   return v;
}

template <class F>
void run(benchmark::State& state, double a, double b, F f)
{
   const std::vector<double> xs = uniform(a, b);
   std::size_t i = 0;
   for (auto _ : state)
   {
      benchmark::DoNotOptimize(f(xs[i++ & 1023]));
   }
}

}

void FresnelCSeries(benchmark::State& state) { run(state, 0, 1.8, [](double x) { return boost::math::fresnel_c(x); }); }
void FresnelCFractionNearSwitch(benchmark::State& state) { run(state, 1.9, 3, [](double x) { return boost::math::fresnel_c(x); }); }
void FresnelCFraction(benchmark::State& state) { run(state, 3, 100, [](double x) { return boost::math::fresnel_c(x); }); }
void FresnelCAsymptotic(benchmark::State& state) { run(state, 1e8, 1e9, [](double x) { return boost::math::fresnel_c(x); }); }
void FresnelCisSeries(benchmark::State& state) { run(state, 0, 1.8, [](double x) { return boost::math::fresnel_cis(x); }); }

BENCHMARK(FresnelCSeries);
BENCHMARK(FresnelCFractionNearSwitch);
BENCHMARK(FresnelCFraction);
BENCHMARK(FresnelCAsymptotic);
BENCHMARK(FresnelCisSeries);

BENCHMARK_MAIN();
