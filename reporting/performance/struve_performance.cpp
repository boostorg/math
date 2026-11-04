//  (C) Copyright Nick Thompson 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Build: clang++ -std=c++17 -O2 -I../../include -I/path/to/boost struve_performance.cpp
//        -I/opt/homebrew/include -L/opt/homebrew/lib -lbenchmark -o struve_performance
//
// One benchmark per evaluation path of struve_h and struve_l, with cyl_neumann and cyl_bessel_i for scale.

#include <random>
#include <vector>
#include <benchmark/benchmark.h>
#include <boost/math/special_functions/struve.hpp>

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
void run(benchmark::State& state, double nu, double a, double b, F f)
{
   const std::vector<double> xs = uniform(a, b);
   std::size_t i = 0;
   for (auto _ : state)
   {
      benchmark::DoNotOptimize(f(nu, xs[i++ & 1023]));
   }
}

}

// Power series.
void StruveHSeries(benchmark::State& state) { run(state, 1.0, 0.1, 3, [](double nu, double x) { return boost::math::struve_h(nu, x); }); }
// Neumann series in J_{nu+1/2+k}.
void StruveHNeumann(benchmark::State& state) { run(state, 40.25, 30, 45, [](double nu, double x) { return boost::math::struve_h(nu, x); }); }
// K by the exp-sinh integral, plus Y.
void StruveHIntegral(benchmark::State& state) { run(state, 1.0, 5, 30, [](double nu, double x) { return boost::math::struve_h(nu, x); }); }
// K by the asymptotic expansion, plus Y.
void StruveHAsymptotic(benchmark::State& state) { run(state, 1.0, 50, 1000, [](double nu, double x) { return boost::math::struve_h(nu, x); }); }
// Negative order: K by the integral, plus Y.
void StruveHNegativeOrder(benchmark::State& state) { run(state, -7.3, 10, 30, [](double nu, double x) { return boost::math::struve_h(nu, x); }); }
void StruveHDerivative(benchmark::State& state)
{
   run(state, 1.0, 5, 30, [](double nu, double x) { double d; double v = boost::math::struve_h(nu, x, &d); return v + d; });
}
void CylNeumann(benchmark::State& state) { run(state, 1.0, 5, 30, [](double nu, double x) { return boost::math::cyl_neumann(nu, x); }); }

void StruveLSeries(benchmark::State& state) { run(state, 1.0, 0.1, 15, [](double nu, double x) { return boost::math::struve_l(nu, x); }); }
// I + M, M by its asymptotic expansion.
void StruveLAsymptotic(benchmark::State& state) { run(state, 1.0, 20, 700, [](double nu, double x) { return boost::math::struve_l(nu, x); }); }
void CylBesselI(benchmark::State& state) { run(state, 1.0, 20, 700, [](double nu, double x) { return boost::math::cyl_bessel_i(nu, x); }); }

BENCHMARK(StruveHSeries);
BENCHMARK(StruveHNeumann);
BENCHMARK(StruveHIntegral);
BENCHMARK(StruveHAsymptotic);
BENCHMARK(StruveHNegativeOrder);
BENCHMARK(StruveHDerivative);
BENCHMARK(CylNeumann);
BENCHMARK(StruveLSeries);
BENCHMARK(StruveLAsymptotic);
BENCHMARK(CylBesselI);

BENCHMARK_MAIN();
