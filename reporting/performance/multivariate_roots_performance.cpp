// Copyright Nick Thompson 2026.
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0. (See accompanying file
// LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

/*
 * Google Benchmark suite: time to solution of multivariate_newton and multivariate_halley, with the iterations
 * and evaluations of f as counters.
 *
 * Example:
 *   c++ -O2 -DNDEBUG -std=c++14 -Iinclude \
 *     reporting/performance/multivariate_roots_performance.cpp \
 *     -lbenchmark -lpthread -o multivariate_roots_bench
 */

#include <boost/math/tools/multivariate_roots.hpp>
#include <benchmark/benchmark.h>

#include <cstdint>
#include <tuple>
#include <vector>

namespace {

using boost::math::tools::multivariate_halley;
using boost::math::tools::multivariate_newton;
using vector = std::vector<double>;

// Broyden's tridiagonal function, J row-major.
struct broyden_tridiagonal
{
   std::size_t* calls;

   auto operator()(vector const& x) const
   {
      ++*calls;
      std::size_t const n = x.size();
      vector F(n);
      vector J(n * n, 0.0);
      for (std::size_t i = 0; i < n; ++i)
      {
         double const left = i > 0 ? x[i - 1] : 0.0;
         double const right = i + 1 < n ? x[i + 1] : 0.0;
         F[i] = (3 - 2 * x[i]) * x[i] - left - 2 * right + 1;
         J[i * n + i] = 3 - 4 * x[i];
         if (i > 0)
            J[i * n + i - 1] = -1;
         if (i + 1 < n)
            J[i * n + i + 1] = -2;
      }
      auto d2 = [](vector const& v) {
         vector w(v.size());
         for (std::size_t i = 0; i < v.size(); ++i)
            w[i] = -4 * v[i] * v[i];
         return w;
      };
      return std::make_tuple(F, J, d2);
   }
};

// The gradient of the Rosenbrock function, from the classic start (-1.2, 1).
struct rosenbrock_gradient
{
   std::size_t* calls;

   auto operator()(vector const& x) const
   {
      ++*calls;
      double const u = x[0];
      double const w = x[1];
      vector F{-2 * (1 - u) - 400 * u * (w - u * u), 200 * (w - u * u)};
      vector J{2 - 400 * (w - u * u) + 800 * u * u, -400 * u, -400 * u, 200.0};
      auto d2 = [u](vector const& v) { return vector{-800 * v[0] * v[1] + 2400 * u * v[0] * v[0], -400 * v[0] * v[0]}; };
      return std::make_tuple(F, J, d2);
   }
};

template <class System, bool Halley>
void run(benchmark::State& state, vector const& x0)
{
   std::size_t calls = 0;
   std::uintmax_t iterations = 0;
   for (auto _ : state)
   {
      calls = 0;
      iterations = 10000;
      System f{&calls};
      auto x = Halley ? multivariate_halley(f, x0, 53, iterations) : multivariate_newton(f, x0, 53, iterations);
      benchmark::DoNotOptimize(x.data());
   }
   state.counters["iterations"] = static_cast<double>(iterations);
   state.counters["f_calls"] = static_cast<double>(calls);
}

template <bool Halley>
void BM_broyden(benchmark::State& state)
{
   run<broyden_tridiagonal, Halley>(state, vector(static_cast<std::size_t>(state.range(0)), -1.0));
}
BENCHMARK_TEMPLATE(BM_broyden, false)->Arg(10)->Arg(100)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_broyden, true)->Arg(10)->Arg(100)->Unit(benchmark::kMicrosecond);

template <bool Halley>
void BM_rosenbrock_gradient(benchmark::State& state)
{
   run<rosenbrock_gradient, Halley>(state, vector{-1.2, 1.0});
}
BENCHMARK_TEMPLATE(BM_rosenbrock_gradient, false)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_rosenbrock_gradient, true)->Unit(benchmark::kMicrosecond);

} // namespace

BENCHMARK_MAIN();
