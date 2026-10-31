// Copyright Nick Thompson 2026.
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0. (See accompanying file
// LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

/*
 * Google Benchmark suite for the Leja rules: building the points and weights, the n-point rule,
 * the adaptive rule, and the sparse grid; the calls counter is the number of evaluations of f.
 *
 * Example:
 *   c++ -O3 -DNDEBUG -std=c++17 -Iinclude \
 *     reporting/performance/leja_quadrature_performance.cpp \
 *     -lbenchmark -lpthread -o leja_bench
 */

#include <boost/math/quadrature/leja.hpp>
#include <boost/math/quadrature/sparse_grid.hpp>
#include <benchmark/benchmark.h>

#include <array>
#include <cmath>
#include <vector>

namespace {

using boost::math::quadrature::leja_hermite;
using boost::math::quadrature::leja_jacobi;
using boost::math::quadrature::leja_laguerre;
using boost::math::quadrature::leja_legendre;

template <class Rule>
void BM_construct(benchmark::State& state)
{
    for (auto _ : state)
    {
        Rule rule;
        benchmark::DoNotOptimize(rule.weights(static_cast<unsigned>(state.range(0))).data());
    }
}
BENCHMARK_TEMPLATE(BM_construct, leja_legendre<double>)->RangeMultiplier(2)->Range(16, 128)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_construct, leja_hermite<double>)->RangeMultiplier(2)->Range(16, 128)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_construct, leja_laguerre<double>)->RangeMultiplier(2)->Range(16, 128)->Unit(benchmark::kMicrosecond);
void BM_construct_jacobi(benchmark::State& state)
{
    for (auto _ : state)
    {
        leja_jacobi<double> rule(-0.5, 2.0);
        benchmark::DoNotOptimize(rule.weights(static_cast<unsigned>(state.range(0))).data());
    }
}
BENCHMARK(BM_construct_jacobi)->RangeMultiplier(2)->Range(16, 128)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_construct, leja_legendre<float>)->Arg(64)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_construct, leja_legendre<long double>)->Arg(64)->Unit(benchmark::kMicrosecond);

// A rule with the points the benchmarks need.
template <class Rule>
const Rule& prepared()
{
    static const Rule r;
    r.weights(128);
    return r;
}
const auto& legendre = prepared<leja_legendre<double>>;
const auto& hermite = prepared<leja_hermite<double>>;

// The n-point rule, n = range(0), or the adaptive rule; s keeps the integrand from being hoisted out of the loop.
template <bool adaptive, class Rule, class F>
void integrate(benchmark::State& state, const Rule& rule, F f, double tolerance = 0)
{
    double s = 0;
    unsigned calls = 0;
    for (auto _ : state)
    {
        calls = 0;
        benchmark::DoNotOptimize(s);
        if constexpr (adaptive)
        {
            auto Q = rule.integrate([&](double x) { ++calls; return f(x + s); }, tolerance);
            benchmark::DoNotOptimize(Q);
        }
        else
        {
            auto Q = rule.integrate([&](double x) { return f(x + s); }, state.range(0));
            benchmark::DoNotOptimize(Q);
        }
    }
    if constexpr (adaptive)
        state.counters["calls"] = calls;
}

auto cos3 = [](double x) { return std::cos(3 * x); };
auto runge = [](double x) { return 1 / (1 + 16 * x * x); };
auto lorentz = [](double x) { return 1 / (1 + x * x); };
auto three = [](double x) { return std::array<double, 3>{std::exp(x), std::cos(3 * x), 1 / (1 + 16 * x * x)}; };

void BM_legendre_exp(benchmark::State& state) { integrate<false>(state, legendre(), [](double x) { return std::exp(x); }); }
BENCHMARK(BM_legendre_exp)->Arg(8)->Arg(16)->Arg(32)->Arg(64);
void BM_legendre_cos3x(benchmark::State& state) { integrate<false>(state, legendre(), cos3); }
BENCHMARK(BM_legendre_cos3x)->Arg(8)->Arg(16)->Arg(32)->Arg(64);
void BM_hermite_cos3x(benchmark::State& state) { integrate<false>(state, hermite(), cos3); }
BENCHMARK(BM_hermite_cos3x)->Arg(8)->Arg(16)->Arg(32)->Arg(64);
void BM_legendre_vector3(benchmark::State& state) { integrate<false>(state, legendre(), three); }
BENCHMARK(BM_legendre_vector3)->Arg(8)->Arg(16)->Arg(32)->Arg(64);

void BM_adaptive_legendre_cos3x(benchmark::State& state) { integrate<true>(state, legendre(), cos3, 1e-12); }
BENCHMARK(BM_adaptive_legendre_cos3x);
void BM_adaptive_legendre_runge(benchmark::State& state) { integrate<true>(state, legendre(), runge, 1e-6); }
BENCHMARK(BM_adaptive_legendre_runge);
void BM_adaptive_hermite_cos3x(benchmark::State& state) { integrate<true>(state, hermite(), cos3, 1e-12); }
BENCHMARK(BM_adaptive_hermite_cos3x);
void BM_adaptive_hermite_lorentz(benchmark::State& state) { integrate<true>(state, hermite(), lorentz, 1e-6); }
BENCHMARK(BM_adaptive_hermite_lorentz);
void BM_adaptive_legendre_vector3(benchmark::State& state) { integrate<true>(state, legendre(), three, 1e-6); }
BENCHMARK(BM_adaptive_legendre_vector3);

// exp(sum_i x_i / 2^i) and cos(sum_i x_i / (i + 1)) in range(0) dimensions.
template <class Rule, class F>
void sparse(benchmark::State& state, const Rule& rule, F f, double tolerance)
{
    const unsigned d = static_cast<unsigned>(state.range(0));
    unsigned calls = 0;
    for (auto _ : state)
    {
        calls = 0;
        auto Q = boost::math::quadrature::sparse_grid([&](const std::vector<double>& x) { ++calls; return f(x); }, rule, d, 1000000, tolerance);
        benchmark::DoNotOptimize(Q);
    }
    state.counters["calls"] = calls;
}

void BM_sparse_legendre_exp(benchmark::State& state)
{
    sparse(state, legendre(), [](const std::vector<double>& x) {
        double s = 0;
        for (std::size_t i = 0; i < x.size(); ++i)
            s += x[i] / (1 << i);
        return std::exp(s);
    }, 1e-8);
}
BENCHMARK(BM_sparse_legendre_exp)->Arg(4)->Arg(8)->Unit(benchmark::kMillisecond);
void BM_sparse_hermite_cos(benchmark::State& state)
{
    sparse(state, hermite(), [](const std::vector<double>& x) {
        double s = 0;
        for (std::size_t i = 0; i < x.size(); ++i)
            s += x[i] / (i + 1);
        return std::cos(s);
    }, 1e-6);
}
BENCHMARK(BM_sparse_hermite_cos)->Arg(4)->Arg(8)->Unit(benchmark::kMillisecond);

} // namespace

BENCHMARK_MAIN();
