// Copyright Nick Thompson 2026.
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0. (See accompanying file
// LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

/*
 * Google Benchmark suite comparing L-BFGS with trust-region Newton: a whole minimization of the n-dimensional
 * Rosenbrock function from (-1.2, 1, -1.2, 1, ...) until the gradient norm is below 1e-8. The counters are the
 * steps, the objective calls and whether the run converged.
 *
 * Example:
 *   c++ -O2 -DNDEBUG -std=c++14 -Iinclude \
 *     reporting/performance/trust_region_newton_performance.cpp \
 *     -lbenchmark -lpthread -o trust_region_newton_bench
 */

#include <boost/math/differentiation/autodiff_reverse.hpp>
#include <boost/math/optimization/lbfgs.hpp>
#include <boost/math/optimization/minimizer.hpp>
#include <boost/math/optimization/trust_region_newton.hpp>
#include <benchmark/benchmark.h>

#include <cstddef>
#include <vector>

namespace {

namespace rdiff = boost::math::differentiation::reverse_mode;
namespace bopt = boost::math::optimization;

std::size_t calls = 0;

template <class Real>
Real rosenbrock(std::vector<Real> const& x)
{
    ++calls;
    Real f = Real(0);
    for (std::size_t i = 0; i + 1 < x.size(); ++i)
    {
        f += 100 * (x[i + 1] - x[i] * x[i]) * (x[i + 1] - x[i] * x[i]) + (1 - x[i]) * (1 - x[i]);
    }
    return f;
}

std::vector<double> start(std::size_t n)
{
    std::vector<double> x0(n);
    for (std::size_t i = 0; i < n; ++i)
    {
        x0[i] = i % 2 ? 1.0 : -1.2;
    }
    return x0;
}

template <std::size_t Order, class MakeOptimizer>
void run(benchmark::State& state, MakeOptimizer make)
{
    using rvar = rdiff::rvar<double, Order>;
    auto const x0 = start(static_cast<std::size_t>(state.range(0)));
    bopt::optimization_result<double> r;
    for (auto _ : state)
    {
        // Untimed: an order 2 tape's buffers hold leaves of the order 1 tape, so clear the order 1 tape last.
        state.PauseTiming();
        rdiff::get_active_tape<double, Order>().clear();
        rdiff::get_active_tape<double, 1>().clear();
        state.ResumeTiming();
        calls = 0;
        std::vector<rvar> x(x0.begin(), x0.end());
        auto opt = make(x);
        r = bopt::minimize(opt, bopt::unconstrained_policy<std::vector<rvar>>{},
                           bopt::gradient_norm_convergence_policy<double>(1e-8), bopt::max_iter_termination_policy(20000),
                           false);
        benchmark::DoNotOptimize(x[0].item());
    }
    state.counters["steps"] = static_cast<double>(r.num_iter);
    state.counters["calls"] = static_cast<double>(calls);
    state.counters["converged"] = r.converged;
}

void BM_lbfgs_rosenbrock(benchmark::State& state)
{
    run<1>(state, [](std::vector<rdiff::rvar<double, 1>>& x) {
        return bopt::make_lbfgs(&rosenbrock<rdiff::rvar<double, 1>>, x);
    });
}

void BM_trust_region_newton_rosenbrock(benchmark::State& state)
{
    run<2>(state, [](std::vector<rdiff::rvar<double, 2>>& x) {
        return bopt::make_trust_region_newton(&rosenbrock<rdiff::rvar<double, 2>>, x);
    });
}

BENCHMARK(BM_lbfgs_rosenbrock)->Arg(2)->Arg(10)->Arg(100)->Unit(benchmark::kMillisecond);
BENCHMARK(BM_trust_region_newton_rosenbrock)->Arg(2)->Arg(10)->Arg(100)->Unit(benchmark::kMillisecond);

} // namespace

BENCHMARK_MAIN();
