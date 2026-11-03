// Copyright Nick Thompson 2026.
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0. (See accompanying file
// LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

/*
 * Google Benchmark suite for wavelet denoising: the transform on the interval (forward + inverse),
 * wavelet_denoise with a fixed wavelet and rule, the cross-validated wavelet_denoise, and denoise,
 * n = 2^10 .. 2^16, plus the odd and non-power-of-two lengths 1000 and 4097. The data is Doppler (Donoho & Johnstone) plus Gaussian noise at SNR 7.
 * The boundary rows are generated once per wavelet and level, before the timing loops.
 *
 * Example:
 *   c++ -O3 -DNDEBUG -std=c++17 -Iinclude \
 *     reporting/performance/wavelet_denoise_performance.cpp \
 *     -lbenchmark -lpthread -o wavelet_denoise_bench
 */

#include <boost/math/filters/wavelet_denoise.hpp>
#include <benchmark/benchmark.h>

#include <cmath>
#include <random>
#include <vector>

namespace {

namespace bf = boost::math::filters;

template <class Real>
std::vector<Real> noisy_doppler(std::size_t n)
{
    std::mt19937_64 gen(12345);
    std::normal_distribution<Real> noise(0, 1);
    std::vector<Real> y(n);
    for (std::size_t i = 0; i < n; ++i)
    {
        const Real t = Real(i) / Real(n);
        y[i] = 7 * 3.3 * std::sqrt(t * (1 - t)) * std::sin(2 * 3.14159265358979 * 1.05 / (t + 0.05)) + noise(gen);
    }
    return y;
}

// state.range(1): 0..8 = db2..db10, 9..17 = sym2..sym10
template <class Real>
bf::orthogonal_wavelet<Real> wavelet_of(benchmark::State& state)
{
    const unsigned t = static_cast<unsigned>(state.range(1));
    return t < 9 ? bf::wavelet::daubechies<Real>(t + 2) : bf::wavelet::symlet<Real>(t - 7);
}

template <class Real>
void BM_transform(benchmark::State& state)
{
    auto x = noisy_doppler<Real>(static_cast<std::size_t>(state.range(0)));
    const auto w = wavelet_of<Real>(state);
    x = bf::wavelet_transform(x, w);
    x = bf::inverse_wavelet_transform(x, w);
    for (auto _ : state)
    {
        x = bf::wavelet_transform(x, w);
        x = bf::inverse_wavelet_transform(x, w);
        benchmark::DoNotOptimize(x.data());
    }
    state.SetComplexityN(state.range(0));
}
BENCHMARK_TEMPLATE(BM_transform, double)->ArgsProduct({{1000, 1024, 4096, 4097, 16384, 65536}, {2, 15, 8}})->Unit(benchmark::kMicrosecond);

// times f(y) on noisy Doppler of state.range(0) samples, after one untimed call
template <class Real, class F>
void time_denoiser(benchmark::State& state, F f)
{
    const auto y = noisy_doppler<Real>(static_cast<std::size_t>(state.range(0)));
    benchmark::DoNotOptimize(f(y).data());
    for (auto _ : state)
        benchmark::DoNotOptimize(f(y).data());
    state.SetComplexityN(state.range(0));
}

template <class Real>
void BM_wavelet_denoise_fixed(benchmark::State& state)
{
    const auto w = wavelet_of<Real>(state);
    const auto rule = static_cast<bf::threshold_rule>(state.range(2));
    time_denoiser<Real>(state, [&](const std::vector<Real>& y) { return bf::wavelet_denoise(y, w, rule); });
}
// sym8 with every rule: wiener, universal_soft, sure, eb_posterior_median
BENCHMARK_TEMPLATE(BM_wavelet_denoise_fixed, double)->ArgsProduct({{1000, 1024, 4096, 4097, 16384, 65536}, {15}, {0, 1, 2, 3}})->Unit(benchmark::kMicrosecond);

template <class Real>
void BM_wavelet_denoise_cv(benchmark::State& state)
{
    time_denoiser<Real>(state, [](const std::vector<Real>& y) { return bf::wavelet_denoise(y); });
}
BENCHMARK_TEMPLATE(BM_wavelet_denoise_cv, float)->RangeMultiplier(4)->Range(1024, 65536)->Unit(benchmark::kMillisecond)->Complexity();
BENCHMARK_TEMPLATE(BM_wavelet_denoise_cv, double)->Arg(1000)->Arg(4097)->Unit(benchmark::kMillisecond);
BENCHMARK_TEMPLATE(BM_wavelet_denoise_cv, double)->RangeMultiplier(4)->Range(1024, 65536)->Unit(benchmark::kMillisecond)->Complexity();

template <class Real>
void BM_denoise(benchmark::State& state)
{
    time_denoiser<Real>(state, [](const std::vector<Real>& y) { return bf::denoise(y); });
}
BENCHMARK_TEMPLATE(BM_denoise, double)->RangeMultiplier(4)->Range(1024, 16384)->Unit(benchmark::kMillisecond)->Complexity();
BENCHMARK_TEMPLATE(BM_denoise, double)->Arg(1000)->Arg(4097)->Unit(benchmark::kMillisecond);

// generating the boundary rows of a fresh wavelet for the levels n needs (both edges)
template <class Real>
void BM_boundary_rows(benchmark::State& state)
{
    const std::size_t n = static_cast<std::size_t>(state.range(0));
    const auto h = bf::detail::wavelet_scaling_filter(wavelet_of<Real>(state));
    for (auto _ : state)
    {
        const auto w = bf::detail::make_orthogonal_wavelet<Real>(h);
        std::vector<Real> x(n, Real(1));
        x = bf::wavelet_transform(x, w);
        benchmark::DoNotOptimize(x.data());
    }
}
BENCHMARK_TEMPLATE(BM_boundary_rows, double)->ArgsProduct({{4096, 4097, 1 << 20, (1 << 20) - 1}, {2, 15, 8}})->Unit(benchmark::kMillisecond);

} // namespace

BENCHMARK_MAIN();
