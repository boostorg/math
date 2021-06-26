#include <benchmark/benchmark.h>
#include <random>
#include <complex>
#include <boost/math/fft.hpp>
#include <boost/math/fft/bsl_backend.hpp>
#include <boost/math/fft/fftw_backend.hpp>
#include <boost/math/fft/gsl_backend.hpp>
#include <boost/math/constants/constants.hpp>

#include "fft_test_helpers.hpp"

std::default_random_engine gen;
std::uniform_real_distribution<double> distribution;


typedef std::complex<double> cd;
std::vector<cd> random_vec(size_t N)
{
    std::vector<cd> V(N);
    for (auto& x : V)
      x = distribution(gen);
    return V;
}

void bench_bsl_dit(benchmark::State& state)
{
    auto A = random_vec(state.range(0));
    for (auto _ : state)
    {
        using boost::math::fft::dft_forward;
        dft_forward<boost::math::fft::test_dft_power2_dit>(A.data(),A.data()+A.size(),A.data());
    }
    state.SetComplexityN(state.range(0));
}
void bench_bsl_dif(benchmark::State& state)
{
    auto A = random_vec(state.range(0));
    for (auto _ : state)
    {
        using boost::math::fft::dft_forward;
        dft_forward<boost::math::fft::test_dft_power2_dif>(A.data(),A.data()+A.size(),A.data());
    }
    state.SetComplexityN(state.range(0));
}
void bench_gsl(benchmark::State& state)
{
    auto A = random_vec(state.range(0));
    for (auto _ : state)
    {
        using boost::math::fft::dft_forward;
        dft_forward<boost::math::fft::gsl_dft>(A.data(),A.data()+A.size(),A.data());
    }
    state.SetComplexityN(state.range(0));
}
void bench_fftw(benchmark::State& state)
{
    auto A = random_vec(state.range(0));
    for (auto _ : state)
    {
        using boost::math::fft::dft_forward;
        dft_forward<boost::math::fft::fftw_dft>(A.data(),A.data()+A.size(),A.data());
    }
    state.SetComplexityN(state.range(0));
}

BENCHMARK(bench_bsl_dit)
    ->RangeMultiplier(8)
    ->Range(1 << 5, 1 << 20)
    ->Complexity(benchmark::oNLogN);

BENCHMARK(bench_bsl_dif)
    ->RangeMultiplier(8)
    ->Range(1 << 5, 1 << 20)
    ->Complexity(benchmark::oNLogN);

BENCHMARK(bench_gsl)
    ->RangeMultiplier(8)
    ->Range(1 << 5, 1 << 20)
    ->Complexity(benchmark::oNLogN);

BENCHMARK(bench_fftw)
    ->RangeMultiplier(8)
    ->Range(1 << 5, 1 << 20)
    ->Complexity(benchmark::oNLogN);

BENCHMARK_MAIN();
