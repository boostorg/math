//  (C) Copyright Sayan Samanta 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Build: c++ -std=c++17 -O2 -I../../include autodiff_performance.cpp -lbenchmark -o autodiff_performance

#include <cmath>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>
#include <benchmark/benchmark.h>
#include <boost/math/constants/constants.hpp>
#include <boost/math/differentiation/autodiff.hpp>
#include <boost/math/tools/random_vector.hpp>
#include <boost/math/tools/rational.hpp>

using boost::math::constants::one_div_root_two;
using boost::math::differentiation::make_ftuple;
using boost::math::differentiation::make_fvar;
using boost::math::generate_random_uniform_vector;
using boost::math::tools::evaluate_polynomial;

template<class Real>
Real transcendental(const Real& x)
{
    using std::exp;
    using std::log;
    using std::sin;
    using std::sqrt;
    return exp(sin(x)) * log(1.0 + x * x) / sqrt(x);
}

template<class X, class Y>
auto two_variables(const X& x, const Y& y)
{
    return exp(sin(x) * y) * log(1.0 + x * x + y * y) / sqrt(x * y);
}

template<class X, class Y, class Z>
auto three_variables(const X& x, const Y& y, const Z& z)
{
    return exp(sin(x) * y) * log(1.0 + x * x + z * z) / sqrt(y * z) + cos(x * z) * erfc(y);
}

// From example/autodiff_black_scholes.cpp.
template<class Price, class Sigma, class Tau, class Rate>
auto black_scholes(const Price& S, const Sigma& sigma, const Tau& tau, const Rate& r)
{
    using std::exp;
    using std::sqrt;
    constexpr double K = 100.0;
    const auto sst = sigma * sqrt(tau);
    const auto d1 = (log(S / K) + (r + sigma * sigma / 2) * tau) / sst;
    const auto d2 = d1 - sst;
    auto Phi = [](const auto& z) { return 0.5 * erfc(-one_div_root_two<double>() * z); };
    return S * Phi(d1) - exp(-r * tau) * K * Phi(d2);
}

// Excess Gibbs energy of a mixed-solvent electrolyte (MSE model) as a function of the
// pressure P in MPa: a Pitzer-Debye-Huckel long-range term, through the solvent density
// and dielectric constant, and a middle-range ionic term. dG/dP is the excess volume.
template<class Real>
Real mse_gibbs_energy(const Real& P, const std::vector<double>& x, const std::vector<double>& z,
                      const std::vector<double>& b0, const std::vector<double>& b1)
{
    using std::exp;
    using std::log;
    using std::pow;
    using std::sqrt;
    const std::size_t n = x.size();
    double I = 0;
    for (std::size_t i = 0; i < n; ++i)
    {
        I += 0.5 * z[i] * z[i] * x[i];
    }
    constexpr double rho = 14.9;
    const Real density = exp(4.5e-4 * (P - 0.1));
    const Real dielectric = 1.0 + 0.0216 * log((300.0 + P) / 300.1);
    const Real Ax = 2.917 * sqrt(density) * pow(dielectric, -1.5);
    Real g = -4.0 * I / rho * std::log(1.0 + rho * std::sqrt(I)) * Ax;
    const double e = std::exp(-std::sqrt(I + 0.01));
    for (std::size_t i = 0; i < n; ++i)
    {
        for (std::size_t j = 0; j < n; ++j)
        {
            g += x[i] * x[j] * e * (b0[i * n + j] + b1[i * n + j] * P);
        }
    }
    return g;
}

template<std::size_t Order, class F>
void run(benchmark::State& state, double x, F f)
{
    const auto v = make_fvar<double, Order>(x);
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(f(v));
    }
}

// Runs f on a variable of the order given by state.range(0); the order is the N of the complexity fit.
template<class F, std::size_t... Orders>
void run_at_order(benchmark::State& state, double x, F f, std::index_sequence<Orders...>)
{
    const auto order = static_cast<std::size_t>(state.range(0));
    ((order == Orders + 1 ? run<Orders + 1>(state, x, f) : void()), ...);
    state.SetComplexityN(state.range(0));
}

template<class F>
void run_at_order(benchmark::State& state, double x, F f)
{
    run_at_order(state, x, f, std::make_index_sequence<32>());
}

// Calls f(std::integral_constant<std::size_t, N>()), N = state.range(0): the order in each variable.
template<class F, std::size_t... Orders>
void run_at_order(benchmark::State& state, F f, std::index_sequence<Orders...>)
{
    const auto order = static_cast<std::size_t>(state.range(0));
    ((order == Orders + 1 ? f(std::integral_constant<std::size_t, Orders + 1>()) : void()), ...);
}

void Horner(benchmark::State& state)
{
    const auto c = generate_random_uniform_vector<double>(64, 12345, -1, 1);
    run_at_order(state, 0.5, [&](const auto& x) { return evaluate_polynomial(c.data(), x, c.size()); });
}

void Transcendental(benchmark::State& state)
{
    run_at_order(state, 2.5, [](const auto& x) { return transcendental(x); });
}

void BlackScholes(benchmark::State& state)
{
    run_at_order(state, 105.0, [](const auto& S) { return black_scholes(S, 0.25, 0.5, 0.01); });
}

void TwoVariables(benchmark::State& state)
{
    auto f = [&](auto order) {
        constexpr std::size_t N = decltype(order)::value;
        const auto v = make_ftuple<double, N, N>(1.3, 0.7);
        for (auto _ : state)
        {
            benchmark::DoNotOptimize(two_variables(std::get<0>(v), std::get<1>(v)));
        }
    };
    run_at_order(state, f, std::make_index_sequence<12>());
}

void ThreeVariables(benchmark::State& state)
{
    auto f = [&](auto order) {
        constexpr std::size_t N = decltype(order)::value;
        const auto v = make_ftuple<double, N, N, N>(1.3, 0.7, 2.1);
        for (auto _ : state)
        {
            benchmark::DoNotOptimize(three_variables(std::get<0>(v), std::get<1>(v), std::get<2>(v)));
        }
    };
    run_at_order(state, f, std::make_index_sequence<6>());
}

// The price as a function of (S, sigma, tau, r).
void BlackScholesFourVariables(benchmark::State& state)
{
    auto f = [&](auto order) {
        constexpr std::size_t N = decltype(order)::value;
        const auto v = make_ftuple<double, N, N, N, N>(105.0, 0.25, 0.5, 0.01);
        for (auto _ : state)
        {
            benchmark::DoNotOptimize(
                black_scholes(std::get<0>(v), std::get<1>(v), std::get<2>(v), std::get<3>(v)));
        }
    };
    run_at_order(state, f, std::make_index_sequence<4>());
}

// 16 species: a solvent, alternating cations and anions, and a neutral solute.
void MseGibbsEnergy(benchmark::State& state)
{
    constexpr std::size_t n = 16;
    std::vector<double> x(n, 0.2 / (n - 1));
    std::vector<double> z(n, 0.0);
    x[0] = 0.8;
    for (std::size_t i = 1; i + 1 < n; ++i)
    {
        z[i] = (i % 2) ? 1 : -1;
    }
    const auto b0 = generate_random_uniform_vector<double>(n * n, 12345, -0.5, 0.5);
    const auto b1 = generate_random_uniform_vector<double>(n * n, 12345, -5e-4, 5e-4);
    run_at_order(state, 0.1, [&](const auto& P) { return mse_gibbs_energy(P, x, z, b0, b1); });
}

BENCHMARK(Horner)->DenseRange(1, 32)->Complexity();
BENCHMARK(Transcendental)->DenseRange(1, 32)->Complexity();
BENCHMARK(BlackScholes)->DenseRange(1, 32)->Complexity();
BENCHMARK(MseGibbsEnergy)->DenseRange(1, 32)->Complexity();
BENCHMARK(TwoVariables)->DenseRange(1, 12);
BENCHMARK(ThreeVariables)->DenseRange(1, 6);
BENCHMARK(BlackScholesFourVariables)->DenseRange(1, 4);

BENCHMARK_MAIN();
