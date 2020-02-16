/*
 * Copyright Nick Thompson, 2020
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */

#include <cmath>
#include <random>
#include <benchmark/benchmark.h>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/math/special_functions/daubechies_scaling.hpp>
#include <boost/math/interpolators/detail/cardinal_cubic_hermite_detail.hpp>
#include <boost/math/interpolators/detail/cardinal_quintic_hermite_detail.hpp>
#include <boost/math/interpolators/detail/septic_hermite_detail.hpp>

static void UnitStep(benchmark::internal::Benchmark* b)
{
    for (int i = 7; i <= 25; ++i)
    {
        b->Args({i});
    }
}

double exponential(benchmark::IterationCount j)
{
    return std::pow(2, j);
}


template<typename Real, int p>
void DyadicGrid(benchmark::State & state)
{
    int j = state.range(0);
    size_t s = 0;
    for (auto _ : state)
    {
        auto v = boost::math::detail::dyadic_grid<Real, 4, 0>(j);
        benchmark::DoNotOptimize(v[0]);
        s = v.size();
    }

    state.counters["RAM"] = s*sizeof(Real);
    state.SetComplexityN(state.range(0));
}

BENCHMARK_TEMPLATE(DyadicGrid, double, 4)->Apply(UnitStep)->Unit(benchmark::kMillisecond)->Complexity(exponential);
//BENCHMARK_TEMPLATE(DyadicGrid, double, 8)->Apply(UnitStep)->Unit(benchmark::kMillisecond)->Complexity(exponential);
//BENCHMARK_TEMPLATE(DyadicGrid, double, 11)->Apply(UnitStep)->Unit(benchmark::kMillisecond)->Complexity(exponential);


template<typename Real, int p>
void ScalingEvaluation(benchmark::State & state)
{
    auto phi = boost::math::daubechies_scaling<Real, p>();
    Real x = 0;
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(phi(x));
        x += std::numeric_limits<Real>::epsilon();
    }
}


BENCHMARK_TEMPLATE(ScalingEvaluation, double, 2);
BENCHMARK_TEMPLATE(ScalingEvaluation, double, 3);
BENCHMARK_TEMPLATE(ScalingEvaluation, double, 4);
BENCHMARK_TEMPLATE(ScalingEvaluation, double, 5);
BENCHMARK_TEMPLATE(ScalingEvaluation, double, 6);
BENCHMARK_TEMPLATE(ScalingEvaluation, double, 7);
BENCHMARK_TEMPLATE(ScalingEvaluation, double, 8);
BENCHMARK_TEMPLATE(ScalingEvaluation, double, 9);
BENCHMARK_TEMPLATE(ScalingEvaluation, double, 10);
BENCHMARK_TEMPLATE(ScalingEvaluation, double, 11);


template<typename Real, int p>
void ScalingConstructor(benchmark::State & state)
{
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(boost::math::daubechies_scaling<Real, p>());
    }
}

BENCHMARK_TEMPLATE(ScalingConstructor, float, 2)->Unit(benchmark::kMillisecond);
BENCHMARK_TEMPLATE(ScalingConstructor, double, 2)->Unit(benchmark::kMillisecond);
BENCHMARK_TEMPLATE(ScalingConstructor, long double, 2)->Unit(benchmark::kMillisecond);

BENCHMARK_TEMPLATE(ScalingConstructor, float, 3)->Unit(benchmark::kMillisecond);
BENCHMARK_TEMPLATE(ScalingConstructor, double, 3)->Unit(benchmark::kMillisecond);
BENCHMARK_TEMPLATE(ScalingConstructor, long double, 3)->Unit(benchmark::kMillisecond);

BENCHMARK_TEMPLATE(ScalingConstructor, float, 4)->Unit(benchmark::kMillisecond);
BENCHMARK_TEMPLATE(ScalingConstructor, double, 4)->Unit(benchmark::kMillisecond);
BENCHMARK_TEMPLATE(ScalingConstructor, long double, 4)->Unit(benchmark::kMillisecond);

BENCHMARK_TEMPLATE(ScalingConstructor, float, 5)->Unit(benchmark::kMillisecond);
BENCHMARK_TEMPLATE(ScalingConstructor, double, 5)->Unit(benchmark::kMillisecond);
BENCHMARK_TEMPLATE(ScalingConstructor, long double, 5)->Unit(benchmark::kMillisecond);

BENCHMARK_TEMPLATE(ScalingConstructor, float, 11)->Unit(benchmark::kMillisecond);
BENCHMARK_TEMPLATE(ScalingConstructor, double, 11)->Unit(benchmark::kMillisecond);
BENCHMARK_TEMPLATE(ScalingConstructor, long double, 11)->Unit(benchmark::kMillisecond);

template<typename Real>
void CardinalCubicHermite(benchmark::State & state)
{
    using boost::math::interpolators::detail::cardinal_cubic_hermite_detail;
    auto n = state.range(0);
    std::vector<Real> y(n);
    std::vector<Real> dydx(n);
    std::random_device rd;
    boost::random::uniform_real_distribution<Real> dis(Real(0), Real(1));
    for (size_t i = 0; i < y.size(); ++i)
    {
        y[i] = dis(rd);
        dydx[i] = dis(rd);
    }

    Real dx = Real(1)/Real(8);
    Real x0 = 0;
    Real xf = x0 + (y.size()-1)*dx;

    auto qh = cardinal_cubic_hermite_detail(std::move(y), std::move(dydx), x0, dx);
    Real x = 0;
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(qh.unchecked_evaluation(x));
        x += xf/128;
        if (x >= xf)
        {
            x = x0;
        }
    }
}

BENCHMARK_TEMPLATE(CardinalCubicHermite, double)->RangeMultiplier(2)->Range(1<<8, 1<<20)->Complexity();

template<typename Real>
void CardinalCubicHermiteAOS(benchmark::State & state)
{
    auto n = state.range(0);
    std::vector<std::array<Real, 2>> dat(n);
    std::random_device rd;
    boost::random::uniform_real_distribution<Real> dis(Real(0), Real(1));
    for (size_t i = 0; i < dat.size(); ++i)
    {
        dat[i][0] = dis(rd);
        dat[i][1] = dis(rd);
    }

    using boost::math::interpolators::detail::cardinal_cubic_hermite_detail_aos;
    Real dx = Real(1)/Real(8);
    Real x0 = 0;
    Real xf = x0 + (dat.size()-1)*dx;
    auto qh = cardinal_cubic_hermite_detail_aos(std::move(dat), x0, dx);
    Real x = 0;
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(qh.unchecked_evaluation(x));
        x += xf/128;
        if (x >= xf)
        {
            x = x0;
        }
    }
}

BENCHMARK_TEMPLATE(CardinalCubicHermiteAOS, double)->RangeMultiplier(2)->Range(1<<8, 1<<20)->Complexity();

template<class Real>
void SineEvaluation(benchmark::State& state)
{
    std::default_random_engine gen;
    std::uniform_real_distribution<Real> x_dis(0, 3.14159);

    Real x = x_dis(gen);
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(std::sin(x));
        x += std::numeric_limits<Real>::epsilon();
    }
}

BENCHMARK_TEMPLATE(SineEvaluation, float);
BENCHMARK_TEMPLATE(SineEvaluation, double);
BENCHMARK_TEMPLATE(SineEvaluation, long double);

template<class Real>
void ExpEvaluation(benchmark::State& state)
{
    std::default_random_engine gen;
    std::uniform_real_distribution<Real> x_dis(0, 3.14159);

    Real x = x_dis(gen);
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(std::exp(x));
        x += std::numeric_limits<Real>::epsilon();
    }
}

BENCHMARK_TEMPLATE(ExpEvaluation, float);
BENCHMARK_TEMPLATE(ExpEvaluation, double);
BENCHMARK_TEMPLATE(ExpEvaluation, long double);

template<class Real>
void PowEvaluation(benchmark::State& state)
{
    std::default_random_engine gen;
    std::uniform_real_distribution<Real> x_dis(0, 3.14159);

    Real x = x_dis(gen);
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(std::pow(x, x+1));
        x += std::numeric_limits<Real>::epsilon();
    }
}

BENCHMARK_TEMPLATE(PowEvaluation, float);
BENCHMARK_TEMPLATE(PowEvaluation, double);
BENCHMARK_TEMPLATE(PowEvaluation, long double);


template<typename Real>
void CardinalQuinticHermite(benchmark::State & state)
{
    using boost::math::interpolators::detail::cardinal_quintic_hermite_detail;
    auto n = state.range(0);
    std::vector<Real> y(n);
    std::vector<Real> dydx(n);
    std::vector<Real> d2ydx2(n);
    std::random_device rd;
    boost::random::uniform_real_distribution<Real> dis(Real(0), Real(1));
    for (size_t i = 0; i < y.size(); ++i)
    {
        y[i] = dis(rd);
        dydx[i] = dis(rd);
        d2ydx2[i] = dis(rd);
    }

    Real dx = Real(1)/Real(8);
    Real x0 = 0;
    Real xf = x0 + (y.size()-1)*dx;

    auto qh = cardinal_quintic_hermite_detail(std::move(y), std::move(dydx), std::move(d2ydx2), x0, dx);
    Real x = 0;
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(qh.unchecked_evaluation(x));
        x += xf/128;
        if (x >= xf)
        {
            x = x0;
        }
    }
}

BENCHMARK_TEMPLATE(CardinalQuinticHermite, float)->RangeMultiplier(2)->Range(1<<8, 1<<22)->Complexity();
BENCHMARK_TEMPLATE(CardinalQuinticHermite, double)->RangeMultiplier(2)->Range(1<<8, 1<<22)->Complexity();
BENCHMARK_TEMPLATE(CardinalQuinticHermite, long double)->RangeMultiplier(2)->Range(1<<8, 1<<22)->Complexity();


template<typename Real>
void CardinalQuinticHermiteAOS(benchmark::State & state)
{
    auto n = state.range(0);
    std::vector<std::array<Real, 3>> dat(n);
    std::random_device rd;
    boost::random::uniform_real_distribution<Real> dis(Real(0), Real(1));
    for (size_t i = 0; i < dat.size(); ++i)
    {
        dat[i][0] = dis(rd);
        dat[i][1] = dis(rd);
        dat[i][2] = dis(rd);
    }

    using boost::math::interpolators::detail::cardinal_quintic_hermite_detail_aos;
    Real dx = Real(1)/Real(8);
    Real x0 = 0;
    Real xf = x0 + (dat.size()-1)*dx;
    auto qh = cardinal_quintic_hermite_detail_aos(std::move(dat), x0, dx);
    Real x = 0;
    for (auto _ : state) {
        benchmark::DoNotOptimize(qh.unchecked_evaluation(x));
        x += xf/128;
        if (x >= xf)
        {
            x = x0;
        }
    }
}

BENCHMARK_TEMPLATE(CardinalQuinticHermiteAOS, float)->RangeMultiplier(2)->Range(1<<8, 1<<22)->Complexity();
BENCHMARK_TEMPLATE(CardinalQuinticHermiteAOS, double)->RangeMultiplier(2)->Range(1<<8, 1<<22)->Complexity();
BENCHMARK_TEMPLATE(CardinalQuinticHermiteAOS, long double)->RangeMultiplier(2)->Range(1<<8, 1<<22)->Complexity();


template<typename Real>
void CardinalSepticHermite(benchmark::State & state)
{
    using boost::math::interpolators::detail::cardinal_septic_hermite_detail;
    auto n = state.range(0);
    std::vector<Real> y(n);
    std::vector<Real> dydx(n);
    std::vector<Real> d2ydx2(n);
    std::vector<Real> d3ydx3(n);
    std::random_device rd;
    boost::random::uniform_real_distribution<Real> dis(Real(0), Real(1));
    for (size_t i = 0; i < y.size(); ++i)
    {
        y[i] = dis(rd);
        dydx[i] = dis(rd);
        d2ydx2[i] = dis(rd);
        d3ydx3[i] = dis(rd);
    }

    Real dx = Real(1)/Real(8);
    Real x0 = 0;
    Real xf = x0 + (y.size()-1)*dx;

    auto sh = cardinal_septic_hermite_detail(std::move(y), std::move(dydx), std::move(d2ydx2), std::move(d3ydx3), x0, dx);
    Real x = 0;
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(sh.unchecked_evaluation(x));
        x += xf/128;
        if (x >= xf)
        {
            x = x0;
        }
    }
}

BENCHMARK_TEMPLATE(CardinalSepticHermite, double)->RangeMultiplier(2)->Range(1<<8, 1<<22)->Complexity();


BENCHMARK_MAIN();