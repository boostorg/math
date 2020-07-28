//  (C) Copyright Nick Thompson 2020.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <benchmark/benchmark.h>
#include <boost/math/constants/constants.hpp>
#include <boost/multiprecision/mpfr.hpp>

using namespace boost::math::constants;
using boost::multiprecision::mpfr_float;
#include <boost/math/special_functions/expint.hpp>

using namespace boost::math::constants;
using boost::multiprecision::mpfr_float;
// Ei(1) = 1.8951178163559367554665209343316342690170605817327075916462284318825138345338041535489007101261\
             3895697181109531794465374258814916416306468808818668253882866963233854509522755525848139221216\
             6459936359948543306285455761625228166868118802856637846665686888646424297019090790472890309099\
             3380190917469997918302494807585852088867837050413737878407416460258387628362162273159751506105\
             1925474727680703251963132680675670740814722824224269686684164320795874672671226770927559887547\
             77656819174952725603954608115985285493707
void Ei1Boost(benchmark::State& state)
{
    
    using boost::math::expint;
    mpfr_float::default_precision(state.range(0));
    mpfr_float x{1};
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(expint(x));
    }
    //std::cout <<std::setprecision(state.range(0)) << "Expint(1) = " << expint(x) << "\n";
    state.SetComplexityN(state.range(0));
}

BENCHMARK(Ei1Boost)->RangeMultiplier(2)->Range(32, 65536)->Complexity()->Unit(benchmark::kMicrosecond);

template<typename Real>
Real Ei1() {
    Real sum = euler<Real>();
    Real k = 1;
    Real fact = 1;
    Real term = 1/fact;
    while (term > std::numeric_limits<Real>::epsilon()) {
        sum += term;
        ++k;
        fact *= k;
        term = 1/(k*fact);
    }
    return sum;
}

void Ei1Nick(benchmark::State& state)
{   
    using boost::math::expint;
    mpfr_float::default_precision(state.range(0));
    mpfr_float x{1};
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(Ei1<mpfr_float>());
    }
    //std::cout << std::setprecision(state.range(0)) << "Expint(1) = " << Ei1<mpfr_float>() << "\n";
    state.SetComplexityN(state.range(0));
}

BENCHMARK(Ei1Nick)->RangeMultiplier(2)->Range(32, 65536)->Complexity()->Unit(benchmark::kMicrosecond);

template<typename Real>
Real Ei1Ramanujan()
{
    using std::pow;
    using std::abs;
    Real term = std::numeric_limits<Real>::max();
    size_t n = 1;
    Real ksum = 0;
    Real fact = 1;
    Real sum = 0;
    while (abs(term) > std::numeric_limits<Real>::epsilon()) {
        fact *= n;
        Real denom = fact;
        // Use shifts if possible:
        /*if (n < 64)
        {
            denom *= (1uLL << n - 1);
        }
        else
        {
            denom *= pow(Real(2), n - 1);
        }
        Real scale = -Real(1)/denom;*/
        Real scale = -Real(1)/(fact*(1uLL << n - 1));
        if (n & 1)
        {
            ksum += 1/Real(n);
            scale = -scale;
        }
        term = scale*ksum;
        sum += term;
        ++n;
    }
    return euler<Real>() + root_e<Real>()*sum;
}

void Ei1Ramanujan(benchmark::State& state)
{
    mpfr_float::default_precision(state.range(0));
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(Ei1Ramanujan<mpfr_float>());
    }
    //std::cout << std::setprecision(state.range(0)) << "Expint(1) = " << Ei1Ramanujan<mpfr_float>() << "\n";
    state.SetComplexityN(state.range(0));
}

BENCHMARK(Ei1Ramanujan)->RangeMultiplier(2)->Range(32, 65536)->Complexity()->Unit(benchmark::kMicrosecond);



void LaplaceLimit(benchmark::State& state)
{
    mpfr_float::default_precision(state.range(0));
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(laplace_limit<mpfr_float>());
    }
    state.SetComplexityN(state.range(0));
}

BENCHMARK(LaplaceLimit)->RangeMultiplier(2)->Range(128, 1<<20)->Complexity()->Unit(benchmark::kMicrosecond);

void Dottie(benchmark::State& state)
{
    mpfr_float::default_precision(state.range(0));
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(dottie<mpfr_float>());
    }
    state.SetComplexityN(state.range(0));
}

BENCHMARK(Dottie)->RangeMultiplier(2)->Range(512, 1<<20)->Complexity()->Unit(benchmark::kMicrosecond);

void ReciprocalFibonacci(benchmark::State& state)
{
    mpfr_float::default_precision(state.range(0));
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(reciprocal_fibonacci<mpfr_float>());
    }
    state.SetComplexityN(state.range(0));
}

BENCHMARK(ReciprocalFibonacci)->RangeMultiplier(2)->Range(512, 1<<20)->Complexity()->Unit(benchmark::kMicrosecond);


void Pi(benchmark::State& state)
{
    mpfr_float::default_precision(state.range(0));
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(pi<mpfr_float>());
    }
    state.SetComplexityN(state.range(0));
}

BENCHMARK(Pi)->RangeMultiplier(2)->Range(512, 1<<20)->Complexity()->Unit(benchmark::kMicrosecond);

void Gauss(benchmark::State& state)
{
    mpfr_float::default_precision(state.range(0));
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(gauss<mpfr_float>());
    }
    state.SetComplexityN(state.range(0));
}

BENCHMARK(Gauss)->RangeMultiplier(2)->Range(512, 1<<20)->Complexity()->Unit(benchmark::kMicrosecond);

void Exp1(benchmark::State& state)
{
    mpfr_float::default_precision(state.range(0));
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(e<mpfr_float>());
    }
    state.SetComplexityN(state.range(0));
}

BENCHMARK(Exp1)->RangeMultiplier(2)->Range(512, 1<<20)->Complexity()->Unit(benchmark::kMicrosecond);

void Catalan(benchmark::State& state)
{
    mpfr_float::default_precision(state.range(0));
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(catalan<mpfr_float>());
    }
    state.SetComplexityN(state.range(0));
}

BENCHMARK(Catalan)->RangeMultiplier(2)->Range(512, 1<<20)->Complexity()->Unit(benchmark::kMicrosecond);

void Plastic(benchmark::State& state)
{
    mpfr_float::default_precision(state.range(0));
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(plastic<mpfr_float>());
    }
    state.SetComplexityN(state.range(0));
}

BENCHMARK(Plastic)->RangeMultiplier(2)->Range(512, 1<<20)->Complexity()->Unit(benchmark::kMicrosecond);

void RootTwo(benchmark::State& state)
{
    mpfr_float::default_precision(state.range(0));
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(root_two<mpfr_float>());
    }
    state.SetComplexityN(state.range(0));
}

BENCHMARK(RootTwo)->RangeMultiplier(2)->Range(512, 1<<20)->Complexity()->Unit(benchmark::kMicrosecond);

void ZetaThree(benchmark::State& state)
{
    mpfr_float::default_precision(state.range(0));
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(zeta_three<mpfr_float>());
    }
    state.SetComplexityN(state.range(0));
}

BENCHMARK(ZetaThree)->RangeMultiplier(2)->Range(512, 1<<20)->Complexity()->Unit(benchmark::kMicrosecond);


void Euler(benchmark::State& state)
{
    mpfr_float::default_precision(state.range(0));
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(euler<mpfr_float>());
    }
    state.SetComplexityN(state.range(0));
}

BENCHMARK(Euler)->RangeMultiplier(2)->Range(512, 1<<20)->Complexity()->Unit(benchmark::kMicrosecond);


void LnTwo(benchmark::State& state)
{
    mpfr_float::default_precision(state.range(0));
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(ln_two<mpfr_float>());
    }
    state.SetComplexityN(state.range(0));
}

BENCHMARK(LnTwo)->RangeMultiplier(2)->Range(512, 1<<20)->Complexity()->Unit(benchmark::kMicrosecond);

void Glaisher(benchmark::State& state)
{
    mpfr_float::default_precision(state.range(0));
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(glaisher<mpfr_float>());
    }
    state.SetComplexityN(state.range(0));
}

BENCHMARK(Glaisher)->RangeMultiplier(2)->Range(512, 4096)->Complexity()->Unit(benchmark::kMicrosecond);


void Khinchin(benchmark::State& state)
{
    mpfr_float::default_precision(state.range(0));
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(khinchin<mpfr_float>());
    }
    state.SetComplexityN(state.range(0));
}

// There is a performance bug in the Khinchin constant:
BENCHMARK(Khinchin)->RangeMultiplier(2)->Range(512, 512)->Complexity()->Unit(benchmark::kMicrosecond);

BENCHMARK_MAIN();
