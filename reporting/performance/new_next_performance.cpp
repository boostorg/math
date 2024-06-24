//  (C) Copyright Matt Borland 2022.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/math/special_functions/next.hpp>
#include <benchmark/benchmark.h>

template <typename T>
void float_distance(benchmark::State& state)
{
    const auto difference = static_cast<int>(state.range(0));
    T left = 2;
    T right = boost::math::float_advance(left, difference);

    for (auto _ : state)
    {
        benchmark::DoNotOptimize(boost::math::float_distance(left, right));
    }
    state.SetComplexityN(state.range(0));
}

BENCHMARK_TEMPLATE(float_distance, float)->RangeMultiplier(2)->Range(1 << 1, 1 << 14)->Complexity()->UseRealTime();
BENCHMARK_TEMPLATE(float_distance, double)->RangeMultiplier(2)->Range(1 << 1, 1 << 14)->Complexity()->UseRealTime();

BENCHMARK_MAIN();

/*
Run on Apple M1 Pro Arch using Apple Clang 14.0.0 (15OCT22)

Original performance (Boost 1.80.0):

Unable to determine clock rate from sysctl: hw.cpufrequency: No such file or directory
This does not affect benchmark measurements, only the metadata output.
2022-10-15T15:24:07-07:00
Running ./new_next_performance
Run on (10 X 24.0916 MHz CPU s)
CPU Caches:
  L1 Data 64 KiB
  L1 Instruction 128 KiB
  L2 Unified 4096 KiB (x10)
Load Average: 1.86, 2.53, 5.83
---------------------------------------------------------------------------------
Benchmark                                       Time             CPU   Iterations
---------------------------------------------------------------------------------
float_distance<float>/2/real_time            61.4 ns         61.4 ns      9074469
float_distance<float>/4/real_time            61.7 ns         61.7 ns     11384150
float_distance<float>/8/real_time            61.4 ns         61.4 ns     10814604
float_distance<float>/16/real_time           61.7 ns         61.7 ns     11348376
float_distance<float>/32/real_time           61.4 ns         61.4 ns     11387167
float_distance<float>/64/real_time           61.6 ns         61.6 ns     11131932
float_distance<float>/128/real_time          61.4 ns         61.4 ns     11382029
float_distance<float>/256/real_time          61.4 ns         61.4 ns     11307649
float_distance<float>/512/real_time          61.4 ns         61.4 ns     11376048
float_distance<float>/1024/real_time         61.4 ns         61.4 ns     11355748
float_distance<float>/2048/real_time         61.8 ns         61.8 ns     11373776
float_distance<float>/4096/real_time         61.4 ns         61.4 ns     11382368
float_distance<float>/8192/real_time         61.4 ns         61.4 ns     11353453
float_distance<float>/16384/real_time        61.4 ns         61.4 ns     11378298
float_distance<float>/real_time_BigO        61.48 (1)       61.47 (1)
float_distance<float>/real_time_RMS             0 %             0 %
float_distance<double>/2/real_time           55.6 ns         55.6 ns     12580218
float_distance<double>/4/real_time           55.6 ns         55.6 ns     12577835
float_distance<double>/8/real_time           55.6 ns         55.6 ns     12564909
float_distance<double>/16/real_time          56.2 ns         56.2 ns     12554909
float_distance<double>/32/real_time          56.0 ns         56.0 ns     12544381
float_distance<double>/64/real_time          55.6 ns         55.6 ns     12566488
float_distance<double>/128/real_time         55.6 ns         55.6 ns     12499581
float_distance<double>/256/real_time         55.6 ns         55.6 ns     12565661
float_distance<double>/512/real_time         56.1 ns         56.1 ns     12550023
float_distance<double>/1024/real_time        55.8 ns         55.8 ns     12568603
float_distance<double>/2048/real_time        55.6 ns         55.6 ns     12546049
float_distance<double>/4096/real_time        55.6 ns         55.6 ns     12528525
float_distance<double>/8192/real_time        55.9 ns         55.9 ns     12563030
float_distance<double>/16384/real_time       56.0 ns         56.0 ns     12447644
float_distance<double>/real_time_BigO       55.78 (1)       55.78 (1)
float_distance<double>/real_time_RMS            0 %             0 %

New performance:

Unable to determine clock rate from sysctl: hw.cpufrequency: No such file or directory
This does not affect benchmark measurements, only the metadata output.
2022-10-15T15:31:37-07:00
Running ./new_next_performance
Run on (10 X 24.122 MHz CPU s)
CPU Caches:
  L1 Data 64 KiB
  L1 Instruction 128 KiB
  L2 Unified 4096 KiB (x10)
Load Average: 2.12, 2.17, 4.26
---------------------------------------------------------------------------------
Benchmark                                       Time             CPU   Iterations
---------------------------------------------------------------------------------
float_distance<float>/2/real_time            15.8 ns         15.8 ns     42162717
float_distance<float>/4/real_time            15.9 ns         15.9 ns     44213877
float_distance<float>/8/real_time            15.8 ns         15.8 ns     43972542
float_distance<float>/16/real_time           15.8 ns         15.8 ns     44209456
float_distance<float>/32/real_time           15.8 ns         15.8 ns     44200244
float_distance<float>/64/real_time           15.8 ns         15.8 ns     44239293
float_distance<float>/128/real_time          15.8 ns         15.8 ns     44171202
float_distance<float>/256/real_time          15.8 ns         15.8 ns     44241507
float_distance<float>/512/real_time          15.9 ns         15.8 ns     44230034
float_distance<float>/1024/real_time         15.8 ns         15.8 ns     44241554
float_distance<float>/2048/real_time         15.8 ns         15.8 ns     44220802
float_distance<float>/4096/real_time         15.8 ns         15.8 ns     44220441
float_distance<float>/8192/real_time         15.9 ns         15.9 ns     44213994
float_distance<float>/16384/real_time        15.8 ns         15.8 ns     44215413
float_distance<float>/real_time_BigO        15.83 (1)       15.83 (1)
float_distance<float>/real_time_RMS             0 %             0 %
float_distance<double>/2/real_time           15.5 ns         15.5 ns     45098165
float_distance<double>/4/real_time           15.6 ns         15.6 ns     45065465
float_distance<double>/8/real_time           15.5 ns         15.5 ns     45058733
float_distance<double>/16/real_time          15.8 ns         15.7 ns     45078404
float_distance<double>/32/real_time          15.5 ns         15.5 ns     44832734
float_distance<double>/64/real_time          15.5 ns         15.5 ns     45077303
float_distance<double>/128/real_time         15.5 ns         15.5 ns     45067255
float_distance<double>/256/real_time         15.5 ns         15.5 ns     45073844
float_distance<double>/512/real_time         15.6 ns         15.6 ns     45109342
float_distance<double>/1024/real_time        15.5 ns         15.5 ns     44845180
float_distance<double>/2048/real_time        15.5 ns         15.5 ns     45051846
float_distance<double>/4096/real_time        15.5 ns         15.5 ns     45064317
float_distance<double>/8192/real_time        15.5 ns         15.5 ns     45115653
float_distance<double>/16384/real_time       15.5 ns         15.5 ns     45067642
float_distance<double>/real_time_BigO       15.54 (1)       15.54 (1)
float_distance<double>/real_time_RMS            0 %             0 %
*/
