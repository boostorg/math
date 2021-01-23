//  (C) Copyright Nick Thompson 2018.
//  (C) Copyright Matt Borland 2021.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/math/statistics/set.hpp>
#include <boost/core/lightweight_test.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <random>
#include <algorithm>
#include <iostream>
#include <cstddef>
#include <cmath>

#if (__cplusplus > 201700 || _MSVC_LANG > 201700) && (__GNUC__ > 9 || (__clang_major__ > 9 && defined __GLIBCXX__)  || _MSC_VER > 1927)
#include <execution>
#endif

using boost::multiprecision::cpp_bin_float_50;

// To stress test, set global_seed = 0, global_size = huge.
static constexpr std::size_t global_seed = 0;
static constexpr std::size_t global_size = 128;

template<typename T>
std::vector<T> generate_random_vector(std::size_t size, std::size_t seed)
{
    if (seed == 0)
    {
        std::random_device rd;
        seed = rd();
    }
    std::vector<T> v(size);

    std::mt19937 gen(seed);

    if constexpr (std::is_floating_point<T>::value)
    {
        std::normal_distribution<T> dis(0, 1);
        for (std::size_t i = 0; i < v.size(); ++i)
        {
         v[i] = dis(gen);
        }
        return v;
    }
    else if constexpr (std::is_integral<T>::value)
    {
        // Rescaling by larger than 2 is UB!
        std::uniform_int_distribution<T> dis(std::numeric_limits<T>::lowest()/2, (std::numeric_limits<T>::max)()/2);
        for (std::size_t i = 0; i < v.size(); ++i)
        {
         v[i] = dis(gen);
        }
        return v;
    }
    else if constexpr (boost::multiprecision::number_category<T>::value == boost::multiprecision::number_kind_floating_point)
    {
        std::normal_distribution<long double> dis(0, 1);
        for (std::size_t i = 0; i < v.size(); ++i)
        {
            v[i] = dis(gen);
        }
        return v;
    }
    else
    {
        BOOST_ASSERT_MSG(false, "Could not identify type for random vector generation.");
        return v;
    }
}

template<typename Real, typename ExecutionPolicy>
void test(ExecutionPolicy&& exec)
{
    using boost::math::statistics::metrics;
    using std::sqrt;
    
    const Real tol = 10*std::numeric_limits<Real>::epsilon();

    std::vector<Real> test {generate_random_vector<Real>(global_size, global_seed)};
    // Add a pair of values to guarantee there is a mode
    test.emplace_back(Real(2));
    test.emplace_back(Real(2));
    
    boost::math::statistics::stats set(exec, test.begin(), test.end());
    set.calc(metrics::all);

    BOOST_TEST(abs(set.mean() - boost::math::statistics::mean(exec, test.begin(), test.end())) < tol);
    BOOST_TEST(abs(set.median() - boost::math::statistics::median(exec, test.begin(), test.end())) < tol);
    BOOST_TEST(set.mode().size() != 0);
    
    const Real var = boost::math::statistics::variance(exec, test.begin(), test.end());
    BOOST_TEST(abs(set.stddev() - sqrt(var)) < tol);
    BOOST_TEST(abs(set.variance() - var) < tol);
    BOOST_TEST(abs(std::get<0>(set.first_four_moments()) - std::get<0>(boost::math::statistics::first_four_moments(exec, test.begin(), test.end()))) < tol);
}

template<typename Z, typename ExecutionPolicy>
void test_integer(ExecutionPolicy&& exec)
{
    using boost::math::statistics::metrics;
    using std::sqrt;
    
    const double tol = 10*std::numeric_limits<double>::epsilon();

    std::vector<double> test {generate_random_vector<double>(global_size, global_seed)};
    // Add a pair of values to guarantee there is a mode
    test.emplace_back(2.0);
    test.emplace_back(2.0);

    boost::math::statistics::stats set(exec, test.begin(), test.end());
    set.calc(metrics::all);

    BOOST_TEST(abs(set.mean() - boost::math::statistics::mean(exec, test.begin(), test.end())) < tol);
    BOOST_TEST(abs(set.median() - boost::math::statistics::median(exec, test.begin(), test.end())) < tol);
    BOOST_TEST(set.mode().size() != 0);
    
    const double var = boost::math::statistics::variance(exec, test.begin(), test.end());
    BOOST_TEST(abs(set.stddev() - sqrt(var)) < tol);
    BOOST_TEST(abs(set.variance() - var) < tol);
    BOOST_TEST(abs(std::get<0>(set.first_four_moments()) - std::get<0>(boost::math::statistics::first_four_moments(exec, test.begin(), test.end()))) < tol);
}

int main(void)
{   
    // Support compilers with P0024R2 implemented without linking TBB
    // https://en.cppreference.com/w/cpp/compiler_support
    #if (__cplusplus > 201700 || _MSVC_LANG > 201700) && (__GNUC__ > 9 || (__clang_major__ > 9 && defined __GLIBCXX__)  || _MSC_VER > 1927)
    
    test<float>(std::execution::seq);
    test<float>(std::execution::par);
    test<double>(std::execution::seq);
    test<double>(std::execution::par);
    test<long double>(std::execution::seq);
    test<long double>(std::execution::par);
    // Expensive using CI:
    //test<cpp_bin_float_50>(std::execution::seq);
    //test<cpp_bin_float_50>(std::execution::par);

    //test<int>(std::execution::seq);
    
    #endif
    return boost::report_errors();
}
