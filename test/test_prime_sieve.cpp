// Copyright 2020 Matt Borland
//
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/math/special_functions/prime_sieve.hpp>
#include <boost/math/special_functions/prime_approximation.hpp>
#include <boost/math/special_functions/detail/linear_prime_sieve.hpp>
#include <boost/math/special_functions/detail/interval_prime_sieve.hpp>
#include <boost/math/special_functions/detail/small_primes.hpp>
#include <boost/math/tools/stopwatch.hpp>
#include <boost/core/lightweight_test.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/gmp.hpp>
#include <list>
#include <deque>
#include <array>
#include <vector>
#include <iostream>
#include <chrono>
#include <algorithm>
#include <execution>

template<typename Integer>
void test_prime_sieve()
{
    constexpr std::size_t array_size {100'000};
    std::array<Integer, array_size> primes;
    std::fill(primes.begin(), primes.end(), 0);

    // 1'000
    boost::math::prime_sieve(std::execution::par, static_cast<Integer>(1'000), primes.begin());
    BOOST_TEST_EQ(array_size - std::count(primes.cbegin(), primes.cend(), 0), 168);

    // 10'000
    std::fill(primes.begin(), primes.end(), 0);
    boost::math::prime_sieve(std::execution::par, static_cast<Integer>(10'000), primes.begin());
    BOOST_TEST_EQ(array_size - std::count(primes.cbegin(), primes.cend(), 0), 1'229);

    // 100'000
    std::fill(primes.begin(), primes.end(), 0);
    boost::math::prime_sieve(std::execution::par, static_cast<Integer>(100'000), primes.begin());
    BOOST_TEST_EQ(array_size - std::count(primes.cbegin(), primes.cend(), 0), 9'592);

    // 1'000'000
    std::fill(primes.begin(), primes.end(), 0);
    boost::math::prime_sieve(std::execution::par, static_cast<Integer>(1'000'000), primes.begin());
    BOOST_TEST_EQ(array_size - std::count(primes.cbegin(), primes.cend(), 0), 78'498);
}

template<typename Integer>
void test_sequential_prime_sieve()
{
    constexpr std::size_t array_size {100'000};
    std::array<Integer, array_size> primes;
    std::fill(primes.begin(), primes.end(), 0);

    // 1'000
    boost::math::prime_sieve(static_cast<Integer>(1'000), primes.begin());
    BOOST_TEST_EQ(array_size - std::count(primes.cbegin(), primes.cend(), 0), 168);

    // 10'000
    std::fill(primes.begin(), primes.end(), 0);
    boost::math::prime_sieve(static_cast<Integer>(10'000), primes.begin());
    BOOST_TEST_EQ(array_size - std::count(primes.cbegin(), primes.cend(), 0), 1'229);

    // 100'000
    std::fill(primes.begin(), primes.end(), 0);
    boost::math::prime_sieve(static_cast<Integer>(100'000), primes.begin());
    BOOST_TEST_EQ(array_size - std::count(primes.cbegin(), primes.cend(), 0), 9'592);

    // 1'000'000
    std::fill(primes.begin(), primes.end(), 0);
    boost::math::prime_sieve(static_cast<Integer>(1'000'000), primes.begin());
    BOOST_TEST_EQ(array_size - std::count(primes.cbegin(), primes.cend(), 0), 78'498);
}

template<typename Integer>
void test_prime_sieve_container_interface()
{
    constexpr std::size_t array_size {100'000};
    std::array<Integer, array_size> primes;
    std::fill(primes.begin(), primes.end(), 0);

    // 1'000
    boost::math::prime_sieve(std::execution::par, static_cast<Integer>(1'000), &primes);
    BOOST_TEST_EQ(array_size - std::count(primes.cbegin(), primes.cend(), 0), 168);

    // 10'000
    std::fill(primes.begin(), primes.end(), 0);
    boost::math::prime_sieve(std::execution::par, static_cast<Integer>(10'000), &primes);
    BOOST_TEST_EQ(array_size - std::count(primes.cbegin(), primes.cend(), 0), 1'229);

    // 100'000
    std::fill(primes.begin(), primes.end(), 0);
    boost::math::prime_sieve(std::execution::par, static_cast<Integer>(100'000), &primes);
    BOOST_TEST_EQ(array_size - std::count(primes.cbegin(), primes.cend(), 0), 9'592);

    // 1'000'000
    std::fill(primes.begin(), primes.end(), 0);
    boost::math::prime_sieve(std::execution::par, static_cast<Integer>(1'000'000), &primes);
    BOOST_TEST_EQ(array_size - std::count(primes.cbegin(), primes.cend(), 0), 78'498);
}

template<typename Integer>
void test_interval_sieve()
{
    const std::size_t array_size {70'000};
    std::array<Integer, array_size> primes;
    std::fill(primes.begin(), primes.end(), 0);

    boost::math::detail::prime_sieve::IntervalSieve sieve(static_cast<Integer>(1'000), static_cast<Integer>(10'000), primes.begin());
    BOOST_TEST_EQ(array_size - std::count(primes.cbegin(), primes.cend(), 0), 1'061);

    std::fill(primes.begin(), primes.end(), 0);
    sieve.NewRange(static_cast<Integer>(10'000), static_cast<Integer>(100'000), primes.begin());
    BOOST_TEST_EQ(array_size - std::count(primes.cbegin(), primes.cend(), 0), 8'363);

    std::fill(primes.begin(), primes.end(), 0);
    sieve.NewRange(static_cast<Integer>(100'000), static_cast<Integer>(1'000'000), primes.begin());
    BOOST_TEST_EQ(array_size - std::count(primes.cbegin(), primes.cend(), 0), 68'906);
}

template<typename Integer>
void test_linear_sieve()
{
    constexpr std::size_t array_size {10'000};
    std::array<Integer, array_size> primes;
    std::fill(primes.begin(), primes.end(), 0);

    boost::math::detail::prime_sieve::linear_sieve(static_cast<Integer>(1'000), primes.begin());
    BOOST_TEST_EQ(array_size - std::count(primes.cbegin(), primes.cend(), 0), 168);

    std::fill(primes.begin(), primes.end(), 0);
    boost::math::detail::prime_sieve::linear_sieve(static_cast<Integer>(10'000), primes.begin());
    BOOST_TEST_EQ(array_size - std::count(primes.cbegin(), primes.cend(), 0), 1'229);

    std::fill(primes.begin(), primes.end(), 0);
    boost::math::detail::prime_sieve::linear_sieve(static_cast<Integer>(100'000), primes.begin());
    BOOST_TEST_EQ(array_size - std::count(primes.cbegin(), primes.cend(), 0), 9'592);

    std::vector<Integer> primes_v;
    //boost::math::prime_reserve(static_cast<Integer>(100'000), primes_v); Prime reserve does not work with OI. Need to use resize.
    primes_v.resize(10'000, 0);
    boost::math::detail::prime_sieve::linear_sieve(static_cast<Integer>(100'000), primes_v.begin());
    BOOST_TEST_EQ(array_size - std::count(primes_v.cbegin(), primes_v.cend(), 0), 9'592);

}

template<typename Integer>
void test_small_primes()
{
    constexpr std::size_t array_size {293};
    std::array<Integer, array_size> primes;
    std::fill(primes.begin(), primes.end(), 0);

    boost::math::detail::prime_sieve::small_primes(static_cast<Integer>(1'000), primes.begin());
    BOOST_TEST_EQ(array_size - std::count(primes.cbegin(), primes.cend(), 0), 168);

    boost::math::detail::prime_sieve::small_primes(static_cast<Integer>(1'920), primes.begin());
    BOOST_TEST_EQ(array_size - std::count(primes.cbegin(), primes.cend(), 0), 293);
}

int main()
{
    // Test prime approximation for constexpr
    static_assert(boost::math::prime_approximation(100) != 0, "log and/or floor is/are not constexpr");

    // Test SFINAE
    std::vector<int> test;
    auto test_ref = &test;
    static_assert(boost::math::detail::prime_sieve::is_container_t<std::remove_pointer_t<decltype(test_ref)>> == 1, "INOP");

    boost::math::tools::stopwatch w;

    // Individual Algorithms
    test_linear_sieve<int>();
    test_linear_sieve<int32_t>();
    test_linear_sieve<int64_t>();
    test_linear_sieve<uint32_t>();
    test_linear_sieve<boost::multiprecision::cpp_int>();
    test_linear_sieve<boost::multiprecision::mpz_int>();
    
    test_interval_sieve<int>();
    test_interval_sieve<int32_t>();
    test_interval_sieve<int64_t>();
    test_interval_sieve<uint32_t>();
    test_interval_sieve<boost::multiprecision::cpp_int>();
    test_interval_sieve<boost::multiprecision::mpz_int>();
    
    test_small_primes<int>();
    test_small_primes<int32_t>();
    test_small_primes<int64_t>();
    test_small_primes<uint32_t>();
    test_small_primes<boost::multiprecision::cpp_int>();
    test_small_primes<boost::multiprecision::mpz_int>();

    // Composite
    test_prime_sieve<int>();
    test_prime_sieve<int32_t>();
    test_prime_sieve<int64_t>();
    test_prime_sieve<uint32_t>();
    test_prime_sieve<boost::multiprecision::cpp_int>();
    test_prime_sieve<boost::multiprecision::mpz_int>();

    test_sequential_prime_sieve<int>();
    test_sequential_prime_sieve<int32_t>();
    test_sequential_prime_sieve<int64_t>();
    test_sequential_prime_sieve<uint32_t>();
    test_sequential_prime_sieve<boost::multiprecision::cpp_int>();
    test_sequential_prime_sieve<boost::multiprecision::mpz_int>();

    test_prime_sieve_container_interface<int>();
    test_prime_sieve_container_interface<int32_t>();
    test_prime_sieve_container_interface<int64_t>();
    test_prime_sieve_container_interface<uint32_t>();
    test_prime_sieve_container_interface<boost::multiprecision::cpp_int>();
    test_prime_sieve_container_interface<boost::multiprecision::mpz_int>();

    boost::math::set_l1d_size(100'000);
    BOOST_ASSERT_MSG(boost::math::detail::prime_sieve::L1D_SIZE == 100'000, "L1 Size not set");

    boost::report_errors();
}
