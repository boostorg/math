// Copyright 2020 Matt Borland
//
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/math/special_functions/prime_sieve.hpp>
#include <boost/core/lightweight_test.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/math/special_functions/interval_sieve.hpp>
#include <list>
#include <deque>
#include <array>

template<typename Integer>
void test_prime_sieve()
{
    std::vector<Integer> primes;
    Integer ref {168}; // Calculated with wolfram-alpha

    // Does the function work with a vector
    boost::math::prime_sieve(1000, primes);
    BOOST_TEST_EQ(primes.size(), ref);

    // Does the sequential policy work
    primes.clear();
    boost::math::prime_sieve(std::execution::seq, 1000, primes);
    BOOST_TEST_EQ(primes.size(), ref);

    // Tests for correctness
    // 2
    primes.clear();
    boost::math::prime_sieve(2, primes);
    BOOST_TEST_EQ(primes.size(), 0);

    // 100
    primes.clear();
    boost::math::prime_sieve(100, primes);
    BOOST_TEST_EQ(primes.size(), 25);

    // 10'000
    primes.clear();
    boost::math::prime_sieve(10000, primes);
    BOOST_TEST_EQ(primes.size(), 1229);

    // 100'000
    primes.clear();
    boost::math::prime_sieve(100000, primes);
    BOOST_TEST_EQ(primes.size(), 9592);

    // 1'000'000
    primes.clear();
    boost::math::prime_sieve(1000000, primes);
    BOOST_TEST_EQ(primes.size(), 78498);

    // Does the function work with a list?
    std::list<Integer> l_primes;
    boost::math::prime_sieve(1000, l_primes);
    BOOST_TEST_EQ(l_primes.size(), ref);

    // Does the function work with a deque?
    std::deque<Integer> d_primes;
    boost::math::prime_sieve(1000, d_primes);
    BOOST_TEST_EQ(d_primes.size(), ref);
}

template<typename Integer>
void test_prime_range()
{
    std::vector<Integer> primes;
    Integer ref {168}; // Calculated with wolfram-alpha

    // Does the upper and lower bound call work
    boost::math::prime_range(static_cast<Integer>(2), static_cast<Integer>(1000), primes);
    BOOST_TEST_EQ(primes.size(), ref);

    // Does the upper bound call work
    primes.clear();
    boost::math::prime_range(static_cast<Integer>(2), static_cast<Integer>(1000), primes);
    BOOST_TEST_EQ(primes.size(), ref);

    // Does it work with a deque?
    std::deque<Integer> d_primes;
    boost::math::prime_range(static_cast<Integer>(2), static_cast<Integer>(1000), d_primes);
    BOOST_TEST_EQ(d_primes.size(), ref);

    // Does it work with a list?
    std::list<Integer> l_primes;
    boost::math::prime_range(static_cast<Integer>(2), static_cast<Integer>(1000), l_primes);
    BOOST_TEST_EQ(l_primes.size(), ref);

    // Does the lower bound change the results?
    ref = 143; // Calculated with wolfram-alpha
    primes.clear();
    boost::math::prime_range(static_cast<Integer>(100), static_cast<Integer>(1000), primes);
    BOOST_TEST_EQ(primes.size(), ref);

    // Will it call the sieve for large input
    ref = 78498; // Calculated with wolfram-alpha
    primes.clear();
    boost::math::prime_range(static_cast<Integer>(2), static_cast<Integer>(1000000), primes);
    BOOST_TEST_EQ(primes.size(), ref);
}

template<typename Integer>
void test_prime_sieve_overflow()
{
    std::vector<Integer> primes;

    // Should die with call to BOOST_ASSERT
    boost::math::prime_sieve(static_cast<Integer>(2), static_cast<Integer>(std::numeric_limits<Integer>::max()), primes);
}

template<typename Integer>
void test_par_prime_sieve_large()
{
    std::vector<Integer> primes;
    Integer ref {1077871}; // Calculated with wolfram-alpha

    // Force the sieve into the multi-threading section
    boost::math::prime_sieve(static_cast<Integer>(16777217), primes);
    BOOST_TEST_EQ(primes.size(), ref);
}

int main()
{
    test_prime_sieve<int>();
    test_prime_sieve<int32_t>();
    test_prime_sieve<int64_t>();
    test_prime_sieve<uint32_t>();
    
    test_prime_range<int>();
    test_prime_range<int32_t>();
    test_prime_range<int64_t>();
    test_prime_range<uint32_t>();

    //test_prime_sieve_overflow<int16_t>();

    test_prime_sieve<boost::multiprecision::cpp_int>();

    //test_par_prime_sieve_large<int64_t>();
    
    boost::report_errors();
}
