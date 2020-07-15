// Copyright 2020 Matt Borland
//
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/math/special_functions/prime_sieve.hpp>
#include <boost/core/lightweight_test.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <list>
#include <deque>
#include <array>

template<typename Z>
void test_prime_sieve()
{
    std::vector<Z> primes;
    Z ref {168}; // Calculated with wolfram-alpha

    // Does the function work with a vector
    boost::math::prime_sieve(1000, std::back_inserter(primes));
    BOOST_TEST_EQ(primes.size(), ref);

    // Tests for correctness
    // 100
    primes.clear();
    boost::math::prime_sieve(100, std::back_inserter(primes));
    BOOST_TEST_EQ(primes.size(), 25);

    // 10'000
    primes.clear();
    boost::math::prime_sieve(10000, std::back_inserter(primes));
    BOOST_TEST_EQ(primes.size(), 1229);

    // 100'000
    primes.clear();
    boost::math::prime_sieve(100000, std::back_inserter(primes));
    BOOST_TEST_EQ(primes.size(), 9592);

    // 1'000'000
    primes.clear();
    boost::math::prime_sieve(1000000, std::back_inserter(primes));
    BOOST_TEST_EQ(primes.size(), 78498);

    // Does the function work with a list?
    std::list<Z> l_primes;
    boost::math::prime_sieve(1000, std::back_inserter(l_primes));
    BOOST_TEST_EQ(l_primes.size(), ref);

    // Does the function work with a deque?
    std::deque<Z> d_primes;
    boost::math::prime_sieve(1000, std::back_inserter(d_primes));
    BOOST_TEST_EQ(d_primes.size(), ref);
}

template<typename Z>
void test_prime_range()
{
    std::vector<Z> primes;
    Z ref {168}; // Calculated with wolfram-alpha

    // Does the upper and lower bound call work
    boost::math::prime_range(2, 1000, std::back_inserter(primes));
    BOOST_TEST_EQ(primes.size(), ref);

    // Does the upper bound call work
    primes.clear();
    boost::math::prime_range(2, 1000, std::back_inserter(primes));
    BOOST_TEST_EQ(primes.size(), ref);

    // Does it work with a deque?
    std::deque<Z> d_primes;
    boost::math::prime_range(2, 1000, std::back_inserter(d_primes));
    BOOST_TEST_EQ(d_primes.size(), ref);

    // Does it work with a list?
    std::list<Z> l_primes;
    boost::math::prime_range(2, 1000, std::front_inserter(l_primes));
    BOOST_TEST_EQ(l_primes.size(), ref);

    // Does the lower bound change the results?
    ref = 143; // Calculated with wolfram-alpha
    primes.clear();
    boost::math::prime_range(100, 1000, std::back_inserter(primes));
    BOOST_TEST_EQ(primes.size(), ref);

    // Does it work with 0 difference?
    primes.clear();
    boost::math::prime_range(2, 2, std::back_inserter(primes));
    BOOST_TEST_EQ(primes.size(), 1);

    // Will it call the sieve for large input
    ref = 78498; // Calculated with wolfram-alpha
    primes.clear();
    boost::math::prime_range(2, 1000000, std::back_inserter(primes));
    BOOST_TEST_EQ(primes.size(), ref);
}

template<typename Z>
void test_prime_sieve_overflow()
{
    std::vector<Z> primes;

    // Should die with call to BOOST_ASSERT
    boost::math::prime_sieve(static_cast<Z>(2), static_cast<Z>(std::numeric_limits<Z>::max()),
                             std::back_inserter(primes));
}

#if __cplusplus >= 201703
template<typename Z>
void test_par_prime_sieve()
{
    std::vector<Z> primes;
    Z ref {168}; // Calculated with wolfram-alpha

    // Does the function work with a vector
    boost::math::prime_sieve(std::execution::par, 1000, std::back_inserter(primes));
    BOOST_TEST_EQ(primes.size(), ref);

    // Tests for correctness
    // 100
    primes.clear();
    boost::math::prime_sieve(std::execution::par, 100, std::back_inserter(primes));
    BOOST_TEST_EQ(primes.size(), 25);

    // 10'000
    primes.clear();
    boost::math::prime_sieve(std::execution::par, 10000, std::back_inserter(primes));
    BOOST_TEST_EQ(primes.size(), 1229);

    // 100'000
    primes.clear();
    boost::math::prime_sieve(std::execution::par, 100000, std::back_inserter(primes));
    BOOST_TEST_EQ(primes.size(), 9592);

    // 1'000'000
    primes.clear();
    boost::math::prime_sieve(std::execution::par, 1000000, std::back_inserter(primes));
    BOOST_TEST_EQ(primes.size(), 78498);

    // Does the function work with a list?
    std::list<Z> l_primes;
    boost::math::prime_sieve(std::execution::par, 1000, std::back_inserter(l_primes));
    BOOST_TEST_EQ(l_primes.size(), ref);

    // Does the function work with a deque?
    std::deque<Z> d_primes;
    boost::math::prime_sieve(std::execution::par, 1000, std::back_inserter(d_primes));
    BOOST_TEST_EQ(d_primes.size(), ref);
}

template<typename Z>
void test_par_prime_sieve_large()
{
    std::vector<Z> primes;
    Z ref {1077871}; // Calculated with wolfram-alpha

    // Force the sieve into the multi-threading section
    boost::math::prime_sieve(std::execution::par, 16777217, std::back_inserter(primes));
    BOOST_TEST_EQ(primes.size(), ref);
}
#endif //__cplusplus >= 201703

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

    #if __cplusplus >= 201703
    test_par_prime_sieve<int>();
    test_par_prime_sieve<int32_t>();
    test_par_prime_sieve<int64_t>();
    test_par_prime_sieve<uint32_t>();

    //test_par_prime_sieve_large<int64_t>();
    #endif

    boost::report_errors();
}
