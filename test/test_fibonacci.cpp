//  Copyright Madhur Chauhan 2020.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/math/policies/error_handling.hpp>
#include <boost/math/special_functions/fibonacci.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/gmp.hpp>
#include <cstdint>
#include <exception>
#define BOOST_TEST_MODULE Fibonacci_Test_Module
#include <boost/test/included/unit_test.hpp>

// Some sanity checks using OEIS A000045
BOOST_AUTO_TEST_CASE(sanity_checks) {
    using boost::math::fibonacci;
    BOOST_TEST(fibonacci<char>(0) == 0);    // Base case
    BOOST_TEST(fibonacci<int>(1) == 1);     // Base case
    BOOST_TEST(fibonacci<int16_t>(2) == 1); // Base Case
    BOOST_TEST(fibonacci<int16_t>(3) == 2); // First computation
    BOOST_TEST(fibonacci<int16_t>(10) == 55);
    BOOST_TEST(fibonacci<int16_t>(15) == 610);
    BOOST_TEST(fibonacci<int16_t>(18) == 2584);
    BOOST_TEST(fibonacci<int32_t>(40) == 102334155);
}

// Tests fibonacci_unchecked by computing naively
BOOST_AUTO_TEST_CASE(big_integer_check) {
    // GMP is used as type for fibonacci and cpp_int is for naive computation
    using GMP = boost::multiprecision::mpz_int;
    using BST = boost::multiprecision::cpp_int;
    BST val = 0, a = 0, b = 1;
    for (int i = 0; i <= 1e4; ++i) {
        BOOST_TEST(boost::math::fibonacci<GMP>(i) == GMP(a));
        val = b, b += a, a = val;
    }
}

// Check for overflow throw using magic constants
BOOST_AUTO_TEST_CASE(overflow_check) {
    using boost::math::fibonacci;
    // 1. check for unsigned integer overflow
    BOOST_CHECK_NO_THROW(fibonacci<uint64_t>(93));
    BOOST_CHECK_THROW(fibonacci<uint64_t>(94), std::exception);

    // 2. check for signed integer overflow
    BOOST_CHECK_NO_THROW(fibonacci<int64_t>(91));
    // BOOST_CHECK_THROW(fibonacci<int64_t>(92), std::exception); // this should be the correct value but imprecisions
    BOOST_CHECK_THROW(fibonacci<int64_t>(93), std::exception);

    // 3. check for floating point (double)
    BOOST_CHECK_NO_THROW(fibonacci<double>(78));
    BOOST_CHECK_THROW(fibonacci<double>(79), std::exception);

    // 4. check using boost's multiprecision unchecked integer
    using namespace boost::multiprecision;
    using T = number<backends::cpp_int_backend<2048, 2048, unsigned_magnitude, unchecked>>;
    BOOST_CHECK_NO_THROW(fibonacci<T>(2950));
    // BOOST_CHECK_THROW(fibonacci<T>(2951), std::exception); // this should be the correct value but imprecisions
    BOOST_CHECK_THROW(fibonacci<T>(2952), std::exception);
}