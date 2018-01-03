/*
 *  (C) Copyright Nick Thompson 2017.
 *  Use, modification and distribution are subject to the
 *  Boost Software License, Version 1.0. (See accompanying file
 *  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 *
 */

#define BOOST_TEST_MODULE sieve_of_eratosthenes_test
#include <boost/test/included/unit_test.hpp>
#include <boost/math/tools/sieve_of_eratosthenes.hpp>
#include <boost/math/special_functions/prime.hpp>
#include <boost/multiprecision/cpp_int.hpp>


using boost::multiprecision::uint256_t;
using boost::math::sieve_of_eratosthenes;

template<class Z>
void test_sieve()
{
    sieve_of_eratosthenes<Z> sieve(boost::math::prime(boost::math::max_prime - 1));
    for (unsigned i = 0; i < sieve.prime_count(); ++i)
    {
        BOOST_CHECK_EQUAL(sieve.prime(i), boost::math::prime((unsigned)i));
    }
}

BOOST_AUTO_TEST_CASE(sieve_of_eratosthenes_test)
{
    test_sieve<int>();
    test_sieve<unsigned>();
    test_sieve<unsigned long long>();
    test_sieve<uint256_t>();
}
