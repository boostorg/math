/*
 *  (C) Copyright Nick Thompson 2017.
 *  Use, modification and distribution are subject to the
 *  Boost Software License, Version 1.0. (See accompanying file
 *  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 *
 */

#define BOOST_TEST_MODULE floor_sqrt_test
#include <boost/test/included/unit_test.hpp>
#include <boost/math/tools/floor_sqrt.hpp>

using boost::math::floor_sqrt;

template<class Z>
void test_floor_sqrt()
{
    Z max_i = 500000;
    for (Z i = 1; i < max_i; ++i)
    {
        Z j = floor_sqrt(i*i);
        BOOST_CHECK_EQUAL(j, i);
    }

    for (Z i = 1; i < max_i; ++i)
    {
        Z  j = floor_sqrt(i);
        BOOST_CHECK_LE(j*j, i);
    }

    for (Z i = 2; i < 100; ++i)
    {
        for (Z k = (i-1)*(i-1)+1; k < i*i; ++k)
        {
            Z j = floor_sqrt(k);
            BOOST_CHECK_GE(j, i-1);
            BOOST_CHECK(j < i);
        }
    }
}

BOOST_AUTO_TEST_CASE(floor_sqrt_test)
{
    test_floor_sqrt<unsigned long long>();
}
