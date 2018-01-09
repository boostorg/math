/*
 *  (C) Copyright Nick Thompson 2018.
 *  Use, modification and distribution are subject to the
 *  Boost Software License, Version 1.0. (See accompanying file
 *  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 *
 */

#define BOOST_TEST_MODULE quasi_monte_carlo_test
#include <iostream>
#include <boost/test/included/unit_test.hpp>
#include <boost/math/quadrature/quasi_monte_carlo.hpp>


using boost::math::quasi_monte_carlo;

template<class Real>
void test_qmc()
{
}

BOOST_AUTO_TEST_CASE(quasi_monte_carlo_test)
{
    test_qmc<double>();
}
