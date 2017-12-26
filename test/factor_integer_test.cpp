/*
 * Copyright Nick Thompson, 2017
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */

#define BOOST_TEST_MODULE factor_test
#include <boost/type_index.hpp>
#include <boost/test/included/unit_test.hpp>
#include <boost/math/tools/factor_integer.hpp>

using boost::math::pollard_rho;

template<class Z>
void test_small_integers()
{
    std::cout << "Testing small integers " << boost::typeindex::type_id<Z>().pretty_name()  << "\n";
    Z max_int = std::numeric_limits<Z>::max();
    Z mprime = boost::math::prime(boost::math::max_prime -1);
    std::cout << "Max prime: " << mprime << std::endl;
    std::cout << "Max provably factored: " << mprime*mprime << std::endl;
    std::cout << "Max for type           " << max_int << std::endl;

    std::ostream cnull(0);
    for(Z i = 2; i < mprime*mprime; ++i)
    {
        auto factors = boost::math::trial_division(i);
        auto p = 1;
        for (auto f : factors)
        {
            p *= std::pow(f.first, f.second);
        }
        BOOST_CHECK_EQUAL(p, i);
        if (i % 100000 == 0)
        {
            std::cout << "All numbers up to i = " << i << " have been tested\n";
        }
    }

}

template<class Z>
void test_pollard_rho()
{
    std::cout << "Testing Pollard rho on type " << boost::typeindex::type_id<Z>().pretty_name()  << "\n";
    Z i = 25279;
    Z factor1 = pollard_rho(i);
    Z factor2;
    if (factor1 == 0)
    {
        std::cout << "Pollard rho failed on i = " << i <<  "\n";
    }
    else
    {
        factor2 = i/factor1;
        BOOST_CHECK_EQUAL(factor1*factor2, i);
    }

    i = 10403;
    factor1 = pollard_rho(i);
    if (factor1 == 0)
    {
        std::cout << "Pollard rho failed on i = " << i <<  "\n";
    }
    else
    {
        factor2 = i/factor1;
        BOOST_CHECK_EQUAL(factor1*factor2, i);
    }


    i = 8051;
    factor1 = pollard_rho(i);
    if (factor1 == 0)
    {
        std::cout << "Pollard rho failed on i = " << i <<  "\n";
    }
    else
    {
        factor2 = i/factor1;
        BOOST_CHECK_EQUAL(factor1*factor2, i);
    }

    i = 13290059;
    factor1 = pollard_rho(i);
    if (factor1 == 0)
    {
        std::cout << "Pollard rho failed on i = " << i <<  "\n";
    }
    else
    {
        factor2 = i/factor1;
        BOOST_CHECK_EQUAL(factor1*factor2, i);
    }

    i = 28028821;
    factor1 = pollard_rho(i);
    if (factor1 == 0)
    {
        std::cout << "Pollard rho failed on i = " << i <<  "\n";
    }
    else
    {
        factor2 = i/factor1;
        BOOST_CHECK_EQUAL(factor1*factor2, i);
    }
}


BOOST_AUTO_TEST_CASE(factor_test)
{
    //test_small_integers<uint32_t>();
    test_pollard_rho<uint64_t>();
}
