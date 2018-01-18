/*
 *  (C) Copyright Nick Thompson 2018.
 *  Use, modification and distribution are subject to the
 *  Boost Software License, Version 1.0. (See accompanying file
 *  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */

#define BOOST_TEST_MODULE atanassov_sequence_test
#include <boost/test/included/unit_test.hpp>
#include <boost/math/tools/atanassov_sequence.hpp>

using boost::math::atanassov_sequence;

template<class Real>
void test_atanassov_sequence()
{
    atanassov_sequence<Real> as(6);
    int i = 0;
    std::cout << std::fixed;
    std::cout << std::setprecision(std::numeric_limits<Real>::digits10);
    std::vector<Real> r(6);
    while(i++ < 500)
    {
        as(r.begin(), r.end());
        for (size_t j = 0; j < r.size();++j)
        {
            std::cout << r[j] << "\t";
        }
        std::cout << std::endl;
    }
}

BOOST_AUTO_TEST_CASE(atanassov_sequence_test)
{
    test_atanassov_sequence<double>();
}
