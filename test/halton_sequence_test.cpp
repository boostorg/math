/*
 *  (C) Copyright Nick Thompson 2018.
 *  Use, modification and distribution are subject to the
 *  Boost Software License, Version 1.0. (See accompanying file
 *  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 *
 */

#define BOOST_TEST_MODULE halton_sequence_test
#include <iostream>
#include <iomanip>
#include <boost/test/included/unit_test.hpp>
#include <boost/math/tools/halton_sequence.hpp>


using boost::math::halton_sequence;
template<class Real>
void test_halton()
{
    halton_sequence<Real> hs(6, true);
    int i = 0;
    std::cout << std::fixed;
    std::cout << std::setprecision(std::numeric_limits<Real>::digits10);
    std::vector<Real> r(6);
    while(i++ < 500)
    {
        hs(r.begin(), r.end());
        for (size_t j = 0; j < r.size();++j)
        {
            std::cout << r[j] << "\t";
        }
        std::cout << std::endl;
    }
}

BOOST_AUTO_TEST_CASE(halton_sequence_test)
{
    test_halton<double>();
}
