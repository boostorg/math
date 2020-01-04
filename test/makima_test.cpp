/*
 * Copyright Nick Thompson, 2020
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */

#include "math_unit_test.hpp"
#include <numeric>
#include <utility>
#include <random>
#include <boost/math/interpolators/makima.hpp>
#ifdef BOOST_HAS_FLOAT128
#include <boost/multiprecision/float128.hpp>
using boost::multiprecision::float128;
#endif


using boost::math::interpolators::makima;

template<typename Real>
void test_linear()
{
    std::vector<Real> x{1,2,3};
    std::vector<Real> y{1,2,3};

    auto x_copy = x;
    auto y_copy = y;
    auto akima = makima(std::move(x_copy), std::move(y_copy));

    
    CHECK_ULP_CLOSE(y[0], akima(x[0]), 0);
    CHECK_ULP_CLOSE(y[1], akima(x[1]), 0);
    CHECK_ULP_CLOSE(y[2], akima(x[2]), 0);
}

int main()
{
    test_linear<double>();
    return boost::math::test::report_errors();
}
