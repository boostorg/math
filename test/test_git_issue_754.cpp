// Copyright Matt Borland 2022.
//
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/math/distributions/chi_squared.hpp>
#include "math_unit_test.hpp"

int main()
{
    int dof = 2;
    boost::math::chi_squared_distribution chiSq_int(dof);
    boost::math::chi_squared_distribution<double> chiSq(dof);

    CHECK_EQUAL(boost::math::cdf(chiSq_int, 3.1415), boost::math::cdf(chiSq, 3.1415));
}
