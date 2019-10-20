/*
 * Copyright Nick Thompson, 2019
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */

#ifndef BOOST_MATH_STATISTICS_LINEAR_REGRESSION_HPP
#define BOOST_MATH_STATISTICS_LINEAR_REGRESSION_HPP

#include <cmath>
#include <algorithm>
#include <utility>

namespace boost::math::statistics {


template<class RandomAccessContainer>
auto ordinary_least_squares(RandomAccessContainer const & x,
                            RandomAccessContainer const & y)
{
    using Real = typename RandomAccessContainer::value_type;
    using std::sqrt;
    using std::abs;
    if (x.size() <= 1)
    {
        throw std::domain_error("At least 2 samples are required to perform a linear regression.");
    }

    if (x.size() != y.size())
    {
        throw std::domain_error("The same number of samples must be in the independent and dependent variable.");
    }
    typedef boost::math::policies::policy<
          boost::math::policies::promote_float<false>,
          boost::math::policies::promote_double<false> >
          no_promote_policy;


    Real c0 = std::numeric_limits<Real>::quiet_NaN();
    Real c1 = std::numeric_limits<Real>::quiet_NaN();
    return std::make_pair(c0, c1);
}

}
#endif
