/*
 * Copyright Nick Thompson, 2019
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */

#ifndef BOOST_MATH_STATISTICS_RUNS_TEST_HPP
#define BOOST_MATH_STATISTICS_RUNS_TEST_HPP

#include <cmath>
#include <algorithm>
#include <boost/math/statistics/univariate_statistics.hpp>

namespace boost::math::statistics {

template<class RandomAccessContainer>
auto runs_above_threshold(RandomAccessContainer const & v,
                          typename RandomAccessContainer::value_type threshold)
{
    typedef boost::math::policies::policy<
          boost::math::policies::promote_float<false>,
          boost::math::policies::promote_double<false> >
          no_promote_policy;

    decltype(v.size()) nabove = 0;
    decltype(v.size()) nbelow = 0;

    for (auto const & x : v) {
        if (x > threshold) {
            ++nabove;
        }
        if (x < threshold) {
            ++nbelow;
        }
    }
    Real n = nabove + nbelow;

    Real mu = Real(1) + Real(2*nabove*nbelow)/Real(n);
    Real var = 2*nabove*nbelow*(2*nabove*nbelow-n)/Real(n*n*(n-1));
}

template<class RandomAccessContainer>
auto runs_above_median(RandomAccessContainer const & v)
{
    using Real = typename RandomAccessContainer::value_type;
    using std::log;
    using std::sqrt;

    Real median;

    {
        // We have to memcpy v because the median does a partial sort,
        // and that would be catastrophic for the runs test.
        auto w = v;
        median = boost::math::statistics::median(w);
    }


    return runs_above_threshold(v, median);
}

}
#endif
