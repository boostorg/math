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
#include <boost/math/distributions/normal.hpp>

namespace boost::math::statistics {

template<class RandomAccessContainer>
auto runs_above_threshold(RandomAccessContainer const & v,
                          typename RandomAccessContainer::value_type threshold)
{
    using Real = typename RandomAccessContainer::value_type;
    using std::sqrt;
    using std::abs;
    if (v.size() == 0)
    {
        throw std::domain_error("Need at least one sample to get number of runs.");
    }
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

    Real expected_runs = Real(1) + Real(2*nabove*nbelow)/Real(n);
    Real var = 2*nabove*nbelow*(2*nabove*nbelow-n)/Real(n*n*(n-1));

    bool run_sign = (v[0] > threshold);
    decltype(v.size()) runs = 1;
    for (decltype(v.size()) i = 1; i < v.size(); ++i) {
      if (v[i] == threshold) {
        // skip values precisely equal to threshold.
        continue;
      }
      if (run_sign == (v[i] > threshold)) {
        continue;
      }
      else {
        run_sign = (v[i] > threshold);
        runs++;
      }
    }
    Real sd = sqrt(var);
    Real statistic = (runs - expected_runs)/sd;
    auto normal = boost::math::normal_distribution<Real, no_promote_policy>(0,1);
    Real pvalue = 2*boost::math::cdf(normal, -abs(statistic));
    return std::make_pair(statistic, pvalue);
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
    std::cout << "Median = " << median << "\n";


    return runs_above_threshold(v, median);
}

}
#endif
