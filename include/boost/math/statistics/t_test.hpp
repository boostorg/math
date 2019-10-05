//  (C) Copyright Nick Thompson 2019.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_STATISTICS_T_TEST_HPP
#define BOOST_MATH_STATISTICS_T_TEST_HPP

#include <boost/math/statistics/univariate_statistics.hpp>

namespace boost::math::statistics {

template<class RandomAccessContainer>
auto one_sample_t_test_statistic(RandomAccessContainer const & v, typename RandomAccessContainer::value_type assumed_mean) {
    using Real = typename RandomAccessContainer::value_type;
    auto [mu, s_sq] = mean_and_sample_variance(v.begin(), v.end());
    return (mu - assumed_mean)/sqrt(s_sq/v.size());
}

}
#endif
