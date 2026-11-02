//  Copyright 2026 Shen-Ta Hsieh
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
//  See: https://github.com/boostorg/math/issues/1621

#include <boost/math/statistics/univariate_statistics.hpp>
#include <boost/core/lightweight_test.hpp>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <vector>

int main()
{
    std::vector<double> v {1.0, 2.0, 2.0, std::numeric_limits<double>::quiet_NaN()};
    std::vector<double> modes;

    // Without exceptions, data that contain a NaN give no mode, and boost::throw_exception is not called
#ifndef BOOST_MATH_NO_EXCEPTIONS
    BOOST_TEST_THROWS(boost::math::statistics::mode(v.begin(), v.end(), std::back_inserter(modes)), std::domain_error);
#else
    boost::math::statistics::mode(v.begin(), v.end(), std::back_inserter(modes));
#endif
    BOOST_TEST(modes.empty());

    return boost::report_errors();
}
