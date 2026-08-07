//  Copyright Matt Borland 2026
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
//  See: https://github.com/boostorg/math/issues/1421

#include <boost/math/statistics/univariate_statistics.hpp>
#include <iterator>
#include <sstream>

int main()
{
    std::istringstream iss {"1 2 3 4 5 6 7 8 9"};
    std::istream_iterator<double> first {iss};
    const std::istream_iterator<double> last {};

    // std::istream_iterator is an input iterator, so this must not compile.
    const auto m = boost::math::statistics::mean(first, last);

    return static_cast<int>(m);
}
