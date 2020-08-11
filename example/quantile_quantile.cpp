//  Copyright Nick Thompson, 2020.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/math/statistics/qq_plot.hpp>
#include <boost/math/distributions/normal.hpp>

using boost::math::statistics::qq_plot;
int main() {
    using Real = double;
    auto f = [](Real x)->Real { return x/2; };
    auto g = [](Real x)->Real { return x*x; };
    auto qq = qq_plot<Real>(f, g);
    qq.write("foo.svg");
    return 0;
}