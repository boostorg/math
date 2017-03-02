// Copyright Nicholas Thompson 2017.
// Copyright Paul A. Bristow 2017.
// Copyright John Maddock 2017.

// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or
// copy at http://www.boost.org/LICENSE_1_0.txt).

#include <iostream>
#include <limits>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <iterator>
#include <cmath>
#include <random>
#include <boost/random/uniform_real_distribution.hpp>

//[cubic_b_spline_example

/*`This example demonstrates how to use the cubic b spline interpolator for regularly spaced data.
*/
#include <boost/math/interpolators/cubic_b_spline.hpp>

//] [/airy_zeros_example_1]

int main()
{
    // We begin with an array of samples:
    std::vector<double> v(500);
    // And decide on a stepsize:
    double step = 0.01;

    // Initialize the vector with a function we'd like to interpolate:
    for (size_t i = 0; i < v.size(); ++i)
    {
        v[i] = sin(i*step);
    }
    // We could define an arbitrary start time, but for now we'll just use 0:
    boost::math::cubic_b_spline<double> spline(v.data(), v.size(), 0 /* start time */, step);

    // Now we can evaluate the spline wherever we please.
    std::mt19937 gen;
    boost::random::uniform_real_distribution<double> absissa(0, v.size()*step);
    for (size_t i = 0; i < 10; ++i)
    {
        double x = absissa(gen);
        std::cout << "sin(" << x << ") = " << sin(x) << ", spline interpolation gives " << spline(x) << std::endl;
        std::cout << "cos(" << x << ") = " << cos(x) << ", spline derivative interpolation gives " << spline.prime(x) << std::endl;
    }
}
