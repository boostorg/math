// Copyright Nick Thompson, 2017

// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or
// copy at http://www.boost.org/LICENSE_1_0.txt).

#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <boost/math/interpolators/catmull_rom.hpp>

using boost::math::catmull_rom;

int main()
{
    std::cout << "This shows how to use Boost's Catmull-Rom spline to create an Archimedean spiral.\n";

    // The Archimedean spiral is given by r = a*theta. We have set a = 1.
    std::vector<std::array<double, 2>> spiral_points(500);
    double theta_max = M_PI;
    for (size_t i = 0; i < spiral_points.size(); ++i)
    {
        double theta = i/theta_max;
        spiral_points[i] = {theta*cos(theta), theta*sin(theta)};
    }

    auto archimedean = catmull_rom<double, std::array<double,2>, 2>(spiral_points.data(), spiral_points.size());
    double max_s = archimedean.max_parameter();
    for (double s = 0; s < max_s; s += 0.1)
    {
        auto p = archimedean(s);
        double x = p[0];
        double y = p[1];
        double r = sqrt(x*x + y*y);
        double theta = atan2(y/r, x/r);
        if (theta < 0)
        {
           theta += 2*M_PI;
        }
        std::cout << "r = " << r << ", theta = " << theta << std::endl;
    }

    return 0;
}
