/*
 * Copyright Nick Thompson, 2017
 *
 * This example shows to to numerically integrate a periodic function using the adaptive_trapezoidal routine provided by boost.
 */

#include <iostream>
#include <cmath>
#include <limits>
#include <boost/math/quadrature/adaptive_trapezoidal.hpp>

int main()
{
    using boost::math::constants::two_pi;
    using boost::math::constants::third;
    // This function has an analytic form for its integral over a period: 2pi/3.
    auto f = [](double x) { return 1/(5 - 4*cos(x)); };

    double Q = boost::math::adaptive_trapezoidal(f, (double) 0, two_pi<double>());

    std::cout << std::setprecision(std::numeric_limits<double>::digits10);
    std::cout << "The adaptive trapezoidal rule gives the integral of our function as " << Q << "\n";
    std::cout << "The exact result is                                                 " << two_pi<double>()*third<double>() << "\n";

}
