// Copyright Matt Borland, 2021
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <vector>
#include <boost/math/interpolators/2d_akima.hpp>

int main(void)
{
    // Data from table 1
    std::vector<double> x {11.16, 24.20, 19.85, 10.35, 19.72,  0.00, 20.87, 19.99, 10.28,  4.51,  0.00, 16.70,  6.08, 25.00, 14.90,  
                            0.00,  9.66,  5.22, 11.77, 15.10, 25.00, 25.00, 14.59, 15.20,  5.23,  2.14,  0.51, 25.00, 21.67,  3.31};
    
    std::vector<double> y { 1.24, 16.23, 10.72,  4.11,  1.39, 20.00, 20.00,  4.62, 15.16, 20.00,  4.48, 19.65,  4.58, 11.87,  3.12,
                            0.00, 20.00, 14.66, 10.47, 17.19,  3.87,  0.00,  8.71,  0.00, 10.72, 15.03,  8.37, 20.00, 14.36,  0.13};
    
    std::vector<double> z {22.15,  2.83,  7.97, 22.33, 16.83, 34.60,  5.74, 14.72, 21.59, 15.61, 61.77,  6.31, 35.74,  4.40, 21.70,
                           58.20,  4.73, 40.36, 13.62, 12.57,  8.74, 12.00, 14.81, 21.60, 26.50, 53.10, 49.43,  0.60,  5.52, 44.08};

    std::vector<std::vector<double>> xyz {x, y, z};

    boost::math::interpolators::akima test_1(x, y, z, 10);
    boost::math::interpolators::akima test_2(xyz[0], xyz[1], xyz[2], 10);
}
