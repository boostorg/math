//  (C) Copyright Matt Borland 2022.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cmath>
#include <array>
#include <numeric>
#include <limits>
#include <iostream>
#include <matplot/matplot.h>
#include <boost/math/statistics/univariate_statistics.hpp>
#include <boost/math/constants/constants.hpp>
#include "../../include/boost/math/ccmath/sin.hpp"

constexpr std::array<double, 1000> array_init()
{
    constexpr double pi = boost::math::constants::pi<double>();
    
    std::array<double, 1000> temp;
    for (std::size_t i {}; i < temp.size(); ++i)
    {
        temp[i] = boost::math::ccmath::sin(pi * (static_cast<double>(i)/1000));
    }

    return temp;
}

int main(void)
{
    using namespace matplot;
    
    constexpr std::array<double, 2> y_range {-std::numeric_limits<double>::epsilon(), 
                                              std::numeric_limits<double>::epsilon()};

    std::array<double, 1000> x_vals;
    for (std::size_t i {}; i < x_vals.size(); ++i)
    {
        x_vals[i] = static_cast<double>(i) / 1000;
    }

    constexpr std::array<double, 1000> constexpr_sin {array_init()};
    
    // In a non-constexpr contect boost::math::ccmath::sin defaults to std::sin
    std::array<double, 1000> cmath_sin {array_init()};

    std::array<double, 1000> delta {};
    for (std::size_t i {}; i < delta.size(); ++i)
    {
        delta[i] = constexpr_sin[i] - cmath_sin[i];
    }

    plot(x_vals, delta);
    title("boost::math::ccmath::sin - std::sin (double precision)");
    xlabel("Increments of PI");
    ylabel("-1 to 1 epsilon");
    save("sin_plot.jpg");

    const double avg_delta = boost::math::statistics::mean(delta) / y_range.back();
    std::cout << "Mean difference: " << avg_delta << std::endl;

    return 0;
}
