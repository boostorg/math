//           Copyright Matthew Pulver 2018 - 2019.
// Distributed under the Boost Software License, Version 1.0.
//      (See accompanying file LICENSE_1_0.txt or copy at
//           https://www.boost.org/LICENSE_1_0.txt)

#include <boost/math/differentiation/autodiff.hpp>
#include <iostream>

template<typename T>
T fourth_power(T x)
{
    x *= x;
    return x *= x;
}

int main()
{
    using namespace boost::math::differentiation;

    constexpr int Order=5; // The highest order derivative to be calculated.
    const autodiff_fvar<double,Order> x = make_fvar<double,Order>(2.0); // Find derivatives at x=2.
    const autodiff_fvar<double,Order> y = fourth_power(x);
    for (int i=0 ; i<=Order ; ++i)
        std::cout << "y.derivative("<<i<<") = " << y.derivative(i) << std::endl;
    return 0;
}
/*
Output:
y.derivative(0) = 16
y.derivative(1) = 32
y.derivative(2) = 48
y.derivative(3) = 48
y.derivative(4) = 24
y.derivative(5) = 0
**/
