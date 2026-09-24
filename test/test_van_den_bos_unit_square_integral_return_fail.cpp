//  Copyright 2026 Nicholas Thompson and Matt Borland
//  Distributed under the Boost Software License, Version 1.0.

#include <boost/math/quadrature/van_den_bos_unit_square.hpp>

int main()
{
    return boost::math::quadrature::van_den_bos_unit_square(
        [](double, double)
        {
            return 1;
        },
        1e-12);
}
