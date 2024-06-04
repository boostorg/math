// Copyright John Maddock 2022.

// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <iomanip>
#include <limits>

import boost.math.constants;

int main()
{
   std::cout << std::setprecision(std::numeric_limits<float>::max_digits10) << boost::math::constants::pi<float>() << std::endl;
   std::cout << std::setprecision(std::numeric_limits<double>::max_digits10) << boost::math::constants::pi<double>() << std::endl;
   std::cout << std::setprecision(std::numeric_limits<long double>::max_digits10) << boost::math::constants::pi<long double>() << std::endl;

   std::cout << std::setprecision(std::numeric_limits<float>::max_digits10) << boost::math::float_constants::pi << std::endl;
   std::cout << std::setprecision(std::numeric_limits<double>::max_digits10) << boost::math::double_constants::pi << std::endl;
   std::cout << std::setprecision(std::numeric_limits<long double>::max_digits10) << boost::math::long_double_constants::pi << std::endl;
}
