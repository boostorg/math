// Copyright Paul A. Bristow 2018

// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/math/special_functions/lambert_w.hpp> // For lambert_w0 and wm1 functions.

#include <iostream>
#include <iomanip>

int main()
{
  double z = 2.0;
  double w0 = boost::math::lambert_w0(z);
  std::cout.setf(std::ios_base::showpoint); // Include any trailing zeros.
  // std::cout.precision(std::numeric_limits<double>::max_digits10); // Show all possibly significant digits.
  // Avoid using max_digfigs10 so as many old compilers can run the most basic lambert_w0 test.
  std::cout << " lambert_w0(" << z << ") = " << w0 << std::endl; // lambert_w0(2.00000) = 0.852606
  z = -0.2;
  double wm1 = boost::math::lambert_wm1(z);
  std::cout << " lambert_wm1(" << z << ") = " << wm1 << std::endl; // lambert_wm1(-0.200000) = -2.54264
  return 0;
} // int main()
