// Copyright John Maddock, 2022
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)


#include <iostream>
#include <sstream>
#include "math_unit_test.hpp"
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/math/quadrature/tanh_sinh.hpp>

using boost::math::quadrature::tanh_sinh;

#include <iostream>

void issue893() {
   typedef boost::multiprecision::cpp_bin_float_100 real;

   auto fun = [](real /*x*/) -> real {
      return 1.0;
   };

   tanh_sinh<real> integrator;
   const real a = 0.0;
   const real b = -0.9999995515592481132478776023609116290187750667053638330158486516399489191171270344610533516275817076;
   real y = integrator.integrate(fun, -b, a);
   
   CHECK_TRUE(y < 1);
}

int main()
{
   issue893();
   return boost::math::test::report_errors();
}
