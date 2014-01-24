//  Copyright John Maddock 2014
//  Copyright Christopher Kormanyos 2014.
//  Copyright Paul A. Bristow 2014.

//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

/* 

/*[cstdfloat_example_1
`This example shows use of a specified-width floating-point typedef
to evaluate a moderately complex math function using
[@http://en.wikipedia.org/wiki/Double-precision_floating-point_format IEEE754 64-bit double-precision].
about 15 decimal digits `std::numeric_limits<boost::float64_t>::digits10`,
(but never exceeding 17 decimal digits `std::numeric_limits<boost::float64_t>::max_digits10`).

The Jahnke Emden lambda function is described at

Weisstein, Eric W. "Lambda Function." From MathWorld--A Wolfram Web Resource. 
http://mathworld.wolfram.com/LambdaFunction.html

E. Jahnke and F. Emde, "Tables of Functions with Formulae and Curves,"
Dover, New York, 4th ed., (1945), pp. 180-188.

//] [/cstdfloat_example_1]
*/

#include <cmath>
#include <iostream>

#include <boost/cstdfloat.hpp>
#include <boost/math/special_functions.hpp>

boost::float64_t jahnke_emden_lambda(boost::float64_t v, boost::float64_t x)
{
  const boost::float64_t gamma_v_plus_one = boost::math::tgamma(v + 1);
  const boost::float64_t x_half_pow_v     = std::pow(x /2, v);

  return gamma_v_plus_one * boost::math::cyl_bessel_j(x, v) / x_half_pow_v;
}

int main()
{
/*[cstdfloat_example_2
`To ensure that all possibly significant digits including trailing zeros are shown.
*/
  std::cout.precision(std::numeric_limits<boost::float64_t>::max_digits10);
  std::cout.setf(std::ios::showpoint); // show all significant trailing zeros.
  try
  { // Always use try'n'catch blocks to ensure any error messages are available.

    // Evaluate a sample of the lambda function:
  boost::float64_t v = 1.;
  boost::float64_t x = 1.;
  
  std::cout << jahnke_emden_lambda(v, x) << std::endl; // 0.88010117148986700
//] [/cstdfloat_example_2]

  for (int i = 0; i < 10; i++)
  {
    std::cout << std::setprecision(2) << boost::float64_t(i) << ' ' 
      << std::setprecision(std::numeric_limits<boost::float64_t>::max_digits10)
      << jahnke_emden_lambda(boost::float64_t(i), v) << std::endl; // 
  }

  }
  catch (std::exception ex)
  {
    std::cout << "Thrown exception " << ex.what() << std::endl;
  }

} // int main()

/*
[cstdfloat_output

  1.0000000000000000
  0.88010117148986700
  0.00 0.00000000000000000
  1.0 0.88010117148986700
  2.0 4.6137984620549854
  3.0 16.274830009244948
  4.0 -25.360637961042869
  5.0 -1257.9038883512258
  6.0 -12749.592182518225
  7.0 -3020.9830849309396
  8.0 2421897.6013183575
  9.0 45577595.449204892

] [/cstdfloat_output]


*/
