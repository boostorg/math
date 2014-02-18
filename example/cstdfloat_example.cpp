//  Copyright John Maddock 2014
//  Copyright Christopher Kormanyos 2014.
//  Copyright Paul A. Bristow 2014.

//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

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

/*!
Function max_digits10
Returns maximum number of possibly significant decimal digits for a floating-point type FPT,
even for older compilers/standard libraries that
lack support for std::numeric_limits<FPT>::max_digits10,
when the Kahan formula is used instead.
Also provides the correct result for Visual Studio 2010 where the value provided for float is wrong.
*/
namespace boost
{
template <typename FPT>
int max_digits10()
{
// Since max_digits10 is not defined (or wrong) on older systems, define a local max_digits10.
  // Usage:   int m = max_digits10<boost::float64_t>();
  const int m =
#if (defined BOOST_NO_CXX11_NUMERIC_LIMITS) || (_MSC_VER == 1600) // is wrongly 8 not 9 for VS2010.
  2 + std::numeric_limits<FPT>::digits * 3010/10000;
#else
  std::numeric_limits<FPT>::max_digits10;
#endif
  return m;
}
} // namespace boost


boost::float64_t jahnke_emden_lambda(boost::float64_t v, boost::float64_t x)
{
  const boost::float64_t gamma_v_plus_one = boost::math::tgamma(v + 1);
  const boost::float64_t x_half_pow_v     = std::pow(x /2, v);

  return gamma_v_plus_one * boost::math::cyl_bessel_j(x, v) / x_half_pow_v;
}

int main()
{
  std::cout.setf(std::ios::showpoint); // show all significant trailing zeros.

    long double p = 1.L;
  //std::cout.precision(std::numeric_limits<long double>::digits10);

  std::cout << "pi = "  << p << std::endl;

/*[cstdfloat_example_2
`To ensure that all possibly significant digits including trailing zeros are shown.
*/
  std::cout.precision(std::numeric_limits<boost::float64_t>::max_digits10);
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

  std::cout << "Floating-point type long double is available." << std::endl;
  std::cout << "  std::numeric_limits<long double>::digits10 == "
    << std::numeric_limits<long double>::digits10 << std::endl;
  std::cout << "  std::numeric_limits<long double>::max_digits10 == "
    << std::numeric_limits<long double>::max_digits10 << std::endl; // 21

  long double p = boost::math::constants::pi<double>();
  std::cout.precision(std::numeric_limits<long double>::digits10);

  std::cout << "pi = "  << p << std::endl;



  // http://www.boost.org/doc/libs/1_54_0/libs/multiprecision/doc/html/boost_multiprecision/tut/floats/float128.html


#ifdef BOOST_FLOAT16_C
  std::cout << "Floating-point type boost::float16_t is available." << std::endl;
#else
  std::cout << "Floating-point type boost::float16_t is NOT available." << std::endl;
#endif
#ifdef BOOST_FLOAT32_C
  std::cout << "Floating-point type boost::float32_t is available." << std::endl;
  std::cout << "  std::numeric_limits<boost::float32_t>::digits10 == "
    << std::numeric_limits<boost::float32_t>::digits10 << std::endl;
  std::cout << "  std::numeric_limits<boost::float32_t>::max_digits10 == "
    << std::numeric_limits<boost::float32_t>::max_digits10 << std::endl;
#else
  std::cout << "Floating-point type boost::float32_t is NOT available." << std::endl;
#endif
#ifdef BOOST_FLOAT64_C
  std::cout << "Floating-point type boost::float64_t is available." << std::endl;
    std::cout << "  std::numeric_limits<boost::float64_t>::digits10 == "
    << std::numeric_limits<boost::float64_t>::digits10 << std::endl;
  std::cout << "  std::numeric_limits<boost::float64_t>::max_digits10 == "
    << std::numeric_limits<boost::float64_t>::max_digits10 << std::endl;

#else
  std::cout << "Floating-point type boost::float64_t is NOT available." << std::endl;
#endif
#ifdef BOOST_FLOAT80_C
  std::cout << "Floating-point type boost::float80_t is available." << std::endl;
  std::cout << "  std::numeric_limits<boost::float80_t>::digits10 == "
    << std::numeric_limits<boost::float80_t>::digits10 << std::endl;
  std::cout << "  std::numeric_limits<boost::float80_t>::max_digits10 == "
    << std::numeric_limits<boost::float80_t>::max_digits10 << std::endl;

#else
  std::cout << "Floating-point type boost::float80_t is NOT available." << std::endl;

#endif
#ifdef BOOST_FLOAT128_C
  std::cout << "Floating-point type boost::float128_t is available." << std::endl;
    std::cout << "  std::numeric_limits<boost::float128_t>::digits10 == "
    << std::numeric_limits<boost::float128_t>::digits10 << std::endl;
  std::cout << "  std::numeric_limits<boost::float128_t>::max_digits10 == "
    << std::numeric_limits<boost::float128_t>::max_digits10 << std::endl;
#else
  std::cout << "Floating-point type boost::float128_t is NOT available." << std::endl;
#endif

// Show some constants that depend on the available type.
#ifdef BOOST_FLOAT16_C
  std::cout.precision(boost::max_digits10<boost::float16_t>()); // Show all significant decimal digits,
  std::cout.setf(std::ios::showpoint); // including all significant trailing zeros.

  std::cout << "BOOST_FLOAT16_C(123.456789012345678901234567890) = "
    << BOOST_FLOAT16_C(123.456789012345678901234567890) << std::endl;
  // BOOST_FLOAT16_C(123.456789012345678901234567890) = 123.45678901234568
#endif

#ifdef BOOST_FLOAT32_C
  std::cout.precision(boost::max_digits10<boost::float32_t>()); // Show all significant decimal digits,
  std::cout.setf(std::ios::showpoint); // including all significant trailing zeros.


  std::cout << "BOOST_FLOAT32_C(123.456789012345678901234567890) = "
    << BOOST_FLOAT32_C(123.456789012345678901234567890) << std::endl;
  //   BOOST_FLOAT32_C(123.456789012345678901234567890) = 123.456787

#endif

#ifdef BOOST_FLOAT64_C
  std::cout.precision(boost::max_digits10<boost::float64_t>()); // Show all significant decimal digits,
  std::cout.setf(std::ios::showpoint); // including all significant trailing zeros.

  std::cout << "BOOST_FLOAT64_C(123.456789012345678901234567890) = "
    << BOOST_FLOAT64_C(123.456789012345678901234567890) << std::endl;
  // BOOST_FLOAT64_C(123.456789012345678901234567890) = 123.45678901234568
#endif

#ifdef BOOST_FLOAT80_C
  std::cout.precision(boost::max_digits10<boost::float80_t>()); // Show all significant decimal digits,
  std::cout.setf(std::ios::showpoint); // including all significant trailing zeros.

  std::cout << "BOOST_FLOAT80_C(123.456789012345678901234567890) = "
    << BOOST_FLOAT80_C(123.456789012345678901234567890) << std::endl;
  // BOOST_FLOAT80_C(123.456789012345678901234567890) = 123.456789012345678901234567890000000000000000
#endif

#ifdef BOOST_FLOAT128_C
  std::cout.precision(boost::max_digits10<boost::float128_t>()); // Show all significant decimal digits,
  std::cout.setf(std::ios::showpoint); // including all significant trailing zeros.

  std::cout << "BOOST_FLOAT128_C(123.456789012345678901234567890) = "
    << BOOST_FLOAT128_C(123.456789012345678901234567890) << std::endl;
  // BOOST_FLOAT128_C(123.456789012345678901234567890) = 123.456789012345678901234567890000000000000000000
#endif
  }
  catch (std::exception ex)
  {
    std::cout << "Thrown exception " << ex.what() << std::endl;
  }

} // int main()

/*
[cstdfloat_output

0.88010117148986700
0.00 0.00000000000000000
1.0 0.88010117148986700
2.0 4.6137984620549872
3.0 16.274830009244951
4.0 -25.360637961042869
5.0 -1257.9038883512264
6.0 -12749.592182518225
7.0 -3020.9830849309437
8.0 2421897.6013183584
9.0 45577595.449204877
Floating-point type boost::float16_t is NOT available.
Floating-point type boost::float32_t is available.
  std::numeric_limits<boost::float32_t>::digits10 == 6
  std::numeric_limits<boost::float32_t>::max_digits10 == 9
Floating-point type boost::float64_t is available.
  std::numeric_limits<boost::float64_t>::digits10 == 15
  std::numeric_limits<boost::float64_t>::max_digits10 == 17
Floating-point type boost::float80_t is NOT available.
Floating-point type boost::float128_t is NOT available.
BOOST_FLOAT32_C(123.456789012345678901234567890) = 123.456787
BOOST_FLOAT64_C(123.456789012345678901234567890) = 123.45678901234568


] [/cstdfloat_output]


*/
