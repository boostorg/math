//  Copyright Paul A. Bristow 2014.
//  Copyright Nakhar Agrawal 2013.
//  Copyright John Maddock 2013.
//  Copyright Christopher Kormanyos 2013.

//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma warning (disable : 4100) // unreferenced formal parameter.
#pragma warning (disable : 4127) // conditional expression is constant.

//#define BOOST_MATH_OVERFLOW_ERROR_POLICY ignore_error

#include <boost/math/concepts/real_concept.hpp>
#include <boost/math/special_functions/bernoulli.hpp>

#include <iostream>

/* First 50 from 2 to 100 inclusive: */
/* TABLE[N[BernoulliB[n], 200], {n,2,100,2}] */

//SC_(0.1666666666666666666666666666666666666666), 
//SC_(-0.0333333333333333333333333333333333333333), 
//SC_(0.0238095238095238095238095238095238095238), 
//SC_(-0.0333333333333333333333333333333333333333), 
//SC_(0.0757575757575757575757575757575757575757), 
//SC_(-0.2531135531135531135531135531135531135531), 
//SC_(1.1666666666666666666666666666666666666666), 
//SC_(-7.0921568627450980392156862745098039215686), 
//SC_(54.9711779448621553884711779448621553884711), 


int main()
{
//[bernoulli_example_1
  try
  { // It is always wise to use try'n'catch blocks around Boost.Math function
    // so that informative error messages will be displayed.


  /*`A simple example computes the value of `Bernoulli(2)` where the return type
is `double`.

[tip All odd Bernoulli numbers (> 1) are zero, so the parameter value 2 computes B[sub 4]. ]

*/
  std::cout
    << std::setprecision(std::numeric_limits<double>::digits10)
    << boost::math::bernoulli_b2n<double>(2) << std::endl;

/*` So B[sub 4] == -1/30 == -0.0333333333333333 

If we use Boost.Multiprecision and its 50 decimal digit floating-point type `cpp_dec_float_50`,
we can calculate the value of much larger numbers like `Bernoulli(100)`
and also obtain much higher precision.
*/

  std::cout
    << std::setprecision(std::numeric_limits<cpp_dec_float_50>::digits10)
    << boost::math::bernoulli_b2n<cpp_dec_float_50>(100) << std::endl;
 
//] //[/bernoulli_example_1]

//[bernoulli_example_2
/*`We can compute and save all the float-precision Bernoulli numbers from one call.
*/
  std::vector<float> bn(32); // Space for all the 32-bit `float` precision Bernoulli numbers.

  std::vector<float>::iterator it = bn.begin();

  // Start with Bernoulli number 0.
  boost::math::bernoulli_b2n<float>(0, bn.size(), it); // Fill vector with even Bernoulli numbers.

  for(size_t i = 0; i < bn.size(); i++)
  { // Show vector of even Bernoulli numbers, showing all significant decimal digits.
      std::cout << std::setprecision(std::numeric_limits<float>::digits10)
          << i*2 << ' '           
          << bn[i]
          << std::endl;
  }
//] //[/bernoulli_example_2]

  }
  catch (std::exception ex)
  {
    std::cout << "Thrown Exception caught: " << ex.what() << std::endl;
  }



//[bernoulli_example_3    
  /*`Of course, for any floating-point type, there is a maximum Bernoulli number than can be computed
  before it overflows the exponent.
  If we try to compute too high a Bernoulli number then an __overflow_error will be raised: the default
  behaviour for which is to throw a `std::overflow_error` exception.
  */
  try
  {

    std::cout
    << std::setprecision(std::numeric_limits<float>::digits10)
    << "Bernoulli number " << 33 * 2 <<std::endl;

     std::cout << boost::math::bernoulli_b2n<float>(33) << std::endl;

/*`
and (provided 'try'n'catch' blocks are used) we will get a helpful error message.
*/

//] //[/bernoulli_example_3]

  }
  catch (std::exception ex)
  {
    std::cout << "Thrown Exception caught: " << ex.what() << std::endl;
    }

} // int main()

/*

//[bernoulli_output_1
  -3.6470772645191354362138308865549944904868234686191e+215
//] //[/bernoulli_output_1]


//[bernoulli_output_2

  0 1
  2 0.166667
  4 -0.0333333
  6 0.0238095
  8 -0.0333333
  10 0.0757576
  12 -0.253114
  14 1.16667
  16 -7.09216
  18 54.9712
  20 -529.124
  22 6192.12
  24 -86580.3
  26 1.42552e+006
  28 -2.72982e+007
  30 6.01581e+008
  32 -1.51163e+010
  34 4.29615e+011
  36 -1.37117e+013
  38 4.88332e+014
  40 -1.92966e+016
  42 8.41693e+017
  44 -4.03381e+019
  46 2.11507e+021
  48 -1.20866e+023
  50 7.50087e+024
  52 -5.03878e+026
  54 3.65288e+028
  56 -2.84988e+030
  58 2.38654e+032
  60 -2.14e+034
  62 2.0501e+036
//] //[/bernoulli_output_2]

//[bernoulli_output_3
  Bernoulli number 66
 Thrown Exception caught: Error in function boost::math::bernoulli<float>: Overflow error while calculating tangent number 2

//] //[/bernoulli_output_3]

*/


