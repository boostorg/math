// Copyright Paul A. Bristow 2016, 2018.

// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or
//  copy at http ://www.boost.org/LICENSE_1_0.txt).

//! Lambert W examples of controlling precision 

// #define BOOST_MATH_INSTRUMENT_LAMBERT_W  // #define only for (much) diagnostic output.

#include <boost/config.hpp> // for BOOST_PLATFORM, BOOST_COMPILER,  BOOST_STDLIB ...
#include <boost/version.hpp>   // for BOOST_MSVC versions.
#include <boost/cstdint.hpp>
#include <boost/exception/exception.hpp>  // boost::exception
#include <boost/math/constants/constants.hpp> // For exp_minus_one == 3.67879441171442321595523770161460867e-01.
#include <boost/math/policies/policy.hpp>

// Built-in/fundamental GCC float128 or Intel Quad 128-bit type, if available.
#ifdef __GNUC__
#include <boost/multiprecision/float128.hpp> // Not available for MSVC.
// sets BOOST_MP_USE_FLOAT128 for GCC
using boost::multiprecision::float128;
#endif //# NOT _MSC_VER

#include <boost/multiprecision/cpp_dec_float.hpp> // boost::multiprecision::cpp_dec_float_50
using boost::multiprecision::cpp_dec_float_50; // 50 decimal digits type.
using boost::multiprecision::cpp_dec_float_100; // 100 decimal digits type.

#include <boost/multiprecision/cpp_bin_float.hpp>
using boost::multiprecision::cpp_bin_float_double_extended;
using boost::multiprecision::cpp_bin_float_double;
using boost::multiprecision::cpp_bin_float_quad;
// For lambert_w function.
#include <boost/math/special_functions/lambert_w.hpp>
// using boost::math::lambert_w0;
// using boost::math::lambert_wm1;

#include <iostream>
#include <exception>
#include <stdexcept>
#include <string>
#include <limits>  // For std::numeric_limits.

int main()
{
  try
  {
    std::cout << "Lambert W examples of precision control." << std::endl;
    std::cout.precision(std::numeric_limits<double>::max_digits10);
    std::cout << std::showpoint << std::endl; // Show any trailing zeros.

    using boost::math::constants::exp_minus_one;

    using boost::math::lambert_w0;
    using boost::math::lambert_wm1;

    // Error handling policy examples.
    using namespace boost::math::policies;
    using boost::math::policies::make_policy;
    using boost::math::policies::policy;
    using boost::math::policies::evaluation_error;
    using boost::math::policies::domain_error;
    using boost::math::policies::overflow_error;
    using boost::math::policies::domain_error;
    using boost::math::policies::throw_on_error;

//[lambert_w_precision_0
    std::cout.precision(std::numeric_limits<float>::max_digits10); // Show all potentially significant decimal digits,
    std::cout << std::showpoint << std::endl; // and show any significant trailing zeros too.

    double x = 10.;
    std::cout << "Lambert W (" << x << ") = " << lambert_w0(x) << std::endl; 
    // 



//] [/lambert_w_precision_0]

/*
//[lambert_w_precision_output_0

//] [/lambert_w_precision_output_0]
*/

//[lambert_w_precision_1

    using boost::math::lambert_w_detail::lambert_w_halley_step;

    double z = 10.;
    double w = lambert_w0(z);
    std::cout << "Lambert W (" << x << ") = " << lambert_w0(z) << std::endl; 
    double ww = lambert_w_halley_step(lambert_w0(z), z);
    std::cout << "Lambert W (" << x << ") = " << lambert_w_halley_step(lambert_w0(z), z) << std::endl; 

    std::cout << "difference from Halley step " << w - ww << std::endl;
//] [/lambert_w_precision_1]

      /*
//[lambert_w_precision_output_1

//] [/lambert_w_precision_output_1]

//     z = 10.F needs Halley to get the last bit or two correct.
//[lambert_w_precision_output_1a

//] [/lambert_w_precision_output_1a]
    */
  
  // Similar example using cpp_bin_float_quad (128-bit floating-point types).

  cpp_bin_float_quad zq = 10.;
  std::cout << "\nTest evaluation of cpp_bin_float_quad Lambert W(" << zq << ")"
    << std::endl;
  std::cout << std::setprecision(3) << "std::numeric_limits<cpp_bin_float_quad>::digits = " << std::numeric_limits<cpp_bin_float_quad>::digits << std::endl;
  std::cout << std::setprecision(3) << "std::numeric_limits<cpp_bin_float_quad>::epsilon() = " << std::numeric_limits<cpp_bin_float_quad>::epsilon() << std::endl;
  std::cout << std::setprecision(3) << "std::numeric_limits<cpp_bin_float_quad>::max_digits10 = " << std::numeric_limits<cpp_bin_float_quad>::max_digits10 << std::endl;
  std::cout << std::setprecision(3) << "std::numeric_limits<cpp_bin_float_quad>::digits10 = " << std::numeric_limits<cpp_bin_float_quad>::digits10 << std::endl;
  std::cout.precision(std::numeric_limits<cpp_bin_float_quad>::max_digits10);
  // All are same precision because double precision first approximation used before Halley.

  /*

      */

  { // Reference value for lambert_w0(10)
    cpp_dec_float_50 z("10");
    cpp_dec_float_50 r;
    std::cout.precision(std::numeric_limits<cpp_dec_float_50>::digits10);
    r = lambert_w0(z); // Default policy.
    std::cout << "lambert_w0(z) cpp_dec_float_50                                   = " << r << std::endl; //  0.56714329040978387299996866221035554975381578718651
    std::cout.precision(std::numeric_limits<cpp_bin_float_quad>::max_digits10);

    std::cout << "lambert_w0(z) cpp_dec_float_50 cast to quad (max_digits10(" << std::numeric_limits<cpp_bin_float_quad>::max_digits10 <<
      " )   = " << static_cast<cpp_bin_float_quad>(r) << std::endl;         // 1.7455280027406993830743012648753899115352881290809
    std::cout.precision(std::numeric_limits<cpp_bin_float_quad>::digits10); // 1.745528002740699383074301264875389837
    std::cout << "lambert_w0(z) cpp_dec_float_50 cast to quad (digits10(" << std::numeric_limits<cpp_bin_float_quad>::digits10 <<
      " )       = " << static_cast<cpp_bin_float_quad>(r) << std::endl;    // 1.74552800274069938307430126487539
    std::cout.precision(std::numeric_limits<cpp_bin_float_quad>::digits10 + 1); // 

    std::cout << "lambert_w0(z) cpp_dec_float_50 cast to quad (digits10(" << std::numeric_limits<cpp_bin_float_quad>::digits10 <<
      " )       = " << static_cast<cpp_bin_float_quad>(r) << std::endl;    // 1.74552800274069938307430126487539

    // [N[productlog[10], 50]] == 1.7455280027406993830743012648753899115352881290809

    // [N[productlog[10], 37]] == 1.745528002740699383074301264875389912
    // [N[productlog[10], 34]] == 1.745528002740699383074301264875390
    // [N[productlog[10], 33]] == 1.74552800274069938307430126487539

    // lambert_w0(z) cpp_dec_float_50 cast to quad = 1.745528002740699383074301264875389837

    // lambert_w0(z) cpp_dec_float_50 = 1.7455280027406993830743012648753899115352881290809
    // lambert_w0(z) cpp_dec_float_50 cast to quad = 1.745528002740699383074301264875389837
    // lambert_w0(z) cpp_dec_float_50 cast to quad = 1.74552800274069938307430126487539
    }

   // std::cout << "Lambert W (" << xd << ") = " << lambert_w0(xd) << std::endl; //
  }
  catch (std::exception& ex)
  {
    std::cout << ex.what() << std::endl;
  }
}  // int main()

   /*
  

  

   */


