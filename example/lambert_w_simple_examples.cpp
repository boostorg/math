// Copyright Paul A. Bristow 2016, 2017.

// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or
//  copy at http ://www.boost.org/LICENSE_1_0.txt).

// Build and run a simple examples of Lambert W function.

// Some macros that will show some(or much) diagnostic values if #defined.
//#define-able macros
//#define BOOST_MATH_INSTRUMENT_LAMBERT_W0 // W0 branch diagnostics.
//#define BOOST_MATH_INSTRUMENT_LAMBERT_Wm1 // W1 branch diagnostics.
//#define BOOST_MATH_INSTRUMENT_LAMBERT_W_HALLEY // Halley refinement diagnostics.
//#define BOOST_MATH_INSTRUMENT_LAMBERT_W_SCHROEDER // Schroeder refinement diagnostics.
//#define BOOST_MATH_INSTRUMENT_LAMBERT_W_TERMS // Number of terms used for near-singularity series.
//#define BOOST_MATH_INSTRUMENT_LAMBERT_W0_NOT_BUILTIN // higher than built-in precision types approximation and refinement.
//#define BOOST_MATH_INSTRUMENT_LAMBERT_W0_BISECTION // Show bisection only estimate.
//#define BOOST_MATH_INSTRUMENT_LAMBERT_W_SINGULARITY_SERIES // Show evaluation of series near branch singularity.
//#define BOOST_MATH_INSTRUMENT_LAMBERT_W_SMALL_Z_SERIES_ITERATIONS  // Show evaluation of series for small z.
//#define BOOST_MATH_INSTRUMENT_LAMBERT_W0_LOOKUP // Show results from lookup table.

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
#endif //# __GNUC__ and NOT _MSC_VER

#include <boost/multiprecision/cpp_dec_float.hpp> // boost::multiprecision::cpp_dec_float_50
using boost::multiprecision::cpp_dec_float_50; // 50 decimal digits type.
using boost::multiprecision::cpp_dec_float_100; // 100 decimal digits type.

#include <boost/multiprecision/cpp_bin_float.hpp>
using boost::multiprecision::cpp_bin_float_double; // == double
using boost::multiprecision::cpp_bin_float_double_extended; // 80-bit long double emulation.
using boost::multiprecision::cpp_bin_float_quad; // 128-bit quad precision.
                                                 
//[lambert_w_simple_examples_includes
#include <boost/math/special_functions/lambert_w.hpp> // For lambert_w function.

using boost::math::lambert_w0;
using boost::math::lambert_wm1;
//] //[/lambert_w_simple_examples_includes]

#include <iostream>
// using std::cout;
// using std::endl;
#include <exception>
#include <stdexcept>
#include <string>
#include <limits>  // For std::numeric_limits.

//! Show value of z to the full possibly-significant max_digits10 precision of type T.
template<typename T>
void show_value(T z)
{
  std::streamsize precision = std::cout.precision(std::numeric_limits<T>::max_digits10);  // Save.
  std::cout.precision(std::numeric_limits<T>::max_digits10); // Show all posssibly significant digits.
  std::ios::fmtflags flags(std::cout.flags());
  std::cout.setf(std::ios_base::showpoint); // Include any trailing zeros.
  std::cout << z;
  // Restore:
  std::cout.precision(precision); 
  std::cout.flags(flags);
} // template<typename T> void show_value(T z)

int main()
{
  try
  {
    std::cout << "Lambert W simple examples." << std::endl;

    using boost::math::constants::exp_minus_one; // The branch point, a singularity.
    
    // using statements needed to change precision policy.
    using boost::math::policies::policy;
    using boost::math::policies::make_policy;
    using boost::math::policies::precision;
    using boost::math::policies::digits2;
    using boost::math::policies::digits10;

    // using statements needed for changing error handling policy.
    using boost::math::policies::evaluation_error;
    using boost::math::policies::domain_error;
    using boost::math::policies::overflow_error;
    using boost::math::policies::ignore_error;
    using boost::math::policies::throw_on_error;

    //[lambert_w_simple_examples_0
    std::cout.precision(std::numeric_limits<double>::max_digits10); // Show all potentially significant decimal digits.
    std::cout << std::showpoint << std::endl; // Show trailing zeros too.
    {
      double z = 10.;
      double r;
      r = lambert_w0(z); // Default policy digits10 = 17, digits2 = 53
      std::cout << "lambert_w0(z) = " << r << std::endl; // lambert_w0(z) = 1.7455280027406992
    //] [/lambert_w_simple_examples_0]

    }
   {
  //[lambert_w_simple_examples_1
    // Other floating-point types can be used too.
    // It is convenient to use function `show_value`
    // to display all potentially significant decimal digits for the type.
    float z = 10.F;
    float r;
    r = lambert_w0(z); // Default policy digits10 = 7, digits2 = 24
    std::cout << "lambert_w0(";
    show_value(z);
    std::cout << ") = "; show_value(r);
    std::cout << std::endl; // lambert_w0(10.0000000) = 1.74552810
   //] //[/lambert_w_simple_examples_1]
    }
   { 
//[lambert_w_simple_examples_2
     // Example of an integer argument to lambert_w,
     // showing that an integer is correctly promoted to a double.
     std::cout.precision(std::numeric_limits<double>::max_digits10);
     double r = lambert_w0(10); // int argument "10" that should be promoted to double argument.
     std::cout << "lambert_w0(10) = " << r << std::endl; // lambert_w0(10) = 1.7455280027406992
     double rp = lambert_w0(10, policy<>()); // int argument "10" that should be promoted to double argument,
     // for simplicity using (default policy.
     std::cout << "lambert_w0(10, policy<>()) = " << rp << std::endl;
     // lambert_w0(10) = 1.7455280027406992
//] //[/lambert_w_simple_examples_2]
    }
   {
     //[lambert_w_simple_examples_3
     // Using multiprecision types to get much higher precision is painless.
     cpp_dec_float_50 z("10"); // Note construction using a decimal digit string, not a double literal.
     cpp_dec_float_50 r;
     r = lambert_w0(z); 
     std::cout << "lambert_w0(";
     show_value(z);
     std::cout << ") = "; show_value(r);
     std::cout << std::endl;
     // lambert_w0(10.000000000000000000000000000000000000000000000000000000000000000000000000000000) = 
     //   1.7455280027406993830743012648753899115352881290809413313533156788409111570000000
     //] //[/lambert_w_simple_examples_3]
   }
   {
     //[lambert_w_simple_examples_4
     // Using multiprecision types to get multiprecision precision wrong!
     cpp_dec_float_50 z(0.9); // Will convert to double, losing precision beyond decimal place 17!
     cpp_dec_float_50 r;
     r = lambert_w0(z);
     std::cout << "lambert_w0(";
     show_value(z);
     std::cout << ") = "; show_value(r);
     std::cout << std::endl;
     // lambert_w0(0.90000000000000002220446049250313080847263336181640625000000000000000000000000000) 
     // = 0.52983296563343442067798493880332411646762461354815017917592839379740431483218821
     std::cout << "lambert_w0(0.9) = " << lambert_w0(static_cast<double>(z))  //  0.52983296563343441
       << std::endl;
     //] //[/lambert_w_simple_examples_4]
   }
   {
     //[lambert_w_simple_examples_4a
     // Using multiprecision types to get multiprecision precision right!
     cpp_dec_float_50 z("0.9"); // Construct from decimal digit string.
     cpp_dec_float_50 r;
     r = lambert_w0(z);
     std::cout << "lambert_w0(";
     show_value(z);
     std::cout << ") = "; show_value(r);
     std::cout << std::endl;
     // 0.90000000000000000000000000000000000000000000000000000000000000000000000000000000)
     // = 0.52983296563343441213336643954546304857788132269804249284012528304239956432480864
     //] //[/lambert_w_simple_examples_4a]
   }
   {
//[lambert_w_simple_examples_precision_policies
      double z = 5.;
    // Define a new, non-default, policy to calculate to precision of approximately 3 decimal digits.
      std::cout << "lambert_w0(" << z << ", policy<digits10<3>>()) = " << lambert_w0(z, policy<digits10<3>>())
        << std::endl; // lambert_w0(0.29999999999999999, policy<digits10<3>>()) = 0.23675531078855935
//] //[/lambert_w_simple_examples_precision_policies]

//[lambert_w_simple_examples_error_policies
      // Define an error handling policy:
      typedef policy<
        domain_error<throw_on_error>,
        overflow_error<ignore_error>
      > throw_policy;

      std::cout << "Lambert W (" << z << ") = " << lambert_w0(z) << std::endl; // 0.23675531078855930
      std::cout << "\nLambert W (" << z << ", throw_policy()) = " << lambert_w0(z, throw_policy()) << std::endl;
    //] //[/lambert_w_simple_example_error_policies]
    }
    {
 //[lambert_w_simple_examples_out_of_range
      // Show error reporting if pass a value to lambert_w0 that is out of range,
      // and probably was meant to be passed to lambert_m1 instead.
      double z = -1.;
      double r = lambert_w0(z);
      std::cout << "lambert_w0(-1.) = " << r << std::endl;
      // Error in function boost::math::lambert_w0<RealType>(<RealType>):
      //  Argument z = %1 out of range (-1/e <= z < (std::numeric_limits<T>::max)())
      //  for Lambert W0 branch (use W-1 branch?).
//] [/lambert_w_simple_examples_out_of_range]
    }
  }
  catch (std::exception& ex)
  {
    std::cout << ex.what() << std::endl;
  }
}  // int main()
 
   /*

   Output:
//[lambert_w_simple_examples_error_message_1
   Error in function boost::math::lambert_w0<RealType>(<RealType>): 
   Argument z = %1 out of range (-1/e <= z < (std::numeric_limits<T>::max)()) 
   for Lambert W0 branch 
   (use W-1 branch?).
//] [/lambert_w_simple_examples_error_message_1]

   /*


   */


