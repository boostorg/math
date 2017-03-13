// Copyright Paul A. Bristow 2016.

// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or
//  copy at http ://www.boost.org/LICENSE_1_0.txt).

// Test that can build and run a simple example of Lambert W function,
// using algorithm of Thomas Luu.
// https://svn.boost.org/trac/boost/ticket/11027

#include <boost/config.hpp> // for BOOST_PLATFORM, BOOST_COMPILER,  BOOST_STDLIB ...
#include <boost/version.hpp>   // for BOOST_MSVC versions.
#include <boost/cstdint.hpp>
#include <boost/exception/exception.hpp>  // boost::exception
#include <boost/math/constants/constants.hpp> // For exp_minus_one == 3.67879441171442321595523770161460867e-01.


// #define BOOST_MATH_INSTRUMENT_LAMBERT_W 

// For lambert_w function.
#include <boost/math/special_functions/lambert_w.hpp>

#include <iostream>
// using std::cout;
// using std::endl;
#include <exception>
#include <stdexcept>
#include <string>
#include <limits>  // For std::numeric_limits.

//! Show information about build, architecture, address model, platform, ...
std::string show_versions()
{
  std::ostringstream message;

  message << "Program: " << __FILE__ << "\n";
#ifdef __TIMESTAMP__
  message << __TIMESTAMP__;
#endif
  message << "\nBuildInfo:\n" "  Platform " << BOOST_PLATFORM;
  // http://stackoverflow.com/questions/1505582/determining-32-vs-64-bit-in-c
#if defined(__LP64__) || defined(_WIN64) || (defined(__x86_64__) && !defined(__ILP32__) ) || defined(_M_X64) || defined(__ia64) || defined (_M_IA64) || defined(__aarch64__) || defined(__powerpc64__)
#define IS64BIT 1
  message << ", 64-bit.";
#else
#define IS32BIT 1
  message << ", 32-bit.";
#endif

  message << "\n  Compiler " BOOST_COMPILER;
#ifdef BOOST_MSC_VER
#ifdef _MSC_FULL_VER
  message << "\n  MSVC version " << BOOST_STRINGIZE(_MSC_FULL_VER) << ".";
#endif
#ifdef __WIN64
  mess age << "\n WIN64" << std::endl;
#endif // __WIN64
#ifdef _WIN32
  message << "\n WIN32" << std::endl;
#endif  // __WIN32
#endif
#ifdef __GNUC__
  //PRINT_MACRO(__GNUC__);
  //PRINT_MACRO(__GNUC_MINOR__);
  //PRINT_MACRO(__GNUC_PATCH__);
  std::cout << "GCC " << __VERSION__ << std::endl;
  //PRINT_MACRO(LONG_MAX);
#endif // __GNUC__

  message << "\n  STL " << BOOST_STDLIB;

  message << "\n  Boost version " << BOOST_VERSION / 100000 << "." << BOOST_VERSION / 100 % 1000 << "." << BOOST_VERSION % 100;

#ifdef BOOST_HAS_FLOAT128
  message << ",  BOOST_HAS_FLOAT128" << std::endl;
#endif
  message << std::endl;
  return message.str();
} // std::string versions()


int main()
{
  try
  {
    //std::cout << "Lambert W example basic!" << std::endl;
    //std::cout << show_versions() << std::endl;
     
    //std::cout << exp(1) << std::endl; // 2.71828
    //std::cout << exp(-1) << std::endl; // 0.367879
    //std::cout << std::numeric_limits<double>::epsilon() / 2 << std::endl; // 1.11022e-16

    using namespace boost::math;
    using boost::math::constants::exp_minus_one;
    double x = 1.;

    double W1 = productlog(1.);
    // Note, NOT integer X, for example: productlog(1); or will get message like
    // error C2338: Must be floating-point, not integer type, for example W(1.), not W(1)!
    //

    std::cout << "Lambert W (" << x << ") = " << productlog(x) << std::endl; // 0.567143
    // This 'golden ratio' for exponentials is http://mathworld.wolfram.com/OmegaConstant.html
    // since exp[-W(1)] = W(1)
    // A030178		Decimal expansion of LambertW(1): the solution to x*exp(x)
    // = 0.5671432904097838729999686622103555497538157871865125081351310792230457930866
      // http://oeis.org/A030178

    double expplogone = exp(-productlog(1.));
    if (expplogone != W1)
    {
      std::cout << expplogone << " " << W1 << std::endl; //
    }


//[lambert_w_example_1
    /*
    x = 0.01;
    std::cout << "Lambert W (" << x << ") = " << productlog(x) << std::endl; // 0.00990147
//] [/lambert_w_example_1]
    x = -0.01;
    std::cout << "Lambert W (" << x << ") = " << productlog(x) << std::endl; // -0.0101015
    x = -0.1;
    std::cout << "Lambert W (" << x << ") = " << productlog(x) << std::endl; //
*/

    /*
    // Near singularity -exp(1) == -0.367879441171442321595523770161460867445811131031767834
    x = -0.367879; // Near -exp(1) = -0.367879
    std::cout << "Lambert W (" << x << ") = " << productlog(x) << std::endl; // -0.99845210378080340
    // N[productlog(-0.367879), 50] = -0.99845210378072725931829498030640227316856835774851
    x = -exp_minus_one<double>(); // At exactly -exp(1). 
    std::cout << "Lambert W (" << x << ") = " << productlog(x) << std::endl; // -1.0000000000000

    x += std::numeric_limits<double>::epsilon() ; // At nearly -exp(1). 
    std::cout << "Lambert W (" << x << ") = " << productlog(x) << std::endl; // 0.0000000000000
  
    x = -exp_minus_one<double>() + 100 * std::numeric_limits<double>::epsilon() ; // At nearly -exp(1). 
    std::cout << "Lambert W (" << x << ") = " << productlog(x) << std::endl; // 0.0000000000000
 */
     // http://www.wolframalpha.com/input/?i=N%5Bproductlog%5B-0.367879%5D,17%5D test value N[productlog[-0.367879],17]
  //  -0.367879441171442321595523770161460867445811131031767834
    x = -0.367879; // < -exp(1) = -0.367879
    std::cout << "Lambert W (" << x << ") = " << productlog(x) << std::endl; // Lambert W (-0.36787900000000001) = -0.99845210378080340
    //  -0.99845210378080340
    //  -0.99845210378072726  N[productlog[-0.367879],17] wolfram  so very close.

    x = -0.3678794; // expect -0.99952696660756813
    std::cout << "Lambert W (" << x << ") = " << productlog(x) << std::endl; // 0.0
    x = -0.36787944; // expect -0.99992019848408340
    std::cout << "Lambert W (" << x << ") = " << productlog(x) << std::endl; // 0.0
    x = -0.367879441; // -0.99996947070054883
    std::cout << "Lambert W (" << x << ") = " << productlog(x) << std::endl; // 0.0
    x = -0.36787944117; // -0.99999719977527159
    std::cout << "Lambert W (" << x << ") = " << productlog(x) << std::endl; // 0.0
    x = -0.367879441171; // -0.99999844928821992
    std::cout << "Lambert W (" << x << ") = " << productlog(x) << std::endl; // 0.0
    x = -0.3678794411714; // N[productlog[-0.3678794411714],17] == -0.99999952032930614  
    std::cout << "Lambert W (" << x << ") = " << productlog(x) << std::endl; // 0.0
    /* */
  }
  catch (std::exception& ex)
  {
    std::cout << ex.what() << std::endl;
  }


}  // int main()

   /*

//[lambert_w_output_1
   Output:

  1>  example_basic.cpp
1>  Generating code
1>  All 237 functions were compiled because no usable IPDB/IOBJ from previous compilation was found.
1>  Finished generating code
1>  LambertW.vcxproj -> J:\Cpp\Misc\x64\Release\LambertW.exe
1>  LambertW.vcxproj -> J:\Cpp\Misc\x64\Release\LambertW.pdb (Full PDB)
1>  Lambert W example basic!
1>  Platform: Win32
1>  Compiler: Microsoft Visual C++ version 14.0
1>  STL     : Dinkumware standard library version 650
1>  Boost   : 1.63.0
1>  _MSC_FULL_VER = 190024123
1>  Win32
1>  x64
1>   (x64)
1>  Iteration #0, w0 0.577547206058041, w1 = 0.567143616915443, difference = 0.0289944962755619, relative 0.018343835374856
1>  Iteration #1, w0 0.567143616915443, w1 = 0.567143290409784, difference = 9.02208135089566e-07, relative 5.75702234328901e-07
1>  Final 0.567143290409784 after 2 iterations, difference = 0
1>  Iteration #0, w0 0.577547206058041, w1 = 0.567143616915443, difference = 0.0289944962755619, relative 0.018343835374856
1>  Iteration #1, w0 0.567143616915443, w1 = 0.567143290409784, difference = 9.02208135089566e-07, relative 5.75702234328901e-07
1>  Final 0.567143290409784 after 2 iterations, difference = 0
1>  Lambert W (1) = 0.567143290409784
1>  Iteration #0, w0 0.577547206058041, w1 = 0.567143616915443, difference = 0.0289944962755619, relative 0.018343835374856
1>  Iteration #1, w0 0.567143616915443, w1 = 0.567143290409784, difference = 9.02208135089566e-07, relative 5.75702234328901e-07
1>  Final 0.567143290409784 after 2 iterations, difference = 0
1>  Iteration #0, w0 0.0099072820916067, w1 = 0.00990147384359511, difference = 5.92416060777624e-06, relative 0.000586604388734591
1>  Final 0.00990147384359511 after 1 iterations, difference = 0
1>  Lambert W (0.01) = 0.00990147384359511
1>  Iteration #0, w0 -0.0101016472705154, w1 = -0.0101015271985388, difference = -1.17664437923951e-07, relative 1.18865171889748e-05
1>  Final -0.0101015271985388 after 1 iterations, difference = 0
1>  Lambert W (-0.01) = -0.0101015271985388
1>  Iteration #0, w0 -0.111843322610692, w1 = -0.111832559158964, difference = -8.54817065376601e-06, relative 9.62461362694622e-05
1>  Iteration #1, w0 -0.111832559158964, w1 = -0.111832559158963, difference = -5.68989300120393e-16, relative 6.43929354282591e-15
1>  Final -0.111832559158963 after 2 iterations, difference = 0
1>  Lambert W (-0.1) = -0.111832559158963
1>  Iteration #0, w0 -0.998452103785573, w1 = -0.998452103780803, difference = -2.72004641033163e-15, relative 4.77662354114727e-12
1>  Final -0.998452103780803 after 1 iterations, difference = 0
1>  Lambert W (-0.367879) = -0.998452103780803

//] [/lambert_w_output_1]
   */
