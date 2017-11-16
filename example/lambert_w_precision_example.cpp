// Copyright Paul A. Bristow 2016.

// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or
//  copy at http ://www.boost.org/LICENSE_1_0.txt).

//! Lambert W examples of controlling precision using Boost.Math policies
//! (and thus also trading less precision for shorter run time).

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
// using std::cout;
// using std::endl;
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

    // Define an error handling policy:
    //    typedef policy<
    //      domain_error<throw_on_error>,
    //      overflow_error<throw_on_error>
    //    > throw_policy;

    //float x = 0.3F;
    // std::cout << "Lambert W (" << x << ") = " << lambert_w(x) << std::endl; // 0.00990147
    // std::cout << "\nLambert W (" << x << ") = " << lambert_w0(x, throw_policy()) << std::endl;
    // std::cout << "Lambert W (" << x << ") = " << lambert_w(x, throw_policy()) << std::endl;

    // Examples of controlling precision using policies.
    using boost::math::policies::policy;
    using boost::math::policies::precision;
    using boost::math::policies::digits10;
    using boost::math::policies::digits2;

    // Creating a policy using an integer literal.
    // typedef policy<digits10<3> > my_pol_3_dec; // Define a new, non-default, policy to calculate to accuracy of approximately 3 decimal digits.

    // Creating a policy using an integer variable (must be const integral, usually const int !).
    // "a variable with non-static storage duration cannot be used as a non-type argument"

    // typedef policy<digits10<decimal_digits> > my_pol_3_dec;
    // Define a new, non-default, policy to calculate to accuracy of approximately 3 decimal digits.
    const int decimal_digits = 3;
    typedef policy<digits10<decimal_digits> > my_pol_3_dec; // Define a custom, non-default, policy to calculate to accuracy of approximately 3 decimal digits.

    std::cout << std::setprecision(3) << "std::numeric_limits<double>::epsilon() = " << std::numeric_limits<double>::epsilon() << std::endl;

    std::cout.precision(std::numeric_limits<double>::max_digits10);
    double x = 1.;
    std::cout << "lambert_w0(x, my_pol_3_dec()) = " << lambert_w0(x, my_pol_3_dec()) <<"\n" << std::endl; // = 0.56714329040978395

    // Can also specify precision using bits, for example:
    typedef policy<digits2<11> > my_prec_11_bits;  // Policy to require only approximately 11 bits of precision.

    std::cout << "lambert_w0(x, my_prec_11_bits()) = " << lambert_w0(x, my_prec_11_bits()) <<"\n" << std::endl; //  0.56714329040978395

    { // Show various precisions that show varying improvements.
      // Specify precision as n bits, base-2, policy<digits2<n> >()

      std::cout << std::setprecision(3)
        << "\nstd::numeric_limits<float>::epsilon() = " << std::numeric_limits<float>::epsilon()
        << "\nstd::numeric_limits<float>::digits = " << std::numeric_limits<float>::digits
        << "\nstd::numeric_limits<float>::digits10 = " << std::numeric_limits<float>::digits10
        << "\nstd::numeric_limits<float>::max_digits10 = " << std::numeric_limits<float>::max_digits10
        << std::endl;

////[lambert_w_precision_0
//      float z(10.F);
//      float r = 0;
//      std::cout.precision(std::numeric_limits<float>::max_digits10);  // Show all possibly significant decimal digits.
//      std::cout << std::showpoint << std::endl; // Show any trailing zeros.
//
//      r = lambert_w0(z, policy<digits2<9> >()); //
//      std::cout << "lambert_w0(z, policy<digits2<9> >())  = " << r << std::endl; //  0.000000000 (nearest using lookup).
//
//      r = lambert_w0(z, policy<digits2<10> >()); // // Just using bisection.
//      std::cout << "lambert_w0(z, policy<digits2<10> >()) = " << r << std::endl; // 0.566406250 (nearest using bisection).
//
//      r = lambert_w0(z, policy<digits2<21> >()); // digits10 = 6, digits2 = 22 - so just using Schroeder refinement.
//      std::cout << "lambert_w0(z, policy<digits2<21> >()) = " << r << std::endl; // 0.567143261
//
//      // Can also specify precision as decimal base-10 digits using policy<digits10<5> >().
//      r = lambert_w0(z, policy<digits10<6> >()); // digits10 = 5, digits2 = 22  Using Schroeder.
//      std::cout << "lambert_w0(z, policy<digits10<5> >()) = " << r << std::endl; // 0.567143261
//
//      r = lambert_w0(z, policy<digits2<23> >()); // digits10 = 7, digits2 = 23 - so using Halley as well (but no improvement).
//      std::cout << "lambert_w0(z, policy<digits2<23> >()) = " << r << std::endl;  // 0.567143261
//
//      r = lambert_w0(z, policy<>()); // Default policy (using lookup, bisection,  Schroeder, and Halley).
//      std::cout << "lambert_w0(z, policy<>()) = " << r << std::endl;  // 0.567143261
//
////] [/lambert_w_precision_0]

/*
//[lambert_w_precision_output_0
      lambert_w0(z, policy<digits2<9> >()) = 0.000000000   // Nearby integer is zero from lookup table.
      lambert_w0(z, policy<digits2<10> >()) = 0.566406250  // Just using bisection from lookup integral values.
      lambert_w0(z, policy<digits2<21> >()) = 0.567143261  // digits10 = 6, digits2 = 22 - so using Schroeder refinement.
      lambert_w0(z, policy<digits10<5> >()) = 0.567143261  // digits10 = 5, digits2 = 19  Using Schroeder.
      lambert_w0(z, policy<digits2<23> >()) = 0.567143261  // digits10 = 7, digits2 = 23 - so using Halley as well (but no improvement).
      lambert_w0(z, policy<>()) = 0.567143261
//] [/lambert_w_precision_output_0]
*/

    }

    { 
      // Test evaluation of Lambert W at varying precisions specified using a policy with digits10 increasing from 1 to 9.
      // For float this covers the full precision range as max_digits10 = 9.
      float z = 10.F;
      std::cout << "\nTest evaluation of float Lambert W (" << z << ") at varying precisions"
        "\nspecified using a policy with digits10 increasing from 1 to 9. "
        << std::endl;
//[lambert_w_precision_1
      std::cout << std::setprecision(3) << "std::numeric_limits<float>::epsilon() = "
        << std::numeric_limits<float>::epsilon() << std::endl;
      std::cout.precision(std::numeric_limits<float>::max_digits10);

      std::cout << "Lambert W (" << z << ", digits10<0>) = " << lambert_w0(z, policy<digits10<0> >()) << std::endl;
      std::cout << "Lambert W (" << z << ", digits10<1>) = " << lambert_w0(z, policy<digits10<1> >()) << std::endl;
      std::cout << "Lambert W (" << z << ", digits10<2>) = " << lambert_w0(z, policy<digits10<2> >()) << std::endl;
      std::cout << "Lambert W (" << z << ", digits10<3>) = " << lambert_w0(z, policy<digits10<3> >()) << std::endl;
      std::cout << "Lambert W (" << z << ", digits10<4>) = " << lambert_w0(z, policy<digits10<4> >()) << std::endl;
      std::cout << "Lambert W (" << z << ", digits10<5>) = " << lambert_w0(z, policy<digits10<5> >()) << std::endl;
      std::cout << "Lambert W (" << z << ", digits10<6>) = " << lambert_w0(z, policy<digits10<6> >()) << std::endl;
      std::cout << "Lambert W (" << z << ", digits10<7>) = " << lambert_w0(z, policy<digits10<7> >()) << std::endl;
      std::cout << "Lambert W (" << z << ", digits10<8>) = " << lambert_w0(z, policy<digits10<8> >()) << std::endl;
      std::cout << "Lambert W (" << z << ", digits10<9>) = " << lambert_w0(z, policy<digits10<9> >()) << std::endl;
      std::cout << "Lambert W (" << z << ") default      = " << lambert_w0(z) << std::endl;
//] [/lambert_w_precision_1]

      /*
//[lambert_w_precision_output_1
      Test evaluation of float Lambert W at varying precisions
      specified using a policy with digits10 increasing from 1 to 9.
      std::numeric_limits<float>::epsilon() = 1.19e-07
      z = 1.F - note that gets to 'exact' result immediately using Schroeder refinement,
      so no advantage using Halley.

      Lambert W (1.00000000, digits10<1>) = 0.000000000  << using lookup.
      Lambert W (1.00000000, digits10<2>) = 0.000000000  << using lookup and bisection.
      Lambert W (1.00000000, digits10<3>) = 0.567143261  << Schroeder refinement, and is sufficient.
      Lambert W (1.00000000, digits10<4>) = 0.567143261
      Lambert W (1.00000000, digits10<5>) = 0.567143261 << Halley refinement, but no improvement.
      Lambert W (1.00000000, digits10<6>) = 0.567143261
      Lambert W (1.00000000, digits10<7>) = 0.567143261
      Lambert W (1.00000000, digits10<8>) = 0.567143261
      Lambert W (1.00000000, digits10<9>) = 0.567143261
      Lambert W (1.00000000) = 0.567143261

//] [/lambert_w_precision_output_1]

//     z = 10.F needs Halley to get the last bit or two correct.
//[lambert_w_precision_output_1a

      z = 10.F needs Halley to get the last bit or two correct.
      Lambert W (10.0000000, digits10<1>) = 1.35914087 << bisection.
      Lambert W (10.0000000, digits10<2>) = 1.35914087
      Lambert W (10.0000000, digits10<3>) = 1.74552798 << Schroeder refinement.
      Lambert W (10.0000000, digits10<4>) = 1.74552798
      Lambert W (10.0000000, digits10<5>) = 1.74552798
      Lambert W (10.0000000, digits10<6>) = 1.74552810 << Halley refinement.
      Lambert W (10.0000000, digits10<7>) = 1.74552810
      Lambert W (10.0000000, digits10<8>) = 1.74552810
      Lambert W (10.0000000, digits10<9>) = 1.74552810
      Lambert W (10.0000000)              = 1.74552810
//] [/lambert_w_precision_output_1a]
    */
    }
    { // similar example using float but specifying with digits2

      float z = 10.;
      std::cout << "\nTest evaluation of float Lambert W (" << z << ") at varying precisions"
        "\nspecified using a policy with (base_2) digits increasing from 1 to 24. "
        << std::endl;
      std::cout << std::setprecision(3) << "std::numeric_limits<float>::digits = " << std::numeric_limits<float>::digits << std::endl;
      std::cout << std::setprecision(3) << "std::numeric_limits<float>::epsilon() = " << std::numeric_limits<float>::epsilon() << std::endl;
      std::cout.precision(std::numeric_limits<float>::max_digits10);
      std::cout << "Lambert W (" << z << ", digits2<0> ) = " << lambert_w0(z, policy<digits2<0> >()) << std::endl;
      std::cout << "Lambert W (" << z << ", digits2<1> ) = " << lambert_w0(z, policy<digits2<1> >()) << std::endl;
      std::cout << "Lambert W (" << z << ", digits2<2> ) = " << lambert_w0(z, policy<digits2<2> >()) << std::endl;
      std::cout << "Lambert W (" << z << ", digits2<3> ) = " << lambert_w0(z, policy<digits2<3> >()) << std::endl;
      std::cout << "Lambert W (" << z << ", digits2<4> ) = " << lambert_w0(z, policy<digits2<4> >()) << std::endl;
      std::cout << "Lambert W (" << z << ", digits2<5> ) = " << lambert_w0(z, policy<digits2<5> >()) << std::endl;
      std::cout << "Lambert W (" << z << ", digits2<6> ) = " << lambert_w0(z, policy<digits2<6> >()) << std::endl;
      std::cout << "Lambert W (" << z << ", digits2<7> ) = " << lambert_w0(z, policy<digits2<7> >()) << std::endl;
      std::cout << "Lambert W (" << z << ", digits2<8> ) = " << lambert_w0(z, policy<digits2<8> >()) << std::endl;
      std::cout << "Lambert W (" << z << ", digits2<9> ) = " << lambert_w0(z, policy<digits2<9> >()) << std::endl;
      std::cout << "Lambert W (" << z << ", digits2<10>) = " << lambert_w0(z, policy<digits2<10> >()) << std::endl;
      std::cout << "Lambert W (" << z << ", digits2<11>) = " << lambert_w0(z, policy<digits2<11> >()) << std::endl;
      std::cout << "Lambert W (" << z << ", digits2<12>) = " << lambert_w0(z, policy<digits2<12> >()) << std::endl;
      std::cout << "Lambert W (" << z << ", digits2<13>) = " << lambert_w0(z, policy<digits2<13> >()) << std::endl;
      std::cout << "Lambert W (" << z << ", digits2<14>) = " << lambert_w0(z, policy<digits2<14> >()) << std::endl;
      std::cout << "Lambert W (" << z << ", digits2<15>) = " << lambert_w0(z, policy<digits2<15> >()) << std::endl;
      std::cout << "Lambert W (" << z << ", digits2<16>) = " << lambert_w0(z, policy<digits2<16> >()) << std::endl;
      std::cout << "Lambert W (" << z << ", digits2<17>) = " << lambert_w0(z, policy<digits2<17> >()) << std::endl;
      std::cout << "Lambert W (" << z << ", digits2<18>) = " << lambert_w0(z, policy<digits2<18> >()) << std::endl;
      std::cout << "Lambert W (" << z << ", digits2<19>) = " << lambert_w0(z, policy<digits2<19> >()) << std::endl;
      std::cout << "Lambert W (" << z << ", digits2<20>) = " << lambert_w0(z, policy<digits2<20> >()) << std::endl;
      std::cout << "Lambert W (" << z << ", digits2<21>) = " << lambert_w0(z, policy<digits2<21> >()) << std::endl;
      std::cout << "Lambert W (" << z << ", digits2<22>) = " << lambert_w0(z, policy<digits2<22> >()) << std::endl;
      std::cout << "Lambert W (" << z << ", digits2<23>) = " << lambert_w0(z, policy<digits2<23> >()) << std::endl;
      std::cout << "Lambert W (" << z << ", digits2<24>) = " << lambert_w0(z, policy<digits2<24> >()) << std::endl;

      std::cout << "Lambert W (" << z << ") default      = " << lambert_w0(z) << std::endl; // Default policy.

      /* This example shows a tiny improvement using Halley after Schroeder.
        Lambert W (10.0000000, digits2<1>) = 1.35914087 
        Lambert W (10.0000000, digits2<2>) = 1.35914087
        Lambert W (10.0000000, digits2<3>) = 1.35914087
        Lambert W (10.0000000, digits2<4>) = 1.35914087
        Lambert W (10.0000000, digits2<5>) = 1.35914087
        Lambert W (10.0000000, digits2<6>) = 1.35914087
        Lambert W (10.0000000, digits2<7>) = 1.35914087
        Lambert W (10.0000000, digits2<8>) = 1.35914087
        Lambert W (10.0000000, digits2<9>) = 1.35914087
        Lambert W (10.0000000, digits2<10>) = 1.74414063 << Bisection.
        Lambert W (10.0000000, digits2<11>) = 1.74552798 << Schroeder.
        Lambert W (10.0000000, digits2<12>) = 1.74552798
        Lambert W (10.0000000, digits2<13>) = 1.74552798
        Lambert W (10.0000000, digits2<14>) = 1.74552798
        Lambert W (10.0000000, digits2<15>) = 1.74552798
        Lambert W (10.0000000, digits2<16>) = 1.74552798
        Lambert W (10.0000000, digits2<17>) = 1.74552798
        Lambert W (10.0000000, digits2<18>) = 1.74552798
        Lambert W (10.0000000, digits2<19>) = 1.74552798
        Lambert W (10.0000000, digits2<20>) = 1.74552798
        Lambert W (10.0000000, digits2<21>) = 1.74552798
        Lambert W (10.0000000, digits2<22>) = 1.74552810 << Halley.
        Lambert W (10.0000000, digits2<23>) = 1.74552810
        Lambert W (10.0000000, digits2<24>) = 1.74552810
        Lambert W (10.0000000) = 1.74552810 << default policy precision.
  */
    }

    { // Similar example using double but specifying with digits2.
      double z = 10.;
      std::cout << "\nTest evaluation of double Lambert W (" << z << ") at varying precisions"
        "\nspecified using a policy with (base_2) digits increasing from 1 to 54. "
        << std::endl;
      std::cout << std::setprecision(3) << "std::numeric_limits<double>::epsilon() = " << std::numeric_limits<double>::epsilon() << std::endl;
      std::cout << std::setprecision(3) << "std::numeric_limits<double>::digits = " << std::numeric_limits<double>::digits << std::endl;

      std::cout.precision(std::numeric_limits<double>::max_digits10);
      std::cout << "Lambert W (" << z << ", digits2<1> ) = " << lambert_w0(z, policy<digits2<1> >()) << std::endl;
      std::cout << "Lambert W (" << z << ", digits2<2> ) = " << lambert_w0(z, policy<digits2<2> >()) << std::endl;
      std::cout << "Lambert W (" << z << ", digits2<3> ) = " << lambert_w0(z, policy<digits2<3> >()) << std::endl;
      std::cout << "Lambert W (" << z << ", digits2<4> ) = " << lambert_w0(z, policy<digits2<4> >()) << std::endl;
      std::cout << "Lambert W (" << z << ", digits2<5> ) = " << lambert_w0(z, policy<digits2<5> >()) << std::endl;
      std::cout << "Lambert W (" << z << ", digits2<6> ) = " << lambert_w0(z, policy<digits2<6> >()) << std::endl;
      std::cout << "Lambert W (" << z << ", digits2<7> ) = " << lambert_w0(z, policy<digits2<7> >()) << std::endl;
      std::cout << "Lambert W (" << z << ", digits2<8> ) = " << lambert_w0(z, policy<digits2<8> >()) << std::endl;
      std::cout << "Lambert W (" << z << ", digits2<9> ) = " << lambert_w0(z, policy<digits2<9> >()) << std::endl;
      std::cout << "Lambert W (" << z << ", digits2<10>) = " << lambert_w0(z, policy<digits2<10> >()) << std::endl;
      std::cout << "Lambert W (" << z << ", digits2<11>) = " << lambert_w0(z, policy<digits2<11> >()) << std::endl;
      std::cout << "Lambert W (" << z << ", digits2<12>) = " << lambert_w0(z, policy<digits2<12> >()) << std::endl;
      std::cout << "Lambert W (" << z << ", digits2<13>) = " << lambert_w0(z, policy<digits2<13> >()) << std::endl;
      std::cout << "Lambert W (" << z << ", digits2<14>) = " << lambert_w0(z, policy<digits2<14> >()) << std::endl;
      std::cout << "Lambert W (" << z << ", digits2<15>) = " << lambert_w0(z, policy<digits2<15> >()) << std::endl;
      std::cout << "Lambert W (" << z << ", digits2<16>) = " << lambert_w0(z, policy<digits2<16> >()) << std::endl;
      std::cout << "Lambert W (" << z << ", digits2<17>) = " << lambert_w0(z, policy<digits2<17> >()) << std::endl;
      std::cout << "Lambert W (" << z << ", digits2<18>) = " << lambert_w0(z, policy<digits2<18> >()) << std::endl;
      std::cout << "Lambert W (" << z << ", digits2<19>) = " << lambert_w0(z, policy<digits2<19> >()) << std::endl;
      std::cout << "Lambert W (" << z << ", digits2<20>) = " << lambert_w0(z, policy<digits2<20> >()) << std::endl;
      std::cout << "Lambert W (" << z << ", digits2<21>) = " << lambert_w0(z, policy<digits2<21> >()) << std::endl;
      std::cout << "Lambert W (" << z << ", digits2<22>) = " << lambert_w0(z, policy<digits2<22> >()) << std::endl;
      std::cout << "Lambert W (" << z << ", digits2<23>) = " << lambert_w0(z, policy<digits2<23> >()) << std::endl;
      std::cout << "Lambert W (" << z << ", digits2<24>) = " << lambert_w0(z, policy<digits2<24> >()) << std::endl;
      std::cout << "Lambert W (" << z << ", digits2<25>) = " << lambert_w0(z, policy<digits2<25> >()) << std::endl;
      std::cout << "Lambert W (" << z << ", digits2<26>) = " << lambert_w0(z, policy<digits2<26> >()) << std::endl;
      std::cout << "Lambert W (" << z << ", digits2<27>) = " << lambert_w0(z, policy<digits2<27> >()) << std::endl;
      std::cout << "Lambert W (" << z << ", digits2<28>) = " << lambert_w0(z, policy<digits2<28> >()) << std::endl;
      std::cout << "Lambert W (" << z << ", digits2<29>) = " << lambert_w0(z, policy<digits2<29> >()) << std::endl;
      std::cout << "Lambert W (" << z << ", digits2<30>) = " << lambert_w0(z, policy<digits2<30> >()) << std::endl;
      std::cout << "Lambert W (" << z << ", digits2<31>) = " << lambert_w0(z, policy<digits2<31> >()) << std::endl;
      std::cout << "Lambert W (" << z << ", digits2<32>) = " << lambert_w0(z, policy<digits2<32> >()) << std::endl;
      std::cout << "Lambert W (" << z << ", digits2<33>) = " << lambert_w0(z, policy<digits2<33> >()) << std::endl;
      std::cout << "Lambert W (" << z << ", digits2<34>) = " << lambert_w0(z, policy<digits2<34> >()) << std::endl;
      std::cout << "Lambert W (" << z << ", digits2<35>) = " << lambert_w0(z, policy<digits2<35> >()) << std::endl;
      std::cout << "Lambert W (" << z << ", digits2<36>) = " << lambert_w0(z, policy<digits2<36> >()) << std::endl;
      std::cout << "Lambert W (" << z << ", digits2<37>) = " << lambert_w0(z, policy<digits2<37> >()) << std::endl;
      std::cout << "Lambert W (" << z << ", digits2<38>) = " << lambert_w0(z, policy<digits2<38> >()) << std::endl;
      std::cout << "Lambert W (" << z << ", digits2<39>) = " << lambert_w0(z, policy<digits2<39> >()) << std::endl;
      std::cout << "Lambert W (" << z << ", digits2<40>) = " << lambert_w0(z, policy<digits2<40> >()) << std::endl;
      std::cout << "Lambert W (" << z << ", digits2<41>) = " << lambert_w0(z, policy<digits2<41> >()) << std::endl;
      std::cout << "Lambert W (" << z << ", digits2<42>) = " << lambert_w0(z, policy<digits2<42> >()) << std::endl;
      std::cout << "Lambert W (" << z << ", digits2<43>) = " << lambert_w0(z, policy<digits2<43> >()) << std::endl;
      std::cout << "Lambert W (" << z << ", digits2<44>) = " << lambert_w0(z, policy<digits2<44> >()) << std::endl;
      std::cout << "Lambert W (" << z << ", digits2<45>) = " << lambert_w0(z, policy<digits2<45> >()) << std::endl;
      std::cout << "Lambert W (" << z << ", digits2<46>) = " << lambert_w0(z, policy<digits2<46> >()) << std::endl;
      std::cout << "Lambert W (" << z << ", digits2<47>) = " << lambert_w0(z, policy<digits2<47> >()) << std::endl;
      std::cout << "Lambert W (" << z << ", digits2<48>) = " << lambert_w0(z, policy<digits2<48> >()) << std::endl;
      std::cout << "Lambert W (" << z << ", digits2<49>) = " << lambert_w0(z, policy<digits2<49> >()) << std::endl;
      std::cout << "Lambert W (" << z << ", digits2<50>) = " << lambert_w0(z, policy<digits2<50> >()) << std::endl;
      std::cout << "Lambert W (" << z << ", digits2<51>) = " << lambert_w0(z, policy<digits2<51> >()) << std::endl;
      std::cout << "Lambert W (" << z << ", digits2<52>) = " << lambert_w0(z, policy<digits2<52> >()) << std::endl;
      std::cout << "Lambert W (" << z << ", digits2<53>) = " << lambert_w0(z, policy<digits2<53> >()) << std::endl;
      std::cout << "Lambert W (" << z << ", digits2<54>) = " << lambert_w0(z, policy<digits2<54> >()) << std::endl;

      std::cout << "Lambert W (" << z << ") default      = " << lambert_w0(z) << std::endl; // Default policy.

      /*
      Test evaluation of double Lambert W at varying precisions
      specified using a policy with (base_2) digits increasing from 1 to 54.
      std::numeric_limits<double>::epsilon() = 2.22e-16
      Lambert W (10.000000000000000, digits2<1>) = 1.3591409142295225
      Lambert W (10.000000000000000, digits2<2>) = 1.3591409142295225
      Lambert W (10.000000000000000, digits2<3>) = 1.3591409142295225
      Lambert W (10.000000000000000, digits2<4>) = 1.3591409142295225
      Lambert W (10.000000000000000, digits2<5>) = 1.3591409142295225
      Lambert W (10.000000000000000, digits2<6>) = 1.3591409142295225
      Lambert W (10.000000000000000, digits2<7>) = 1.3591409142295225
      Lambert W (10.000000000000000, digits2<8>) = 1.3591409142295225
      Lambert W (10.000000000000000, digits2<9>) = 1.3591409142295225
      Lambert W (10.000000000000000, digits2<10>) = 1.7441406250000000
      Lambert W (10.000000000000000, digits2<11>) = 1.7455280027406996
      Lambert W (10.000000000000000, digits2<12>) = 1.7455280027406996
      Lambert W (10.000000000000000, digits2<13>) = 1.7455280027406996
      Lambert W (10.000000000000000, digits2<14>) = 1.7455280027406996
      Lambert W (10.000000000000000, digits2<15>) = 1.7455280027406996
      Lambert W (10.000000000000000, digits2<16>) = 1.7455280027406996
      Lambert W (10.000000000000000, digits2<17>) = 1.7455280027406996
      Lambert W (10.000000000000000, digits2<18>) = 1.7455280027406996
      Lambert W (10.000000000000000, digits2<19>) = 1.7455280027406996
      Lambert W (10.000000000000000, digits2<20>) = 1.7455280027406996
      Lambert W (10.000000000000000, digits2<21>) = 1.7455280027406996
      Lambert W (10.000000000000000, digits2<22>) = 1.7455280027406996
      Lambert W (10.000000000000000, digits2<23>) = 1.7455280027406996
      Lambert W (10.000000000000000, digits2<24>) = 1.7455280027406996
      Lambert W (10.000000000000000, digits2<25>) = 1.7455280027406996
      Lambert W (10.000000000000000, digits2<26>) = 1.7455280027406996
      Lambert W (10.000000000000000, digits2<27>) = 1.7455280027406996
      Lambert W (10.000000000000000, digits2<28>) = 1.7455280027406996
      Lambert W (10.000000000000000, digits2<29>) = 1.7455280027406996
      Lambert W (10.000000000000000, digits2<30>) = 1.7455280027406996
      Lambert W (10.000000000000000, digits2<31>) = 1.7455280027406996
      Lambert W (10.000000000000000, digits2<32>) = 1.7455280027406996
      Lambert W (10.000000000000000, digits2<33>) = 1.7455280027406996
      Lambert W (10.000000000000000, digits2<34>) = 1.7455280027406996
      Lambert W (10.000000000000000, digits2<35>) = 1.7455280027406996
      Lambert W (10.000000000000000, digits2<36>) = 1.7455280027406996
      Lambert W (10.000000000000000, digits2<37>) = 1.7455280027406996
      Lambert W (10.000000000000000, digits2<38>) = 1.7455280027406996
      Lambert W (10.000000000000000, digits2<39>) = 1.7455280027406996
      Lambert W (10.000000000000000, digits2<40>) = 1.7455280027406996
      Lambert W (10.000000000000000, digits2<41>) = 1.7455280027406996
      Lambert W (10.000000000000000, digits2<42>) = 1.7455280027406996
      Lambert W (10.000000000000000, digits2<43>) = 1.7455280027406996
      Lambert W (10.000000000000000, digits2<44>) = 1.7455280027406996
      Lambert W (10.000000000000000, digits2<45>) = 1.7455280027406996
      Lambert W (10.000000000000000, digits2<46>) = 1.7455280027406996
      Lambert W (10.000000000000000, digits2<47>) = 1.7455280027406996
      Lambert W (10.000000000000000, digits2<48>) = 1.7455280027406996
      Lambert W (10.000000000000000, digits2<49>) = 1.7455280027406996
      Lambert W (10.000000000000000, digits2<50>) = 1.7455280027406996
      Lambert W (10.000000000000000, digits2<51>) = 1.7455280027406992
      Lambert W (10.000000000000000, digits2<52>) = 1.7455280027406992
      Lambert W (10.000000000000000, digits2<53>) = 1.7455280027406992
      Lambert W (10.000000000000000, digits2<54>) = 1.7455280027406992
      Lambert W (10.000000000000000) default      = 1.7455280027406992

      // cpp_dec_float_50 reference value:
      lambert_w0(10.000000000000000)              = 1.7455280027406993830743012648753899115352881290809

      */
    }
    { // Similar example using cpp_bin_float_quad (128-bit floating-point types).
      // multiprecision type to show use of Schroeder double approximation before Halley refinement.
      // Specifying precision with digits2 not digits10.

      cpp_bin_float_quad z = 10.;
      std::cout << "\nTest evaluation of cpp_bin_float_quad Lambert W(" << z << ") at varying precisions"
        "\nspecified using a policy with (base_2) digits increasing from 1 to 54. "
        << std::endl;
      std::cout << std::setprecision(3) << "std::numeric_limits<cpp_bin_float_quad>::digits = " << std::numeric_limits<cpp_bin_float_quad>::digits << std::endl;
      std::cout << std::setprecision(3) << "std::numeric_limits<cpp_bin_float_quad>::epsilon() = " << std::numeric_limits<cpp_bin_float_quad>::epsilon() << std::endl;
      std::cout << std::setprecision(3) << "std::numeric_limits<cpp_bin_float_quad>::max_digits10 = " << std::numeric_limits<cpp_bin_float_quad>::max_digits10 << std::endl;
      std::cout << std::setprecision(3) << "std::numeric_limits<cpp_bin_float_quad>::digits10 = " << std::numeric_limits<cpp_bin_float_quad>::digits10 << std::endl;
      std::cout.precision(std::numeric_limits<cpp_bin_float_quad>::max_digits10);
      std::cout << "Lambert W (" << z << ", digits2<1>)  = " << lambert_w0(z, policy<digits2<1> >()) << std::endl;
      std::cout << "Lambert W (" << z << ", digits2<2>)  = " << lambert_w0(z, policy<digits2<2> >()) << std::endl;
      std::cout << "..." << std::endl;
      std::cout << "Lambert W (" << z << ", digits2<52>) = " << lambert_w0(z, policy<digits2<52> >()) << std::endl;
      std::cout << "Lambert W (" << z << ", digits2<53>) = " << lambert_w0(z, policy<digits2<53> >()) << std::endl;
      std::cout << "Lambert W (" << z << ", digits2<54>) = " << lambert_w0(z, policy<digits2<54> >()) << std::endl;
      std::cout << "..." << std::endl;
      std::cout << "Lambert W (" << z << ", digits2<113>) = " << lambert_w0(z, policy<digits2<113> >()) << std::endl;
      std::cout << "Lambert W (" << z << ") default       = " << lambert_w0(z) << std::endl; // Default policy.
      // All are same precision because double precision first approximation used before Halley.

      /*

      Test evaluation of cpp_bin_float_quad Lambert W at varying precisions
      specified using a policy with (base_2) digits increasing from 1 to 113.
      std::numeric_limits<cpp_bin_float_quad>::epsilon() = 1.93e-34
      Lambert W (10.00000000000000000000000000000000000, digits2<1>) = 1.745528002740699383074301264875390030
      Lambert W (10.00000000000000000000000000000000000, digits2<2>) = 1.745528002740699383074301264875390030
      ...
      Lambert W (10.00000000000000000000000000000000000) default      = 1.745528002740699383074301264875390030
      lambert_w0(z) cpp_dec_float_50                                  = 1.7455280027406993830743012648753899115352881290809 // For comparison.

      */
    }

    { // Reference value for lambert_w(1)
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
      std::cout.precision(std::numeric_limits<cpp_bin_float_quad>::digits10 +1); // 

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


    // Examples of values near singularity when only a final Halley refinement is available.

    //using boost::math::detail::lambert_w_singularity_series;
    //std::cout << series(z) << std::endl;

    //std::cout << lambert_w_singularity_series(z) << std::endl;

    /*
    1>  Lambert W example single spot tests!
    1>  epsilon 2.2204460492503131e-16
    1>  digits2 53
    1>  digits10 15
    1>  digits <RealType, policies::policy<>() 53
    1>  Lambert W (0.29999999999999999) = 0.23675531078855933
    1>  epsilon 3.814697265625e-06
    1>  digits2 19
    1>  digits10 5
    1>  digits <RealType, policies::policy<>() 53
    1>  Lambert W (0.29999999999999999) = 0.2367553107885593


    for (double xd = 0.001; xd < 1e20; xd *= 10)
    {

    // 1.  0.56714329040978387
    //     0.56714329040978384

    // 10 1.7455280027406994
    //    1.7455280027406994

    // 100 3.3856301402900502
    //     3.3856301402900502
    // 1000 5.2496028524015959
    //      5.249602852401596227126056319697306282521472386059592844451465483991362228320942832739693150854347718

    // 1e19 40.058769161984308
    //      40.05876916198431163898797971203180915622644925765346546858291325452428038208071849105889199253335063
    std::cout << "Lambert W (" << xd << ") = " << lambert_w(xd) << std::endl; //


    1>  Iteration #0, w0 0.0099072820916067030, w1 = 0.0099014738435951096, difference = -5.8082480115934088e-06, relative 0.00058660438873459064, float distance = 3415046080725.0000
    1>  f'(x) = 5.8082144415180436e-06, f''(x) = -5.8082480115936214e-06
    1>  Iteration #1, w0 0.0099014738435951096, w1 = 0.0099014738435950107, difference = -9.8879238130678004e-17, relative 9.9920072216264089e-15, float distance = 58.000000000000000
    1>  f'(x) = 9.8645911007260505e-17, f''(x) = -9.8645911007260505e-17
    1>  Iteration #2, w0 0.0099014738435950107, w1 = 0.0099014738435950125, difference = 1.7347234759768071e-18, relative -2.2204460492503131e-16, float distance = 1.0000000000000000
    1>  f'(x) = -1.7007915690906991e-18, f''(x) = 1.7007915690906991e-18
    1>
    1>  Return refined 0.0099014738435950125 after 3 iterations, difference = -2.2204460492503131e-16, Float distance = 1.0000000000000000
    1>  Lambert W (0.010000000000000000) = 0.0099014738435950125


    1>  Iteration #0, w0 0.091722597168998055, w1 = 0.091276527200431806, difference = -0.00044606996856624836, relative 0.0048870173115453941, float distance = 38445375960087.000
    1>  f'(x) = 0.00044587943029702348, f''(x) = -0.00044606996856624685
    1>  Iteration #1, w0 0.091276527200431806, w1 = 0.091276527160862264, difference = -3.9569542087392051e-11, relative 4.3351278122827352e-10, float distance = 3408918.0000000000
    1>  f'(x) = 3.9569548190331561e-11, f''(x) = -3.9569548191831827e-11
    1>  Iteration #2, w0 0.091276527160862264, w1 = 0.091276527160862264, exact.
    1>  f'(x) = 0.00000000000000000, f''(x) = -0.00000000000000000
    1>


    1>  Iteration #0, w0 0.87095897919507137, w1 = 0.85260701286390628, difference = -0.018351966331165093, relative 0.021524531295515459, float distance = 182161642476426.00
    1>  f'(x) = 0.018097151008287501, f''(x) = -0.018351966331165141
    1>  Iteration #1, w0 0.85260701286390628, w1 = 0.85260550201372554, difference = -1.5108501807414854e-06, relative 1.7720389760000899e-06, float distance = 14784835581.000000
    1>  f'(x) = 1.5108484233744529e-06, f''(x) = -1.5108501807757508e-06
    1>  Iteration #2, w0 0.85260550201372554, w1 = 0.85260550201372554, exact.
    1>  f'(x) = 0.00000000000000000, f''(x) = -0.00000000000000000
    1>
    1>  Return refined 0.85260550201372554 after 3 iterations, exact.
    1>  Lambert W (2.0000000000000000) = 0.85260550201372554

    1>  Iteration #0, w0 0.87095897919507137, w1 = 0.85260701286390628, difference = -0.018351966331165093, relative 0.021524531295515459, float distance = 182161642476426.00
    1>  f'(x) = 0.018097151008287501, f''(x) = -0.018351966331165141
    1>  Iteration #1, w0 0.85260701286390628, w1 = 0.85260550201372554, difference = -1.5108501807414854e-06, relative 1.7720389760000899e-06, float distance = 14784835581.000000
    1>  f'(x) = 1.5108484233744529e-06, f''(x) = -1.5108501807757508e-06
    1>  Iteration #2, w0 0.85260550201372554, w1 = 0.85260550201372554, exact.
    1>  f'(x) = 0.00000000000000000, f''(x) = -0.00000000000000000
    1>
    1>  Return refined 0.85260550201372554 after 3 iterations, exact.
    1>  Lambert W (2.0000000000000000) = 0.85260550201372554

    }
    */

  }
  catch (std::exception& ex)
  {
    std::cout << ex.what() << std::endl;
  }
}  // int main()

   /*
   BOOST_CHECK_CLOSE_FRACTION(lambert_w(BOOST_MATH_TEST_VALUE(RealType, 0.01)),
   BOOST_MATH_TEST_VALUE(RealType, 0.0099014738435950118853363268165701079536277464949174),
   // Output from https://www.wolframalpha.com/input/ N[lambert_w[0.01],50])
   tolerance * 25);  // <<< Needs a much bigger tolerance???
   // 0.0099014738435951096 this test max_digits10
   // 0.00990147384359511  digits10
   // 0.0099014738435950118  wolfram
   // 0.00990147384359501  wolfram  digits10
   // 0.0099014738435950119 N[lambert_w[0.01],17]
   // 0.00990147384359501   N[lambert_w[0.01],15] which really is more different than expected.

   // 0.00990728209160670  approx
   // 0.00990147384359511  previous

   after allowing iterations, get the Wolfram result

   x = 0.0100000000000000, 1st approximation  = 0.00990728209160670
   1>  Iteration #0, w0 0.00990728209160670, w1 = 0.00990147384359511, difference = 0.000586604388734591, relative 0.000586604388734591
   1>  f'(x) = 5.80821444151804e-06, f''(x) = -5.80824801159362e-06
   1>  Iteration #1, w0 0.00990147384359511, w1 = 0.00990147384359501, difference = 9.99200722162641e-15, relative 9.99200722162641e-15
   1>  f'(x) = 9.86459110072605e-17, f''(x) = -9.86459110072605e-17
   1>  Iteration #2, w0 0.00990147384359501, w1 = 0.00990147384359501, difference = -2.22044604925031e-16, relative -2.22044604925031e-16
   1>  f'(x) = -1.70079156909070e-18, f''(x) = 1.70079156909070e-18
   1>  Estimate 0.00990147384359501 refined after 2 iterations.
   1>  Exact.
   1>  Lambert W (0.0100000000000000) = 0.00990147384359501


   using float_distance < 2

   x = 0.0100000000000000, 1st approximation  = 0.00990728209160670
   1>  Iteration #0, w0 0.00990728209160670, w1 = 0.00990147384359511, difference = 0.000586604388734591, relative 0.000586604388734591
   1>  f'(x) = 5.80821444151804e-06, f''(x) = -5.80824801159362e-06
   1>  Iteration #1, w0 0.00990147384359511, w1 = 0.00990147384359501, difference = 9.99200722162641e-15, relative 9.99200722162641e-15
   1>  f'(x) = 9.86459110072605e-17, f''(x) = -9.86459110072605e-17
   1>  Estimate 0.00990147384359501 refined after 1 iterations.
   1>  diff = -1.73472347597681e-18, relative 0.000000000000000
   1>  Lambert W (0.0100000000000000) = 0.00990147384359501
   */

   /*

   Output:

   1>  Lambert W example single spot tests!
   1>  epsilon 0.03125
   1>  Lambert W (0.300000012, digits10<1>) = 0.239290372

   1>  epsilon 0.00390625
   1>  Iteration #1, w0 0.239290372, w1 = 0.236755311, difference = -0.00253506005, relative 0.0107074976, float distance = 133610
   1>  f'(x) = 0.00252926024, f''(x) = -0.00253505306

   1>  Lambert W (0.300000012, digits10<2>) = 0.236755311
   1>  epsilon 0.000244140625
   1>  Iteration #1, w0 0.239290372, w1 = 0.236755311, difference = -0.00253506005, relative 0.0107074976, float distance = 133610
   1>  f'(x) = 0.00252926024, f''(x) = -0.00253505306

   1>  Lambert W (0.300000012, digits10<3>) = 0.236755311
   1>  epsilon 3.05175781e-05
   1>  Iteration #1, w0 0.239290372, w1 = 0.236755311, difference = -0.00253506005, relative 0.0107074976, float distance = 133610
   1>  f'(x) = 0.00252926024, f''(x) = -0.00253505306

   1>  Lambert W (0.300000012, digits10<4>) = 0.236755311
   1>  epsilon 3.81469727e-06
   1>  Iteration #1, w0 0.239290372, w1 = 0.236755311, difference = -0.00253506005, relative 0.0107074976, float distance = 133610
   1>  f'(x) = 0.00252926024, f''(x) = -0.00253505306
   1>  Iteration #2, w0 0.236755311, w1 = 0.236755326, difference = 1.49011612e-08, relative -5.96046448e-08, float distance = 1
   1>  f'(x) = -1.90171221e-08, f''(x) = 1.90171221e-08
   1>
   1>  Return refined 0.236755326 after 2 iterations, difference = -5.96046448e-08, Float distance = 1

   1>  Lambert W (0.300000012, digits10<5>) = 0.236755326
   1>  epsilon 2.38418579e-07
   1>  Iteration #1, w0 0.239290372, w1 = 0.236755311, difference = -0.00253506005, relative 0.0107074976, float distance = 133610
   1>  f'(x) = 0.00252926024, f''(x) = -0.00253505306
   1>  Iteration #2, w0 0.236755311, w1 = 0.236755326, difference = 1.49011612e-08, relative -5.96046448e-08, float distance = 1
   1>  f'(x) = -1.90171221e-08, f''(x) = 1.90171221e-08
   1>
   1>  Return refined 0.236755326 after 2 iterations, difference = -5.96046448e-08, Float distance = 1

   1>  Lambert W (0.300000012, digits10<6>) = 0.236755326
   1>  epsilon 1.1920929e-07

   1>  Lambert W (0.300000012, digits10<7>) = 0.239290372
   1>  epsilon 1.1920929e-07

   1>  Lambert W (0.300000012, digits10<8>) = 0.239290372
   1>  epsilon 1.1920929e-07

   1>  Lambert W (0.300000012, digits10<9>) = 0.239290372
   */


   /*
   1>  std::numeric_limits<float>::epsilon() = 1.19209290e-07
   1>  Lambert_w argument = 10.0000000
   1>  Argument Type = float
   1>  digits10 = 1, digits2 = 5
   1>  Result Integer W between g[0] = 0.000000000 < 10.0000000 < 2.71828175 = g[1], bisect 1.35914087
   1>
   1>  Lambert_w argument = 10.0000000
   1>  Argument Type = float
   1>  digits10 = 3, digits2 = 10
   1>  Result Bisection = 1.74414063
   1>
   1>  Lambert_w argument = 10.0000000
   1>  Argument Type = float
   1>  digits10 = 6, digits2 = 22
   1>  Result Schroeder refinement = 1.74552798
   1>
   1>  Lambert_w argument = 10.0000000
   1>  Argument Type = float
   1>  digits10 = 6, digits2 = 23
   1>  Result Halley refinement = 1.74552810
   1>
   1>  Lambert_w argument = 10.0000000
   1>  Argument Type = float
   1>  digits10 = 7, digits2 = 24
   1>  Result Halley refinement = 1.74552810

   1>
   1>  std::numeric_limits<float>::epsilon() = 1.19209290e-07
   1>  Lambert_w argument = 10.000000000000000
   1>  Argument Type = double
   1>  digits10 = 1, digits2 = 5
   1>  Result Integer W bisection between g[0] = 0.00000000000000000 < 10.000000000000000 < 2.7182818284590451 = g[1], bisect 1.3591409142295225
   1>
   1>  Lambert_w argument = 10.000000000000000
   1>  Argument Type = double
   1>  digits10 = 3, digits2 = 10
   1>  Result Bisection = 1.7441406250000000
   1>
   1>  Lambert_w argument = 10.000000000000000
   1>  Argument Type = double
   1>  digits10 = 6, digits2 = 22
   1>  Result Schroeder refinement = 1.7455280027406996
   1>
   1>  Lambert_w argument = 10.000000000000000
   1>  Argument Type = double
   1>  digits10 = 6, digits2 = 23
   1>  Result Halley refinement = 1.7455280027406994
   1>
   1>  Lambert_w argument = 10.000000000000000
   1>  Argument Type = double
   1>  digits10 = 15, digits2 = 53
   1>  Result Halley refinement = 1.7455280027406994

   Codeblocks long double

   Lambert W example single spot tests!

   std::numeric_limits<float>::epsilon() = 1.19209290e-007
   Lambert_w argument = 10.0000000000000000000
   Argument Type = e, max_digits10 = 21, epsilon = 1.08420217248550443401e-019
   digits10 = 1, digits2 = 5
   Result Integer W bisection between g[0] = 0.00000000000000000000 < 10.0000000000000000000 < 2.71828182845904523543 = g[1],
   bisect =                      1.35914091422952261771

   Lambert_w argument = 10.0000000000000000000
   Argument Type = e, max_digits10 = 21, epsilon = 1.08420217248550443401e-019
   digits10 = 3, digits2 = 10
   Result Bisection =            1.74414062500000000000

   Lambert_w argument = 10.0000000000000000000
   Argument Type = e, max_digits10 = 21, epsilon = 1.08420217248550443401e-019
   digits10 = 6, digits2 = 22
   Result Schroeder refinement = 1.74552800274069937864

   Lambert_w argument = 10.0000000000000000000
   Argument Type = e, max_digits10 = 21, epsilon = 1.08420217248550443401e-019
   digits10 = 6, digits2 = 23
   Result Halley refinement =    1.74552800274069938309

   Lambert_w argument = 10.0000000000000000000
   Argument Type = e, max_digits10 = 21, epsilon = 1.08420217248550443401e-019
   digits10 = 19, digits2 = 64
   Result Halley refinement =    1.74552800274069938309


   Process returned 0 (0x0)   execution time : 0.079 s
   Press any key to continue.


   After Halley reordering

   Lambert W example single spot tests!

   std::numeric_limits<float>::epsilon() = 1.19209290e-007
   Lambert_w argument = 10.0000000000000000000
   Argument Type = e, max_digits10 = 21, epsilon = 1.08420217248550443401e-019
   digits10 = 2, digits2 = 9
   Result Integer W bisection between g[0] = 0.00000000000000000000 < 10.0000000000000000000 < 2.71828182845904523543 = g[1],
   bisect mean =                 1.35914091422952261771

   Lambert_w argument = 10.0000000000000000000
   Argument Type = e, max_digits10 = 21, epsilon = 1.08420217248550443401e-019
   digits10 = 3, digits2 = 10
   Result Bisection =            1.74414062500000000000

   Lambert_w argument = 10.0000000000000000000
   Argument Type = e, max_digits10 = 21, epsilon = 1.08420217248550443401e-019
   digits10 = 6, digits2 = 21
   Result Schroeder refinement = 1.74552800274069937864

   Lambert_w argument = 10.0000000000000000000
   Argument Type = e, max_digits10 = 21, epsilon = 1.08420217248550443401e-019
   digits10 = 5, digits2 = 19
   Result Schroeder refinement = 1.74552800274069937864

   Lambert_w argument = 10.0000000000000000000
   Argument Type = e, max_digits10 = 21, epsilon = 1.08420217248550443401e-019
   digits10 = 6, digits2 = 23
   Result Schroeder refinement = 1.74552800274069937864

   Lambert_w argument = 10.0000000000000000000
   Argument Type = e, max_digits10 = 21, epsilon = 1.08420217248550443401e-019
   digits10 = 19, digits2 = 64
   w = 1.74552800274069937864, z = 10.0000000000000000000, exp(w) = 5.72892556538694148333, diff = -6.93889390390722837765e-017
   1 Halley iterations.
   Distance = 41
   Result Halley refinement =    1.74552800274069938309


   Process returned 0 (0x0)   execution time : 0.070 s
   Press any key to continue.


   14 Aug 2017 using merged #include <boost/math/special_functions/lambert_w.hpp>

   lambert_w_pb_precision.vcxproj -> J:\Cpp\Misc\x64\Release\lambert_w_pb_precision.exe

   lambert_w_pb_precision.vcxproj -> J:\Cpp\Misc\x64\Release\lambert_w_pb_precision.exe
   lambert_w_pb_precision.vcxproj -> J:\Cpp\Misc\x64\Release\lambert_w_pb_precision.pdb (Full PDB)
   Lambert W examples of precision control.

   std::numeric_limits<double>::epsilon() = 2.22e-16
   lambert_w0(x, my_pol_3_dec()) = 0.567

   lambert_w0(x, my_prec_10_bits()) = 0.566

   lambert_w0(z, policy<digits2<9> >())  = 0.000
   lambert_w0(z, policy<digits2<10> >()) = 0.566
   lambert_w0(z, policy<digits2<21> >()) = 0.567
   lambert_w0(z, policy<digits10<5> >()) = 0.567
   lambert_w0(z, policy<digits2<23> >()) = 0.567

   Test evaluation of float Lambert W at varying precisions specified using a policy with digits10 increasing from 1 to 9.
   std::numeric_limits<float>::epsilon() = 1.19e-07
   Lambert W (10.0000000, digits10<1>) = 1.35914087
   Lambert W (10.0000000, digits10<2>) = 1.35914087
   Lambert W (10.0000000, digits10<3>) = 1.74552798
   Lambert W (10.0000000, digits10<4>) = 1.74552798
   Lambert W (10.0000000, digits10<5>) = 1.74552798
   Lambert W (10.0000000, digits10<6>) = 1.74552810
   Lambert W (10.0000000, digits10<7>) = 1.74552810
   Lambert W (10.0000000, digits10<8>) = 1.74552810
   Lambert W (10.0000000, digits10<9>) = 1.74552810
   Lambert W (10.0000000) = 1.74552810

   Test evaluation of double Lambert W at varying precisions specified using a policy with digits10 increasing from 1 to 24.
   std::numeric_limits<double>::epsilon() = 2.22e-16
   Lambert W (10.000000000000000, digits2<1>) = 1.3591409142295225
   Lambert W (10.000000000000000, digits2<2>) = 1.3591409142295225
   Lambert W (10.000000000000000, digits2<3>) = 1.3591409142295225
   Lambert W (10.000000000000000, digits2<4>) = 1.3591409142295225
   Lambert W (10.000000000000000, digits2<5>) = 1.3591409142295225
   Lambert W (10.000000000000000, digits2<6>) = 1.3591409142295225
   Lambert W (10.000000000000000, digits2<7>) = 1.3591409142295225
   Lambert W (10.000000000000000, digits2<8>) = 1.3591409142295225
   Lambert W (10.000000000000000, digits2<9>) = 1.3591409142295225
   Lambert W (10.000000000000000, digits2<10>) = 1.7441406250000000
   Lambert W (10.000000000000000, digits2<11>) = 1.7455280027406996
   Lambert W (10.000000000000000, digits2<12>) = 1.7455280027406996
   Lambert W (10.000000000000000, digits2<13>) = 1.7455280027406996
   Lambert W (10.000000000000000, digits2<14>) = 1.7455280027406996
   Lambert W (10.000000000000000, digits2<15>) = 1.7455280027406996
   Lambert W (10.000000000000000, digits2<16>) = 1.7455280027406996
   Lambert W (10.000000000000000, digits2<17>) = 1.7455280027406996
   Lambert W (10.000000000000000, digits2<18>) = 1.7455280027406996
   Lambert W (10.000000000000000, digits2<19>) = 1.7455280027406996
   Lambert W (10.000000000000000, digits2<20>) = 1.7455280027406996
   Lambert W (10.000000000000000, digits2<21>) = 1.7455280027406996
   Lambert W (10.000000000000000, digits2<22>) = 1.7455280027406996
   Lambert W (10.000000000000000, digits2<23>) = 1.7455280027406996
   Lambert W (10.000000000000000, digits2<24>) = 1.7455280027406996
   Lambert W (10.000000000000000) = 1.7455280027406992
   lambert_w0(z) = 1.7455280027406993830743012648753899115352881290809


   30 Aug 2017

   lambert_w_pb_precision.cpp
   lambert_w_pb_precision.vcxproj -> J:\Cpp\Misc\x64\Debug\lambert_w_pb_precision.exe
   lambert_w_pb_precision.vcxproj -> J:\Cpp\Misc\x64\Debug\lambert_w_pb_precision.pdb (Partial PDB)
   Lambert W examples of precision control.

   std::numeric_limits<double>::epsilon() = 2.22e-16
   lambert_w0(x, my_pol_3_dec()) = 0.56714329040978395

   lambert_w0(x, my_prec_11_bits()) = 0.56714329040978395


   std::numeric_limits<float>::epsilon() = 1.19e-07
   std::numeric_limits<float>::digits = 24
   std::numeric_limits<float>::digits10 = 6
   std::numeric_limits<float>::max_digits10 = 9
   lambert_w0(z, policy<digits2<9> >())  = 0.000000000
   lambert_w0(z, policy<digits2<10> >()) = 0.566406250
   lambert_w0(z, policy<digits2<21> >()) = 0.567143261
   lambert_w0(z, policy<digits10<5> >()) = 0.567143261
   lambert_w0(z, policy<digits2<23> >()) = 0.567143261
   lambert_w0(z, policy<>()) = 0.567143261

   Test evaluation of float Lambert W at varying precisions
   specified using a policy with digits10 increasing from 1 to 9.
   std::numeric_limits<float>::epsilon() = 1.19e-07
   Lambert W (1.00000000, digits10<1>) = 0.000000000
   Lambert W (1.00000000, digits10<2>) = 0.000000000
   Lambert W (1.00000000, digits10<3>) = 0.567143261
   Lambert W (1.00000000, digits10<4>) = 0.567143261
   Lambert W (1.00000000, digits10<5>) = 0.567143261
   Lambert W (1.00000000, digits10<6>) = 0.567143261
   Lambert W (1.00000000, digits10<7>) = 0.567143261
   Lambert W (1.00000000, digits10<8>) = 0.567143261
   Lambert W (1.00000000, digits10<9>) = 0.567143261
   Lambert W (1.00000000) default      = 0.567143261

   Test evaluation of float Lambert W at varying precisions
   specified using a policy with (base_2) digits increasing from 1 to 24.
   std::numeric_limits<float>::epsilon() = 1.19e-07
   Lambert W (10.0000000, digits2<1> ) = 1.35914087
   Lambert W (10.0000000, digits2<2> ) = 1.35914087
   Lambert W (10.0000000, digits2<3> ) = 1.35914087
   Lambert W (10.0000000, digits2<4> ) = 1.35914087
   Lambert W (10.0000000, digits2<5> ) = 1.35914087
   Lambert W (10.0000000, digits2<6> ) = 1.35914087
   Lambert W (10.0000000, digits2<7> ) = 1.35914087
   Lambert W (10.0000000, digits2<8> ) = 1.35914087
   Lambert W (10.0000000, digits2<9> ) = 1.35914087
   Lambert W (10.0000000, digits2<10>) = 1.74414063
   Lambert W (10.0000000, digits2<11>) = 1.74552798
   Lambert W (10.0000000, digits2<12>) = 1.74552798
   Lambert W (10.0000000, digits2<13>) = 1.74552798
   Lambert W (10.0000000, digits2<14>) = 1.74552798
   Lambert W (10.0000000, digits2<15>) = 1.74552798
   Lambert W (10.0000000, digits2<16>) = 1.74552798
   Lambert W (10.0000000, digits2<17>) = 1.74552798
   Lambert W (10.0000000, digits2<18>) = 1.74552798
   Lambert W (10.0000000, digits2<19>) = 1.74552798
   Lambert W (10.0000000, digits2<20>) = 1.74552798
   Lambert W (10.0000000, digits2<21>) = 1.74552798
   Lambert W (10.0000000, digits2<22>) = 1.74552810
   Lambert W (10.0000000, digits2<23>) = 1.74552810
   Lambert W (10.0000000, digits2<24>) = 1.74552810
   Lambert W (10.0000000) default      = 1.74552810

   Test evaluation of double Lambert W at varying precisions
   specified using a policy with (base_2) digits increasing from 1 to 54.
   std::numeric_limits<double>::epsilon() = 2.22e-16
   Lambert W (10.000000000000000, digits2<1> ) = 1.3591409142295225
   Lambert W (10.000000000000000, digits2<2> ) = 1.3591409142295225
   Lambert W (10.000000000000000, digits2<3> ) = 1.3591409142295225
   Lambert W (10.000000000000000, digits2<4> ) = 1.3591409142295225
   Lambert W (10.000000000000000, digits2<5> ) = 1.3591409142295225
   Lambert W (10.000000000000000, digits2<6> ) = 1.3591409142295225
   Lambert W (10.000000000000000, digits2<7> ) = 1.3591409142295225
   Lambert W (10.000000000000000, digits2<8> ) = 1.3591409142295225
   Lambert W (10.000000000000000, digits2<9> ) = 1.3591409142295225
   Lambert W (10.000000000000000, digits2<10>) = 1.7441406250000000
   Lambert W (10.000000000000000, digits2<11>) = 1.7455280027406996
   Lambert W (10.000000000000000, digits2<12>) = 1.7455280027406996
   Lambert W (10.000000000000000, digits2<13>) = 1.7455280027406996
   Lambert W (10.000000000000000, digits2<14>) = 1.7455280027406996
   Lambert W (10.000000000000000, digits2<15>) = 1.7455280027406996
   Lambert W (10.000000000000000, digits2<16>) = 1.7455280027406996
   Lambert W (10.000000000000000, digits2<17>) = 1.7455280027406996
   Lambert W (10.000000000000000, digits2<18>) = 1.7455280027406996
   Lambert W (10.000000000000000, digits2<19>) = 1.7455280027406996
   Lambert W (10.000000000000000, digits2<20>) = 1.7455280027406996
   Lambert W (10.000000000000000, digits2<21>) = 1.7455280027406996
   Lambert W (10.000000000000000, digits2<22>) = 1.7455280027406996
   Lambert W (10.000000000000000, digits2<23>) = 1.7455280027406996
   Lambert W (10.000000000000000, digits2<24>) = 1.7455280027406996
   Lambert W (10.000000000000000, digits2<25>) = 1.7455280027406996
   Lambert W (10.000000000000000, digits2<26>) = 1.7455280027406996
   Lambert W (10.000000000000000, digits2<27>) = 1.7455280027406996
   Lambert W (10.000000000000000, digits2<28>) = 1.7455280027406996
   Lambert W (10.000000000000000, digits2<29>) = 1.7455280027406996
   Lambert W (10.000000000000000, digits2<30>) = 1.7455280027406996
   Lambert W (10.000000000000000, digits2<31>) = 1.7455280027406996
   Lambert W (10.000000000000000, digits2<32>) = 1.7455280027406996
   Lambert W (10.000000000000000, digits2<33>) = 1.7455280027406996
   Lambert W (10.000000000000000, digits2<34>) = 1.7455280027406996
   Lambert W (10.000000000000000, digits2<35>) = 1.7455280027406996
   Lambert W (10.000000000000000, digits2<36>) = 1.7455280027406996
   Lambert W (10.000000000000000, digits2<37>) = 1.7455280027406996
   Lambert W (10.000000000000000, digits2<38>) = 1.7455280027406996
   Lambert W (10.000000000000000, digits2<39>) = 1.7455280027406996
   Lambert W (10.000000000000000, digits2<40>) = 1.7455280027406996
   Lambert W (10.000000000000000, digits2<41>) = 1.7455280027406996
   Lambert W (10.000000000000000, digits2<42>) = 1.7455280027406996
   Lambert W (10.000000000000000, digits2<43>) = 1.7455280027406996
   Lambert W (10.000000000000000, digits2<44>) = 1.7455280027406996
   Lambert W (10.000000000000000, digits2<45>) = 1.7455280027406996
   Lambert W (10.000000000000000, digits2<46>) = 1.7455280027406996
   Lambert W (10.000000000000000, digits2<47>) = 1.7455280027406996
   Lambert W (10.000000000000000, digits2<48>) = 1.7455280027406996
   Lambert W (10.000000000000000, digits2<49>) = 1.7455280027406996
   Lambert W (10.000000000000000, digits2<50>) = 1.7455280027406996
   Lambert W (10.000000000000000, digits2<51>) = 1.7455280027406994
   Lambert W (10.000000000000000, digits2<52>) = 1.7455280027406994
   Lambert W (10.000000000000000, digits2<53>) = 1.7455280027406992
   Lambert W (10.000000000000000, digits2<54>) = 1.7455280027406992
   Lambert W (10.000000000000000) default      = 1.7455280027406992

   Test evaluation of cpp_bin_float_quad Lambert W at varying precisions
   specified using a policy with (base_2) digits increasing from 1 to 54.
   std::numeric_limits<cpp_bin_float_quad>::epsilon() = 1.93e-34
   Lambert W (10.00000000000000000000000000000000000, digits2<1> ) = 1.745528002740699383074301264875389837
   Lambert W (10.00000000000000000000000000000000000, digits2<2> ) = 1.745528002740699383074301264875389837
   ...
   Lambert W (10.00000000000000000000000000000000000, digits2<52>) = 1.745528002740699383074301264875389837
   Lambert W (10.00000000000000000000000000000000000, digits2<53>) = 1.745528002740699383074301264875389837
   Lambert W (10.00000000000000000000000000000000000, digits2<54>) = 1.745528002740699383074301264875390030
   Lambert W (10.00000000000000000000000000000000000) default      = 1.745528002740699383074301264875390030
   lambert_w0(z) cpp_dec_float_50       = 1.7455280027406993830743012648753899115352881290809




   */


