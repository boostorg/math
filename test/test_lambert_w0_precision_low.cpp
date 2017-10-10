// Copyright Paul A. Bristow 2017.

// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or
//  copy at http ://www.boost.org/LICENSE_1_0.txt).

// Test evaluations of Lambert W function and comparison with reference values computed at 100 decimal digits precision.

//#define BOOST_MATH_LAMBERT_W_TEST_OUTPUT // Define to output all values of argument z and lambert W and differences.

#include <boost/cstdfloat.hpp> // For float_64_t, float128_t. Must be first include!
#include <boost/config.hpp> // for BOOST_PLATFORM, BOOST_COMPILER,  BOOST_STDLIB ...Fw =
#include <boost/version.hpp>   // for BOOST_MSVC versions.
#include <boost/cstdint.hpp>
#include <boost/exception/exception.hpp>  // boost::exception
#include <boost/math/constants/constants.hpp> // For exp_minus_one == 3.67879441171442321595523770161460867e-01.
#include <boost/math/special_functions/next.hpp>
//#include "test_value.hpp"  // for create_test_value and macro BOOST_MATH_TEST_VALUE.

// Built-in/fundamental Quad 128-bit type, if available.
#ifdef BOOST_HAS_FLOAT128
#include <boost/multiprecision/float128.hpp> // Not available for MSVC.
//using boost::multiprecision::float128;
#endif

#include <boost/multiprecision/cpp_dec_float.hpp> // boost::multiprecision::cpp_dec_float_50
using boost::multiprecision::cpp_dec_float_50; // 50 decimal digits type.
using boost::multiprecision::cpp_dec_float_100; // 100 decimal digits type.

#include <boost/multiprecision/cpp_bin_float.hpp>
using boost::multiprecision::cpp_bin_float_50; // 50 decimal digits type.

// Multiprecision types:
//test_spots(static_cast<boost::multiprecision::cpp_dec_float_50>(0));
//test_spots(static_cast<boost::multiprecision::cpp_bin_float_double_extended>(0));
//test_spots(static_cast<boost::multiprecision::cpp_bin_float_quad>(0));

// #define BOOST_MATH_INSTRUMENT_LAMBERT_W  // #define only for (much) diagnostic output.

// For lambert_w function.
#include <boost/math/special_functions/lambert_w.hpp>
using boost::multiprecision::cpp_bin_float_double_extended;
using boost::multiprecision::cpp_bin_float_double;
using boost::multiprecision::cpp_bin_float_quad;

#include <boost/type_traits/is_constructible.hpp>
#include <boost/array.hpp>
#include <boost/lexical_cast.hpp>
using boost::lexical_cast;
#include <boost/type_traits/is_constructible.hpp>

#include <iostream>
#include <exception>
#include <stdexcept>
#include <string>
#include <limits>  // For std::numeric_limits.
#include <cmath>
#include <typeinfo>

// Include 'Known good' Lambert w values using N[productlog(-0.3), 100]
// evaluated to full precision of RealType (up to 100 decimal digits).
// Checked for a few z values using WolframAlpha command: N[productlog(-0.3), 100]
#include "J:\Cpp\Misc\lambert_w_1000_test_values\lambert_w_mp_values.ipp"
// #include "lambert_w_mp_values.ipp"

// Optional functions to list z arguments and reference lambert w values.
template<typename RealType>
void list_zws()
{
  std::cout.precision(std::numeric_limits<RealType>::max_digits10);
  std::cout << noof_tests << " Test values of type: " << typeid(RealType).name() << std::endl;

  init_ws<RealType>();
  init_zs<RealType>();

  for (size_t i = 0; i != noof_tests; i++)
  {
    std::cout << i << " " << zs<RealType>[i] << " " << ws<RealType>[i] << std::endl;
  }
} // list_zws()template<typename RealType>

//! Test difference between 'known good' reference values (100 decimal digit precision)
//! of Lambert W and values computed for the RealType by lambert_w0 (or lambert_wm1) function.
template<typename RealType>
float test_zws()
{
  std::cout << "Tests with type: " << typeid(RealType).name() << std::endl;
  // std::cout << noof_tests << " test values computed and compared." << std::endl;
  // Better to show just once at start of all tests?
  std::cout << "digits10 = " << std::numeric_limits<RealType>::digits10
    << ", max_digits10 = " << std::numeric_limits<RealType>::max_digits10
    << std:: setprecision(3) << ", epsilon = " << std::numeric_limits<RealType>::epsilon()
    << std::endl;

  std::cout.precision(std::numeric_limits<RealType>::max_digits10); // Show all possibly significant digits.
  init_zs<RealType>(); // Test z arguments.
  init_ws<RealType>(); // Lambert W reference values.

  using boost::math::lambert_w0;
  using boost::math::lambert_wm1;
  using boost::math::float_distance;

  //! Biggest float distance between reference and computed Lambert W.
  long long int max_distance = 0;
  //! Total of all float distances for all the values tested.
  long long int total_distance = 0;
  // Mean distance = total_distance / noof_tests.
  float mean_distance = 0.F;

  for (size_t i = 0; i != noof_tests; i++)
  {
    RealType z = zs<RealType>[i];
    RealType w = lambert_w0<RealType>(z); // Only for float or double!
    RealType kw = ws<RealType>[i];  // 'Known good' 100 decimal digits.

    RealType fd = float_distance<RealType>(w, kw);
    RealType diff = (w - kw);

    long long int distance = static_cast<long long int>(fd);
    long long int abs_distance = abs(distance);
    total_distance += abs_distance;
    if (abs_distance > max_distance)
    {
      max_distance = abs_distance;
    }
    mean_distance = (float)total_distance / noof_tests;

  #ifdef BOOST_MATH_LAMBERT_W_TEST_OUTPUT
    std::cout << z << " " << w << " " << ws<RealType>[i] // << ", Distance = " << distance  not useful for multiprecision types?
      << std::setprecision(3)
      << ", diff abs " << diff
      << ", diff rel " << (w - kw) / kw
      << ", epsilons = " << diff / std::numeric_limits<RealType>::epsilon()
      << std::setprecision(std::numeric_limits<RealType>::max_digits10)
      << std::endl;
  #endif
  } // for
  std::streamsize precision = std::cout.precision(3);
  std::cout << "max " << max_distance << ", total distance = " << total_distance << " bits, mean " << mean_distance << std::endl;
  std::cout.precision(precision);
  return mean_distance;
} // void test_zws()

//! Compare a single evaluation of lambert_w0 with reference value (100 decimal digit precision).
//! \returns distance between values in bits.
//! This falls foul of conversion using multiprecision.
//! And of loss of precision using multiprecision.
//! using macro BOOST_MATH_TEST_VALUE avoids this pitfall.
template<typename RealType>
int test_spot(RealType z)
{
  using boost::math::lambert_w0;
  using boost::math::float_distance;

  RealType lw = lambert_w0<RealType>(z);
  RealType rlw = static_cast<RealType>(lambert_w0<cpp_dec_float_50>(z));

  std::cout.precision(std::numeric_limits<RealType>::max_digits10); // or
   //std::cout.precision(std::numeric_limits<RealType>::digits10);
 std::cout << lw << "\n" << rlw << std::endl;

  int distance = static_cast<int>(float_distance(rlw, lw));
  std::cout << "distance " << distance << std::endl;
  return distance;
} // int test_spot(RealType z)

//! Compare computed Lambert W(z) with 'known good' reference values.
//! Both as the same RealType.
template<typename RealType>
int compare_spot(RealType z, RealType w)
{
  using boost::math::lambert_w0;
  using boost::math::lambert_wm1;
  using boost::math::nextafter;
  using boost::math::float_next;
  using boost::math::float_distance;

  init_zs<RealType>(); // Test z arguments.
  init_ws<RealType>(); // Lambert W reference values.

  std::cout << std::boolalpha << "Test with type: " << typeid(RealType).name() << std::endl;

  std::cout  << "digits10 = " << std::numeric_limits<RealType>::digits10
    << ", max_digits10 = " << std::numeric_limits<RealType>::max_digits10
    << ", epsilon = " << std::numeric_limits<RealType>::epsilon()
    << std::endl;

  // std::cout.precision(std::numeric_limits<RealType>::max_digits10);
  std::cout.precision(std::numeric_limits<RealType>::digits10);

  RealType lw = lambert_w0<RealType>(z);
  RealType flw = w;

  std::cout << lw << "\n" << flw << std::endl;

  int distance = static_cast<int>(float_distance(flw, lw));

  std::cout << "distance " << distance << std::endl;

  return distance;
} // template<typename RealType> int compare_spot(RealType z, RealType w)

//! THis is no longer useful?
//template<typename RealType>
//int test_spots(RealType)
//{
//  std::cout << std::boolalpha << "Test with type: " << typeid(RealType).name() << std::endl;
//  //std::cout.precision(std::numeric_limits<float>::max_digits10); // or digits10?
//  std::cout.precision(std::numeric_limits<float>::digits10);
//  std::size_t n = noof_tests;
//  std::cout << n << " test values. " << std::endl;
//  int distances = 0;
//  for (std::size_t i = 0; i != n; i++)
//  {
//    // RealType z = BOOST_MATH_TEST_VALUE(RealType, 0.5);
//    //RealType w = BOOST_MATH_TEST_VALUE(RealType, "0.351733711249195826024909300929951065171464215517111804046");
//    //RealType z = lexical_cast<RealType>(zs[i]);  // lexical_cast doesn't do more than double precision!
//    //RealType w = lexical_cast<RealType>(ws[i]);
//    // So can only do this safely using macro BOOST_MATH_TEST_VALUE.
//    //RealType z = BOOST_MATH_TEST_VALUE(RealType, test_zs[i]);
//    //RealType w = BOOST_MATH_TEST_VALUE(RealType, test_ws[i]);
//    // error Missing suffix L
//
//    RealType z = zs<RealType>[i];
//    RealType w = ws<RealType>[i];
//    std::cout << z << " " << w << std::endl;
//    distances += compare_spot(z, w);
//  }
//  return distances;
//} //  template<typename RealType>int test_spots(RealType)

int main()
{
  try
  {
    std::cout << "Lambert W tests,\n"
      << noof_tests <<
      " single spot tests comparing with 100 decimal digits precision reference values." << std::endl;
    std::cout.precision(std::numeric_limits<double>::max_digits10);
    std::cout << std::showpoint << std::endl; // Trailing zeros.

    using boost::math::constants::exp_minus_one;

    using boost::math::lambert_w0;
    using boost::math::lambert_wm1;
    using boost::math::nextafter;
    using boost::math::float_next;

    using boost::math::policies::make_policy;
    using boost::math::policies::digits10;
    using boost::math::policies::policy;

    using boost::math::policies::evaluation_error;
    using boost::math::policies::domain_error;
    using boost::math::policies::overflow_error;
    using boost::math::policies::domain_error;
    using boost::math::policies::throw_on_error;

    // Define a policy:
    typedef policy<
      domain_error<throw_on_error>,
      overflow_error<throw_on_error>
    > throw_policy;

    // Examples of some single tests show pitfalls using mutliprecision.

    // Single spot value test of lambert_w0(0.5).
    //using boost::multiprecision::cpp_bin_float_quad;
    //std::cout.precision(std::numeric_limits<cpp_bin_float_quad>::max_digits10);
    //cpp_bin_float_quad zq(0.5);
    //std::cout << "\nLambert W (" << zq << ") = " << LambertW0<cpp_bin_float_quad>(zq) << std::endl; //
    // 0.351733711249195826024909300929951065171464215517111804046  expected from Wolfram (and Luu triple Halley).
    // 0.3517337112491958262237012910577188474 // cpp_bin_float_quad
    // 0.35173371124919584 // double
    
    test_spot<float>(0.5F); // 0.351733714
    test_spot<double>(0.5); //  0.35173371124919584
    // Multiprecision won't convert at full precision.
    //test_spot<cpp_bin_float_quad>(0.5);
    //test_spot<cpp_dec_float_50>(0.5);
    //test_spot<cpp_bin_float_double_extended>(0.5);

    //std::cout <<"float total distances = " << test_spots<float>(0.1f) << std::endl;
    //std::cout <<"double total distances = " << test_spots<double>(0.1) << std::endl;
    //std::cout <<"cpp_bin_float_quad total distances = " << test_spots<cpp_bin_float_quad>(0.1) << std::endl;

    // Compare double spot value:
    compare_spot(0.5, 0.351733711249195826024909300929951065171464215517111804046); // 0.351733711249195826024909300929951

    //cpp_bin_float_quad hq(0.5); // Care - only convert 17 decimal digits precision for double.
    // but this value 0.5 is exactly representable, so exact, and so OK.
    cpp_bin_float_quad hq("0.5"); // 
    // For multiprecision types, need to use a decimal digit string, not a floating-point literal:
    cpp_bin_float_quad rhq("0.351733711249195826024909300929951065171464215517111804046");
    compare_spot(hq, rhq); // 0.351733711249195826024909300929951 and distance == 0
                           // Useful to check a particular z argument.
    init_zs<double>();
    init_ws<double>();

    compare_spot(zs<double>[1], ws<double>[1]); // -0.653694501269090 -0.653694501269089, distance = -1

    //  Simple group of test values, including near zero, and near singularity.

    // Optionally list the z arguments and reference test lambert_w values.
    //list_zws<float>();
    //list_zws<double>();
    //list_zws<long double>();
    //list_zws<cpp_bin_float_double>();
    //list_zws<cpp_bin_float_double_extended>();
    //list_zws<__float128>(); // OK, but fails test below.
//#ifdef BOOST_HAS_FLOAT128
//    list_zws<float128>();
//#endif
    // list_zws<cpp_bin_float_quad>();
    //list_zws<cpp_dec_float_50>();

    // Test several floating-point types:
    test_zws<float>();
    test_zws<double>();

    //    test_zws<cpp_bin_float_double>(); expect exactly as double.
    if (std::numeric_limits<long double>::digits == std::numeric_limits<double>::digits)
    { 
      // Emulate 80-bit long double.
      test_zws<cpp_bin_float_double_extended>();
    }
    else
    { // True 80-bit (or 128-bit) long double.
      test_zws<long double>();
    }
    test_zws<cpp_bin_float_quad>();
   //  test_zws<cpp_dec_float_50>();
  }
  catch (std::exception& ex)
  {
    std::cout << ex.what() << std::endl;
  }
}  // int main()

/*
test_lambert_w0_precision_low.cpp
Generating code
235 of 3718 functions ( 6.3%) were compiled, the rest were copied from previous compilation.
0 functions were new in current compilation
618 functions had inline decision re-evaluated but remain unchanged
Finished generating code
Lambert_w_pb2_tests.vcxproj -> J:\Cpp\Misc\x64\Release\Lambert_w_pb2_tests.exe
Lambert_w_pb2_tests.vcxproj -> J:\Cpp\Misc\x64\Release\Lambert_w_pb2_tests.pdb (Full PDB)
Lambert W tests,
100 single spot tests comparing with 100 decimal digits precision reference values.

0.351733714
0.351733714
distance 0
0.35173371124919584
0.35173371124919584
distance 0
Test with type: double
digits10 = 15, max_digits10 = 17, epsilon = 2.2204460492503131e-16
0.351733711249196
0.351733711249196
distance 0
Test with type: class boost::multiprecision::number<class boost::multiprecision::backends::cpp_bin_float<113,2,void,short,-16382,16383>,0>
digits10 = 33, max_digits10 = 37, epsilon = 1.92592994438724e-34
0.351733711249195826024909300929951
0.351733711249195826024909300929951
distance 0
Test with type: double
digits10 = 15, max_digits10 = 17, epsilon = 2.22044604925031308084726333618164e-16
-0.653694501269090
-0.653694501269089
distance -1

Tests with type: float
digits10 = 6, max_digits10 = 9, epsilon = 1.19e-07
max 2, total distance = 37 bits, mean 0.370
Tests with type: double
digits10 = 15, max_digits10 = 17, epsilon = 2.22e-16
max 9, total distance = 49 bits, mean 0.490
Tests with type: class boost::multiprecision::number<class boost::multiprecision::backends::cpp_bin_float<64,2,void,short,-16382,16383>,0>
digits10 = 18, max_digits10 = 22, epsilon = 1.08e-19
max 4, total distance = 47 bits, mean 0.470
Tests with type: class boost::multiprecision::number<class boost::multiprecision::backends::cpp_bin_float<113,2,void,short,-16382,16383>,0>
digits10 = 33, max_digits10 = 37, epsilon = 1.93e-34
max 4, total distance = 54 bits, mean 0.540

*/