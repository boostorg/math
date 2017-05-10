// Copyright Paul A. Bristow 2016, 2017.
// Copyright John Maddock 2016.

// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

// test_lambertw.cpp
//! \brief Basic sanity tests for Lambert W function plog or lambert_w using algorithm from Thomas Luu.

#include <boost/math/concepts/real_concept.hpp> // for real_concept
#define BOOST_TEST_MAIN
#define BOOST_LIB_DIAGNOSTIC

#include <boost/test/unit_test.hpp> // Boost.Test
#include <boost/test/floating_point_comparison.hpp>
#include "test_value.hpp"  // for create_test_value

#include <boost/multiprecision/cpp_dec_float.hpp> // boost::multiprecision::cpp_dec_float_50
using boost::multiprecision::cpp_dec_float_50;

#include <boost/multiprecision/cpp_bin_float.hpp> 

#ifdef BOOST_HAS_FLOAT128
#include <boost/multiprecision/float128.hpp> // Not available for MSVC.
#endif

#include <boost/fixed_point/fixed_point.hpp>

//#define BOOST_MATH_INSTRUMENT_LAMBERT_W  // #define only for Lambert_w diagnostic output.
#include <boost/math/special_functions/lambert_w.hpp> // For Lambert W lambert_w function.
using boost::math::lambert_w;
#include <boost/math/special_functions/fpclassify.hpp> // isnan, ifinite.

#include <boost/array.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/type_traits/is_constructible.hpp>

#include <limits>
#include <cmath>
#include <typeinfo>
#include <iostream>
#include <type_traits>
#include <exception>

static const unsigned int noof_tests = 2;

static const boost::array<const char*, noof_tests> test_data =
{{
  "1",
  "6"
}}; // array test_data

//static const boost::array<const char*, noof_tests> test_expected =
//{ {
//    BOOST_MATH_TEST_VALUE(T, 0.56714329040978387299996866221035554975381578718651), // Output from https://www.wolframalpha.com/input/?i=lambert_w(1)
//    BOOST_MATH_TEST_VALUE(T, 1.432404775898300311234078007212058694786434608804302025655) // Output from https://www.wolframalpha.com/input/?i=lambert_w(6)

//  } }; // array test_data

//! Build a message of information about build, architecture, address model, platform, ...
std::string show_versions()
{
  // Some of this information can also be obtained from running with a Custom Post-build step
  // adding the option --build_info=yes
    // "$(TargetDir)$(TargetName).exe" --build_info=yes

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

template <class RealType>
void test_spots(RealType)
{
  // (Unused Parameter value, arbitrarily zero, only communicates the floating point type).
  // test_spots(0.F); test_spots(0.); test_spots(0.L);

  using boost::math::lambert_w;
  using boost::math::constants::exp_minus_one;
  using boost::math::policies::policy;

  typedef policy <
    boost::math::policies::domain_error<boost::math::policies::ignore_error>,
    boost::math::policies::overflow_error<boost::math::policies::ignore_error>,
    boost::math::policies::underflow_error<boost::math::policies::ignore_error>,
    boost::math::policies::denorm_error<boost::math::policies::ignore_error>,
    boost::math::policies::pole_error<boost::math::policies::ignore_error>,
    boost::math::policies::evaluation_error<boost::math::policies::ignore_error>
  > ignore_all_policy;

//   // Test some bad parameters to the function, with default policy and with ignore_all policy.
#ifndef BOOST_NO_EXCEPTIONS
  BOOST_CHECK_THROW(boost::math::lambert_w<RealType>(-1.), std::domain_error);
  BOOST_CHECK_THROW(lambert_w<RealType>(std::numeric_limits<RealType>::quiet_NaN()), std::domain_error); // Would be NaN.
  BOOST_CHECK_THROW(lambert_w<RealType>(std::numeric_limits<RealType>::infinity()), std::domain_error); // Would be infinite.
  BOOST_CHECK_THROW(lambert_w<RealType>(-static_cast<RealType>(0.4)), std::domain_error); // Would be complex.

#else
  BOOST_MATH_CHECK_EQUAL(boost::math::lambert_w<RealType>(std::numeric_limits<RealType>::quiet_NaN(), ignore_all_policy), std::numeric_limits<RealType::quiet_NaN()); // nan
  BOOST_MATH_CHECK_EQUAL(boost::math::lambert_w<RealType>(std::numeric_limits<RealType>::infinity(), ignore_all_policy), std::numeric_limits<RealType::quiet_NaN()); // infinity
#endif

  std::cout << "\nTesting type " << typeid(RealType).name() << std::endl;
  int epsilons = 2;
  RealType tolerance = boost::math::tools::epsilon<RealType>() * epsilons; // 2 eps as a fraction.
  std::cout << "Tolerance " << epsilons << " * epsilon == " << tolerance << std::endl;
  std::cout << "Precision " << std::numeric_limits<RealType>::digits10 << " decimal digits." << std::endl;
  // std::cout.precision(std::numeric_limits<RealType>::digits10);
  std::cout.precision(std::numeric_limits <RealType>::max_digits10);
  std::cout.setf(std::ios_base::showpoint);  // show trailing significant zeros.
  std::cout << "-exp(-1) = " << -exp_minus_one<RealType>() << std::endl;

  // Test at singularity.  Three tests because some failed previously - bug now gone?
  //RealType test_value = BOOST_MATH_TEST_VALUE(RealType, -0.36787944117144232159552377016146086744581113103176783450783680169746149574489980335714727434591964374662732527);
  RealType singular_value = -exp_minus_one<RealType>();
  // -exp(-1) = -0.36787944117144232159552377016146086744581113103176783450783680169746149574489980335714727434591964374662732527
  // lambert_w[-0.367879441171442321595523770161460867445811131031767834] == -1
  //           -0.36787945032119751 
  RealType minus_one_value = BOOST_MATH_TEST_VALUE(RealType, -1.);
  //std::cout << "singular_value " << singular_value << ", expected Lambert W = " << minus_one_value << std::endl;

  BOOST_CHECK_CLOSE_FRACTION( // Check -exp(-1) = -0.367879450 = -1
    lambert_w(singular_value),
    minus_one_value,
    tolerance);  // OK

  BOOST_CHECK_CLOSE_FRACTION(  // Check -exp(-1) ~= -0.367879450 == -1
    lambert_w(BOOST_MATH_TEST_VALUE(RealType, -0.36787944117144232159552377016146086744581113103176783450783680169746149574489980335714727434591964374662732527)),
    BOOST_MATH_TEST_VALUE(RealType, -1.),
    tolerance);

  BOOST_CHECK_CLOSE_FRACTION(  // Check -exp(-1) ~= -0.367879450 == -1
    lambert_w<RealType>(-exp_minus_one<RealType>()),
    BOOST_MATH_TEST_VALUE(RealType, -1.),
    tolerance);

  // Tests with some spot values computed using
  // https://www.wolframalpha.com/input
  // For example: N[lambert_w[1], 50] outputs:
  // 0.56714329040978387299996866221035554975381578718651

  BOOST_CHECK_CLOSE_FRACTION(lambert_w(BOOST_MATH_TEST_VALUE(RealType, 0.1)),
    BOOST_MATH_TEST_VALUE(RealType, 0.091276527160862264299895721423179568653119224051472),
    // Output from https://www.wolframalpha.com/input/?i=lambert_w(0.2)
    tolerance);

  BOOST_CHECK_CLOSE_FRACTION(lambert_w(BOOST_MATH_TEST_VALUE(RealType, 0.2)),
    BOOST_MATH_TEST_VALUE(RealType, 0.16891597349910956511647490370581839872844691351073),
    // Output from https://www.wolframalpha.com/input/?i=lambert_w(0.2)
    tolerance);

  BOOST_CHECK_CLOSE_FRACTION(lambert_w(BOOST_MATH_TEST_VALUE(RealType, 0.5)),
    BOOST_MATH_TEST_VALUE(RealType, 0.351733711249195826024909300929951065171464215517111804046),
    // Output from https://www.wolframalpha.com/input/?i=lambert_w(0.5)
    tolerance);

  BOOST_CHECK_CLOSE_FRACTION(
    lambert_w(BOOST_MATH_TEST_VALUE(RealType, 1.)),
    BOOST_MATH_TEST_VALUE(RealType, 0.56714329040978387299996866221035554975381578718651),
   // Output from https://www.wolframalpha.com/input/?i=lambert_w(1)
   tolerance);

  BOOST_CHECK_CLOSE_FRACTION(lambert_w(BOOST_MATH_TEST_VALUE(RealType, 2.)),
    BOOST_MATH_TEST_VALUE(RealType, 0.852605502013725491346472414695317466898453300151403508772),
    // Output from https://www.wolframalpha.com/input/?i=lambert_w(2.)
    tolerance);

  BOOST_CHECK_CLOSE_FRACTION(lambert_w(BOOST_MATH_TEST_VALUE(RealType, 3.)),
    BOOST_MATH_TEST_VALUE(RealType, 1.049908894964039959988697070552897904589466943706341452932),
    // Output from https://www.wolframalpha.com/input/?i=lambert_w(3.)
    tolerance);

  BOOST_CHECK_CLOSE_FRACTION(lambert_w(BOOST_MATH_TEST_VALUE(RealType, 5.)),
    BOOST_MATH_TEST_VALUE(RealType, 1.326724665242200223635099297758079660128793554638047479789),
    // Output from https://www.wolframalpha.com/input/?i=lambert_w(0.5)
    tolerance);

  BOOST_CHECK_CLOSE_FRACTION(lambert_w(BOOST_MATH_TEST_VALUE(RealType, 6.)),
    BOOST_MATH_TEST_VALUE(RealType, 1.432404775898300311234078007212058694786434608804302025655),
    // Output from https://www.wolframalpha.com/input/?i=lambert_w(6)
    tolerance);

  BOOST_CHECK_CLOSE_FRACTION(lambert_w(BOOST_MATH_TEST_VALUE(RealType, 100.)),
    BOOST_MATH_TEST_VALUE(RealType, 3.3856301402900501848882443645297268674916941701578),
    // Output from https://www.wolframalpha.com/input/?i=lambert_w(100)
    tolerance);

  // This fails for fixed_point type used for other tests because out of range?
    //BOOST_CHECK_CLOSE_FRACTION(lambert_w(BOOST_MATH_TEST_VALUE(RealType, 1.0e6)),
    //BOOST_MATH_TEST_VALUE(RealType, 11.383358086140052622000156781585004289033774706019),
    //// Output from https://www.wolframalpha.com/input/?i=lambert_w(1e6)
    //// tolerance * 1000); // fails for fixed_point type exceeds 0.00015258789063
    //  // 15.258789063
    //  // 11.383346558
    // tolerance * 100000);

  // So need to use some spot tests for specific types, or use a bigger fixed_point type.

  // Check zero.
  BOOST_CHECK_CLOSE_FRACTION(lambert_w(BOOST_MATH_TEST_VALUE(RealType, 0.0)),
    BOOST_MATH_TEST_VALUE(RealType, 0.0),
    tolerance);
  // these fail for cpp_dec_float_50
  // 'boost::multiprecision::detail::expression<boost::multiprecision::detail::negate,boost::multiprecision::number<boost::multiprecision::backends::cpp_dec_float<50,int32_t,void>,boost::multiprecision::et_on>,void,void,void>'
  // : no appropriate default constructor available
  // TODO ???????????
/* 
  if (std::numeric_limits<RealType>::is_specialized)
  { // Check +/- values near to zero.
    BOOST_CHECK_CLOSE_FRACTION(lambert_w(std::numeric_limits<RealType>::epsilon()),
      BOOST_MATH_TEST_VALUE(RealType, 0.0),
      tolerance);

      BOOST_CHECK_CLOSE_FRACTION(lambert_w(-std::numeric_limits<RealType>::epsilon()),
        BOOST_MATH_TEST_VALUE(RealType, 0.0),
        tolerance);
  } // is_specialized
  
  // Check infinite if possible.
  if (std::numeric_limits<RealType>::has_infinity)
  {
    BOOST_CHECK_CLOSE_FRACTION(lambert_w(+std::numeric_limits<RealType>::infinity()),
      BOOST_MATH_TEST_VALUE(RealType, 0.0),
      tolerance);

    BOOST_CHECK_CLOSE_FRACTION(lambert_w(-std::numeric_limits<RealType>::infinity()),
      BOOST_MATH_TEST_VALUE(RealType, 0.0),
      tolerance);
  }

  // Check NaN if possible.
  if (std::numeric_limits<RealType>::has_quiet_NaN)
  {
    BOOST_CHECK_CLOSE_FRACTION(lambert_w(+std::numeric_limits<RealType>::quiet_NaN()),
    BOOST_MATH_TEST_VALUE(RealType, 0.0),
    tolerance);

  BOOST_CHECK_CLOSE_FRACTION(lambert_w(-std::numeric_limits<RealType>::quiet_NaN()),
    BOOST_MATH_TEST_VALUE(RealType, 0.0),
    tolerance);
  }

   */



   // Checks on input that should throw.

 /* This should throw when implemented.
  BOOST_CHECK_CLOSE_FRACTION(lambert_w(-0.5),
  BOOST_MATH_TEST_VALUE(RealType, ),
    // Output from https://www.wolframalpha.com/input/?i=lambert_w(-0.5)
    tolerance);
   */
} //template <class RealType>void test_spots(RealType)

BOOST_AUTO_TEST_CASE(test_types)
{
  BOOST_MATH_CONTROL_FP;
  // BOOST_TEST_MESSAGE( output only appears if command line has --log_level="message"
  // or call set_threshold_level function:
  boost::unit_test_framework::unit_test_log.set_threshold_level(boost::unit_test_framework::log_messages);
  BOOST_TEST_MESSAGE("Test Lambert W function for several types.");
  BOOST_TEST_MESSAGE(show_versions());  // Full version of Boost, STL and compiler info.

  // Fundamental built-in types:
  test_spots(0.0F); // float
  test_spots(0.0); // double
  //if (sizeof(long double) > sizeof (double))
  //{ // Avoid pointless re-testing if double and long double are identical (for example, MSVC).
  //  test_spots(0.0L); // long double
  //}
  // Built-in/fundamental Quad 128-bit type, if available.
  #ifdef BOOST_HAS_FLOAT128
  test_spots(static_cast<boost::multiprecision::float128>(0));
#endif

  // Multiprecision types:
  test_spots(static_cast<boost::multiprecision::cpp_dec_float_50>(0));
  test_spots(static_cast<boost::multiprecision::cpp_bin_float_double_extended>(0));
  test_spots(static_cast<boost::multiprecision::cpp_bin_float_quad>(0));

  // Fixed-point types:

  // Some fail 0.1 to 1.0 ??? 
  //test_spots(static_cast<boost::fixed_point::negatable<15,-16> >(0));
  // fixed_point::negatable<15,-16> has range 1.52587891e-05 to 32768, epsilon 3.05175781e-05

  //test_spots(boost::math::concepts::real_concept(0.1));  // "real_concept" - was OK.

} // BOOST_AUTO_TEST_CASE( test_types )

BOOST_AUTO_TEST_CASE(test_range_of_values)
{
  using boost::math::lambert_w;
  using boost::math::constants::exp_minus_one;

  BOOST_TEST_MESSAGE("Test Lambert W function type double for range of values.");
  // Want to test almost largest value.
  // test_value = (std::numeric_limits<RealType>::max)() / 4;
  // std::cout << std::setprecision(std::numeric_limits<RealType>::max_digits10) << "Max value = " << test_value << std::endl;
  // Can't use a test like this for all types because max_value depends on RealType and thus the expected result of lambert_w does too.
  //BOOST_CHECK_CLOSE_FRACTION(lambert_w<RealType>(test_value),
  //  BOOST_MATH_TEST_VALUE(RealType, ???),
  //  tolerance);
  // So this section just tests a single type, say IEEE 64-bit double, for a range of spot values.

  typedef double RealType;
  
  int epsilons = 2;
  RealType tolerance = boost::math::tools::epsilon<RealType>() * epsilons; // 2 eps as a fraction.
  std::cout << "Tolerance " << epsilons  << " * epsilon == " << tolerance << std::endl;

  BOOST_CHECK_CLOSE_FRACTION(lambert_w(BOOST_MATH_TEST_VALUE(RealType, 1.0e-6)),
    BOOST_MATH_TEST_VALUE(RealType, 9.9999900000149999733333854165586669000967020964243e-7),
    // Output from https://www.wolframalpha.com/input/ N[lambert_w[1e-6],50])
    tolerance);  
  BOOST_CHECK_CLOSE_FRACTION(lambert_w(BOOST_MATH_TEST_VALUE(RealType, 0.0001)),
    BOOST_MATH_TEST_VALUE(RealType, 0.000099990001499733385405869000452213835767629477903460),
    // Output from https://www.wolframalpha.com/input/ N[lambert_w[0.001],50])
    tolerance);
  BOOST_CHECK_CLOSE_FRACTION(lambert_w(BOOST_MATH_TEST_VALUE(RealType, 0.001)),
    BOOST_MATH_TEST_VALUE(RealType, 0.00099900149733853088995782787410778559957065467928884),
    // Output from https://www.wolframalpha.com/input/ N[lambert_w[0.001],50])
    tolerance);
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

  BOOST_CHECK_CLOSE_FRACTION(lambert_w(BOOST_MATH_TEST_VALUE(RealType, 0.05)),
    BOOST_MATH_TEST_VALUE(RealType, 0.047672308600129374726388900514160870747062965933891),
    // Output from https://www.wolframalpha.com/input/ N[lambert_w[0.01],50])
    tolerance);

  BOOST_CHECK_CLOSE_FRACTION(lambert_w(BOOST_MATH_TEST_VALUE(RealType, 0.1)),
    BOOST_MATH_TEST_VALUE(RealType, 0.091276527160862264299895721423179568653119224051472),
    // Output from https://www.wolframalpha.com/input/ N[lambert_w[1],50])
    tolerance);

  BOOST_CHECK_CLOSE_FRACTION(lambert_w(BOOST_MATH_TEST_VALUE(RealType, 1.)),
    BOOST_MATH_TEST_VALUE(RealType, 0.56714329040978387299996866221035554975381578718651),
    // Output from https://www.wolframalpha.com/input/ N[lambert_w[1],50])
    tolerance);

  BOOST_CHECK_CLOSE_FRACTION(lambert_w(BOOST_MATH_TEST_VALUE(RealType, 2.)),
    BOOST_MATH_TEST_VALUE(RealType, 0.852605502013725491346472414695317466898453300151403508772),
    // Output from https://www.wolframalpha.com/input/?i=lambert_w(2.)
    tolerance);

  BOOST_CHECK_CLOSE_FRACTION(lambert_w(BOOST_MATH_TEST_VALUE(RealType, 3.)),
    BOOST_MATH_TEST_VALUE(RealType, 1.049908894964039959988697070552897904589466943706341452932),
    // Output from https://www.wolframalpha.com/input/?i=lambert_w(3.)
    tolerance);

  BOOST_CHECK_CLOSE_FRACTION(lambert_w(BOOST_MATH_TEST_VALUE(RealType, 5.)),
    BOOST_MATH_TEST_VALUE(RealType, 1.326724665242200223635099297758079660128793554638047479789),
    // Output from https://www.wolframalpha.com/input/?i=lambert_w(0.5)
    tolerance);

  BOOST_CHECK_CLOSE_FRACTION(lambert_w(BOOST_MATH_TEST_VALUE(RealType, 6.)),
    BOOST_MATH_TEST_VALUE(RealType, 1.432404775898300311234078007212058694786434608804302025655),
    // Output from https://www.wolframalpha.com/input/?i=lambert_w(6)
    tolerance);

  BOOST_CHECK_CLOSE_FRACTION(lambert_w(BOOST_MATH_TEST_VALUE(RealType, 10.)),
    BOOST_MATH_TEST_VALUE(RealType, 1.7455280027406993830743012648753899115352881290809),
    // Output from https://www.wolframalpha.com/input/ N[lambert_w[10],50])
    tolerance);

    BOOST_CHECK_CLOSE_FRACTION(lambert_w(BOOST_MATH_TEST_VALUE(RealType, 100.)),
    BOOST_MATH_TEST_VALUE(RealType, 3.3856301402900501848882443645297268674916941701578),
    // Output from https://www.wolframalpha.com/input/?i=lambert_w(100)
    tolerance);

  BOOST_CHECK_CLOSE_FRACTION(lambert_w(BOOST_MATH_TEST_VALUE(RealType, 1000.)),
    BOOST_MATH_TEST_VALUE(RealType, 5.2496028524015962271260563196973062825214723860596),
    // Output from https://www.wolframalpha.com/input/?i=lambert_w(1000)
    tolerance);

  // This fails for fixed_point type used for other tests because out of range of the type?
  BOOST_CHECK_CLOSE_FRACTION(lambert_w(BOOST_MATH_TEST_VALUE(RealType, 1.0e6)),
    BOOST_MATH_TEST_VALUE(RealType, 11.383358086140052622000156781585004289033774706019),
    // Output from https://www.wolframalpha.com/input/?i=lambert_w(1e6)
    tolerance); //
 
  // Tests for double only near the max and the singularity where Lambert_w estimates are less precise.
  // 
  if (std::numeric_limits<RealType>::is_specialized)
  {
    // Check near std::numeric_limits<>::max() for type.

    //std::cout << std::setprecision(std::numeric_limits<RealType>::max_digits10) 
    //  << (std::numeric_limits<double>::max)() // == 1.7976931348623157e+308
    //  << " " << (std::numeric_limits<double>::max)() / 4  // == 4.4942328371557893e+307
    //  << std::endl;

    BOOST_CHECK_CLOSE_FRACTION(boost::math::lambert_w(4.4942328371557893e+307), // max_value/4 for IEEE 64-bit double.
      static_cast<double>(702.02379914670587),
      // N[lambert_w[4.4942328371557893e+307], 35]  == 701.84270921429200143342782556643059
      // as a   double == 701.83341468208209
      // Lambert computed 702.02379914670587
     // std::numeric_limits<double>::epsilon() * 256);
      0.0003); // Much less precise near the max edge.

    BOOST_CHECK_CLOSE_FRACTION(boost::math::lambert_w(4.4942328371557893e+307/2), // near max_value for IEEE 64-bit double.
      static_cast<double>(701.15054872492476914094824907722937),
      // N[lambert_w[4.4942328371557893e+307], 35]  == 701.84270921429200143342782556643059
      // as a   double == 701.83341468208209
      // Lambert computed 702.02379914670587
      std::numeric_limits<double>::epsilon());

    BOOST_CHECK_CLOSE_FRACTION(boost::math::lambert_w(4.4942328371557893e+307/4), // near max_value for IEEE 64-bit double.
      static_cast<double>(700.45838920868939857588606393559517),
      // N[lambert_w[4.4942328371557893e+307], 35]  == 701.84270921429200143342782556643059
      // as a   double == 701.83341468208209
      // Lambert computed 702.02379914670587
      std::numeric_limits<double>::epsilon());

    // This test value is one epsilon close to the singularity at -exp(1) * x
    // (below which the result has a non-zero imaginary part).
    RealType test_value = -exp_minus_one<RealType>();
    test_value += (std::numeric_limits<RealType>::epsilon() * 1);
    BOOST_CHECK_CLOSE_FRACTION(lambert_w(test_value),
      BOOST_MATH_TEST_VALUE(RealType, -0.99999996349975895),
      tolerance * 1000000000);
    // -0.99999996788201051
    // -0.99999996349975895
    // Would not expect to get a result closer than sqrt(epsilon)?
  } //  if (std::numeric_limits<RealType>::is_specialized)

} // BOOST_AUTO_TEST_CASE( test_main )

  /*

  Output:


  */




