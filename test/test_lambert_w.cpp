// Copyright Paul A. Bristow 2016, 2017.
// Copyright John Maddock 2016.

// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

// test_lambert_w.cpp
//! \brief Basic sanity tests for Lambert W function plog productlog
//! or lambert_w using algorithm informed by Thomas Luu, Veberic and Tosio Fukushima.

#include <boost/config.hpp>   // for BOOST_MSVC definition.
#include <boost/version.hpp>   // for BOOST_MSVC versions.
#include <boost/math/concepts/real_concept.hpp> // for real_concept

// Boost macros
#define BOOST_TEST_MAIN
#define BOOST_LIB_DIAGNOSTIC "on" // Report library file details.
#include <boost/test/unit_test.hpp> // Boost.Test
#include <boost/test/floating_point_comparison.hpp>

#include <boost/array.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/type_traits/is_constructible.hpp>

#include <boost/multiprecision/cpp_dec_float.hpp> // boost::multiprecision::cpp_dec_float_50
using boost::multiprecision::cpp_dec_float_50;

#include <boost/multiprecision/cpp_bin_float.hpp>
//#include <boost/fixed_point/fixed_point.hpp> // If available.
#include <boost/math/special_functions/fpclassify.hpp> // isnan, ifinite.
#include <boost/math/special_functions/next.hpp> // float_next, float_prior
using boost::math::float_next;
using boost::math::float_prior;

#include <boost/math/tools/test_value.hpp>  // for create_test_value and macro BOOST_MATH_TEST_VALUE.
#include <boost/math/policies/policy.hpp>
using boost::math::policies::digits2;
#include <boost/math/special_functions/lambert_w.hpp> // For Lambert W lambert_w function.

// #include "C:\Users\Paul\Desktop\lambert_w0 - Copy before amalgamation.hpp"  // JM latest providing jm_lambert_w0 and lambert_wm1
//#include "C:\Users\Paul\Desktop\lambert_w0.hpp"  // JM latest providing jm_lambert_w0 and lambert_wm1
using boost::math::lambert_wm1;
using boost::math::lambert_w0; // Use jm version instead.

#include <limits>
#include <cmath>
#include <typeinfo>
#include <iostream>
#include <type_traits>
#include <exception>

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

#ifdef __MINGW64__
std::cout << "MINGW64 " << __MINGW32_MAJOR_VERSION << __MINGW32_MINOR_VERSION << std::endl;
//
//  << __MINGW64_MAJOR_VERSION << __MINGW64_MINOR_VERSION << std::endl; not declared in this scope???
#endif // __MINGW64__

#ifdef __MINGW32__
std::cout << "MINGW64 " << __MINGW32_MAJOR_VERSION << __MINGW32_MINOR_VERSION << std::endl;
#endif // __MINGW32__

  message << "\n  STL " << BOOST_STDLIB;
  message << "\n  Boost version " << BOOST_VERSION / 100000 << "." << BOOST_VERSION / 100 % 1000 << "." << BOOST_VERSION % 100;

#ifdef BOOST_HAS_FLOAT
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

  using boost::math::lambert_w0;
  using boost::math::lambert_wm1;
  using boost::math::constants::exp_minus_one;
  using boost::math::constants::e;
  using boost::math::policies::policy;

  /*  Example of an exception-free 'ignore_all' policy (possibly ill-advised?).
  typedef policy <
    boost::math::policies::domain_error<boost::math::policies::ignore_error>,
    boost::math::policies::overflow_error<boost::math::policies::ignore_error>,
    boost::math::policies::underflow_error<boost::math::policies::ignore_error>,
    boost::math::policies::denorm_error<boost::math::policies::ignore_error>,
    boost::math::policies::pole_error<boost::math::policies::ignore_error>,
    boost::math::policies::evaluation_error<boost::math::policies::ignore_error>
  > ignore_all_policy;
  */

//  Test some bad parameters to the function, with default policy and with ignore_all policy.
#ifndef BOOST_NO_EXCEPTIONS
  BOOST_CHECK_THROW(lambert_w0<RealType>(-1.), std::domain_error);
  BOOST_CHECK_THROW(lambert_wm1<RealType>(-1.), std::domain_error);
  BOOST_CHECK_THROW(lambert_w0<RealType>(std::numeric_limits<RealType>::quiet_NaN()), std::domain_error); // Would be NaN.
  // BOOST_CHECK_THROW(lambert_w0<RealType>(std::numeric_limits<RealType>::infinity()), std::domain_error); // Was if infinity should throw, now infinity.
  BOOST_CHECK_THROW(lambert_w0<RealType>(-static_cast<RealType>(0.4)), std::domain_error); // Would be complex.
#else // No exceptions so set policy to ignore and check result is NaN
  BOOST_MATH_CHECK_EQUAL(boost::math::lambert_w0<RealType>(std::numeric_limits<RealType>::quiet_NaN(), ignore_all_policy), std::numeric_limits<RealType::quiet_NaN()); // NaN.
  BOOST_MATH_CHECK_EQUAL(boost::math::lambert_w0<RealType>(std::numeric_limits<RealType>::infinity(), ignore_all_policy), std::numeric_limits<RealType::infinity()); // infinity.
  BOOST_MATH_CHECK_EQUAL(boost::math::lambert_w0<RealType>(std::numeric_limits<RealType>::infinity(), ignore_all_policy), std::numeric_limits<RealType::infinity()); // infinity.
#endif

  std::cout << "\nTesting type " << typeid(RealType).name() << std::endl;
  int epsilons = 1;
  if (std::numeric_limits<RealType>::digits > 53)
  { // Multiprecision types.
    epsilons *= 8; // (Perhaps needed because need slightly longer (55) reference values?).
  }
  RealType tolerance = boost::math::tools::epsilon<RealType>() * epsilons; // 2 eps as a fraction.
  std::cout << "Tolerance " << epsilons << " * epsilon == " << tolerance << std::endl;
  std::cout << "Precision " << std::numeric_limits<RealType>::digits10 << " decimal digits." << std::endl;
  // std::cout.precision(std::numeric_limits<RealType>::digits10);
  std::cout.precision(std::numeric_limits <RealType>::max_digits10);
  std::cout.setf(std::ios_base::showpoint);  // show trailing significant zeros.
  std::cout << "-exp(-1) = " << -exp_minus_one<RealType>() << std::endl;

  // Test at singularity.
  //RealType test_value = BOOST_MATH_TEST_VALUE(RealType, -0.36787944117144232159552377016146086744581113103176783450783680169746149574489980335714727434591964374662732527);
  RealType singular_value = -exp_minus_one<RealType>();
  // -exp(-1) = -0.36787944117144232159552377016146086744581113103176783450783680169746149574489980335714727434591964374662732527
  // lambert_w0[-0.367879441171442321595523770161460867445811131031767834] == -1
  //           -0.36787945032119751
  RealType minus_one_value = BOOST_MATH_TEST_VALUE(RealType, -1.);
  //std::cout << "singular_value " << singular_value << ", expected Lambert W = " << minus_one_value << std::endl;

  BOOST_CHECK_CLOSE_FRACTION( // Check -exp(-1) = -0.367879450 = -1
    lambert_w0(singular_value),
    minus_one_value,
    tolerance);  // OK

  BOOST_CHECK_CLOSE_FRACTION(  // Check -exp(-1) ~= -0.367879450 == -1
    lambert_w0(BOOST_MATH_TEST_VALUE(RealType, -0.36787944117144232159552377016146086744581113103176783450783680169746149574489980335714727434591964374662732527)),
    BOOST_MATH_TEST_VALUE(RealType, -1.),
    tolerance);

  BOOST_CHECK_CLOSE_FRACTION(  // Check -exp(-1) ~= -0.367879450 == -1
    lambert_w0<RealType>(-exp_minus_one<RealType>()),
    BOOST_MATH_TEST_VALUE(RealType, -1.),
    tolerance);

  // Tests with some spot values computed using
  // https://www.wolframalpha.com/input
  // For example: N[lambert_w[1], 50] outputs:
  // 0.56714329040978387299996866221035554975381578718651

  // At branch junction singularity.
  BOOST_CHECK_CLOSE_FRACTION(  // Check -exp(-1) ~= -0.367879450 == -1
    lambert_w0(BOOST_MATH_TEST_VALUE(RealType, -0.36787944117144232159552377016146086744581113103176783450783680169746149574489980335714727434591964374662732527)),
    BOOST_MATH_TEST_VALUE(RealType, -1.),
    tolerance);

    BOOST_CHECK_CLOSE_FRACTION(lambert_w0(BOOST_MATH_TEST_VALUE(RealType, 0.1)),
    BOOST_MATH_TEST_VALUE(RealType, 0.091276527160862264299895721423179568653119224051472),
    // Output from https://www.wolframalpha.com/input/?i=lambert_w0(0.2)
    tolerance);

  BOOST_CHECK_CLOSE_FRACTION(lambert_w0(BOOST_MATH_TEST_VALUE(RealType, 0.2)),
    BOOST_MATH_TEST_VALUE(RealType, 0.16891597349910956511647490370581839872844691351073),
    // Output from https://www.wolframalpha.com/input/?i=lambert_w0(0.2)
    tolerance);

  BOOST_CHECK_CLOSE_FRACTION(lambert_w0(BOOST_MATH_TEST_VALUE(RealType, 0.5)),
    BOOST_MATH_TEST_VALUE(RealType, 0.351733711249195826024909300929951065171464215517111804046),
    // Output from https://www.wolframalpha.com/input/?i=lambert_w0(0.5)
    tolerance);

  BOOST_CHECK_CLOSE_FRACTION(
    lambert_w0(BOOST_MATH_TEST_VALUE(RealType, 1.)),
    BOOST_MATH_TEST_VALUE(RealType, 0.56714329040978387299996866221035554975381578718651),
   // Output from https://www.wolframalpha.com/input/?i=lambert_w0(1)
   tolerance);

  BOOST_CHECK_CLOSE_FRACTION(lambert_w0(BOOST_MATH_TEST_VALUE(RealType, 2.)),
    BOOST_MATH_TEST_VALUE(RealType, 0.852605502013725491346472414695317466898453300151403508772),
    // Output from https://www.wolframalpha.com/input/?i=lambert_w0(2.)
    tolerance);

  BOOST_CHECK_CLOSE_FRACTION(lambert_w0(BOOST_MATH_TEST_VALUE(RealType, 3.)),
    BOOST_MATH_TEST_VALUE(RealType, 1.049908894964039959988697070552897904589466943706341452932),
    // Output from https://www.wolframalpha.com/input/?i=lambert_w0(3.)
    tolerance);

  BOOST_CHECK_CLOSE_FRACTION(lambert_w0(BOOST_MATH_TEST_VALUE(RealType, 5.)),
    BOOST_MATH_TEST_VALUE(RealType, 1.326724665242200223635099297758079660128793554638047479789),
    // Output from https://www.wolframalpha.com/input/?i=lambert_w0(0.5)
    tolerance);

  BOOST_CHECK_CLOSE_FRACTION(lambert_w0(BOOST_MATH_TEST_VALUE(RealType, 6.)),
    BOOST_MATH_TEST_VALUE(RealType, 1.432404775898300311234078007212058694786434608804302025655),
    // Output from https://www.wolframalpha.com/input/?i=lambert_w0(6)
    tolerance);

  BOOST_CHECK_CLOSE_FRACTION(lambert_w0(BOOST_MATH_TEST_VALUE(RealType, 100.)),
    BOOST_MATH_TEST_VALUE(RealType, 3.3856301402900501848882443645297268674916941701578),
    // Output from https://www.wolframalpha.com/input/?i=lambert_w0(100)
    tolerance);

  if (std::numeric_limits<RealType>::has_infinity)
  {
    // BOOST_CHECK_THROW(lambert_w0(std::numeric_limits<RealType>::infinity()), std::domain_error); // If should throw exception.
    //BOOST_CHECK_EQUAL(lambert_w0(std::numeric_limits<RealType>::infinity()), +std::numeric_limits<RealType>::infinity()); // message is:
    // Error in "test_types": class boost::exception_detail::clone_impl<struct boost::exception_detail::error_info_injector<class std::domain_error> > :
    // Error in function boost::math::lambert_w0<RealType>(<RealType>) : Argument z is infinite!
    //BOOST_CHECK_EQUAL(lambert_w0(std::numeric_limits<RealType>::infinity()), +std::numeric_limits<RealType>::infinity()); // If infinity allowed.
    BOOST_CHECK_THROW(lambert_wm1(std::numeric_limits<RealType>::infinity()), std::domain_error); // Infinity NOT allowed.
  }
  if (std::numeric_limits<RealType>::has_quiet_NaN)
  { // Argument Z == NaN is always an throwable error for both branches.
    // BOOST_CHECK_EQUAL(lambert_w0(std::numeric_limits<RealType>::quiet_NaN()), +std::numeric_limits<RealType>::infinity()); // message is:
    // Error in function boost::math::lambert_w0<RealType>(<RealType>): Argument z is NaN!
    BOOST_CHECK_THROW(lambert_w0(std::numeric_limits<RealType>::quiet_NaN()), std::domain_error);
    BOOST_CHECK_THROW(lambert_wm1(std::numeric_limits<RealType>::quiet_NaN()), std::domain_error);
  }

  // denorm - but might be == min or zero?
  if (std::numeric_limits<RealType>::has_denorm == true)
  { // Might also return infinity like z == 0?
    BOOST_CHECK_THROW(lambert_wm1(std::numeric_limits<RealType>::denorm_min()), std::overflow_error);
  }

    // Tests of Lambert W-1 branch.
    BOOST_CHECK_CLOSE_FRACTION(  // Check -exp(-1) ~= -0.367879450 == -1 at the singularity branch point.
    lambert_wm1(BOOST_MATH_TEST_VALUE(RealType, -0.36787944117144232159552377016146086744581113103176783450783680169746149574489980335714727434591964374662732527)),
    BOOST_MATH_TEST_VALUE(RealType, -1.),
    tolerance);

    // Near singularity and using series approximation.
    // N[productlog(-1, -0.36), 50] = -1.2227701339785059531429380734238623131735264411311
    BOOST_CHECK_CLOSE_FRACTION(lambert_wm1(BOOST_MATH_TEST_VALUE(RealType, -0.36)),
      BOOST_MATH_TEST_VALUE(RealType, -1.2227701339785059531429380734238623131735264411311),
    10 *   tolerance); // tolerance OK for quad
    // -1.2227701339785059531429380734238623131735264411311
    // -1.222770133978505953142938073423862313173526441131033

    // Just using series approximation (switch at -0.35).
    // N[productlog(-0.351), 50] = -0.72398644140937651483634596143951001600417138085814
    BOOST_CHECK_CLOSE_FRACTION(lambert_w0(BOOST_MATH_TEST_VALUE(RealType, -0.351)),
      BOOST_MATH_TEST_VALUE(RealType, -0.72398644140937651483634596143951001600417138085814),
      // 2 * tolerance); // Note 2 * tolerance for PB fukushima
    // got -0.723986441409376931150560229265736446 without Halley
    // exp -0.72398644140937651483634596143951001
    // got -0.72398644140937651483634596143951029 with Halley
     10 * tolerance); // expect -0.72398644140937651 float -0.723987103 needs 10 * tolerance
     // 2 * tolerance is fine for double and up.
    // Float is OK

    // Same for W-1 branch
    BOOST_CHECK_CLOSE_FRACTION(lambert_wm1(BOOST_MATH_TEST_VALUE(RealType, -0.351)),
      BOOST_MATH_TEST_VALUE(RealType, -1.3385736984773431852492145715526995809854973408320),
      10 * tolerance); // 2 tolerance OK for quad

    // Near singularity and NOT using series approximation (switch at -0.35)
    // N[productlog(-1, -0.34), 50]
    BOOST_CHECK_CLOSE_FRACTION(lambert_wm1(BOOST_MATH_TEST_VALUE(RealType, -0.34)),
      BOOST_MATH_TEST_VALUE(RealType, -1.4512014851325470735077533710339268100722032730024),
     10 * tolerance); // tolerance OK for quad
    //

    // Decreasing z until near zero (small z) .
    //N[productlog(-1, -0.3), 50] = -1.7813370234216276119741702815127452608215583564545
    BOOST_CHECK_CLOSE_FRACTION(lambert_wm1(BOOST_MATH_TEST_VALUE(RealType, -0.3)),
     BOOST_MATH_TEST_VALUE(RealType, -1.7813370234216276119741702815127452608215583564545),
    2 * tolerance);
    //                               -1.78133702342162761197417028151274526082155835645446

    //N[productlog(-1, -0.2), 50] = -2.5426413577735264242938061566618482901614749075294
  BOOST_CHECK_CLOSE_FRACTION(lambert_wm1(BOOST_MATH_TEST_VALUE(RealType, -0.2)),
    BOOST_MATH_TEST_VALUE(RealType, -2.5426413577735264242938061566618482901614749075294),
   2 * tolerance);

  BOOST_CHECK_CLOSE_FRACTION(lambert_wm1(BOOST_MATH_TEST_VALUE(RealType, -0.1)),
    BOOST_MATH_TEST_VALUE(RealType, -3.577152063957297218409391963511994880401796257793),
    tolerance);

     //N[productlog(-1, -0.01), 50] = -6.4727751243940046947410578927244880371043455902257
  BOOST_CHECK_CLOSE_FRACTION(lambert_wm1(BOOST_MATH_TEST_VALUE(RealType, -0.01)),
    BOOST_MATH_TEST_VALUE(RealType, -6.4727751243940046947410578927244880371043455902257),
    tolerance);

  //  N[productlog(-1, -0.001), 50] = -9.1180064704027401212583371820468142742704349737639
  BOOST_CHECK_CLOSE_FRACTION(lambert_wm1(BOOST_MATH_TEST_VALUE(RealType, -0.001)),
    BOOST_MATH_TEST_VALUE(RealType, -9.1180064704027401212583371820468142742704349737639),
    tolerance);

  //  N[productlog(-1, -0.000001), 50] = -16.626508901372473387706432163984684996461726803805
  BOOST_CHECK_CLOSE_FRACTION(lambert_wm1(BOOST_MATH_TEST_VALUE(RealType, -0.000001)),
    BOOST_MATH_TEST_VALUE(RealType, -16.626508901372473387706432163984684996461726803805),
    tolerance);

  BOOST_CHECK_CLOSE_FRACTION(lambert_wm1(BOOST_MATH_TEST_VALUE(RealType, -1e-12)),
    BOOST_MATH_TEST_VALUE(RealType, -31.067172842017230842039496250208586707880448763222),
    tolerance);

  BOOST_CHECK_CLOSE_FRACTION(lambert_wm1(BOOST_MATH_TEST_VALUE(RealType, -1e-25)),
    BOOST_MATH_TEST_VALUE(RealType, -61.686695602074505366866968627049381352503620377944),
    tolerance);

  // z nearly too small.
  BOOST_CHECK_CLOSE_FRACTION(lambert_wm1(BOOST_MATH_TEST_VALUE(RealType, -2e-26)),
    BOOST_MATH_TEST_VALUE(RealType, -63.322302839923597803393585145387854867226970485197),
    tolerance* 2);

  // z very nearly too small.  G(k=64) g[63] = -1.0264389699511303e-26 to using 1.027e-26
  BOOST_CHECK_CLOSE_FRACTION(lambert_wm1(BOOST_MATH_TEST_VALUE(RealType, -1.027e-26)),
    BOOST_MATH_TEST_VALUE(RealType, -63.999444896732265186957073549916026532499356695343),
    tolerance);
  // So -64 is the most negative value that can be determined using lookup.
  // N[productlog(-1, -1.0264389699511303 * 10^-26 ), 50] -63.999999999999997947255011093606206983577811736472 == -64
  // G[k=64] = g[63] = -1.0264389699511303e-26

  // z too small for G(k=64) g[63] = -1.0264389699511303e-26 to using 1.027e-26
  //  N[productlog(-1, -10 ^ -26), 50]    = -31.067172842017230842039496250208586707880448763222
  BOOST_CHECK_CLOSE_FRACTION(lambert_wm1(BOOST_MATH_TEST_VALUE(RealType, -1e-26)),
    BOOST_MATH_TEST_VALUE(RealType, -64.026509628385889681156090340691637712441162092868),
    tolerance); //                  -64.0265121

  if (std::numeric_limits<RealType>::has_infinity)
  {
    BOOST_CHECK_EQUAL(lambert_wm1(0), -std::numeric_limits<RealType>::infinity());
  }
  if (std::numeric_limits<RealType>::has_quiet_NaN)
  {
    // BOOST_CHECK_EQUAL(lambert_w0(std::numeric_limits<RealType>::quiet_NaN()), +std::numeric_limits<RealType>::infinity()); // message is:
    // Error in function boost::math::lambert_w0<RealType>(<RealType>): Argument z is NaN!
    BOOST_CHECK_THROW(lambert_wm1(std::numeric_limits<RealType>::quiet_NaN()), std::domain_error);
  }

   // W0 Tests for too big and too small to use lookup table.
   // Exactly W = 64, not enough to be OK for lookup.
    BOOST_CHECK_CLOSE_FRACTION(lambert_w0(BOOST_MATH_TEST_VALUE(RealType, 3.9904954117194348050619127737142206366920907815909119e+29)),
    BOOST_MATH_TEST_VALUE(RealType, 64.0),
    tolerance);

    // Just below z for F[64]
    BOOST_CHECK_CLOSE_FRACTION(lambert_w0(BOOST_MATH_TEST_VALUE(RealType, 3.99045411719434e+29)),
     BOOST_MATH_TEST_VALUE(RealType, 63.999989810930513468726486827408823607175844852495), tolerance);
    // Fails for quad_float -1.22277013397850595265
    //                      -1.22277013397850595319

  // Just too big, so using log approx and Halley refinement.
  BOOST_CHECK_CLOSE_FRACTION(lambert_w0(BOOST_MATH_TEST_VALUE(RealType, 4e+29)),
    BOOST_MATH_TEST_VALUE(RealType, 64.002342375637950350970694519073803643686041499677),
    tolerance);

  // Check at reduced precision.
  BOOST_CHECK_CLOSE_FRACTION(lambert_w0(BOOST_MATH_TEST_VALUE(RealType, 4e+29), policy<digits2<11> >()),
    BOOST_MATH_TEST_VALUE(RealType, 64.002342375637950350970694519073803643686041499677),
    0.00002);  // 0.00001 fails.

  // Tests to ensure that all JM rational polynomials are being checked.

  // 1st polynomal if (z < 0.5)   // 0.05 < z < 0.5
  BOOST_CHECK_CLOSE_FRACTION(lambert_w0(BOOST_MATH_TEST_VALUE(RealType, 0.49)),
    BOOST_MATH_TEST_VALUE(RealType, 0.3465058086974944293540338951489158955895910665452626949),
    tolerance);
  BOOST_CHECK_CLOSE_FRACTION(lambert_w0(BOOST_MATH_TEST_VALUE(RealType, 0.051)),
    BOOST_MATH_TEST_VALUE(RealType, 0.04858156174600359264950777241723801201748517590507517888),
    tolerance);



  // 2st polynomal if 0.5 < z < 2
  BOOST_CHECK_CLOSE_FRACTION(lambert_w0(BOOST_MATH_TEST_VALUE(RealType, 0.51)),
    BOOST_MATH_TEST_VALUE(RealType, 0.3569144916935871518694242462560450385494399307379277704),
    tolerance);

  BOOST_CHECK_CLOSE_FRACTION(lambert_w0(BOOST_MATH_TEST_VALUE(RealType, 1.9)),
    BOOST_MATH_TEST_VALUE(RealType, 0.8291763302658400337004358009672187071638421282477162293),
    tolerance);

  // 3rd polynomials 2 < z < 6
  BOOST_CHECK_CLOSE_FRACTION(lambert_w0(BOOST_MATH_TEST_VALUE(RealType, 2.1)),
    BOOST_MATH_TEST_VALUE(RealType, 0.8752187586805470099843211502166029752154384079916131962),
    tolerance);

  BOOST_CHECK_CLOSE_FRACTION(lambert_w0(BOOST_MATH_TEST_VALUE(RealType, 5.9)),
    BOOST_MATH_TEST_VALUE(RealType, 1.422521411785098213935338853943459424120416844150520831),
    tolerance);

  // 4th polynomials 6 < z < 18
  BOOST_CHECK_CLOSE_FRACTION(lambert_w0(BOOST_MATH_TEST_VALUE(RealType, 6.1)),
    BOOST_MATH_TEST_VALUE(RealType, 1.442152194116056579987235881273412088690824214100254315),
    tolerance);

  BOOST_CHECK_CLOSE_FRACTION(lambert_w0(BOOST_MATH_TEST_VALUE(RealType, 17.9)),
    BOOST_MATH_TEST_VALUE(RealType, 2.129100923757568114366514708174691237123820852409339147),
    tolerance);

  // 5th polynomials if (z < 9897.12905874)  // 2.8 < log(z) < 9.2
  BOOST_CHECK_CLOSE_FRACTION(lambert_w0(BOOST_MATH_TEST_VALUE(RealType, 18.1)),
    BOOST_MATH_TEST_VALUE(RealType, 2.136665501382339778305178680563584563343639180897328666),
    tolerance);

  BOOST_CHECK_CLOSE_FRACTION(lambert_w0(BOOST_MATH_TEST_VALUE(RealType, 9897.)),
    BOOST_MATH_TEST_VALUE(RealType, 7.222751047988674263127929506116648714752441161828893633),
    tolerance);

  // 6th polynomials if (z < 7.896296e+13)  // 9.2 < log(z) <= 32
  BOOST_CHECK_CLOSE_FRACTION(lambert_w0(BOOST_MATH_TEST_VALUE(RealType, 9999.)),
    BOOST_MATH_TEST_VALUE(RealType, 7.231758181708737258902175236106030961433080976032516996),
    tolerance);

  BOOST_CHECK_CLOSE_FRACTION(lambert_w0(BOOST_MATH_TEST_VALUE(RealType, 7.7e+13)),
    BOOST_MATH_TEST_VALUE(RealType, 28.62069643025822480911439831021393125282095606713326376),
    tolerance);

  // 7th polynomial // 32 < log(z) < 100
  BOOST_CHECK_CLOSE_FRACTION(lambert_w0(BOOST_MATH_TEST_VALUE(RealType, 8.0e+18)),
    BOOST_MATH_TEST_VALUE(RealType, 39.84107480517853176296156400093560722439428484537515586),
    tolerance);

  // Largest 32-bit float. (Larger values for other types tested using max())
  BOOST_CHECK_CLOSE_FRACTION(lambert_w0(BOOST_MATH_TEST_VALUE(RealType, 1.e38)),
    BOOST_MATH_TEST_VALUE(RealType, 83.07844821316409592720410446942538465411465113447713574),
    tolerance);

  // Small z < 0.05  if (z < -0.051)  -0.27 < z < -0.051

  BOOST_CHECK_CLOSE_FRACTION(lambert_w0(BOOST_MATH_TEST_VALUE(RealType, -0.28)),
    BOOST_MATH_TEST_VALUE(RealType, -0.4307588745271127579165306568413721388196459822705155385),
    tolerance);

  BOOST_CHECK_CLOSE_FRACTION(lambert_w0(BOOST_MATH_TEST_VALUE(RealType, -0.05)),
    BOOST_MATH_TEST_VALUE(RealType, -0.05270598355154634795995650617915721289427674396592395160),
    tolerance);

  // Using z small series function < 0.05
  BOOST_CHECK_CLOSE_FRACTION(lambert_w0(BOOST_MATH_TEST_VALUE(RealType, 0.049)),
    BOOST_MATH_TEST_VALUE(RealType, 0.04676143671340832342497289393737051868103596756298863555),
    tolerance);

  BOOST_CHECK_CLOSE_FRACTION(lambert_w0(BOOST_MATH_TEST_VALUE(RealType, -0.049)),
    BOOST_MATH_TEST_VALUE(RealType, -0.05159448479219405354564920228913331280713177046648170658),
    tolerance);

  BOOST_CHECK_CLOSE_FRACTION(lambert_w0(BOOST_MATH_TEST_VALUE(RealType, 1e-6)),
    BOOST_MATH_TEST_VALUE(RealType, 9.999990000014999973333385416558666900096702096424344715e-7),
    tolerance);

  BOOST_CHECK_CLOSE_FRACTION(lambert_w0(BOOST_MATH_TEST_VALUE(RealType, -1e-6)),
    BOOST_MATH_TEST_VALUE(RealType, -1.000001000001500002666671875010800023343107568372593753e-6),
    tolerance);

  // Near Smallest float.
  BOOST_CHECK_CLOSE_FRACTION(lambert_w0(BOOST_MATH_TEST_VALUE(RealType, 1e-38)),
    BOOST_MATH_TEST_VALUE(RealType, 9.99999999999999999999999999999999999990000000000000000e-39),
    tolerance);

  BOOST_CHECK_CLOSE_FRACTION(lambert_w0(BOOST_MATH_TEST_VALUE(RealType, -1e-38)),
    BOOST_MATH_TEST_VALUE(RealType, -1.000000000000000000000000000000000000010000000000000000e-38),
    tolerance);





  // Similar too near zero tests for W-1 branch
  // lambert_wm1(-1.0264389699511283e-26) = -64.000000000000000
  // Exactly z for W=-64
  BOOST_CHECK_CLOSE_FRACTION(lambert_wm1(BOOST_MATH_TEST_VALUE(RealType, -1.026438969951128225904695701851094643838952857740385870e-26)),
    BOOST_MATH_TEST_VALUE(RealType, -64.000000000000000000000000000000000000),
   2 * tolerance);

  // Just more negative than G[64 max] = wm1zs[63] so can't use lookup table.
  BOOST_CHECK_CLOSE_FRACTION(lambert_wm1(BOOST_MATH_TEST_VALUE(RealType, -1.5e-27)),
    BOOST_MATH_TEST_VALUE(RealType, -65.953279000145077719128800110134854577850889171784),
    tolerance); //                  -65.9532776

  // Just less negative than G[64 max] = wm1zs[63] so can use lookup table.
  BOOST_CHECK_CLOSE_FRACTION(lambert_wm1(BOOST_MATH_TEST_VALUE(RealType, -1.1e-26)),
    BOOST_MATH_TEST_VALUE(RealType, -63.929686062157630858625440758283127600360210072859),
    tolerance);

   // N[productlog(-1, -10 ^ -26), 50]    = -31.067172842017230842039496250208586707880448763222
  BOOST_CHECK_CLOSE_FRACTION(lambert_wm1(BOOST_MATH_TEST_VALUE(RealType, -1e-26)),
    BOOST_MATH_TEST_VALUE(RealType, -64.026509628385889681156090340691637712441162092868),
    tolerance);

  // 1e-28 is too small
  //  N[productlog(-1, -10 ^ -28), 50]    = -31.067172842017230842039496250208586707880448763222
  BOOST_CHECK_CLOSE_FRACTION(lambert_wm1(BOOST_MATH_TEST_VALUE(RealType, -1e-28)),
    BOOST_MATH_TEST_VALUE(RealType, -68.702163291525429160769761667024460023336801014578),
    tolerance);

  // Check for overflow when using a double (including when using for approximate value for refinement for higher precision).

  // N[productlog(-1, -10 ^ -30), 50]    = -73.373110313822976797067478758120874529181611813766
  //BOOST_CHECK_CLOSE_FRACTION(lambert_wm1(BOOST_MATH_TEST_VALUE(RealType, -1e30)),
  //  BOOST_MATH_TEST_VALUE(RealType, -73.373110313822976797067478758120874529181611813766),
  //  tolerance);
  //unknown location : fatal error : in "test_types" :
  //class boost::exception_detail::clone_impl<struct boost::exception_detail::error_info_injector<class std::domain_error> >
  //  : Error in function boost::math::lambert_wm1<RealType>(<RealType>) :
  //  Argument z = -1.00000002e+30 out of range(z < -exp(-1) = -3.6787944) for Lambert W - 1 branch!

  BOOST_CHECK_THROW(lambert_wm1(BOOST_MATH_TEST_VALUE(RealType, -1e30)), std::domain_error);

  // Too negative
  BOOST_CHECK_THROW(lambert_wm1(-0.5), std::domain_error);

  // This fails for fixed_point type used for other tests because out of range?
    //BOOST_CHECK_CLOSE_FRACTION(lambert_w0(BOOST_MATH_TEST_VALUE(RealType, 1.0e6)),
    //BOOST_MATH_TEST_VALUE(RealType, 11.383358086140052622000156781585004289033774706019),
    //// Output from https://www.wolframalpha.com/input/?i=lambert_w0(1e6)
    //// tolerance * 1000); // fails for fixed_point type exceeds 0.00015258789063
    //  // 15.258789063
    //  // 11.383346558
    // tolerance * 100000);
  // So need to use some spot tests for specific types, or use a bigger fixed_point type.

  // Check zero.
  BOOST_CHECK_CLOSE_FRACTION(lambert_w0(BOOST_MATH_TEST_VALUE(RealType, 0.0)),
    BOOST_MATH_TEST_VALUE(RealType, 0.0),
    tolerance);
  // these fail for cpp_dec_float_50
  // 'boost::multiprecision::detail::expression<boost::multiprecision::detail::negate,boost::multiprecision::number<boost::multiprecision::backends::cpp_dec_float<50,int32_t,void>,boost::multiprecision::et_on>,void,void,void>'
  // : no appropriate default constructor available
  // TODO ???????????

 } //template <class RealType>void test_spots(RealType)

BOOST_AUTO_TEST_CASE(test_types)
{
  BOOST_MATH_CONTROL_FP;
  // BOOST_TEST_MESSAGE output only appears if command line has --log_level="message"
  // or call set_threshold_level function:
  boost::unit_test_framework::unit_test_log.set_threshold_level(boost::unit_test_framework::log_messages);
  BOOST_TEST_MESSAGE("\nTest Lambert W function for several types.");
  BOOST_TEST_MESSAGE(show_versions());  // Full version of Boost, STL and compiler info.

  // Fundamental built-in types:
  test_spots(0.0F); // float
  test_spots(0.0); // double
  if (sizeof(long double) > sizeof(double))
  { // Avoid pointless re-testing if double and long double are identical (for example, MSVC).
    test_spots(0.0L); // long double
  }

  // Multiprecision types:
  test_spots(static_cast<boost::multiprecision::cpp_bin_float_double_extended>(0));
  test_spots(static_cast<boost::multiprecision::cpp_bin_float_quad>(0));
  test_spots(static_cast<boost::multiprecision::cpp_bin_float_50>(0));
  //test_spots(static_cast<boost::multiprecision::cpp_dec_float_50>(0));


  // Fixed-point types:
  // Some fail 0.1 to 1.0 ???
  // test_spots(static_cast<boost::fixed_point::negatable<15,-16> >(0));
  // fixed_point::negatable<15,-16> has range 1.52587891e-05 to 32768, epsilon 3.05175781e-05

  //test_spots(boost::math::concepts::real_concept(0.1));  // "real_concept" - was OK.

} // BOOST_AUTO_TEST_CASE( test_types)

BOOST_AUTO_TEST_CASE(test_range_of_double_values)
{
  using boost::math::constants::exp_minus_one;
  using boost::math::lambert_w0;

  BOOST_TEST_MESSAGE("\nTest Lambert W function type double for range of values.");

  // Want to test almost largest value.
  // test_value = (std::numeric_limits<RealType>::max)() / 4;
  // std::cout << std::setprecision(std::numeric_limits<RealType>::max_digits10) << "Max value = " << test_value << std::endl;
  // Can't use a test like this for all types because max_value depends on RealType
  // and thus the expected result of lambert_w0 does too.
  //BOOST_CHECK_CLOSE_FRACTION(lambert_w0<RealType>(test_value),
  //  BOOST_MATH_TEST_VALUE(RealType, ???),
  //  tolerance);
  // So this section just tests a single type, say IEEE 64-bit double, for a range of spot values.

  typedef double RealType; // Some tests assume type is double.

  int epsilons = 1;
  RealType tolerance = boost::math::tools::epsilon<RealType>() * epsilons; // 2 eps as a fraction.
  std::cout << "Tolerance " << epsilons  << " * epsilon == " << tolerance << std::endl;

  BOOST_CHECK_CLOSE_FRACTION(lambert_w0(BOOST_MATH_TEST_VALUE(RealType, 1.0e-6)),
    BOOST_MATH_TEST_VALUE(RealType, 9.9999900000149999733333854165586669000967020964243e-7),
    // Output from https://www.wolframalpha.com/input/ N[lambert_w[1e-6],50])
    tolerance);
  BOOST_CHECK_CLOSE_FRACTION(lambert_w0(BOOST_MATH_TEST_VALUE(RealType, 0.0001)),
    BOOST_MATH_TEST_VALUE(RealType, 0.000099990001499733385405869000452213835767629477903460),
    // Output from https://www.wolframalpha.com/input/ N[lambert_w[0.001],50])
    tolerance);
  BOOST_CHECK_CLOSE_FRACTION(lambert_w0(BOOST_MATH_TEST_VALUE(RealType, 0.001)),
    BOOST_MATH_TEST_VALUE(RealType, 0.00099900149733853088995782787410778559957065467928884),
    // Output from https://www.wolframalpha.com/input/ N[lambert_w[0.001],50])
    tolerance);
  BOOST_CHECK_CLOSE_FRACTION(lambert_w0(BOOST_MATH_TEST_VALUE(RealType, 0.01)),
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

  BOOST_CHECK_CLOSE_FRACTION(lambert_w0(BOOST_MATH_TEST_VALUE(RealType, 0.05)),
    BOOST_MATH_TEST_VALUE(RealType, 0.047672308600129374726388900514160870747062965933891),
    // Output from https://www.wolframalpha.com/input/ N[lambert_w[0.01],50])
    tolerance);

  BOOST_CHECK_CLOSE_FRACTION(lambert_w0(BOOST_MATH_TEST_VALUE(RealType, 0.1)),
    BOOST_MATH_TEST_VALUE(RealType, 0.091276527160862264299895721423179568653119224051472),
    // Output from https://www.wolframalpha.com/input/ N[lambert_w[1],50])
    tolerance);

  BOOST_CHECK_CLOSE_FRACTION(lambert_w0(BOOST_MATH_TEST_VALUE(RealType, 1.)),
    BOOST_MATH_TEST_VALUE(RealType, 0.56714329040978387299996866221035554975381578718651),
    // Output from https://www.wolframalpha.com/input/ N[lambert_w[1],50])
    tolerance);

  BOOST_CHECK_CLOSE_FRACTION(lambert_w0(BOOST_MATH_TEST_VALUE(RealType, 2.)),
    BOOST_MATH_TEST_VALUE(RealType, 0.852605502013725491346472414695317466898453300151403508772),
    // Output from https://www.wolframalpha.com/input/?i=lambert_w0(2.)
    tolerance);

  BOOST_CHECK_CLOSE_FRACTION(lambert_w0(BOOST_MATH_TEST_VALUE(RealType, 3.)),
    BOOST_MATH_TEST_VALUE(RealType, 1.049908894964039959988697070552897904589466943706341452932),
    // Output from https://www.wolframalpha.com/input/?i=lambert_w0(3.)
    tolerance);

  BOOST_CHECK_CLOSE_FRACTION(lambert_w0(BOOST_MATH_TEST_VALUE(RealType, 5.)),
    BOOST_MATH_TEST_VALUE(RealType, 1.326724665242200223635099297758079660128793554638047479789),
    // Output from https://www.wolframalpha.com/input/?i=lambert_w0(0.5)
    tolerance);

  BOOST_CHECK_CLOSE_FRACTION(lambert_w0(BOOST_MATH_TEST_VALUE(RealType, 6.)),
    BOOST_MATH_TEST_VALUE(RealType, 1.432404775898300311234078007212058694786434608804302025655),
    // Output from https://www.wolframalpha.com/input/?i=lambert_w0(6)
    tolerance);

  BOOST_CHECK_CLOSE_FRACTION(lambert_w0(BOOST_MATH_TEST_VALUE(RealType, 10.)),
    BOOST_MATH_TEST_VALUE(RealType, 1.7455280027406993830743012648753899115352881290809),
    // Output from https://www.wolframalpha.com/input/ N[lambert_w[10],50])
    tolerance);

    BOOST_CHECK_CLOSE_FRACTION(lambert_w0(BOOST_MATH_TEST_VALUE(RealType, 100.)),
    BOOST_MATH_TEST_VALUE(RealType, 3.3856301402900501848882443645297268674916941701578),
    // Output from https://www.wolframalpha.com/input/?i=lambert_w0(100)
    tolerance);

  BOOST_CHECK_CLOSE_FRACTION(lambert_w0(BOOST_MATH_TEST_VALUE(RealType, 1000.)),
    BOOST_MATH_TEST_VALUE(RealType, 5.2496028524015962271260563196973062825214723860596),
    // Output from https://www.wolframalpha.com/input/?i=lambert_w0(1000)
    tolerance);

  // This fails for fixed_point type used for other tests because out of range of the type?
  BOOST_CHECK_CLOSE_FRACTION(lambert_w0(BOOST_MATH_TEST_VALUE(RealType, 1.0e6)),
    BOOST_MATH_TEST_VALUE(RealType, 11.383358086140052622000156781585004289033774706019),
    // Output from https://www.wolframalpha.com/input/?i=lambert_w0(1e6)
    tolerance); //

  // Tests for double only near the max and the singularity where Lambert_w estimates are less precise.
  //
  if (std::numeric_limits<RealType>::is_specialized)
  {
    // Check near std::numeric_limits<>::max() for type.
    //std::cout << std::setprecision(std::numeric_limits<RealType>::max_digits10)
    //  << (std::numeric_limits<double>::max)()          // == 1.7976931348623157e+308
    //  << " " << (std::numeric_limits<double>::max)()/4 // == 4.4942328371557893e+307
    //  << std::endl;

    // All these result in faulty error message
    // unknown location : fatal error : in "test_range_of_values": class boost::exception_detail::clone_impl<struct boost::exception_detail::error_info_injector<class std::domain_error> >: Error in function boost::math::lambert_w0<RealType>(<RealType>): Argument z = %1 too large.
    // I:\modular - boost\libs\math\test\test_lambert_w.cpp(456) : last checkpoint

    BOOST_CHECK_CLOSE_FRACTION(lambert_w0(1.7976931348623157e+308 / 2), // max_value/2 for IEEE 64-bit double.
      static_cast<double>(702.53487067487671916110655783739076368512998658347L),
      // N[productlog[0, 1.7976931348623157*10^308 /2],50] == 702.53487067487671916110655783739076368512998658347
      tolerance);

    BOOST_CHECK_CLOSE_FRACTION(lambert_w0(1.7976931348623157e+308 /4), // near max_value/4 for IEEE 64-bit double.
      static_cast<double>(701.8427092142920014223182853764045476L),
      // N[productlog(0, 1.7976931348623157* 10^308 /4 ), 37] =701.8427092142920014223182853764045476
      // N[productlog(0, 0.25 * 1.7976931348623157*10^307), 37]
      tolerance);

    BOOST_CHECK_CLOSE_FRACTION(lambert_w0(4.4942328371557893e+307), // max_value/4 for IEEE 64-bit double.
      static_cast<double>(701.84270921429200143342782556643059L),
      // N[lambert_w[4.4942328371557893e+307], 35]  == 701.8427092142920014334278255664305887
      // as a double == 701.83341468208209
      // Lambert computed 702.02379914670587
      0.000003); // OK Much less precise at the max edge???

    BOOST_CHECK_CLOSE_FRACTION(lambert_w0((std::numeric_limits<double>::max)()), // max_value for IEEE 64-bit double.
      static_cast<double>(703.2270331047701868711791887193075930),
      // N[productlog(0, 1.7976931348623157* 10^308), 37] = 703.2270331047701868711791887193075930
      //                                                    703.22700325995515 lambert W
      //                                                    703.22703310477016  Wolfram
      tolerance * 2e8); // OK but much less accurate near max.

  // Compare precisions very close to the singularity.
    // This test value is one epsilon close to the singularity at -exp(-1) * x
    // (below which the result has a non-zero imaginary part).
    RealType test_value = -exp_minus_one<RealType>();
    test_value += (std::numeric_limits<RealType>::epsilon() * 1);
    BOOST_CHECK_CLOSE_FRACTION(lambert_w0(test_value),
      BOOST_MATH_TEST_VALUE(RealType, -0.99999996349975895),
      tolerance * 1000000000);
    // -0.99999996788201051
    // -0.99999996349975895
    // Would not expect to get a result closer than sqrt(epsilon)?
  } //  if (std::numeric_limits<RealType>::is_specialized)

  // Can only compare float_next for specific type T = double.
  // Comparison with Wolfram N[productlog(0,-0.36787944117144228 ), 17]
  // Note big loss of precision and big tolerance needed to pass.
  BOOST_CHECK_CLOSE_FRACTION(  // Check float_next(-exp(-1) )
    lambert_w0(BOOST_MATH_TEST_VALUE(double, -0.36787944117144228)),
    BOOST_MATH_TEST_VALUE(RealType, -0.99999998496215738),
    1e8 * tolerance); // diff 6.03558e-09 v 2.2204460492503131e-16

   BOOST_CHECK_CLOSE_FRACTION(  // Check  float_next(float_next(-exp(-1) ))
    lambert_w0(BOOST_MATH_TEST_VALUE(double, -0.36787944117144222)),
    BOOST_MATH_TEST_VALUE(RealType, -0.99999997649828679),
    2e7 * tolerance);// diff 2.30785e-09 v 2.2204460492503131e-16

  // Compare with previous PB/FK computations at double precision.
  BOOST_CHECK_CLOSE_FRACTION(  // Check float_next(-exp(-1) )
    lambert_w0(BOOST_MATH_TEST_VALUE(RealType, -0.36787944117144228)),
    BOOST_MATH_TEST_VALUE(RealType, -0.99999997892657588),
    tolerance); // diff 6.03558e-09 v 2.2204460492503131e-16

  BOOST_CHECK_CLOSE_FRACTION(  // Check  float_next(float_next(-exp(-1) ))
    lambert_w0(BOOST_MATH_TEST_VALUE(RealType, -0.36787944117144222)),
    BOOST_MATH_TEST_VALUE(RealType, -0.99999997419043196),
    tolerance);// diff 2.30785e-09 v 2.2204460492503131e-16

  // Comparison with Wolfram N[productlog(0,-0.36787944117144228 ), 17]
  // Using conversion from double to higher precision cpp_bin_float_quad
  using boost::multiprecision::cpp_bin_float_quad;
  BOOST_CHECK_CLOSE_FRACTION(  // Check float_next(-exp(-1) )
    lambert_w0(BOOST_MATH_TEST_VALUE(cpp_bin_float_quad, -0.36787944117144228)),
    BOOST_MATH_TEST_VALUE(RealType, -0.99999998496215738),
    tolerance); // OK

  BOOST_CHECK_CLOSE_FRACTION(  // Check  float_next(float_next(-exp(-1) ))
    lambert_w0(BOOST_MATH_TEST_VALUE(cpp_bin_float_quad, -0.36787944117144222)),
    BOOST_MATH_TEST_VALUE(RealType, -0.99999997649828679),
    tolerance);// OK

  // z increasingly close to singularity.
  BOOST_CHECK_CLOSE_FRACTION(lambert_w0(BOOST_MATH_TEST_VALUE(RealType, -0.36)),
    BOOST_MATH_TEST_VALUE(RealType, -0.8060843159708177782855213616209920019974599683466713016),
    2 * tolerance); // -0.806084335

  BOOST_CHECK_CLOSE_FRACTION(lambert_w0(BOOST_MATH_TEST_VALUE(RealType, -0.365)),
    BOOST_MATH_TEST_VALUE(RealType, -0.8798200914159538111724840007674053239388642469453350954),
    5 * tolerance); // Note 5 * tolerance

  BOOST_CHECK_CLOSE_FRACTION(lambert_w0(BOOST_MATH_TEST_VALUE(RealType, -0.3678)),
    BOOST_MATH_TEST_VALUE(RealType, -0.9793607149578284774761844434886481686055949229547379368),
   15 * tolerance); // Note 15 * tolerance when this close to singularity.

  // Just using series approximation (Fukushima switch at -0.35, but JM at 0.01 of singularity < -0.3679).
  // N[productlog(-0.351), 50] = -0.72398644140937651483634596143951001600417138085814
  // N[productlog(-0.351), 55] = -0.7239864414093765148363459614395100160041713808581379727
  BOOST_CHECK_CLOSE_FRACTION(lambert_w0(BOOST_MATH_TEST_VALUE(RealType, -0.351)),
    BOOST_MATH_TEST_VALUE(RealType, -0.72398644140937651483634596143951001600417138085814),
    10 * tolerance); // Note was 2 * tolerance

   // Check value just not using near_singularity series approximation (and using rational polynomial instead).
    BOOST_CHECK_CLOSE_FRACTION(lambert_wm1(BOOST_MATH_TEST_VALUE(RealType, -0.3)),
      BOOST_MATH_TEST_VALUE(RealType, -1.7813370234216276119741702815127452608215583564545),
      // Output from https://www.wolframalpha.com/input/
     //N[productlog(-1, -0.3), 50] = -1.7813370234216276119741702815127452608215583564545
     tolerance);

    // Using table lookup and schroeder with decreasing z to zero.
    BOOST_CHECK_CLOSE_FRACTION(lambert_wm1(BOOST_MATH_TEST_VALUE(RealType, -0.2)),
      BOOST_MATH_TEST_VALUE(RealType, -2.5426413577735264242938061566618482901614749075294),
    // N[productlog[-1, -0.2],50] -2.5426413577735264242938061566618482901614749075294
      tolerance);

    BOOST_CHECK_CLOSE_FRACTION(lambert_wm1(BOOST_MATH_TEST_VALUE(RealType, -0.1)),
      BOOST_MATH_TEST_VALUE(RealType, -3.5771520639572972184093919635119948804017962577931),
    //N[productlog(-1, -0.1), 50] = -3.5771520639572972184093919635119948804017962577931
      tolerance);

    BOOST_CHECK_CLOSE_FRACTION(lambert_wm1(BOOST_MATH_TEST_VALUE(RealType, -0.001)),
      BOOST_MATH_TEST_VALUE(RealType, -9.1180064704027401212583371820468142742704349737639),
    //  N[productlog(-1, -0.001), 50] = -9.1180064704027401212583371820468142742704349737639
      tolerance);

    BOOST_CHECK_CLOSE_FRACTION(lambert_wm1(BOOST_MATH_TEST_VALUE(RealType, -0.000001)),
      BOOST_MATH_TEST_VALUE(RealType, -16.626508901372473387706432163984684996461726803805),
    //  N[productlog(-1, -0.000001), 50] = -16.626508901372473387706432163984684996461726803805
      tolerance);

    BOOST_CHECK_CLOSE_FRACTION(lambert_wm1(BOOST_MATH_TEST_VALUE(RealType, -1e-6)),
      BOOST_MATH_TEST_VALUE(RealType, -16.626508901372473387706432163984684996461726803805),
    //  N[productlog(-1, -10 ^ -6), 50] = -16.626508901372473387706432163984684996461726803805
      tolerance);

    BOOST_CHECK_CLOSE_FRACTION(lambert_wm1(BOOST_MATH_TEST_VALUE(RealType, -1.0e-26)),
      BOOST_MATH_TEST_VALUE(RealType, -64.026509628385889681156090340691637712441162092868),
      // Output from https://www.wolframalpha.com/input/
      // N[productlog(-1, -1 * 10^-26 ), 50] = -64.026509628385889681156090340691637712441162092868
      tolerance);

    BOOST_CHECK_CLOSE_FRACTION(lambert_wm1(BOOST_MATH_TEST_VALUE(RealType, -2e-26)),
      BOOST_MATH_TEST_VALUE(RealType, -63.322302839923597803393585145387854867226970485197),
      // N[productlog[-1, -2*10^-26],50] = -63.322302839923597803393585145387854867226970485197
      tolerance *2);

    // Smaller than lookup table, so must use approx and Halley refinements.
    BOOST_CHECK_CLOSE_FRACTION(lambert_wm1(BOOST_MATH_TEST_VALUE(RealType, -1e-30)),
      BOOST_MATH_TEST_VALUE(RealType, -73.373110313822976797067478758120874529181611813766),
    //  N[productlog(-1, -10 ^ -30), 50]    = -73.373110313822976797067478758120874529181611813766
      tolerance);

    // std::numeric_limits<RealType>::min
    std::cout.precision(std::numeric_limits<RealType>::max_digits10);
    std::cout << "(std::numeric_limits<RealType>::min)() " << (std::numeric_limits<RealType>::min)() << std::endl;

    BOOST_CHECK_CLOSE_FRACTION(lambert_wm1(BOOST_MATH_TEST_VALUE(RealType, -2.2250738585072014e-308)),
      BOOST_MATH_TEST_VALUE(RealType, -714.96865723796647086868547560654825435542227693935),
    // N[productlog[-1, -2.2250738585072014e-308],50] =  -714.96865723796647086868547560654825435542227693935
      tolerance);

    // For z = 0, W = -infinity
    if (std::numeric_limits<RealType>::has_infinity)
    {
      BOOST_CHECK_EQUAL(lambert_wm1(BOOST_MATH_TEST_VALUE(RealType, 0.)),
        -std::numeric_limits<RealType>::infinity());
    }

} // BOOST_AUTO_TEST_CASE( test_main )

  /*

  Output:


  */




