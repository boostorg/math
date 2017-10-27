// Copyright Paul A. Bristow 2016, 2017.

// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or
//  copy at http ://www.boost.org/LICENSE_1_0.txt).

/*
Implementation of the Fukushima algorithm for the Lambert W real-only function.

This code is based on the algorithm by
Toshio Fukushima, "Precise and fast computation of Lambert W-functions without
transcendental function evaluations", J. Comp. Appl. Math. 244 (2013) 77-89,
and its author's FORTRAN code,
and on a C/C++ version by Darko Veberic, darko.veberic@ijs.si
based on the algorithm and a FORTRAN version of Toshio Fukushima.

Some macros that will show some (or much) diagnostic values if #defined.

#define-able macros
BOOST_MATH_INSTRUMENT_LAMBERT_W_SMALL_Z_SERIES
BOOST_MATH_INSTRUMENT_LAMBERT_W0 // W0 branch diagnostics.
BOOST_MATH_INSTRUMENT_LAMBERT_W0 // W1 branch diagnostics.
BOOST_MATH_INSTRUMENT_LAMBERT_W_HALLEY // Halley refinement diagnostics.
BOOST_MATH_INSTRUMENT_LAMBERT_W_HALLEY_W0 // Halley refinement diagnostics only for W0 branch.
BOOST_MATH_INSTRUMENT_LAMBERT_W_HALLEY_WM1 // Halley refinement diagnostics only for W-1 branch.
BOOST_MATH_INSTRUMENT_LAMBERT_WM1_TINY // K > 64, z > -1.0264389699511303e-26
BOOST_MATH_INSTRUMENT_LAMBERT_W_SCHROEDER // Schroeder refinement diagnostics.
BOOST_MATH_INSTRUMENT_LAMBERT_W_TERMS // Number of terms used for near-singularity series.
BOOST_MATH_INSTRUMENT_LAMBERT_W0_NOT_BUILTIN // higher than built-in precision types approximation and refinement.
BOOST_MATH_INSTRUMENT_LAMBERT_W0_BISECTION // Show bisection only estimate.
BOOST_MATH_INSTRUMENT_LAMBERT_W_SINGULARITY_SERIES // Show evaluation of series near branch singularity.
BOOST_MATH_INSTRUMENT_LAMBERT_W_SMALL_Z_SERIES_ITERATIONS  // Show evaluation of series for small z.
#define BOOST_MATH_INSTRUMENT_LAMBERT_W0_LOOKUP // Show results from lookup table.
*/

#ifndef BOOST_MATH_SF_LAMBERT_W_HPP
#define BOOST_MATH_SF_LAMBERT_W_HPP

#ifdef _MSC_VER
#  pragma warning (disable: 4127) // warning C4127: conditional expression is constant
#endif // _MSC_VER

#include <boost/math/policies/error_handling.hpp>
#include <boost/math/policies/policy.hpp>
#include <boost/math/tools/promotion.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/math/special_functions/log1p.hpp> // for log (1 + x)
#include <boost/math/constants/constants.hpp> // For exp_minus_one == 3.67879441171442321595523770161460867e-01.
#include <boost/math/special_functions/next.hpp> // for float_distance
#include <boost/math/special_functions/pow.hpp> // Link fails for pow without this, even though uses std::pow.
#include <boost/math/tools/series.hpp> // series functor.
#include <boost/math/tools/polynomial.hpp>  // polynomial.
#include <boost/math/tools/rational.hpp>  // evaluate_polynomial.
#include <boost/mpl/int.hpp>
#include <boost/type_traits/is_integral.hpp>
#include <boost/math/tools/precision.hpp> // boost::math::tools::max_value().
#include <boost/math/tools/test_value.hpp> //  Macro BOOST_MATH_TEST_VALUE.

#include <limits>
#include <cmath>
#include <limits>
#include <exception>

// Needed for testing and diagnostics only.
#include <iostream>
#include <typeinfo>
#include <boost/math/special_functions/next.hpp>  // For float_distance.
//#include <boost/math/test/test_value.hpp> // For create_test_value and macro BOOST_MATH_TEST_VALUE.

typedef double lookup_t; // Type for lookup table (double or float, or even long double?)

#include "J:\Cpp\Misc\lambert_w_lookup_table_generator\lambert_w_lookup_table.ipp"
// #include "lambert_w_lookup_table.ipp" // TODO copy final version.

namespace boost
{
namespace math
{      

  namespace detail
  {

  //! Series expansion used near the singularity/branch point z = -exp(-1) = -3.6787944.
  // Some integer constants overflow so use largest integer size available (LL).
  // Wolfram InverseSeries[Series[sqrt[2(p Exp[1 + p] + 1)], {p,-1, 20}]]
  // T. Fukushima / Journal of Computational and Applied Mathematics 244 (2013) Page 85, Table 3.
  // Wolfram command used to obtain 40 series terms at 50 decimal digit precision was
  // N[InverseSeries[Series[Sqrt[2(p Exp[1 + p] + 1)], { p,-1,40 }]], 50]
  // -1+p-p^2/3+(11 p^3)/72-(43 p^4)/540+(769 p^5)/17280-(221 p^6)/8505+(680863 p^7)/43545600 ...
  // Made these constants full precision for any T using original fractions from Table 3.
  // Decimal values of specifications for built-in floating-point types below
  // are at least 21 digits precision == max_digits10 for long double.
  // Longer decimal digits strings are rationals evaluated using Wolfram.

    template<typename T = double>
    T lambert_w_singularity_series(const T p)
    {
#ifdef BOOST_MATH_INSTRUMENT_LAMBERT_W_SINGULARITY_SERIES
      std::size_t saved_precision = std::cout.precision(3);
      std::cout << "Singularity_series Lambert_w p argument = " << p << std::endl;
      std::cout
        //<< "Argument Type = " << typeid(T).name()
        //<< ", max_digits10 = " << std::numeric_limits<T>::max_digits10
        //<< ", epsilon = " << std::numeric_limits<T>::epsilon()
      << std::endl;
      std::cout.precision(saved_precision);
#endif // BOOST_MATH_INSTRUMENT_LAMBERT_W_SINGULARITY_SERIES

      static const T q[] =
      {
        -T(1), // j0
        +T(1), // j1
        -T(1) / 3, // 1/3  j2
        +T(11) / 72, // 0.152777777777777778, // 11/72 j3
        -T(43) / 540, // 0.0796296296296296296, // 43/540 j4
        +T(769) / 17280, // 0.0445023148148148148,  j5
        -T(221) / 8505, // 0.0259847148736037625,  j6
        //+T(0.0156356325323339212L), // j7
        //+T(0.015635632532333921222810111699000587889476778365667L), // j7 from Wolfram N[680863/43545600, 50]
        +T(680863uLL) / 43545600uLL, // +0.0156356325323339212, j7
        //-T(0.00961689202429943171L), // j8
        -T(1963uLL) / 204120uLL, // 0.00961689202429943171, j8
        //-T(0.0096168920242994317068391142465216539290613364687439L), // j8 from Wolfram N[1963/204120, 50]
        +T(226287557uLL) / 37623398400uLL, // 0.00601454325295611786, j9
        -T(5776369uLL) / 1515591000uLL, // 0.00381129803489199923, j10
        //+T(0.00244087799114398267L), j11 0.0024408779911439826658968585286437530215699919795550
        +T(169709463197uLL) / 69528040243200uLL, // j11
        // -T(0.00157693034468678425L), // j12  -0.0015769303446867842539234095399314115973161850314723
        -T(1118511313uLL) / 709296588000uLL, // j12
        +T(667874164916771uLL) / 650782456676352000uLL, // j13
        //+T(0.00102626332050760715L), // j13 0.0010262633205076071544375481533906861056468041465973
        -T(500525573uLL) / 744761417400uLL, // j14
        // -T(0.000672061631156136204L), j14
        //+T(1003663334225097487uLL) / 234281684403486720000uLL, // j15 0.00044247306181462090993020760858473726479232802068800 error C2177: constant too big
        //+T(0.000442473061814620910L, // j15
        BOOST_MATH_TEST_VALUE(T, +0.000442473061814620910), // j15
        // -T(0.000292677224729627445L), // j16
        BOOST_MATH_TEST_VALUE(T, -0.000292677224729627445), // j16
        //+T(0.000194387276054539318L), // j17
        BOOST_MATH_TEST_VALUE(T, 0.000194387276054539318), // j17
        //-T(0.000129574266852748819L), // j18
        BOOST_MATH_TEST_VALUE(T, -0.000129574266852748819), // j18
        //+T(0.0000866503580520812717L), // j19 N[+1150497127780071399782389/13277465363600276402995200000, 50] 0.000086650358052081271660451590462390293190597827783288
        BOOST_MATH_TEST_VALUE(T, +0.0000866503580520812717), // j19
        //-T(0.0000581136075044138168L) // j20  N[2853534237182741069/49102686267859224000000, 50] 0.000058113607504413816772205464778828177256611844221913
        // -T(2853534237182741069uLL) / 49102686267859224000000uLL  // j20 // error C2177: constant too big,
        // so must use BOOST_MATH_TEST_VALUE(T, ) format in hope of using suffix Q for quad or decimal digits string for others.
         //-T(0.000058113607504413816772205464778828177256611844221913L), // j20  N[2853534237182741069/49102686267859224000000, 50] 0.000058113607504413816772205464778828177256611844221913
        BOOST_MATH_TEST_VALUE(T, -0.000058113607504413816772205464778828177256611844221913) // j20  - last used by Fukushima
        // more terms don't seem to give any improvement (worse in fact) and are not use for many z values.
        //BOOST_MATH_TEST_VALUE(T, +0.000039076684867439051635395583044527492132109160553593), // j21
        //BOOST_MATH_TEST_VALUE(T, -0.000026338064747231098738584082718649443078703982217219), // j22
        //BOOST_MATH_TEST_VALUE(T, +0.000017790345805079585400736282075184540383274460464169), // j23
        //BOOST_MATH_TEST_VALUE(T, -0.000012040352739559976942274116578992585158113153190354), // j24
        //BOOST_MATH_TEST_VALUE(T, +8.1635319824966121713827512573558687050675701559448E-6), // j25
        //BOOST_MATH_TEST_VALUE(T, -5.5442032085673591366657251660804575198155559225316E-6) // j26
       // -T(5.5442032085673591366657251660804575198155559225316E-6L) // j26
      // 21 to 26 Added for long double.
      }; // static const T q[]

      /*
      // Temporary copy of original double values for comparison and these are reproduced well.
      static const T q[] =
      {
        -1L,  // j0
        +1L,  // j1
        -0.333333333333333333L, // 1/3 j2
        +0.152777777777777778L, // 11/72 j3
        -0.0796296296296296296L, // 43/540
        +0.0445023148148148148L,
        -0.0259847148736037625L,
        +0.0156356325323339212L,
        -0.00961689202429943171L,
        +0.00601454325295611786L,
        -0.00381129803489199923L,
        +0.00244087799114398267L,
        -0.00157693034468678425L,
        +0.00102626332050760715L,
        -0.000672061631156136204L,
        +0.000442473061814620910L,
        -0.000292677224729627445L,
        +0.000194387276054539318L,
        -0.000129574266852748819L,
        +0.0000866503580520812717L,
        -0.0000581136075044138168L // j20
      };
      */

      // Decide how many series terms to use, increasing as z approaches the singularity,
      // balancing run time and computational noise from round-off.
      // In practice, we truncate the series expansion at a certain order.
      // If the order is too large, not only does the amount of computation increase,
      // but also the round-off errors accumulate.
      // See Fukushima equation 35, page 85 for logic of choice of number of series terms.

      // TODO remove temporary check on different between double and quad values. worst 1e-17
      //std::cout << "q[21]" << std::endl;
      //for (int i = 0; i != 21; i++)
      //{
      //  std::cout << i << " " << q[i] << " " << q[i] - double_q[i] << std::endl;
      //}

      BOOST_MATH_STD_USING // Aid argument dependent lookup (ADL) of abs.

      const T ap = abs(p);

#ifdef BOOST_MATH_INSTRUMENT_LAMBERT_W_TERMS
      {
        int terms = 20; // Default to using all terms.
        if (ap < 0.001150)
        { // Very near singularity.
          terms = 6;
        }
        else if (ap < 0.0766)
        { // Near singularity.
          terms = 10;
        }
        std::streamsize saved_precision = std::cout.precision(3);
        std::cout << "abs(p) = " << ap << ", terms = " << terms << std::endl;
        std::cout.precision(saved_precision);
      }
#endif // BOOST_MATH_INSTRUMENT_LAMBERT_W_TERMS

    if (ap < 0.01159)
    { // Only 6 near-singularity series terms are useful.
      return
        -1 +
        p*(1 +
          p*(q[2] +
            p*(q[3] +
              p*(q[4] +
                p*(q[5] +
                  p*q[6]
                  )))));
    }
    else if (ap < 0.0766) // Use 10 near-singularity series terms.
    { // Use 10 near-singularity series terms.
      return
        -1 +
        p*(1 +
          p*(q[2] +
            p*(q[3] +
              p*(q[4] +
                p*(q[5] +
                  p*(q[6] +
                    p*(q[7] +
                      p*(q[8] +
                        p*(q[9] +
                          p*q[10]
                          )))))))));
    }
    else
    { // Use all 20 near-singularity series terms.
      return
        -1 +
        p*(1 +
          p*(q[2] +
            p*(q[3] +
              p*(q[4] +
                p*(q[5] +
                  p*(q[6] +
                    p*(q[7] +
                      p*(q[8] +
                        p*(q[9] +
                          p*(q[10] +
                            p*(q[11] +
                              p*(q[12] +
                                p*(q[13] +
                                  p*(q[14] +
                                    p*(q[15] +
                                      p*(q[16] +
                                        p*(q[17] +
                                          p*(q[18] +
                                            p*(q[19] +
                                              p*q[20] // Last Fukushima term.
                                               )))))))))))))))))));
      //                                                + // more terms for more precise T: long double ...
      //// but makes almost no difference, so don't use more terms?
      //                                          p*q[21] +
      //                                            p*q[22] +
      //                                              p*q[23] +
      //                                                p*q[24] +
      //                                                 p*q[25]
      //                                         )))))))))))))))))));
    }
  } // template<typename T = double> T lambert_w_singularity_series(const T p)

  /////////////////////////////////////////////////////////////////////////////////////////////

  //! Series expansion used near zero (abs(z) < 0.05).
  // Coefficients of the inverted series expansion of the Lambert W function around z = 0.
  // Tosio Fukushima always uses all 17 terms of a Taylor series computed using Wolfram with
  //   InverseSeries[Series[z Exp[z],{z,0,17}]]
  // Tosio Fukushima / Journal of Computational and Applied Mathematics 244 (2013) page 86.

  // InverseSeries[Series[z Exp[z],{z,0,34}]] also computed.

  // Decimal values of specifications for built-in floating-point types below
  // are 21 digits precision == max_digits10 for long double.
  // TODO Might these overflow some fixed_point?

  // This version is intended to allow use by user-defined types like Boost.Multiprecision quad and cpp_dec_float types.
  // The three specializations below for built-in float, double (and perhaps long double) will be chosen in preference for these types.

  // This version uses rationals computed by Wolfram as far as possible,
  // limited by maximum size of uLL integers.
  // For higher term, uses decimal digit strings computed by Wolfram up to the maximum possible using uLL rationals,
  // and then higher coefficients are computed as necessary using function lambert_w0_small_z_series_term
  // until the precision required by the policy is achieved.

  // Series evaluation for LambertW(z) as z -> 0.
  // See http://functions.wolfram.com/ElementaryFunctions/ProductLog/06/01/01/0003/
  //  http://functions.wolfram.com/ElementaryFunctions/ProductLog/06/01/01/0003/MainEq1.L.gif

  //! \brief  lambert_w0_small_z uses a tag_type to select a variant depending on the size of the type.
  //! The Lambert W is computed by lambert_w0_small_z for small z.
  //! The cutoff for z smallness determined by Tosio Fukushima by trial and error is (abs(z) < 0.05),
  //! but the optimum might be a function of the size of the type of z.

  //! \details
  //! The tag_type selection is based on the value std::numeric_limits<T>::max_digits10.
  //! This allows distinguishing between long double types that commonly vary between 64 and 80-bits,
  //! and also compilers that have a float type using 64 bits and/or long double using 128-bits.
  //! It assumes that max_digits10 is defined correctly or this might fail to make the correct selection.
  //! causing very small differences in computing lambert_w that would be very difficult to detect and diagnose.
  //! Cannot switch on std::numeric_limits<>::max() because comparison values may overflow the compiler limit.
  //! Cannot switch on std::numeric_limits<long double>::max_exponent10() because both 80 and 128 bit float use 11 bits for exponent.
  //! So must rely on std::numeric_limits<long double>::max_digits10.

  //! Specialization of float zero series expansion used for small z (abs(z) < 0.05).
  //! Specializations of lambert_w0_small_z for built-in types.
  //! These specializations should be chosen in preference to T version.
  //! For example: lambert_w0_small_z(0.001F) should use the float version.
  // Policy is not used by built-in types when all terms are used during an inline computation,
  // but for the ta_type selection to work, they all must include Policy the their signature.

  // Forward declaration of variants of lambert_w0_small_z.
  template <class T, class Policy>
  T lambert_w0_small_z(T x, const Policy&, boost::mpl::int_<0> const&);   //  for float (32-bit) type.

  template <class T, class Policy>
  T lambert_w0_small_z(T x, const Policy&, boost::mpl::int_<1> const&);   //  for double (64-bit) type.

  template <class T, class Policy>
  T lambert_w0_small_z(T x, const Policy&, boost::mpl::int_<2> const&);   //  for long double (double extended 80-bit) type.

  template <class T, class Policy>
  T lambert_w0_small_z(T x, const Policy&, boost::mpl::int_<3> const&);   //  for float128

  template <class T, class Policy>
  T lambert_w0_small_z(T x, const Policy&, boost::mpl::int_<4> const&);   //  Generic multiprecision.

  // Set tag_type depending on max_digits10.
  template <class T, class Policy>
  T lambert_w0_small_z(T x, const Policy& pol)
  {
    typedef boost::mpl::int_
      <
      std::numeric_limits<T>::max_digits10 <= 9 ? 0 :  //  for float
      std::numeric_limits<T>::max_digits10 <= 17 ? 1 : //  for double 64-bit
      std::numeric_limits<T>::max_digits10 <= 22 ? 2 :  //  for  80-bit double extended
      std::numeric_limits<T>::max_digits10 < 37 ? 3 : //  for float128  128-bit quad
      4 // Generic multiprecision types.
      > tag_type;
    // std::cout << "tag type = " << tag_type << std::endl; // error C2275: 'tag_type': illegal use of this type as an expression
    return lambert_w0_small_z(x, pol, tag_type());
  } // template <class T> T lambert_w0_small_z(T x)

  //! Specialization of float (32-bit) series expansion used for small z (abs(z) < 0.05).
  // Only 9 Coefficients are computed to 21 decimal digits precision, ample for 32-bit float used by most platforms.
  // Taylor series coefficients used are computed by Wolfram to 50 decimal digits using instruction
  // N[InverseSeries[Series[z Exp[z],{z,0,34}]],50],
  // as proposed by Tosio Fukushima and implemented by Darko Veberic.

  template <class T, class Policy>
  T lambert_w0_small_z(T z, const Policy&, boost::mpl::int_<0> const&)
  {
#ifdef BOOST_MATH_INSTRUMENT_LAMBERT_W_SMALL_Z_SERIES
    std::cout << "float lambert_w0_small_z called with z = " << z << " using " << 9 << " terms of precision "
      << std::numeric_limits<float>::max_digits10 << " decimal digits. " << std::endl;
#endif // BOOST_MATH_INSTRUMENT_LAMBERT_W_SMALL_Z_SERIES
    return
      z*(1 - // j1 z^1 term = 1
        z*(1 -  // j2 z^2 term = -1
          z*(static_cast<float>(3uLL) / 2uLL - // 3/2 // j3 z^3 term = 1.5.
            z*(2.6666666666666666667F -  // 8/3 // j4
              z*(5.2083333333333333333F - // -125/24 // j5
                z*(10.8F - // j6
                  z*(23.343055555555555556F - // j7
                    z*(52.012698412698412698F - // j8
                      z*118.62522321428571429F)))))))); // j9
  } // template <class T>   T lambert_w0_small_z(T x, mpl::int_<0> const&)

    //! Specialization of long double (80-bit double extended) series expansion used for small z (abs(z) < 0.05).
    // 17 Coefficients are computed to 21 decimal digits precision suitable for 64-bit double used by most platforms.
    // Taylor series coefficients used are computed by Wolfram to 50 decimal digits using instruction
    // N[InverseSeries[Series[z Exp[z],{z,0,34}]],50], as proposed by Tosio Fukushima and implemented by Veberic.

  template <class T, class Policy>
  T lambert_w0_small_z(const T z, const Policy&, boost::mpl::int_<1> const&)
  {
#ifdef BOOST_MATH_INSTRUMENT_LAMBERT_W_SMALL_Z_SERIES
    std::cout << "double lambert_w0_small_z called with z = " << z << " using " << 17 << " terms of precision, "
      << std::numeric_limits<double>::max_digits10 << " decimal digits. " << std::endl;
#endif // BOOST_MATH_INSTRUMENT_LAMBERT_W_SMALL_Z_SERIES
    T  result =
      z*(1. - // j1 z^1
        z*(1. -  // j2 z^2
          z*(1.5 - // 3/2 // j3 z^3
            z*(2.6666666666666666667 -  // 8/3 // j4
              z*(5.2083333333333333333 - // -125/24 // j5
                z*(10.8 - // j6
                  z*(23.343055555555555556 - // j7
                    z*(52.012698412698412698 - // j8
                      z*(118.62522321428571429 - // j9
                        z*(275.57319223985890653 - // j10
                          z*(649.78717234347442681 - // j11
                            z*(1551.1605194805194805 - // j12
                              z*(3741.4497029592385495 - // j13
                                z*(9104.5002411580189358 - // j14
                                  z*(22324.308512706601434 - // j15
                                    z*(55103.621972903835338 - // j16
                                      z*136808.86090394293563)))))))))))))))); // j17 z^17
    return result;
  }  // T lambert_w0_small_z(const T z, boost::mpl::int_<1> const&)

     //! Specialization of long double (80-bit double extended) series expansion used for small z (abs(z) < 0.05).
     // 21 Coefficients are computed to 21 decimal digits precision suitable for 80-bit long double used by some
     // platforms including GCC and Clang when generating for Intel X86 floating-point processors with 80-bit operations enabled (the default).
     // (This is NOT used by Microsoft Visual Studio where double and long always both use only 64-bit type.)
  template <class T, class Policy>
  T lambert_w0_small_z(const T z, const Policy&, boost::mpl::int_<2> const&)
  {
#ifdef BOOST_MATH_INSTRUMENT_LAMBERT_W_SMALL_Z_SERIES
    std::cout << "long double (80-bit double extended) lambert_w0_small_z called with z = " << z << " using " << 21 << " terms of precision "
      << std::numeric_limits<long double>::max_digits10 << " decimal digits. " << std::endl;
#endif // BOOST_MATH_INSTRUMENT_LAMBERT_W_SMALL_Z_SERIES
    T  result =
      z*(1.L - // j1 z^1
        z*(1.L -  // j2 z^2
          z*(1.5L - // 3/2 // j3
            z*(2.6666666666666666667L -  // 8/3 // j4
              z*(5.2083333333333333333L - // -125/24 // j5
                z*(10.800000000000000000L - // j6
                  z*(23.343055555555555556L - // j7
                    z*(52.012698412698412698L - // j8
                      z*(118.62522321428571429L - // j9
                        z*(275.57319223985890653L - // j10
                          z*(649.78717234347442681L - // j11
                            z*(1551.1605194805194805L - // j12
                              z*(3741.4497029592385495L - // j13
                                z*(9104.5002411580189358L - // j14
                                  z*(22324.308512706601434L - // j15
                                    z*(55103.621972903835338L - // j16
                                      z*(136808.86090394293563L - // j17 z^17  last term used by Fukushima double.
                                        z*(341422.050665838363317L - // z^18
                                          z*(855992.9659966075514633L - // z^19
                                            z*(2.154990206091088289321e6L - // z^20
                                              z*5.4455529223144624316423e6L   // z^21
                                              ))))))))))))))))))));

    return result;
  }  // T lambert_w0_small_z(const T z, boost::mpl::int_<1> const&)

    //! Specialization of long double (128-bit quad) series expansion used for small z (abs(z) < 0.05).
    // 34 Taylor series coefficients used are computed by Wolfram to 50 decimal digits using instruction
    // N[InverseSeries[Series[z Exp[z],{z,0,34}]],50],
    // and are suffixed by L as they are assumed of type long double.
    // (This is NOT used for 128-bit quad which required a suffix Q
    // nor multiprecision type cpp_bin_float_quad that can only be initialised at full precision of the type
    // using a decimal digit string like "2.6666666666666666666666666666666666666666666666667".)

  template <class T, class Policy>
  T lambert_w0_small_z(const T z, const Policy&, boost::mpl::int_<3> const&)
  {
#ifdef BOOST_MATH_INSTRUMENT_LAMBERT_W_SMALL_Z_SERIES
    std::cout << "long double (128-bit quad) lambert_w0_small_z called with z = " << z << " using " << 34 << " terms of precision "
      << std::numeric_limits<double>::max_digits10 << " decimal digits. " << std::endl;
#endif // BOOST_MATH_INSTRUMENT_LAMBERT_W_SMALL_Z_SERIES
    T  result =
      z*(1.L - // j1
        z*(1.L -  // j2
          z*(1.5L - // 3/2 // j3
            z*(2.6666666666666666666666666666666666L -  // 8/3 // j4
              z*(5.2052083333333333333333333333333333L - // -125/24 // j5
                z*(-10.800000000000000000000000000000000L - // j6
                  z*(23.343055555555555555555555555555555L - // j7
                    z*(52.0126984126984126984126984126984126L - // j8
                      z*(118.625223214285714285714285714285714L - // j9
                        z*(275.57319223985890652557319223985890L - // * z ^ 10 - // j10
                          z*(649.78717234347442680776014109347442680776014109347L - // j11
                            z*(1551.1605194805194805194805194805194805194805194805L - // j12
                              z*(3741.4497029592385495163272941050718828496606274384L - // j13
                                z*(9104.5002411580189357967135744913522691300469078247L - // j14
                                  z*(22324.308512706601434280005708577137148565719994291L - // j15
                                    z*(55103.621972903835337697771560205422639285073147507L - // j16
                                      z*136808.86090394293563342215789305736395683485630576L)))))))))))))))); // j17
    return result;
  }  // T lambert_w0_small_z(const T z, boost::mpl::int_<1> const&)

  //! Series functor to compute series term using pow and factorial.
  //! \details Functor is called after evaluating polynomial with the coefficients as rationals below.
  template <class T>
  struct lambert_w0_small_z_series_term
  {
    typedef T result_type;
    //! \param _z Lambert W argument z.
    //! \param -term  -pow<18>(z) / 6402373705728000uLL
    //! \param _k number of terms == initially 18

    //  Note *after* evaluating N terms, its internal state has k = N and term = (-1)^N z^N.

    lambert_w0_small_z_series_term(T _z, T _term, int _k)
      : k(_k), z(_z), term(_term) { }

    T operator()()
    { // Called by sum_series until needs precision set by factor (policy::get_epsilon).
      using std::pow;
      ++k;
      term *= -z / k;
      //T t = pow(z, k) * pow(T(k), -1 + k) / factorial<T>(k); // (z^k * k(k-1)^k) / k!
      T result = term * pow(T(k), -1 + k); // term * k^(k-1)
                                           // std::cout << " k = " << k << ", term = " << term << ", result = " << result << std::endl;
      return result; //
    }
  private:
    int k;
    T z;
    T term;
  }; // template <class T> struct lambert_w0_small_z_series_term

  //! Generic variant for User-defined types like Boost.Multiprecision.
  template <class T, class Policy>
  inline T lambert_w0_small_z(T z, const Policy& pol, boost::mpl::int_<4> const&)
  {
#ifdef BOOST_MATH_INSTRUMENT_LAMBERT_W_SMALL_Z_SERIES
    std::cout << "Generic lambert_w0_small_z called with z = " << z << " using as many terms needed for precision." << std::endl;
    std::cout << "Argument z is of type " << typeid(T).name() << std::endl;
#endif // BOOST_MATH_INSTRUMENT_LAMBERT_W_SMALL_Z_SERIES

    // First several terms of the series are tabulated and evaluated as a polynomial:
    // this will save us a bunch of expensive calls to pow.
    // Then our series functor is initialized "as if" it had already reached term 18,
    // enough evaluation of built-in 64-bit double and float (and 80-bit long double?) types.

    // Coefficients should be stored such that the coefficients for the x^i terms are in poly[i].
    static const T coeff[] =
    {
      0, // z^0  Care: zeroth term needed by tools::evaluate_polynomial, but not in the Wolfram equation, so indexes are one different!
      1, // z^1 term.
      -1, // z^2 term
      static_cast<T>(3uLL) / 2uLL, // z^3 term.
      -static_cast<T>(8uLL) / 3uLL, // z^4
      static_cast<T>(125uLL) / 24uLL, // z^5
      -static_cast<T>(54uLL) / 5uLL, // z^6
      static_cast<T>(16807uLL) / 720uLL, // z^7
      -static_cast<T>(16384uLL) / 315uLL, // z^8
      static_cast<T>(531441uLL) / 4480uLL, // z^9
      -static_cast<T>(156250uLL) / 567uLL, // z^10
      static_cast<T>(2357947691uLL) / 3628800uLL, // z^11
      -static_cast<T>(2985984uLL) / 1925uLL, // z^12
      static_cast<T>(1792160394037uLL) / 479001600uLL, // z^13
      -static_cast<T>(7909306972uLL) / 868725uLL, // z^14
      static_cast<T>(320361328125uLL) / 14350336uLL, // z^15
      -static_cast<T>(35184372088832uLL) / 638512875uLL, // z^16
      static_cast<T>(2862423051509815793uLL) / 20922789888000uLL, // z^17 term
      -static_cast<T>(5083731656658uLL) / 14889875uLL,
      // z^18 term. = 136808.86090394293563342215789305735851647769682393

      // z^18 is biggest that can be computed as rational using the largest possible uLL integers,
      // so higher terms cannot be potentially compiler-computed as uLL rationals.
      // Wolfram (5083731656658 z ^ 18) / 14889875 or
      // -341422.05066583836331735491399356945575432970390954 z^18

      // See note below calling the functor to compute another term,
      // sufficient for 80-bit long double precision.
      // Wolfram -341422.05066583836331735491399356945575432970390954 z^19 term.
      // (5480386857784802185939 z^19)/6402373705728000
      // But now this variant is not used to compute long double
      // as specializations are provided above.
    }; // static const T coeff[]

       /*
       Table of 19 computed coefficients:

       #0 0
       #1 1
       #2 -1
       #3 1.5
       #4 -2.6666666666666666666666666666666665382713370408509
       #5 5.2083333333333333333333333333333330765426740817019
       #6 -10.800000000000000000000000000000000616297582203915
       #7 23.343055555555555555555555555555555076212991619177
       #8 -52.012698412698412698412698412698412659282693193402
       #9 118.62522321428571428571428571428571146835390992496
       #10 -275.57319223985890652557319223985891400375196748314
       #11 649.7871723434744268077601410934743969785223845882
       #12 -1551.1605194805194805194805194805194947599566007429
       #13 3741.4497029592385495163272941050719510009019331763
       #14 -9104.5002411580189357967135744913524243896052869184
       #15 22324.308512706601434280005708577137322392070452582
       #16 -55103.621972903835337697771560205423203318720697224
       #17 136808.86090394293563342215789305735851647769682393
       136808.86090394293563342215789305735851647769682393   == Exactly same as Wolfram computed value.
       #18 -341422.05066583836331735491399356947486381600607416
       341422.05066583836331735491399356945575432970390954  z^19  Wolfram value differs at 36 decimal digit, as expected.
       */

    using boost::math::policies::get_epsilon;
    using boost::math::tools::sum_series;
    using boost::math::tools::evaluate_polynomial;
    // http://www.boost.org/doc/libs/release/libs/math/doc/html/math_toolkit/roots/rational.html

    // std::streamsize prec = std::cout.precision(std::numeric_limits <T>::max_digits10);

    T result = evaluate_polynomial(coeff, z);
    //  template <std::size_t N, class T, class V>
    //  V evaluate_polynomial(const T(&poly)[N], const V& val);
    // Size of coeff found from N
    //std::cout << "evaluate_polynomial(coeff, z); == " << result << std::endl;
    //std::cout << "result = " << result << std::endl;
    // It's an artefact of the way I wrote the functor: *after* evaluating N
    // terms, its internal state has k = N and term = (-1)^N z^N.  So after
    // evaluating 18 terms, we initialize the functor to the term we've just
    // evaluated, and then when it's called, it increments itself to the next term.
    // So 18!is 6402373705728000, which is where that comes from.

    // The 19th coefficient of the polynomial is actually, 19 ^ 18 / 19!=
    // 104127350297911241532841 / 121645100408832000 which after removing GCDs
    // reduces down to Wolfram rational 5480386857784802185939 / 6402373705728000.
    // Wolfram z^19 term +(5480386857784802185939 z^19) /6402373705728000
    // +855992.96599660755146336302506332246623424823099755 z^19

    //! Evaluate Functor.
    lambert_w0_small_z_series_term<T> s(z, -pow<18>(z) / 6402373705728000uLL, 18);

    // Temporary to list the coefficients.
    //std::cout << " Table of coefficients" << std::endl;
    //std::streamsize saved_precision = std::cout.precision(50);
    //for (size_t i = 0; i != 19; i++)
    //{
    //  std::cout << "#" << i << " " << coeff[i] << std::endl;
    //}
    //std::cout.precision(saved_precision);

    boost::uintmax_t max_iter = policies::get_max_series_iterations<Policy>(); // Max iterations from policy.

                                                                               //std::cout << "max iter from policy = " << max_iter << std::endl; //   max iter from policy = 1000000
                                                                               // max_iter = 0;
    result = sum_series(s, get_epsilon<T, Policy>(), max_iter, result);
    // result == evaluate_polynomial.
    //sum_series(Functor& func, int bits, boost::uintmax_t& max_terms, const U& init_value)
    // std::cout << "sum_series(s, get_epsilon<T, Policy>(), max_iter, result); = " << result << std::endl;

    //T epsilon = get_epsilon<T, Policy>();
    //std::cout << "epilson from policy = " << epsilon << std::endl;
    // epilson from policy = 1.93e-34 for T == quad
    //  5.35e-51 for t = cpp_bin_float_50

    // std::cout << " get eps = " << get_epsilon<T, Policy>() << std::endl; // quad eps = 1.93e-34, bin_float_50 eps = 5.35e-51
    policies::check_series_iterations<T>("boost::math::lambert_w0_small_z<%1%>(%1%)", max_iter, pol);
#ifdef BOOST_MATH_INSTRUMENT_LAMBERT_W_SMALL_Z_SERIES_ITERATIONS
    std::cout << "z = " << z << " needed  " << max_iter << " iterations." << std::endl;
#endif // BOOST_MATH_INSTRUMENT_LAMBERT_W_SMALL_Z_SERIES_ITERATIONS
    //std::cout.precision(prec); // Restore.
    return result;
  } // template <class T, class Policy> inline T lambert_w0_small_z_series(T z, const Policy& pol)

  /////////////////////////////////////////////////////////////////////////////

  // Halley refinement, if required, for precision in policy.
  // Only perform this if the most precise value possible (usually within one epsilon) is desired.
  // The previous Fukushima Schroeder refinement usually gives within several epsilon.
  // Halley refinement needs to perform some evaluations of the exp function, and
  // may nearly double the execution time.
  // However Boost.Math prioritizes accuracy over speed,
  // so this extra step is taken by the default policy.

  // Two versions:
  //  do_while always makes one iteration before testing against a tolerance,
  //  while_do does a test first, perhaps avoiding an iteration if already within tolerance.

  template<typename T = double, class Policy>
  inline
  T halley_update(T w0, const T z, const Policy&)
  {
    // Iterate a few times to refine value using Halley's method.
    // Inline Halley iteration, rather than calling boost::math::tools::halley_iterate
    // since we can simplify the expressions algebraically,
    // and don't need most of the error checking of the boilerplate version
    // as we know in advance that the function is reasonably well-behaved,
    // (and in any case the derivatives require evaluation of Lambert W!)

    BOOST_MATH_STD_USING // Aid argument dependent (ADL) lookup of abs.

    using boost::math::constants::exp_minus_one; // 0.36787944
    using boost::math::tools::max_value;

    T tolerance = boost::math::policies::get_epsilon<T, Policy>();

    // TODO should get_precision here.
    int iterations = 0;
    int iterations_required = 6;
    int max_iterations = 20;

    T w1 = w0; // Refined estimate.
    T previous_diff = boost::math::tools::max_value<T>();

    T expw0 = exp(w0); // Compute z == W * exp(W); from best Lambert W estimate so far.
    // Hope that  w0 * expw0 == z;
    T diff = (w0 * expw0) - z; // Difference from argument z.
    // std::cout << "w = " << w0 << ", z = " << z << ", exp(w) = " << expw0 << ", diff = " << diff << std::endl;

#ifdef BOOST_MATH_INSTRUMENT_LAMBERT_W_HALLEY
    std::streamsize precision = std::cout.precision(std::numeric_limits<T>::max_digits10);  // Save.
    std::cout.precision(std::numeric_limits<T>::max_digits10); // Show all posssibly significant digits.
    std::ios::fmtflags flags(std::cout.flags()); // Save.
    std::cout.setf(std::ios_base::showpoint); // Include any trailing zeros.
    std::cout << "w = " << w0 << ", z = " << z << ", w * exp(w) = " << w0 * expw0 << ", diff = " << diff << std::endl;
    std::cout.precision(precision); // Restore.
    std::cout.flags(flags);
#endif // BOOST_MATH_INSTRUMENT_LAMBERT_W_HALLEY
if (diff == 0)  // Exact result - common.
    {
      return w0;
    }
    // Refine.
    do
    {
      iterations++;
#ifdef BOOST_MATH_INSTRUMENT_LAMBERT_W_HALLEY
      std::cout << iterations << " Halley iterations." << std::endl;
#endif // BOOST_MATH_INSTRUMENT_LAMBERT_W_HALLEY
      previous_diff = diff; // Remember so that can check if any new estimate is better.

      // Halley's method from Luu equation 6.39, line 17.
      // https://en.wikipedia.org/wiki/Halley%27s_method
      // http://www.wolframalpha.com/input/?i=differentiate+productlog(x)
      // differentiate productlog(x)  d/dx(W(x)) = (W(x))/(x W(x) + x) = e^w (1 + w)
      // Wolfram Alpha (d^2 )/(dw^2)(w exp(w) - z) = e^w (w + 2)
      // f'(w) = e^w (1 + w)
      // f''(w) = e^w (2 + w) ,
      // f''(w)/f'(w) = (2+w) / (1+w),  Luu equation 6.38.

      w1 = w0 // Refine a new estimate using Halley's method using Luu equation 6.39.
        - diff / ((expw0 * (w0 + 1) - (w0 + 2) * diff / (w0 + w0 + 2)));
      diff = (w1 * exp(w1)) - z;

#ifdef BOOST_MATH_INSTRUMENT_LAMBERT_W_HALLEY
      T dis = boost::math::float_distance<T>(w0, w1);
      std::cout << "w = " << w0 << ", z = " << z << ", w * exp(w) = " << w0 * expw0 << ", diff = " << diff << std::endl;
      int d = static_cast<int>(dis);
      if (abs(dis) < 2147483647)
      {
        std::cout << "float_distance = " << d << std::endl;
      }
#endif // BOOST_MATH_INSTRUMENT_LAMBERT_W_HALLEY

      if (diff == 0) // Exact.
      {
        return w1; // Exact, or no improvement, so as good as it can be.
      }
      if (abs(diff) >= abs(previous_diff)) // Latest estimate is not better, or worse, so avoid oscillating result (common when near singularity).
      {
        return w1;  // Or mean (w0 + w1)/2 ??
      }
      if (fabs((w0 / w1) - static_cast<T>(1)) < tolerance)
      { // Reached estimate of Lambert W within relative tolerance (usually an epsilon or few).
        return w1; // Or mean (w0 + w1)/2 ??
      }
      w0 = w1;
      expw0 = exp(w0);
    }
    while ((iterations < iterations_required) || (iterations <= max_iterations)); // Absolute limit during testing - is looping if need this.

    return w1;
  } // T halley_update(T w0, const T z)

  template<typename T = double>
  inline
  T halley_update_do_while(T w0, const T z)
  {
    BOOST_MATH_STD_USING // Aid argument dependent lookup of abs.

    using boost::math::constants::exp_minus_one; // 0.36787944
    using boost::math::tools::max_value;
    T tolerance = std::numeric_limits<T>::epsilon();
    int iterations = 0;
    int min_iterations = 1; // Specify a minimum number of iterations to perform.
    int max_iterations = 6; // Should not be needed, but include to avoid risk of looping forever.
#ifdef BOOST_MATH_INSTRUMENT_LAMBERT_W_HALLEY
    std::cout << "At least = " << min_iterations << " Halley iterations." << std::endl;
    std::cout << "(But not more than " << max_iterations << " Halley iterations.)" << std::endl;
#endif // BOOST_MATH_INSTRUMENT_LAMBERT_W_HALLEY

    T w1; // Refined estimate.
    T previous_diff = boost::math::tools::max_value<T>();

    // TODO perform check first with while( diff small) do {refinement}
    do
    { // Iterate a few times to refine value using Halley's method.
      // Inline Halley iteration, rather than calling boost::math::tools::halley_iterate
      // since we can simplify the expressions algebraically,
      // and don't need most of the error checking of the boilerplate version
      // as we know in advance that the function is reasonably well-behaved,
      // and in any case the derivatives require evaluation of Lambert W!

      // z == W * exp(W);
      T expw0 = exp(w0); // Compute z from best Lambert W estimate so far.
                                // Hope that  w0 * expw0 == z;
      T diff = (w0 * expw0) - z; // Difference from argument z.

#ifdef BOOST_MATH_INSTRUMENT_LAMBERT_W_HALLEY
        std::cout << "w = " << w0 << ", z = " << z << ", exp(w) = " << expw0 << ", diff = " << diff << std::endl;
#endif // BOOST_MATH_INSTRUMENT_LAMBERT_W_HALLEY

      iterations++;
      // Halley's method from Luu equation 6.39, line 17.
      // https://en.wikipedia.org/wiki/Halley%27s_method
      // http://www.wolframalpha.com/input/?i=differentiate+productlog(x)
      // differentiate productlog(x)  d/dx(W(x)) = (W(x))/(x W(x) + x) = e^w (1 + w)
      // Wolfram Alpha (d^2 )/(dw^2)(w exp(w) - z) = e^w (w + 2)
      // f'(w) = e^w (1 + w)
      // f''(w) = e^w (2 + w) ,
      // f''(w)/f'(w) = (2+w) / (1+w),  Luu equation 6.38.

      w1 = w0 // Refine a new estimate using Halley's method using Luu equation 6.39.
        - diff / ((expw0 * (w0 + 1) - (w0 + 2) * diff / (w0 + w0 + 2)));
      if ( // Reached estimate of Lambert W within relative tolerance (usually an epsilon or few).
        (fabs((w0 / w1) - static_cast<T>(1)) < tolerance)
        || // Or latest estimate is not better, or worse, so avoid oscillating result (common when near singularity).
        (abs(diff) >= abs(previous_diff))
        )
      {
        return w1; // No improvement, so as good as it can be.
      }

      w0 = w1;
      previous_diff = diff; // Remember so that can check if any new estimate is better.
    }
    while (
      (iterations < min_iterations) || // do minimum iterations.
      (iterations >= max_iterations) // Absolute max limit during testing
           // (unexpected looping if need this).
    );

    return w1; //
  } // Halley do while

  //! Schroeder method, fifth-order update formula, see T. Fukushima page 80-81, and
  // A. Householder, The Numerical Treatment of a Single Nonlinear Equation,
  // McGraw-Hill, New York, 1970, section 4.4.
  // Switch to schroeder_update after a pre-computed jb/jmax bisections,
  // chosen to ensure that the result will be achieve the +/- 10 epsilon target.
  //! \param w Lambert w estimate from bisection.
  //! \param y bracketing value from bisection.
  //! \returns Refined estimate of Lambert w.

  // Schroeder refinement, called unless not required by precision policy.
  template<typename T = double>
  inline
  T schroeder_update(const T w, const T y)
  {
    // Compute derivatives using Schroeder refinement, Fukushima equations 18, page 6.
    // Since this is the final step, it will always use the highest precision type T.
    // Called: result = schroeder_update(w, y);
    // w is estimate of Lambert W (from bisection or series).
    // y is z * e^-w

#ifdef BOOST_MATH_INSTRUMENT_LAMBERT_W_SCHROEDER
    std::streamsize saved_precision = std::cout.precision(std::numeric_limits<T>::max_digits10);
    using boost::math::float_distance;
    T fd = float_distance<T>(w, y);
    std::cout << "Pre-Schroder Distance = " << static_cast<int>(fd) << std::endl;
#endif // BOOST_MATH_INSTRUMENT_LAMBERT_W_SCHROEDER

    const T f0 = w - y; // f0 = w - y.
    const T f1 = 1 + y; // f1 = df/dW
    const T f00 = f0 * f0;
    const T f11 = f1 * f1;
    const T f0y = f0 * y;
    const T result =
     w - 4 * f0 * (6 * f1 * (f11 + f0y)  +  f00 * y) /
      (f11 * (24 * f11 + 36 * f0y) +
        f00 * (6 * y * y  +  8 * f1 * y  +  f0y)); // Fukushima Page 81, equation 21 from equation 20.

#ifdef BOOST_MATH_INSTRUMENT_LAMBERT_W_SCHROEDER
    std::cout << "Schroeder refined " << w << "  " << y << " to  " << result << ", difference " << w - result << std::endl;
    std::cout.precision(saved_precision); // Restore.
#endif // BOOST_MATH_INSTRUMENT_LAMBERT_W_SCHROEDER

    return result;
  } // template<typename T = double> T schroeder_update(const T w, const T y)

  // Temporary to list out the static arrays.
  template<typename T>
  inline void print_collection(const T& coll, std::string name, std::string sep)
  {
    std::streamsize precision = std::cout.precision();
    std::cout.precision(std::numeric_limits<lookup_t>::max_digits10);
    std::cout << name << ": ";

    for (auto it = std::begin(coll); it != std::end(coll); ++it)
    {
      std::cout << *(it) << sep;
    }
    std::cout.precision(precision); // Restore.
    std::cout << std::endl;
  } // print_collection

  /////////////  namespace detail implementations of Lambert W for W0 and W-1 branches.

  //! Compute Lambert W for W0 (or W+) branch.
  template<typename T, class Policy>
  T lambert_w0_imp(const T z, const Policy& pol)
  {
    // Catch providing an integer value as parameter x to lambert_w, for example, lambert_w(1).
    // Need to ensure it is a floating-point type (of the desired type, float 1.F, double 1., or long double 1.L),
    // or static_casted integer, for example:  static_cast<float>(1) or static_cast<cpp_dec_float_50>(1).
    // Want to allow fixed_point types too, so do not just test for floating-point.
    // Integral types should been promoted to double by user Lambert w functions.
    BOOST_STATIC_ASSERT_MSG(!std::is_integral<T>::value, "Must be floating-point or fixed type (not integer type), for example: lambert_w0(1.), not lambert_w0(1)!");

#ifdef BOOST_MATH_INSTRUMENT_LAMBERT_W0
    {
      std::size_t saved_precision = std::cout.precision(std::numeric_limits<T>::max_digits10);
      std::cout << "Lambert_w argument = " << z << std::endl;
      std::cout << "Argument Type = " << typeid(T).name()
        << ", max_digits10 = " << std::numeric_limits<T>::max_digits10
        << std::setprecision(3) << ", epsilon = " << std::numeric_limits<T>::epsilon() << std::endl;
      std::cout.precision(saved_precision); // Restore precision.
    }
#endif // BOOST_MATH_INSTRUMENT_LAMBERT_W0

    const char* function = "boost::math::lambert_w0<RealType>(<RealType>)"; // Used for error messages.

    // Test for edge and corner cases first:
    if (boost::math::isnan(z))
    {
      return policies::raise_domain_error(function,
        "Argument z is NaN!",
        z, pol);
    }
    if (boost::math::isinf(z))
    {
      if (std::numeric_limits<T>::has_infinity)
      {
        return std::numeric_limits<T>::infinity();
      }
      else
      { // Arbitrary precision type.
        return tools::max_value<T>();
      }
      // Or might opt to throw an exception?
      //return policies::raise_domain_error(function,
      //  "Argument z is infinite!",
      //  z, pol);
    }
    if (z > (std::numeric_limits<T>::max)()/4)
    { // TODO not sure this is right - what is largest possible for which can compute Lambert_w0?
      return policies::raise_domain_error(function,
        "Argument z %1% is too big!",
        z, pol);
      // Or could we return a suitable large value, computed from RealType T?
    }

    BOOST_MATH_STD_USING // Aid argument dependent lookup of abs.
    if (abs(z) <= static_cast<T>(0.05))
    { // z is near zero, so use near-zero series expansion.
      // Change to use this series approximation made by Fukushima at 0.05,
      // found by trial and error, and may not be best for higher precision.
      // Chose float, double, or long double, or multiprecision specializations.
      // Higher precision types either need more terms, or change in range of abs(z),
      // or use as an approximation?
      return detail::lambert_w0_small_z(z, pol);
    }
    else if (z == -boost::math::constants::exp_minus_one<T>())
    { // At singularity, so return exactly minus one.
#ifdef BOOST_MATH_INSTRUMENT_LAMBERT_W_SINGULARITY_SERIES
      std::cout << "At singularity. " << std::endl; // within epsilon of branch singularity 0.36.
#endif
      return static_cast<T>(-1);
    }
    else if (z < -boost::math::constants::exp_minus_one<T>()) // z < -1/e so out of range of W0 branch (should using W1 branch instead).
    {
      return policies::raise_domain_error(function,
        "Argument z = %1% out of range (-1/e <= z < (std::numeric_limits<T>::max)()) for Lambert W0 branch (use W-1 branch?).",
        z, pol);
    }
    else if (z < static_cast<T>(-0.35))
    { // Near singularity/branch point at z = -0.3678794411714423215955237701614608727
      const T p2 = 2 * (boost::math::constants::e<T>() * z + 1);
      if (p2 > 0)
      { // Use near-singularity series expansion.
         // http://www.boost.org/doc/libs/release/libs/multiprecision/doc/html/boost_multiprecision/intro.html
         // Taking care to avoid trouble from expression templates.
        T w_series = lambert_w_singularity_series(T(sqrt(p2)));
        //Schroeder does not improve result - differences about 5e-15, so do NOT use refinement step below.
        //int d2 = policies::digits<T, Policy>(); // Precision as digits base 2 from policy.
        //if (d2 > (std::numeric_limits<T>::digits - 6))
        //{ // If more (fullest) accuracy required, then use refinement.
        //  T y = z * exp(-w_series);
        //  T s_result = schroeder_update(w_series, y);
        //}
        // neither does Halley
        //T y = z * exp(-w_series);
        //return halley_update(w_series, y);
        return w_series;
      }
    } //  z < -0.35
 // z < -1/e
    else // Argument z is in the 'normal' range and float or double precision, so use Lookup, Bracket, Bisection and Schroeder (and Halley).
    {

      ///////////////////////////////////////////////////////////
      // TODO take this table out as a templated function,
      // to avoid potential multithreading initialization problems.
      // (Would like to avoid initializing table of static data?
      // but cannot avoid - always computed at load time).
      // This data is not used if z is small or z < -0.35 near branch point,
      // so this code is good at avoiding unnecessary initialisation.
      // See Fukushima section 2.2, page 81.
      // Defining integral parts of Lambert W function values to be used for bisection.
      // e and g have space for 64-bit reals.

      /*
      static const size_t noof_wm1s = 66;
      static T e[noof_wm1s]; // lambert_w[k] for W-1 branch. 2.7182 1. 0.3678 0.1353 0.04978 ... 4.359e-28 1.603e-28
      static const size_t noof_w0s = 65;
      static T g[noof_w0s]; // lambert_w[k] for W0 branch. 0 2.7182 14.77 60.2566 ... 1.445e+29 3.990e+29.
      // lambert_w0(2.77182) == 1.;
      // lambert_w0(14.7700) == 2.;
      // lambert_w0(2.77182) == 3.;
      // ...
      // lambert_w0(3.990e+29) == 64.
      // 3.990e+29 == g[64];

      // Defining lookup table sqrts of 1/e.
      static const size_t noof_sqrts = 12;
      static T a[noof_sqrts]; // 0.6065 0.7788, 0.8824 ... 0.9997, sqrt of previous elements.
      static T b[noof_sqrts]; // 0.5 0.25 0.125, 0.0625 ...  0.0002441, halves of previous elements.

      if (!e[0])
      { // Fill static data on first use.
        // Fukushima calls array F[k] = ke^k  (0 <= k <= 64 for f)
        // and G[k] -ke^-k (1 <= k <= 64 for g)

        T ej = 1; // F[k] for W-1 branch.
        e[0] = boost::math::constants::e<T>(); // e = 2.7182
        e[1] = 1;
        g[0] = 0; // G[k] for W0 branch.
        for (int j = 1, jj = 2; jj != noof_wm1s; ++jj)
        {
          e[jj] = e[j] * boost::math::constants::exp_minus_one<T>(); // For W-1 branch.
          // exp(-1) - 1/e exp_minus_one = 0.36787944.
          ej *= boost::math::constants::e<T>(); // 2.7182
          g[j] = j * ej; // For W0 branch.
          j = jj;
        } // for

        // a[0] = sqrt(e1);
        // TODO compile time constant sqrt(1/e) - add to math constants, and add to documentation too.
        // This would be slightly more accurate - conversion from long double could lose precision,
        // but might not matter if using Schroeder and/or Halley?
        // a[0] = boost::math::constants::sqrt_one_div_e<T>();
        a[0] = static_cast<T>(0.606530659712633423603799534991180453441918135487186955682L);
        b[0] = static_cast<T>(0.5); // Exactly representable.
        for (int j = 0, jj = 1; jj != noof_sqrts; ++jj)
        {
          a[jj] = sqrt(a[j]); //  Sqrt of previous element. sqrt(1/e),sqrt(sqrt(1/e)) ...
        //  b[jj] = b[j] * static_cast<T>(0.5);  // Half previous element.
          b[jj] = b[j] / 2;  // Half previous element (/2 will be optimised better?).
          j = jj;
        }
        print_collection(e, "e", " "); // 2.71828175 1 0.36787945 0.135335296 0.0497870743 0.0183156412 0.00673794793 ... 1.60381359e-28
        print_collection(g, "g", " ");
        print_collection(a, "a", " ");
        print_collection(b, "b", " ");

        // Values saved as Fukushima bisection static arrays.txt
        // These might be stored as constants in data?
        //e: 2.71828 1 0.367879 0.135335 0.0497871 0.0183156 0.00673795 0.00247875 0.000911882 0.000335463 0.00012341 4.53999e-05 1.67017e-05 6.14421e-06 2.26033e-06 8.31529e-07 3.05902e-07 1.12535e-07 4.13994e-08 1.523e-08 5.6028e-09 2.06115e-09 7.58256e-10 2.78947e-10 1.02619e-10 3.77514e-11 1.3888e-11 5.10909e-12 1.87953e-12 6.91441e-13 2.54367e-13 9.35763e-14 3.44248e-14 1.26642e-14 4.65889e-15 1.71391e-15 6.30512e-16 2.31952e-16 8.53305e-17 3.13914e-17 1.15482e-17 4.24836e-18 1.56288e-18 5.74953e-19 2.11513e-19 7.78114e-20 2.86252e-20 1.05306e-20 3.874e-21 1.42517e-21 5.24289e-22 1.92875e-22 7.09548e-23 2.61028e-23 9.60269e-24 3.53263e-24 1.29958e-24 4.7809e-25 1.75879e-25 6.47024e-26 2.38027e-26 8.75652e-27 3.22135e-27 1.18507e-27 4.35962e-28 1.60381e-28
        //g : 0 2.71828 14.7781 60.2566 218.393 742.066 2420.57 7676.43 23847.7 72927.7 220265 658615 1.95306e+06 5.75137e+06 1.68365e+07 4.90352e+07 1.42178e+08 4.10634e+08 1.18188e+09 3.39116e+09 9.7033e+09 2.76951e+10 7.8868e+10 2.2413e+11 6.35738e+11 1.80012e+12 5.08897e+12 1.43653e+13 4.04952e+13 1.14009e+14 3.20594e+14 9.00514e+14 2.52681e+15 7.08323e+15 1.98377e+16 5.55104e+16 1.55204e+17 4.33608e+17 1.21052e+18 3.37714e+18 9.4154e+18 2.62336e+19 7.30495e+19 2.03297e+20 5.6547e+20 1.57204e+21 4.36821e+21 1.21322e+22 3.36803e+22 9.34599e+22 2.59235e+23 7.18767e+23 1.99212e+24 5.51929e+24 1.5286e+25 4.23213e+25 1.17133e+26 3.24086e+26 8.96411e+26 2.47871e+27 6.85203e+27 1.89362e+28 5.23177e+28 1.44508e+29 3.99049e+29
        //a : 0.606531 0.778801 0.882497 0.939413 0.969233 0.984496 0.992218 0.996101 0.998049 0.999024 0.999512 0.999756
        //b : 0.5 0.25 0.125 0.0625 0.03125 0.015625 0.0078125 0.00390625 0.00195313 0.000976563 0.000488281 0.000244141


      }  //  if (!e[0]) end of static data fill.
      */
      ///////////////////////////////////////////////////////////////////////////

      if (std::numeric_limits<T>::digits > 53)
      { // T is more precise than 64-bit double (or long double, or ?),
        // so compute an approximate value using only one Schroeder refinement,
        // (avoiding any double-precision Halley refinement from policy double_digits2<50> 53 - 3 = 50
        // because are next going to use Halley refinement at full/high precision using this as an approximation).
        using boost::math::policies::precision;
        using boost::math::policies::digits10;
        using boost::math::policies::digits2;
        using boost::math::policies::policy;
        // Compute a 50-bit precision approximate W0 in a double (no Halley refinement).
        T double_approx = static_cast<T>(lambert_w0_imp<double>(static_cast<double>(z), policy<digits2<50> >()));
#ifdef  BOOST_MATH_INSTRUMENT_LAMBERT_W0_NOT_BUILTIN
        std::cout << "Argument Type " << typeid(T).name() << " approximation double = " << result << std::endl;
#endif // BOOST_MATH_INSTRUMENT_LAMBERT_W0
        // Perform additional Halley refinement(s) to ensure that
        // get a near as possible to correct result (usually +/- one epsilon).
        T result = halley_update(double_approx, z, pol);
#ifdef  BOOST_MATH_INSTRUMENT_LAMBERT_W0
        std::streamsize saved_precision = std::cout.precision(std::numeric_limits<T>::max_digits10);
        std::cout << "Result " << typeid(T).name() << " precision Halley refinement =    " << result << std::endl;
        std::cout.precision(saved_precision);
#endif // BOOST_MATH_INSTRUMENT_LAMBERT_W0
        return result;
      } // digits > 53  - higher precision than double.
      else // T is double or less precision.
      {
        // Use a lookup table to find the nearest integer part of Lambert W as starting point for Bisection.
        using namespace boost::math::detail::lambert_w_lookup;
        // Test sequence is n is (0, 1, 2, 4, 8, 16, 32, 64) for W0 branch.
        // Since z is probably quite small, start with lowest (n = 0).
        int n; // Indexing W0 lookup table w0s[n] of z.
        for (n = 0; n <= 2; ++n)
        { // Try 1st few.
          if (w0s[n] > z)
          { //
            goto bisect;
          }
        }
        n = 2;
        for (int j = 1; j <= 5; ++j)
        {
          n *= 2; // Try bigger steps.
          if (w0s[n] > z) // z <= w0s[n]
          {
            goto overshot; //
          }
        } // for
        // get here if
        //std::cout << "Too big " << z << std::endl;
        //// else z too large :-(
        //const char* function = "boost::math::lambert_w0<RealType>(<RealType>)";
        return policies::raise_domain_error("boost::math::lambert_w0<RealType>(<RealType>)",
          "Argument z = %1% is too large!",
          z, pol);

      overshot:
        {
          int nh = n / 2;
          for (int j = 1; j <= 5; ++j)
          {
            nh /= 2;
            if (nh <= 0)
            {
              break;
            }
            if (w0s[n - nh] > z)
            {
              n -= nh;
            }
          } // for
        }

      bisect:
        --n; // w0s[n] is nearest below, so n is integer part of W,
        // and w0s[n+1] is next integral part of W.

        // These are used as initial values for bisection.

        // typedef typename policies::precision<T, Policy>::type prec;
        //std::cout << "digits2 = " << digits2 << std::endl;  //  digits = 24 for float.
        // Check precision specified by policy.
        //int digits2 = policies::digits<T, Policy>();
        //std::cout << "digits2 = " << digits2 << std::endl;  //  digits = 24 for float.

        // for example = 1 from lambert_w(x, policy<digits10<1> >()
        int d2 = policies::digits<T, Policy>(); // digits base 2 from policy.
#ifdef BOOST_MATH_INSTRUMENT_LAMBERT_W0
        int d10 = policies::digits_base10<T, Policy>(); // policy template parameter digits10
        std::cout << "digits10 = " << d10 << ", digits2 = " << d2 // For example: digits10 = 1, digits2 = 5
          << std::endl;
#endif // BOOST_MATH_INSTRUMENT_LAMBERT_W0

#ifdef BOOST_MATH_INSTRUMENT_LAMBERT_W0_LOOKUP
          std::cout << "Result lookup W(" << z <<  ") bisection between w0s[" << n - 1 << "] = " << w0s[n - 1] << " and w0s[" << n << "] = " << w0s[n]
            << ", bisect mean = " << (w0s[n - 1] + w0s[n]) / 2 << std::endl;
#endif // BOOST_MATH_INSTRUMENT_LAMBERT_W0_LOOKUP

        if (d2 <= 7)
        { // Only 7 binary digits precision required to hold integer size of w0s[65],
          // so just return the mean of two nearby integral values.
          // This is a *very* approximate estimate, and perhaps not very useful?
          return static_cast<T>((w0s[n - 1] + w0s[n]) / 2);
        }
        // Compute the number of bisections (jmax) that ensure that result is close enough that
        // a single 5th order Schroeder update is sufficient to achieve near double (53-bit) accuracy.
        int jmax = 8; //  Assume minimum number of bisections will be needed (most common case).
        if (z <= -0.36)
        { // Very near singularity, so need several more bisections.
          jmax = noof_sqrts;
        }
        else if (z <= -0.3)
        { // Near singularity, so need a few more bisections.
          jmax = 11;
        }
        else if (n <= 0)
        { // Very small z, so need 3 more bisections.
          jmax = 10;
        }
        else if (n <= 1)
        { // Small z, so need just 1 extra bisection.
          jmax = 9;
        }
        lookup_t dz = static_cast<lookup_t>(z); 
        // Only use lookup_t precision, default double, for bisection.
        // Use Halley refinement for higher precisions. 
        // Perform the jmax fractional bisections for necessary precision.
        // Avoid using exponential function for speed.
        lookup_t w = static_cast<lookup_t>(n);  // Equation 25, Integral Lambert W.
       // lookup_t y = dz * e[n + 1]; // Integral +1 Lambert W.
        lookup_t y = dz * static_cast<lookup_t>(wm1s[n + 1]); // Equation 26, Integral +1 Lambert W.
        for (int j = 0; j < jmax; ++j)
        { // Equation 27.
          const lookup_t wj = w + static_cast<lookup_t>(halves[j]); // Add 1/2, 1/4, 1/8 ...
          const lookup_t yj = y * static_cast<lookup_t>(sqrts[j]); // Multiply by sqrt(1/e), ...
          if (wj < yj)
          {
            w = wj;
            y = yj;
          }
        } // for jmax

        if (d2 <= 10)
        { // Only 10 digits2 required (== digits10 = 3 decimal digits precision),
          // so just return the nearest bisection.
          // (Might make this test depend on size of T?)
#ifdef BOOST_MATH_INSTRUMENT_LAMBERT_W0_BISECTION
          std::cout << "Bisection estimate =            " << w << " " << y << std::endl;
#endif // BOOST_MATH_INSTRUMENT_LAMBERT_W0_BISECTION
          return static_cast<T>(w); // Bisection only.
        }
        else // More than 10 digits2 wanted, so continue with Fukushima's Schroeder refinement.
        {
          T result = static_cast<T>(schroeder_update(w, y));
          if (d2 <= (std::numeric_limits<T>::digits - 3))
          { // Only enough digits2 required, so just return the Schroeder refined value.
            // For example, for float, if (d2 <= 22) then
            // digits-3 returns Schroeder result if up to 3 bits might be wrong.
#ifdef BOOST_MATH_INSTRUMENT_LAMBERT_W0
            std::cout << "Schroeder refinement estimate = " << result << std::endl;
#endif // BOOST_MATH_INSTRUMENT_LAMBERT_W0
            return result; // Schroeder only.
          }
          else // Perform additional Halley refinement(s) to ensure that
          // get a near as possible to correct result (usually +/- epsilon).
          {
            result = halley_update(result, z, pol);
#ifdef BOOST_MATH_INSTRUMENT_LAMBERT_W0_HALLEY
            std::cout.precision(std::numeric_limits<T>::max_digits10);
            std::cout << "Halley refinement estimate =    " << result << std::endl;
            std::cout.precision(saved_precision);
#endif // BOOST_MATH_INSTRUMENT_LAMBERT_W0_HALLEY
            return result; // Halley
          } // schroeder or Schroeder and Halley.
        } // more than 10 digits2
      }
    } // normal range and precision.
    throw std::logic_error("No result from lambert_w0_imp");  // Should never get here.
  } //  T lambert_w0_imp(const T z, const Policy& /* pol */)

  //! Lambert W for W-1 branch, -max(z) < z <= -1/e.
  // TODO is -max(z) allowed?
  // TODO check changes to W0 branch have also been made here.
  template<typename T, class Policy>
  T lambert_wm1_imp(const T z, const Policy&  pol)
  {
    // Catch providing an integer value as parameter x to lambert_w, for example, lambert_w(1).
    // Need to ensure it is a floating-point type (of the desired type, float 1.F, double 1., or long double 1.L),
    // or static_casted integer, for example:  static_cast<float>(1) or static_cast<cpp_dec_float_50>(1).
    // Want to allow fixed_point types too, so do not just test for floating-point.
    // Integral types should be promoted to double by user Lambert w functions.
    // If integral type provided to user function lambert_w0 or lambert_wm1,
    // then should already have been promoted to double.
    BOOST_STATIC_ASSERT_MSG(!std::is_integral<T>::value, "Must be floating-point or fixed type (not integer type), for example: lambert_wm1(1.), not lambert_wm1(1)!");

    BOOST_MATH_STD_USING // Aid argument dependent lookup (ADL) of abs.

    using boost::math::tools::max_value;

    const char* function = "boost::math::lambert_wm1<RealType>(<RealType>)"; // Used for error messages.

    if (z == static_cast<T>(0))
    { // z is exactly zero.return -std::numeric_limits<T>::infinity();
      return -std::numeric_limits<T>::infinity();
    } 
    if (z > static_cast<T>(0))
    { // 
      return policies::raise_domain_error(function,
        "Argument z = %1% is out of range (z > 0) for Lambert W-1 branch! (Try Lambert W0 branch?)",
        z, pol);
    }
    if(z > -(std::numeric_limits<T>::min)())
    { // z is denormalized, so cannot be computed. 
      // -std::numeric_limits<T>::min() is smallest for type T, 
      // for example, for double: lambert_wm1(-2.2250738585072014e-308) = -714.96865723796634
      return policies::raise_domain_error(function,
        "Argument z = %1% is too small (z < -std::numeric_limits<T>::min so denormalized) for Lambert W-1 branch!",
        z, pol);
    }
    if (z == -boost::math::constants::exp_minus_one<T>()) // == singularity/branch point z = -exp(-1) = -3.6787944.
    { // At singularity, so return exactly -1.
      return -static_cast<T>(1);
    }
    // z is too negative for the W-1 (or W0) branch.
    if (z < -boost::math::constants::exp_minus_one<T>()) // > singularity/branch point z = -exp(-1) = -3.6787944.
    {
      return policies::raise_domain_error(function,
        "Argument z = %1% is out of range (z < -exp(-1) = -3.6787944... <= 0) for Lambert W-1 (or W0) branch!",
        z, pol);
    }
    if (z < -0.35)
    { // Close to singularity/branch point but on w-1 branch.  TODO.
      const T p2 = 2 * (boost::math::constants::e<T>() * z + 1);
      if (p2 > 0)
      {
        return lambert_w_singularity_series(-sqrt(p2));
        // TODO Halley refinement options here.
      }
      if (p2 == 0)
      { // At the singularity at branch point.
        return -1;
      }
      return policies::raise_domain_error(function,
        "Argument z = %1% is out of range!",
        z, pol);
    } // if (z < -0.35)

    // Create static arrays of Lambert Wm1 function
    // TODO these are NOT the same size or values as the w0 constant array values?????
    
    // double LambertW0(const double z)
    // static double e[66];
    // static double g[65];
    // static double a[12];
    // static double b[12];

    // Create static arrays of Lambert W function for Wm1
    // with lambert_w values of integers.
    // Used for bisection of W-1 branch.
    static const size_t noof_sqrts = 12; // and halves
    // F[k] is 0 <= k <= 64 so size is 65 (but is 64 below????).
    // G[k] is 1 <= k <= 64 so size is 64
    // So e and F are not the same!

    static const size_t noof_wm1s = 64;

    static lookup_t e[noof_wm1s]; // e[0] to e[63] 
    static lookup_t g[noof_wm1s]; // g[0] == Gk when k = 1,  to g[63] == Gk=64
    static lookup_t a[noof_sqrts];
    static lookup_t b[noof_sqrts];

    if (!e[0])
    {
      const lookup_t e1 = 1 / boost::math::constants::e<lookup_t>();
      lookup_t ej = e1; // 1/e= 0.36787945 
      e[0] = boost::math::constants::e<lookup_t>();
      g[0] = -e1;
      for (int j = 0, jj = 1; jj < 64; ++jj)
      { // jj is j+1
        ej *= e1;
        e[jj] = e[j] * boost::math::constants::e<lookup_t>();
        g[jj] = -(jj + 1) * ej;
        j = jj;
      }
      a[0] = sqrt(boost::math::constants::e<lookup_t>());
      b[0] = static_cast<lookup_t>(0.5);
      for (int j = 0, jj = 1; jj < noof_sqrts; ++jj)
      {
        a[jj] = sqrt(a[j]); 
     //   b[jj] = b[j] * static_cast<lookup_t>(0.5);
        b[jj] = b[j] / 2; // More efficient for multiprecision?
        j = jj;
      }  // for j


//std::cout << " W-1 e[0] = " << e[0] << ", e[1] = " << e[1] << "... e[63] = " << e[63]  << std::endl;
// W-1 e[0] = 2.71828175, e[1] = 7.38905573... e[63] = 6.23513822e+27
// W-1 e[0] = 2.7182818284590451, e[1] = 7.3890560989306495... e[63] = 6.235149080811597e+27, e[64] = 58265.164084420219
//std::cout << " W-1 g[0] = " << g[0] << ", g[1] = " << g[1] << "... g[63] = " << g[63]  << std::endl;
// W-1 g[0] = -0.36787945 k=1, -1e^-1 = -0.3678794 ,
//  g[1] = -0.270670593  ... g[63] = -1.0264407e-26
// W-1 g[0] = -0.36787944117144233, g[1] = -0.2706705664732254... g[63] = -1.0264389699511303e-26
      // TODO temporary output.
      //std::cout << "lambert_wm1 version of arrays - NOT same as w0 version." << std::endl;
     // print_collection(e, "e", " "); // e[0] = 2.71828175, e[1] = 1, e[2] = 0.36787945, 0.135335296 0.0497870743 0.0183156412 0.00673794793 ... e[63] = 2.293783159, e[64]= 6.235149080e+27
     // print_collection(g, "g", " "); // g[0] = -0.3678794, g[1] = -0.2706705...  , g[63] = -1.026438e-26
      //print_collection(a, "a", " ");
      //print_collection(b, "b", " ");
    } // if (!e[0])
// static arrays filled.

    // TODO move up when using precomputed array.
    //if (abs(z) < abs(g[63]))

if (z > g[63])
    {  // z > -1.0264389699511303e-26 (but != 0 and >= std::numeric_limits<T>::min() and so NOT denormalized).
      std::streamsize saved_precision = std::cout.precision(std::numeric_limits<T>::max_digits10);

     // std::cout << "abs(z) " << abs(z) << ",  abs(g[63] = " << abs(g[63]) << std::endl;

     // std::cout << "-std::numeric_limits<float>::min() = " << -(std::numeric_limits<float>::min)() << std::endl;
     // std::cout << "-std::numeric_limits<double>::min() = " << -(std::numeric_limits<double>::min)() << std::endl;
      // -std::numeric_limits<float>::min() = -1.1754943508222875e-38
      // -std::numeric_limits<double>::min() = -2.2250738585072014e-308

      // N[productlog(-1, -1.1754943508222875 * 10^-38 ), 50] = -91.856775324595479509567756730093823993834155027858

      // N[productlog(-1, -2.2250738585072014e-308 * 10^-308 ), 50] = -1424.8544521230553853558132180518404363617968042942

      T guess ; // bisect lowest possible Gk[=64] (for lookup_t type?)

      //int exponent = ilogb(z);
      //std::cout << exponent << std::endl;
      //guess = static_cast<T>(exponent);

     // R.M.Corless, G.H.Gonnet, D.E.G.Hare, D.J.Jeffrey, and D.E.Knuth, “On the Lambert W function, ” Adv.Comput.Math., vol. 5, pp. 329–359, 1996.
      // François Chapeau-Blondeau and Abdelilah Monir
      // Numerical Evaluation of the Lambert W Function
      // IEEE Transactions On Signal Processing, VOL. 50, NO. 9, Sep 2002
      // https://pdfs.semanticscholar.org/7a5a/76a9369586dd0dd34dda156d8f2779d1fd59.pdf
      // Estimate Lambert W using ln(-z)  ...
// This is roughly the power of ten * ln(10) ~= 2.3.   n ~= 10^n  
//  and improve by adding a second term -ln(ln(-z))
      T lz = log(-z);
      T llz = log(-lz);
      guess = lz - llz + (llz/lz); // Chapeau-Blondeau equation 20, page 2162.
#ifdef BOOST_MATH_INSTRUMENT_LAMBERT_WM1_TINY
      std::cout << "z = " << z << ", guess = " << guess << ", ln(-z) = " << lz << ", ln(-ln(-z) = " << llz  << ", llz/lz = " << (llz / lz) << std::endl;
      // z = -1.0000000000000001e-30, guess = -73.312782616731482, ln(-z) = -69.077552789821368, ln(-ln(-z) = 4.2352298269101114, llz/lz = -0.061311231447304194
      // z = -9.9999999999999999e-91, guess = -212.56650048504233, ln(-z) = -207.23265836946410, ln(-ln(-z) = 5.3338421155782205, llz/lz = -0.025738424423764311
      // >z = -2.2250738585072014e-308, guess = -714.95942238244606, ln(-z) = -708.39641853226408, ln(-ln(-z) = 6.5630038501819854, llz/lz = -0.0092645920821846622
      int d10 = policies::digits_base10<T, Policy>(); // policy template parameter digits10
      int d2 = policies::digits<T, Policy>(); // digits base 2 from policy.
      std::cout << "digits10 = " << d10 << ", digits2 = " << d2 // For example: digits10 = 1, digits2 = 5
        << std::endl;
#endif // BOOST_MATH_INSTRUMENT_LAMBERT_WM1_TINY
      if (policies::digits<T, Policy>() < 12)
      { // For the worst case near w = 64, the error in the 'guess' is ~0.008, ratio ~ 0.0001 or 1 in 10,000 digits 10 ~= 4, or digits2 ~= 12.
        return guess;
      }
      T result = halley_update(guess, z, pol);
      std::cout.precision(saved_precision);
      return result;

      //// G[k=64] == g[63] == -1.02643897e-26
      //return policies::raise_domain_error(function,
      //  "Argument z = %1% is too small (< -1.02643897e-26) !",
      //  z, pol);
    }



    // Use a lookup table to find the nearest integer part of Lambert W as starting point for Bisection.
    // Bracketing sequence  n = (2, 4, 8, 16, 32, 64) for W-1 branch.
    // Since z is probably quite small, start with lowest n (=2).
    int n = 2;
    if (g[n - 1] > z)
    {
      goto bisect;
    }
    for (int j = 1; j <= 5; ++j)
    {
      n *= 2;
      if (g[n - 1] > z)
      {
        goto overshot;
      }
    }
    // else z < g[63] == -1.0264389699511303e-26, so Lambert W integer part > 64.

 // Might assign the lower bracket to -1.0264389699511303e-26 and upper to (std::numeric_limits<T>::min)()
// and jump to schroeder to halley? But that would only allow use with lookup_t precision.
    // So try using Halley direct?
    return policies::raise_domain_error(function, 
      "Argument z = %1% is too small (< -1.026439e-26) !",
      z, pol);

  overshot:
    {
      int nh = n / 2; 
      for (int j = 1; j <= 5; ++j)
      {
        nh /= 2; // halve step size.
        if (nh <= 0)
        {
          break; // goto bisect;
        }
        if (g[n - nh - 1] > z)
        {
          n -= nh;
        }
      }
    }
  bisect:
    --n;   // g[n] now holds lambert W of floor integer n and g[n+1] the ceil part.

    // Check precision specified by policy.
    int d2 = policies::digits<T, Policy>();
#ifdef BOOST_MATH_INSTRUMENT_LAMBERT_W1
    std::cout << "digits2 = " << d2 << std::endl;  //  digits = 24 for float.
#endif // BOOST_MATH_INSTRUMENT_LAMBERT_W1
    if (d2 < 5)
    { // Only a very few digits2 required, so just return the floor integral part of Lambert_w.
#ifdef BOOST_MATH_INSTRUMENT_LAMBERT_W1
      std::cout << "between g[" << n << "] = " << g[n] << " < " << z << " < " << g[n + 1] << std::endl;
#endif // BOOST_MATH_INSTRUMENT_LAMBERT_W1
      return static_cast<T>(n);
    }

    // jmax is the number of bisections computed from n,
    // such that a single application of the fifth-order Schroeder update formula
    // after the bisections is enough to evaluate Lambert W-1 with 53-bit accuracy.
    int jmax = 11; //  Assume maximum number of bisections will be needed (most common case).
    if (n >= 8)
    {
      jmax = 8;
    }
    else if (n >= 3)
    {
      jmax = 9;
    }
    else if (n >= 2)
    {
      jmax = 10;
    }
    // Bracketing, Fukushima section 2.3, page 82:
    // (Avoid using exponential function for speed).
    // Only use lookup_t precision, default double, for bisection (for speed).
    lookup_t w = -static_cast<lookup_t>(n); // Equation 25,
    lookup_t y = static_cast<lookup_t>(z * e[n - 1]); // Equation 26,
    for (int j = 0; j < jmax; ++j)
    { // Equation 27.
      lookup_t wj = w - b[j]; // Subtract 1/2, 1/4, 1/8 ...
      lookup_t yj = y * a[j]; // Multiply by sqrt(1/e), ...
      if (wj < yj)
      {
        w = wj;
        y = yj;
      }
    } // for j

    T result = static_cast<T>(schroeder_update(w, y)); // Schroeder 5th order method refinement.
    // Only use Halley refinement (using exp) for higher precisions. 
    result = halley_update(result, z, pol);
    return result;
  } // template<typename T = double> T lambert_wm1_imp(const T z)

} // namespace detail ////////////////////////////////////////////////////////////

  // User Lambert W0 and Lambert W-1 functions.

  //! W0 branch, -1/e <= z < max(z)
  //! W-1 branch -1/e >= z > min(z)

  //! Lambert W0 using User-defined policy.
  template <class T, class Policy>
  inline
  typename tools::promote_args<T>::type
  lambert_w0(T z, const Policy& pol)
  {
    // Promote integer or expression template arguments to double,
    // without doing any other internal promotion like float to double.
    typedef typename tools::promote_args<T>::type result_type;
    return detail::lambert_w0_imp(result_type(z), pol); //
  } // lambert_w0(T z, const Policy& pol)

  //! Lambert W0 using default policy.
  template <class T>
  inline
  typename tools::promote_args<T>::type
  lambert_w0(T z)
  {
    typedef typename tools::promote_args<T>::type result_type;
    return detail::lambert_w0_imp(result_type(z), policies::policy<>());
  } // lambert_w0(T z)


  //! W-1 branch (-max(z) < z <= -1/e).

  //! Lambert W-1 using User-defined policy.
  template <class T, class Policy>
  inline
  typename tools::promote_args<T>::type
  lambert_wm1(T z, const Policy& pol)
  {
    // Promote integer or expression template arguments to double,
    // without doing any other internal promotion like float to double.
    typedef typename tools::promote_args<T>::type result_type;
    return detail::lambert_wm1_imp(result_type(z), pol); //
  }

  //! Lambert W-1 using default policy.
  template <class T>
  inline
  typename tools::promote_args<T>::type
  lambert_wm1(T z)
  {
    typedef typename tools::promote_args<T>::type result_type;
    return detail::lambert_wm1_imp(result_type(z), policies::policy<>());
  } // lambert_wm1(T z)

  } // namespace math
} // namespace boost

#endif //  #ifndef BOOST_MATH_SF_LAMBERT_W_HPP
