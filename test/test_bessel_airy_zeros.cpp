//  Copyright John Maddock 2013
//  Copyright Christopher Kormanyos 2013.
//  Copyright Paul A. Bristow 2013.

//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifdef _MSC_VER
#  pragma warning(disable : 4127) // conditional expression is constant.
#endif

//#include <pch_light.hpp> // commente dout during testing.

#include <boost/math/special_functions/bessel.hpp>

#include <boost/math/concepts/real_concept.hpp> // for real_concept
#include <boost/test/test_exec_monitor.hpp> // Boost.Test
#include <boost/test/floating_point_comparison.hpp>
//
// DESCRIPTION:
// ~~~~~~~~~~~~
//
// This file tests the functions that evaluate zeros (or roots) of Bessel and Airy functions.

// Spot tests which compare our results with selected values computed 
// using the online special function calculator at functions.wolfram.com,
// and values generated with Boost.Multiprecision at about 1000-bit or 100 decimal digits precision.

// Weisstein, Eric W. "Bessel Function Zeros." From MathWorld--A Wolfram Web Resource.
// http://mathworld.wolfram.com/BesselFunctionZeros.html 

// See also http://dlmf.nist.gov/10.21

template <class RealType>
void test_bessel_zeros(RealType)
{
  // Basic sanity checks for finding zeros of Bessel and Airy function.
  // where template parameter RealType can be float, double, long double,
  // or real_concept, a prototype for user-defined floating-point types.

  // Parameter RealType is only used to communicate the RealType, float, double...
  // and is an arbitrary zero for all tests.
   RealType tolerance = (std::max)(
     static_cast<RealType>(boost::math::tools::epsilon<long double>()),
     boost::math::tools::epsilon<RealType>());
   std::cout << "Tolerance for type " << typeid(RealType).name()  << " is " << tolerance << "." << std::endl;

   // http://www.wolframalpha.com/
/*
Table[N[BesselJZero[0, n], 50], {n, 1, 5, 1}]
n | 
1 | 2.4048255576957727686216318793264546431242449091460
2 | 5.5200781102863106495966041128130274252218654787829
3 | 8.6537279129110122169541987126609466855657952312754
4 | 11.791534439014281613743044911925458922022924699695
5 | 14.930917708487785947762593997388682207915850115633
  
Table[N[BesselJZero[1, n], 50], {n, 1, 4, 1}]
n | 
1 | 3.8317059702075123156144358863081607665645452742878
2 | 7.0155866698156187535370499814765247432763115029142
3 | 10.173468135062722077185711776775844069819512500192
4 | 13.323691936314223032393684126947876751216644731358

Table[N[BesselJZero[5, n], 50], {n, 1, 5, 1}]
n | 
1 | 8.7714838159599540191228671334095605629810770148974
2 | 12.338604197466943986082097644459004412683491122239
3 | 15.700174079711671037587715595026422501346662246893
4 | 18.980133875179921120770736748466932306588828411497
5 | 22.217799896561267868824764947529187163096116704354
*/
   using boost::math::cyl_bessel_j_zero; // (nu, j) 
   using boost::math::isnan;

  if (std::numeric_limits<RealType>::has_quiet_NaN)
  {
    BOOST_CHECK(isnan(cyl_bessel_j_zero(static_cast<RealType>(0), 0U))); // yes - returns NaN - is this right?
  }
  BOOST_CHECK_CLOSE_FRACTION(cyl_bessel_j_zero(static_cast<RealType>(0), 1U), static_cast<RealType>(2.4048255576957727686216318793264546431242449091460L), tolerance);
  BOOST_CHECK_CLOSE_FRACTION(cyl_bessel_j_zero(static_cast<RealType>(1), 1U), static_cast<RealType>(3.8317059702075123156144358863081607665645452742878L), tolerance);
  BOOST_CHECK_CLOSE_FRACTION(cyl_bessel_j_zero(static_cast<RealType>(1), 2U), static_cast<RealType>(7.0155866698156187535370499814765247432763115029142L), tolerance);
  BOOST_CHECK_CLOSE_FRACTION(cyl_bessel_j_zero(static_cast<RealType>(5), 5U), static_cast<RealType>(22.217799896561267868824764947529187163096116704354L), tolerance);

  // Some none integral tests.
  BOOST_CHECK_CLOSE_FRACTION(cyl_bessel_j_zero(static_cast<RealType>(3.736842105263157894736842105263157894736842105263157894736842105263157894736842105263157894736842105L), 1U), static_cast<RealType>(7.273175193831648950318569426229076558896319670162279791988152000556091140599946365217211157877052381L), tolerance);
  BOOST_CHECK_CLOSE_FRACTION(cyl_bessel_j_zero(static_cast<RealType>(3.736842105263157894736842105263157894736842105263157894736842105263157894736842105263157894736842105L), 20U), static_cast<RealType>(67.81514561969629092555679137555595116511146058545787883557679231060644931096494584364894743334132014L), tolerance);

  // Some none integral tests in 'tough' regions.
  BOOST_CHECK_CLOSE_FRACTION(cyl_bessel_j_zero(static_cast<RealType>(219)/100, 1U), static_cast<RealType>(5.37568854370623186731066365697341253761466705063679L), tolerance);
  BOOST_CHECK_CLOSE_FRACTION(cyl_bessel_j_zero(static_cast<RealType>(219)/100, 2U), static_cast<RealType>(8.67632060963888122764226633146460596009874991130394L), tolerance);
  BOOST_CHECK_CLOSE_FRACTION(cyl_bessel_j_zero(static_cast<RealType>(219)/100, 20U), static_cast<RealType>(65.4517712237598926858973399895944886397152223643028L), tolerance);

  BOOST_CHECK_CLOSE_FRACTION(cyl_bessel_j_zero(static_cast<RealType>(221)/100, 1U), static_cast<RealType>(5.40084731984998184087380740054933778965260387203942L), tolerance);
  BOOST_CHECK_CLOSE_FRACTION(cyl_bessel_j_zero(static_cast<RealType>(221)/100, 2U), static_cast<RealType>(8.70347906513509618445695740167369153761310106851599L), tolerance);
  BOOST_CHECK_CLOSE_FRACTION(cyl_bessel_j_zero(static_cast<RealType>(221)/100, 20U), static_cast<RealType>(65.4825314862621271716158606625527548818843845600782L), tolerance);

  BOOST_CHECK_CLOSE_FRACTION(cyl_bessel_j_zero(static_cast<RealType>(7001)/19, 1U), static_cast<RealType>(381.922015230244893869172044704348426991540311353476L), tolerance);
  BOOST_CHECK_CLOSE_FRACTION(cyl_bessel_j_zero(static_cast<RealType>(7001)/19, 2U), static_cast<RealType>(392.175086576487375026512998530998525670012392177242L), tolerance);
  BOOST_CHECK_CLOSE_FRACTION(cyl_bessel_j_zero(static_cast<RealType>(7001)/19, 20U), static_cast<RealType>(496.394350379382525575353754985779897202722983108025L), tolerance);

  // Check expected exceptions :
  //BOOST_CHECK_THROW(boost::math::cyl_bessel_j_zero(static_cast<RealType>(0), 0L), std::domain_error);

  //BOOST_CHECK_THROW(boost::math::cyl_bessel_j_zero(static_cast<RealType>(0), -1), std::domain_error);
  //  warning C4245: 'argument' : conversion from 'int' to 'unsigned int', signed/unsigned mismatch

 // BOOST_CHECK_THROW(cyl_bessel_j_zero(static_cast<RealType>(std::numeric_limits<RealType>::quiet_NaN()), 0U), std::domain_error); // No exception?
  //BOOST_CHECK_THROW(boost::math::cyl_bessel_j_zero(static_cast<RealType>(std::numeric_limits<RealType>::quiet_NaN()), -1), std::domain_error); // OK if unsigned.
  if (std::numeric_limits<RealType>::has_quiet_NaN)
  {
    BOOST_CHECK_THROW(cyl_bessel_j_zero(static_cast<RealType>(std::numeric_limits<RealType>::quiet_NaN()), 1U), std::domain_error);
   }
 // BOOST_CHECK_THROW(cyl_bessel_j_zero(static_cast<RealType>(std::numeric_limits<RealType>::infinity()), 0U), std::domain_error);
  //BOOST_CHECK_THROW(boost::math::cyl_bessel_j_zero(static_cast<RealType>(std::numeric_limits<RealType>::infinity()), -1), std::domain_error); // OK if unsigned.
  if (std::numeric_limits<RealType>::has_infinity)
  {
     BOOST_CHECK_THROW(cyl_bessel_j_zero(static_cast<RealType>(std::numeric_limits<RealType>::infinity()), 1U), std::domain_error);
  }

  // BOOST_CHECK_THROW(static_cast<RealType>(0.L), 0L, std::domain_error);

} // template <class RealType> void test_spots(RealType)

int test_main(int, char* [])
{
#ifdef TEST_GSL
   gsl_set_error_handler_off();
#endif
   //expected_results();
   BOOST_MATH_CONTROL_FP;

   test_bessel_zeros(0.1F);
   test_bessel_zeros(0.1);
#ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
   test_bessel_zeros(0.1L);
#ifndef BOOST_MATH_NO_REAL_CONCEPT_TESTS
   test_bessel_zeros(boost::math::concepts::real_concept(0.1));
#endif
#else
   std::cout << "<note>The long double tests have been disabled on this platform "
      "either because the long double overloads of the usual math functions are "
      "not available at all, or because they are too inaccurate for these tests "
      "to pass.</note>" << std::cout;
#endif
   return 0;
} // int test_main(int, char* [])




