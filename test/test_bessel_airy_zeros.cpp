//  Copyright John Maddock 2013
//  Copyright Christopher Kormanyos 2013.
//  Copyright Paul A. Bristow 2013.

//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifdef _MSC_VER
#  pragma warning(disable : 4127) // conditional expression is constant.
#  pragma warning(disable : 4512) // assignment operator could not be generated.
#endif

//#include <pch_light.hpp> // commente dout during testing.

#include <boost/math/special_functions/math_fwd.hpp>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/airy.hpp>

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



// Tests of Airy zeros.


/*  The algorithms use tabulated values for the first 10 zeros,
whereby algorithms are used for rank 11 and higher.
So testing the zeros of Ai and Bi from 1 through 20 handles
this cross-over nicely.

In addition, the algorithms for the estimates of the zeros
become increasingly accurate for larger, negative argument.

On the other hand, the zeros become increasingly close
for large, negative argument. So another nice test
involves testing pairs of zeros for different orders of
magnitude of the zeros, to insure that the program
properly resolves very closely spaced zeros.
*/

  // Test Data for airy_ai
   using boost::math::airy_ai_zero; // 

   using boost::math::isnan;

  if (std::numeric_limits<RealType>::has_quiet_NaN)
  {
    BOOST_CHECK(isnan(airy_ai_zero<RealType>(0)) );
  }

  // WolframAlpha  Table[N[AiryAiZero[n], 51], {n, 1, 20, 1}]

  BOOST_CHECK_CLOSE_FRACTION(airy_ai_zero<RealType>(1U), static_cast<RealType>(-2.33810741045976703848919725244673544063854014567239L), tolerance);
  BOOST_CHECK_CLOSE_FRACTION(airy_ai_zero<RealType>(2U), static_cast<RealType>(-4.08794944413097061663698870145739106022476469910853L), tolerance);
  BOOST_CHECK_CLOSE_FRACTION(airy_ai_zero<RealType>(3U), static_cast<RealType>(-5.52055982809555105912985551293129357379721428061753L), tolerance);
  BOOST_CHECK_CLOSE_FRACTION(airy_ai_zero<RealType>(4U), static_cast<RealType>(-6.78670809007175899878024638449617696605388247739349L), tolerance);
  BOOST_CHECK_CLOSE_FRACTION(airy_ai_zero<RealType>(5U), static_cast<RealType>(-7.94413358712085312313828055579826853214067439697221L), tolerance);
  BOOST_CHECK_CLOSE_FRACTION(airy_ai_zero<RealType>(6U), static_cast<RealType>(-9.02265085334098038015819083988008925652467753515608L), tolerance);
  BOOST_CHECK_CLOSE_FRACTION(airy_ai_zero<RealType>(7U), static_cast<RealType>(-10.0401743415580859305945567373625180940429025691058L), tolerance);
  BOOST_CHECK_CLOSE_FRACTION(airy_ai_zero<RealType>(8U), static_cast<RealType>(-11.0085243037332628932354396495901510167308253815040L), tolerance);
  BOOST_CHECK_CLOSE_FRACTION(airy_ai_zero<RealType>(9U), static_cast<RealType>(-11.9360155632362625170063649029305843155778862321198L), tolerance);
  BOOST_CHECK_CLOSE_FRACTION(airy_ai_zero<RealType>(10U), static_cast<RealType>(-12.8287767528657572004067294072418244773864155995734L), tolerance);
  BOOST_CHECK_CLOSE_FRACTION(airy_ai_zero<RealType>(11U), static_cast<RealType>(-13.6914890352107179282956967794669205416653698092008L), tolerance);
  BOOST_CHECK_CLOSE_FRACTION(airy_ai_zero<RealType>(12U), static_cast<RealType>(-14.5278299517753349820739814429958933787141648698348L), tolerance);
  BOOST_CHECK_CLOSE_FRACTION(airy_ai_zero<RealType>(13U), static_cast<RealType>(-15.3407551359779968571462085134814867051175833202480L), tolerance);
  BOOST_CHECK_CLOSE_FRACTION(airy_ai_zero<RealType>(14U), static_cast<RealType>(-16.1326851569457714393459804472025217905182723970763L), tolerance);
  BOOST_CHECK_CLOSE_FRACTION(airy_ai_zero<RealType>(15U), static_cast<RealType>(-16.9056339974299426270352387706114765990900510950317L), tolerance);
  BOOST_CHECK_CLOSE_FRACTION(airy_ai_zero<RealType>(16U), static_cast<RealType>(-17.6613001056970575092536503040180559521532186681200L), tolerance);
  BOOST_CHECK_CLOSE_FRACTION(airy_ai_zero<RealType>(17U), static_cast<RealType>(-18.4011325992071154158613979295043367545938146060201L), tolerance);
  BOOST_CHECK_CLOSE_FRACTION(airy_ai_zero<RealType>(18U), static_cast<RealType>(-19.1263804742469521441241486897324946890754583847531L), tolerance);
  BOOST_CHECK_CLOSE_FRACTION(airy_ai_zero<RealType>(19U), static_cast<RealType>(-19.8381298917214997009475636160114041983356824945389L), tolerance);
  BOOST_CHECK_CLOSE_FRACTION(airy_ai_zero<RealType>(20U), static_cast<RealType>(-20.5373329076775663599826814113081017453042180147375L), tolerance);

  // Table[N[AiryAiZero[n], 51], {n, 1000, 1001, 1}]

  BOOST_CHECK_CLOSE_FRACTION(airy_ai_zero<RealType>(1000U), static_cast<RealType>(-281.031519612521552835336363963709689055717463965420L), tolerance);
  BOOST_CHECK_CLOSE_FRACTION(airy_ai_zero<RealType>(1001U), static_cast<RealType>(-281.218889579130068414512015874511112547569713693446L), tolerance);

  // Table[N[AiryAiZero[n], 51], {n, 1000000, 1000001, 1}]
  BOOST_CHECK_CLOSE_FRACTION(airy_ai_zero<RealType>(1000000U), static_cast<RealType>(-28107.8319793795834876064419863203282898723750036048L), tolerance);
  BOOST_CHECK_CLOSE_FRACTION(airy_ai_zero<RealType>(1000001U), static_cast<RealType>(-28107.8507179357979542838020057465277368471496446555L), tolerance);


  // Table[N[AiryAiZero[n], 51], {n, 1000000000, 1000000001, 1}]
  BOOST_CHECK_CLOSE_FRACTION(airy_ai_zero<RealType>(1000000000U), static_cast<RealType>(-2.81078366593344513918947921096193426320298300481145E+6L), tolerance);
  BOOST_CHECK_CLOSE_FRACTION(airy_ai_zero<RealType>(1000000001U), static_cast<RealType>(-2.81078366780730091663459728526906320267920607427246E+6L), tolerance);




  // Test Data for airy_bi
  using boost::math::airy_bi_zero;

  if (std::numeric_limits<RealType>::has_quiet_NaN)
  {
    BOOST_CHECK(isnan(airy_bi_zero<RealType>(0)) );
  }

  if (std::numeric_limits<RealType>::has_infinity)
  {
    BOOST_CHECK(isnan(airy_bi_zero<RealType>(std::numeric_limits<RealType>::infinity)) );
  }


  // Table[N[AiryBiZero[n], 51], {n, 1, 20, 1}]
  BOOST_CHECK_CLOSE_FRACTION(airy_bi_zero<RealType>(1U), static_cast<RealType>(-1.17371322270912792491997996247390210454364638917570L), tolerance);
  BOOST_CHECK_CLOSE_FRACTION(airy_bi_zero<RealType>(2U), static_cast<RealType>(-3.27109330283635271568022824016641380630093596910028L), tolerance);
  BOOST_CHECK_CLOSE_FRACTION(airy_bi_zero<RealType>(3U), static_cast<RealType>(-4.83073784166201593266770933990517817696614261732301L), tolerance);
  BOOST_CHECK_CLOSE_FRACTION(airy_bi_zero<RealType>(4U), static_cast<RealType>(-6.16985212831025125983336452055593667996554943427563L), tolerance);
  BOOST_CHECK_CLOSE_FRACTION(airy_bi_zero<RealType>(5U), static_cast<RealType>(-7.37676207936776371359995933044254122209152229939710L), tolerance);
  BOOST_CHECK_CLOSE_FRACTION(airy_bi_zero<RealType>(6U), static_cast<RealType>(-8.49194884650938801344803949280977672860508755505546L), tolerance);
  BOOST_CHECK_CLOSE_FRACTION(airy_bi_zero<RealType>(7U), static_cast<RealType>(-9.53819437934623888663298854515601962083907207638247L), tolerance);
  BOOST_CHECK_CLOSE_FRACTION(airy_bi_zero<RealType>(8U), static_cast<RealType>(-10.5299135067053579244005555984531479995295775946214L), tolerance);
  BOOST_CHECK_CLOSE_FRACTION(airy_bi_zero<RealType>(9U), static_cast<RealType>(-11.4769535512787794379234649247328196719482538148877L), tolerance);
  BOOST_CHECK_CLOSE_FRACTION(airy_bi_zero<RealType>(10U), static_cast<RealType>(-12.3864171385827387455619015028632809482597983846856L), tolerance);
  BOOST_CHECK_CLOSE_FRACTION(airy_bi_zero<RealType>(11U), static_cast<RealType>(-13.2636395229418055541107433243954907752411519609813L), tolerance);
  BOOST_CHECK_CLOSE_FRACTION(airy_bi_zero<RealType>(12U), static_cast<RealType>(-14.1127568090686577915873097822240184716840428285509L), tolerance);
  BOOST_CHECK_CLOSE_FRACTION(airy_bi_zero<RealType>(13U),  static_cast<RealType>(-14.9370574121541640402032143104909046396121763517782L), tolerance);
  BOOST_CHECK_CLOSE_FRACTION(airy_bi_zero<RealType>(14U), static_cast<RealType>(-15.7392103511904827708949784797481833807180162767841L), tolerance);
  BOOST_CHECK_CLOSE_FRACTION(airy_bi_zero<RealType>(15U), static_cast<RealType>(-16.5214195506343790539179499652105457167110310370581L), tolerance);
  BOOST_CHECK_CLOSE_FRACTION(airy_bi_zero<RealType>(16U), static_cast<RealType>(-17.2855316245812425329342366922535392425279753602710L), tolerance);
  BOOST_CHECK_CLOSE_FRACTION(airy_bi_zero<RealType>(17U), static_cast<RealType>(-18.0331132872250015721711125433391920008087291416406L), tolerance);
  BOOST_CHECK_CLOSE_FRACTION(airy_bi_zero<RealType>(18U), static_cast<RealType>(-18.7655082844800810413429789236105128440267189551421L), tolerance);
  BOOST_CHECK_CLOSE_FRACTION(airy_bi_zero<RealType>(19U), static_cast<RealType>(-19.4838801329892340136659986592413575122062977793610L), tolerance);
  BOOST_CHECK_CLOSE_FRACTION(airy_bi_zero<RealType>(20U), static_cast<RealType>(-20.1892447853962024202253232258275360764649783583934L), tolerance);


 // Table[N[AiryBiZero[n], 51], {n, 1000, 1001, 1}]
  BOOST_CHECK_CLOSE_FRACTION(airy_bi_zero<RealType>(1000U), static_cast<RealType>(-280.937811203415240157883427412260300146245056425646L), tolerance);
  BOOST_CHECK_CLOSE_FRACTION(airy_bi_zero<RealType>(1001U), static_cast<RealType>(-281.125212400956392021977771104562061554648675044114L), tolerance);

  // Table[N[AiryBiZero[n], 51], {n, 1000000, 1000001, 1}]
  BOOST_CHECK_CLOSE_FRACTION(airy_bi_zero<RealType>(1000000U), static_cast<RealType>(-28107.8226100991339342855024130953986989636667226163L), tolerance);
  BOOST_CHECK_CLOSE_FRACTION(airy_bi_zero<RealType>(1000001U), static_cast<RealType>(-28107.8413486584714939255315213519230566014624895515L), tolerance);

  //Table[N[AiryBiZero[n], 51], {n, 1000000000, 1000000001, 1}]
  BOOST_CHECK_CLOSE_FRACTION(airy_bi_zero<RealType>(1000000000U), static_cast<RealType>(-2.81078366499651725023268820158218492845371527054171E+6L), tolerance);
  BOOST_CHECK_CLOSE_FRACTION(airy_bi_zero<RealType>(1000000001U), static_cast<RealType>(-2.81078366687037302799011557215619265502627118526716E+6L), tolerance);

  
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




