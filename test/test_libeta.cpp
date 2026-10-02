// Copyright Nick Thompson, 2026
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "math_unit_test.hpp"
#include <cmath>
#include <limits>
#include <stdexcept>
#include <boost/math/concepts/real_concept.hpp>
#include <boost/math/special_functions/beta.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>

using boost::math::libeta;
using boost::math::libetac;
using boost::multiprecision::cpp_bin_float_100;

// log(ibeta(a, b, x)) and log(ibetac(a, b, x)) from Arb's arb_hypgeom_beta_lower, with the precision
// doubled until each logarithm's ball had relative radius below 1e-115. The inputs are dyadic, so exact
// in every type. Several targets underflow double: there log(ibeta) is -infinity, but libeta is not.
struct libeta_case
{
   double a;
   double b;
   double x;
   const char* log_p;
   const char* log_q;
};

const libeta_case cases[] = {
      { 2.5, 7.25, 0.3125, "-3.7167156277665188049290840111055425693593488651567342000012502851004994412078034526975515049850137362945980883e-01", "-1.1698312803088233349651352879295389693292834586979826954540514935872825554887239020131638530480843027129424824e+00" },
      { 2.5, 7.25, 0.875, "-4.4345947122224081046530969749310990093659919691233663114625273621430941801405299283156934670374411191393943876e-06", "-1.2326076547615867168151504108295024750650355692244371499856114021272344401647386276346484452551499377402488966e+01" },
      { 0.5, 0.5, 0.5, "-6.9314718055994530941723212145817656807550013436025525412068000949339362196969471560586332699641868754200148102e-01", "-6.9314718055994530941723212145817656807550013436025525412068000949339362196969471560586332699641868754200148102e-01" },
      { 1.0, 1.0, 0.1875, "-1.6739764335716715462736832489101805676545099796182715647480257043360801946601698955498375531719426613800358950e+00", "-2.0763936477824450161544104426738766749673259268081390006367452731010982055436868302387864454518526139874773199e-01" },
      { 0.5, 0.5, 9.094947017729282e-13, "-1.4314526316488209470620542120468834675861013543469559344594411813858688169108548092137802186276440566259295234e+01", "-6.0712811052568164727952258383323139837273179818800005468192845698497291569843003152216679057253352649397217909e-07" },
      { 9.5367431640625e-7, 3.0, 0.5, "-6.4990279784362527021851305783498086889378815860578000369818737418560387836805289126143048273527262421196071331e-08", "-1.6549028152507438492978740557867797551584908242395123816839275331964007995320601266648979794701183754577299207e+01" },
      { 3.0, 9.5367431640625e-7, 0.5, "-1.6549028152507438492978740557867797551584908242395123816839275331964007995320601266648979794701183754577299207e+01", "-6.4990279784362527021851305783498086889378815860578000369818737418560387836805289126143048273527262421196071331e-08" },
      { 0.0009765625, 0.0009765625, 0.125, "-6.9504617998484206497446276820368611445317909319707817733433332393193935558041868358136995705973397981340663124e-01", "-6.9125178049975242775597449072452416228716139105677147030161509217245257586337866435580760659714648626361431213e-01" },
      { 50.0, 2.5, 0.0625, "-1.3310392296311938936477859354681779646430842972669870304195913133829620589539560964307738984041510250414505628e+02", "-1.5620708837902603944950917943401674846306871183800437733485803592936921101711474787720997003896385889333011021e-58" },
      { 1000.0, 1000.0, 0.375, "-6.7879098782070203571535741057973743782955536957303318080122891535681449203670406857800476852551450175275314762e+01", "-3.3149880133573969922824995974224965637056758043476670891119727642919243631695704636348416155978166270733913420e-30" },
      { 1000.0, 1000.0, 0.625, "-3.3149880133573969922824995974224965637056758043476670891119727642919243631695704636348416155978166270733913420e-30", "-6.7879098782070203571535741057973743782955536957303318080122891535681449203670406857800476852551450175275314762e+01" },
      { 10000.0, 10000.0, 0.25, "-2.8819984220755331098839194770117315448486264030839286705152670660549762353701257939070968841957833171644866865e+03", "-1.3633688684767639881765687928093097744250603074917173616473757806904196132829559298164511082526880657885781775e-613" },
      { 10000.0, 10000.0, 0.75, "-1.3633688684767639881765687928093097744250603074917173616473757806904196132829559298164511082526880657885781775e-613", "-2.8819984220755331098839194770117315448486264030839286705152670660549762353701257939070968841957833171644866865e+03" },
      { 10000.0, 0.5, 0.875, "-1.3394521027626777099748952446997573462746702885509053896612331431976821464436216161037222597806676708845142582e+03", "-1.9201846626087125029027612859359907113707395383457549239532866991395985295359825524988696116301598187398883665e-582" },
      { 1.048576e6, 1024.0, 0.5, "-7.1941557634186688357844532955373761625884349655139799481873092808405534876172700544822832792998463119052012317e+05", "-3.1661843852186122127237313070884491832685258910359624873249688080229877326939015766751471565635099226598406911e-611" },
      { 100.0, 1.048576e6, 0.0009765625, "-1.4486608968995453830084851361306300995014059641294937436625787355612796662632052991551882269877592472658236946e-303", "-6.9731264356684974023614850681051528231841414568944325979176586916712546741354321795568654597356458739661894994e+02" },
      { 37.5, 40.25, 0.5, "-4.7307555088411157725979633902922107691136049613019416337906660727145158645379390792662893477597617498250318147e-01", "-9.7573026146792365826661180434525064594695070782057324369556972488301090797505588648171444949724195815954006735e-01" },
      { 4.0, 3.0, 0.0078125, "-1.6712608432463251991856384283020918341020810329590680592279080948835199946873096273593449850993901626536894712e+01", "-5.5183137805310855384254552540924836333280609811359952446688756907101759825770857836273702787048662110831118110e-08" }
};

template <class Real>
Real from_reference(const cpp_bin_float_100& x)
{
   return static_cast<Real>(x);
}

#ifndef BOOST_MATH_NO_REAL_CONCEPT_TESTS
template <>
boost::math::concepts::real_concept from_reference<boost::math::concepts::real_concept>(const cpp_bin_float_100& x)
{
   return boost::math::concepts::real_concept(static_cast<long double>(x));
}
#endif

// Where |log| < 1 it passes through zero, so relative error is meaningless: check the absolute error.
template <class Real>
void check(const Real& expected, const Real& computed, std::size_t ulps, const Real& absolute_tolerance)
{
   using std::abs;
   if (abs(expected) >= 1)
      CHECK_ULP_CLOSE(expected, computed, ulps);
   else
      CHECK_LE(Real(abs(expected - computed)), absolute_tolerance);
}

template <class Real>
void test_spots(std::size_t ulps, Real absolute_tolerance)
{
   for (const libeta_case& c : cases)
   {
      Real a(c.a), b(c.b), x(c.x);
      check(from_reference<Real>(cpp_bin_float_100(c.log_p)), Real(libeta(a, b, x)), ulps, absolute_tolerance);
      check(from_reference<Real>(cpp_bin_float_100(c.log_q)), Real(libetac(a, b, x)), ulps, absolute_tolerance);
   }

   // The boundaries, which follow ibeta: a == 0 gives P = 1 and b == 0 gives P = 0 whatever x is.
   using ignore_overflow = boost::math::policies::policy<boost::math::policies::overflow_error<boost::math::policies::ignore_error>>;
   CHECK_EQUAL(Real(libeta(Real(2), Real(3), Real(1))), Real(0));
   CHECK_EQUAL(Real(libetac(Real(2), Real(3), Real(0))), Real(0));
   CHECK_EQUAL(Real(libeta(Real(0), Real(3), Real(0))), Real(0));
   CHECK_EQUAL(Real(libetac(Real(3), Real(0), Real(1))), Real(0));
   if (std::numeric_limits<Real>::has_infinity)
   {
      const Real minus_infinity = -std::numeric_limits<Real>::infinity();
      CHECK_EQUAL(Real(libeta(Real(2), Real(3), Real(0), ignore_overflow())), minus_infinity);
      CHECK_EQUAL(Real(libetac(Real(2), Real(3), Real(1), ignore_overflow())), minus_infinity);
      CHECK_EQUAL(Real(libeta(Real(2), Real(0), Real(0.5), ignore_overflow())), minus_infinity);
      CHECK_EQUAL(Real(libetac(Real(0), Real(2), Real(0.5), ignore_overflow())), minus_infinity);
   }
#ifndef BOOST_NO_EXCEPTIONS
   CHECK_THROW(libeta(Real(2), Real(3), Real(0)), std::overflow_error);
   CHECK_THROW(libeta(Real(2), Real(3), Real(-0.5)), std::domain_error);
   CHECK_THROW(libeta(Real(2), Real(3), Real(1.5)), std::domain_error);
   CHECK_THROW(libetac(Real(-2), Real(3), Real(0.5)), std::domain_error);
   CHECK_THROW(libetac(Real(0), Real(0), Real(0.5)), std::domain_error);
#endif
}

// With a = 1, the continued fraction for an underflowing target starts with 0/0, so libeta must use
// P = 1 - (1 - x)^b instead. In double these targets are subnormal; to full precision,
// log(P) = log(b) + log(x) when b x is far below epsilon.
void test_a_equal_to_one_underflow()
{
   using std::log;
   const double ln_two = boost::math::constants::ln_two<double>();
   // P = 1 - (1 - 2^-1070)^(2^17), so log(P) = (17 - 1070) log(2):
   CHECK_ULP_CLOSE(-1053 * ln_two, libeta(1.0, std::ldexp(1.0, 17), std::ldexp(1.0, -1070)), 1);
   // Q = 1 - (1/2)^(2^-1040) after swapping a and b, so log(Q) = -1040 log(2) + log(log(2)):
   CHECK_ULP_CLOSE(-1040 * ln_two + log(ln_two), libetac(std::ldexp(1.0, -1040), 1.0, 0.5), 1);
}

// There is no Lanczos approximation at 100 digits, so libeta takes its generic code path. With the exponent range
// narrowed to about double's, the table's targets below exp(-800) underflow, so only the generic log-space power
// terms produce them. (Elsewhere ibeta itself currently loses digits at this precision, which libeta would inherit.)
void test_generic_underflow()
{
   using narrow_100 = boost::multiprecision::number<boost::multiprecision::cpp_bin_float<100, boost::multiprecision::digit_base_10, void, std::int32_t, -1100, 1100>>;
   static_assert(std::is_same<boost::math::lanczos::lanczos<narrow_100, boost::math::policies::policy<>>::type, boost::math::lanczos::undefined_lanczos>::value,
                 "This test needs a type without a Lanczos approximation.");
   int checked = 0;
   for (const libeta_case& c : cases)
   {
      narrow_100 a(c.a), b(c.b), x(c.x);
      narrow_100 log_p = from_reference<narrow_100>(cpp_bin_float_100(c.log_p));
      narrow_100 log_q = from_reference<narrow_100>(cpp_bin_float_100(c.log_q));
      if (log_p < -800)
      {
         CHECK_ULP_CLOSE(log_p, narrow_100(libeta(a, b, x)), 8);
         ++checked;
      }
      if (log_q < -800)
      {
         CHECK_ULP_CLOSE(log_q, narrow_100(libetac(a, b, x)), 8);
         ++checked;
      }
   }
   CHECK_LE(4, checked);
}

int main()
{
   test_a_equal_to_one_underflow();
   // Near the mean, libeta is log(ibeta), so it inherits ibeta's error: 18.5 epsilon at (37.5, 40.25, 0.5) in double.
   test_spots<float>(2, 2 * std::numeric_limits<float>::epsilon());
   test_spots<double>(8, 32 * std::numeric_limits<double>::epsilon());
#ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
   test_spots<long double>(16, 64 * std::numeric_limits<long double>::epsilon());
#endif
#ifndef BOOST_MATH_NO_REAL_CONCEPT_TESTS
   test_spots<boost::math::concepts::real_concept>(16, boost::math::concepts::real_concept(64 * std::numeric_limits<long double>::epsilon()));
#endif
   using boost::multiprecision::cpp_bin_float_50;
   test_spots<cpp_bin_float_50>(16, 16 * std::numeric_limits<cpp_bin_float_50>::epsilon());
   test_generic_underflow();
   return boost::math::test::report_errors();
}
