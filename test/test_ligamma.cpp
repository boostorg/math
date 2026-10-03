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
#include <boost/math/special_functions/gamma.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>

using boost::math::ligamma;
using boost::math::ligamma_lower;
using boost::multiprecision::cpp_bin_float_100;

// log(tgamma(a, x)) and log(tgamma_lower(a, x)), computed with Arb: arb_hypgeom_gamma_upper for x >= a, and for x < a
// the series for the lower function, each with the other as gamma(a) minus the first. The precision was doubled until
// each logarithm's ball had relative radius below 1e-115. The inputs are dyadic, so exact in every type. Several values
// overflow or underflow double, where log(tgamma(a, x)) is infinite but ligamma is not.
struct ligamma_case
{
   double a;
   double x;
   const char* log_upper;
   const char* log_lower;
};

const ligamma_case cases[] = {
      { 2.5, 1.0, "1.2115759487826881505275305651186870694041754810602241686386550907439660274410122689632242508378186319701322744e-01", "-1.6067535371453137835485637734565884075661272599144183153328939345813381648202869722143154438944755967497311409e+00" },
      { 2.5, 30.0, "-2.4848633033033568585775606150217063231776860531174716723611152050121158069210433588962288107825722546981486249e+01", "2.8468287046076458985523776395111859010124529460512663635052087420960128616776934174356639565068600583951170758e-01" },
      { 0.5, 0.5, "-5.7550952152461810928183727609782309764852436800536248513008065098880573654672044680219886267372404171768343111e-01", "1.9064979662257401479753066816587055020203900191383989139279624179991169101923287236368302011386826091354535869e-01" },
      { 9.5367431640625e-7, 0.25, "4.3329681976764116812962168205761668110531802122293663061271899236280959010116297976855105750784261579707721890e-02", "1.3862942064817819371597083203484523981378986248191528178643086339948011084679626489495691758190580840604079834e+01" },
      { 9.5367431640625e-7, 9.313225746154785e-10, "3.0065235630469892943086953167194409014444406877094897610620142931242359707297922698528910680764788181213409437e+00", "1.3862923780098997617405608217022900716380592969540924868089321252596890761811200226264389634053398086943925588e+01" },
      { 1.0, 0.125, "-1.2500000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000e-01", "-2.1412905847632012003460915474560573950973053346029827451505631309617872378509999247775869790609354267681006076e+00" },
      { 100.0, 7.62939453125e-6, "3.5913420536957539877604401046028690961262171808562972877561279307484079922862431495362900583179335272293183799e+02", "-1.1829553846917510861148888194102544355366943034135567498326091371553742251900583033628108371168959536911942860e+03" },
      { 0.5, 800.0, "-8.0334292989026914277011788712147523027625048445796001213977030571616076706566234902264761738165995628239877020e+02", "5.7236494292470008707171367567652935582364740645765578575681153573606888494241303989181163513774485385100490611e-01" },
      { 37.5, 37.5, "9.6784219013624273538564367039851830199348647345043818460498983911692030778614396820970016069787243324610329474e+01", "9.6871148482431460918087441810862444588322980641847852624659906118444232717140238995838074935912952693191695927e+01" },
      { 150.0, 160.0, "5.9842166600398933896473633740479780587741555270045023808144087418030819202526400004141534851719087907025996299e+02", "5.9978084473911611070673933986747762225349459031614663814392284583053929946405261709380560441291096424584533341e+02" },
      { 1000.0, 1000.0, "5.9045188299725352393303969775776785919325346037738600473424508213569779194944742662149261807627638282089426729e+03", "5.9045356513458909016799290303089150623692120881673413020117569273880583353428974787238677453871076071638650619e+03" },
      { 1000.0, 500.0, "5.9052204232091812118260769123614407898489424097154325900233875198883841333644788925921620543460926574075386671e+03", "5.7083915041483186807926091879564531470086287259517343889685970127418854617508689583124643132768875701028822116e+03" },
      { 1000.0, 4000.0, "4.2860428284350736577552641746493656020777712625262765706663139050781686399320799273897846364221909714378263326e+03", "5.9052204232091812118260769123614407898489424097154325900233875198883841333644788925921620873288206280785351523e+03" },
      { 1.0e6, 1.5e6, "1.2720962543703061159755127230382120390711888790333146381432032604790510263473733787158079665363156769867669015e+07", "1.2815504569147611659976971785017113153687975196214851361632706092068992333829159968797138805737861576248904294e+07" },
      { 1.0e6, 500000.0, "1.2815504569147611659976971785017113153687975196214851361632706092068992333829159968797138805737861576248904294e+07", "1.2622350255038951404361743252556656958733649336038873758330500797737834380859841254882259339362580687585976636e+07" },
      { 1.048576e6, 1.0496e6, "1.3487766106501424078633459785724088776137287211137638456602466106862593851691460461761679317593366554025664762e+07", "1.3487767774769577991040370044266971477801670978480955073684733273999936390880788611365025536637318183046721133e+07" },
      { 3.0, 9.332636185032189e-302, "6.9314718055994530941723212145817656807550013436025525412068000949339362196969471560586332699641868754200148102e-01", "-2.0805401539685040379430916096114522299311478936385885118137747228138183602023027557844635967440697947147924131e+03" },
      { 4.0, 712.0, "-6.9229154983203152364214453894071131752969051036854104094512292173514623019574780148370046658687497631154999836e+02", "1.7917594692280550008124773583807022727229906921830047058553743431308879151883036824794790818101507763299715101e+00" }
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
   using std::ldexp;
   const int max_exponent = boost::math::itrunc(boost::math::tools::log_max_value<Real>() / boost::math::constants::ln_two<Real>()) - 2;
   const int min_exponent = boost::math::itrunc(boost::math::tools::log_min_value<Real>() / boost::math::constants::ln_two<Real>()) + 2;
   for (const ligamma_case& c : cases)
   {
      int ea, ex;
      std::frexp(c.a, &ea);
      std::frexp(c.x, &ex);
      if ((ea > max_exponent) || (ex > max_exponent) || (ex < min_exponent))
         continue;
      Real a(c.a), x(c.x);
      check(from_reference<Real>(cpp_bin_float_100(c.log_upper)), Real(ligamma(a, x)), ulps, absolute_tolerance);
      check(from_reference<Real>(cpp_bin_float_100(c.log_lower)), Real(ligamma_lower(a, x)), ulps, absolute_tolerance);
   }

   // At the ends of the range the functions are lgamma(a) or the logarithm of zero:
   using ignore_overflow = boost::math::policies::policy<boost::math::policies::overflow_error<boost::math::policies::ignore_error>>;
   CHECK_EQUAL(Real(ligamma(Real(2.5), Real(0))), Real(boost::math::lgamma(Real(2.5))));
   if (std::numeric_limits<Real>::has_infinity)
   {
      const Real infinity = std::numeric_limits<Real>::infinity();
      CHECK_EQUAL(Real(ligamma_lower(Real(2.5), infinity)), Real(boost::math::lgamma(Real(2.5))));
      CHECK_EQUAL(Real(ligamma(Real(2.5), infinity, ignore_overflow())), Real(-infinity));
      CHECK_EQUAL(Real(ligamma_lower(Real(2.5), Real(0), ignore_overflow())), Real(-infinity));
   }
   // NaN arguments are domain errors; ignoring the error gives NaN:
   if (std::numeric_limits<Real>::has_quiet_NaN)
   {
      using ignore_domain = boost::math::policies::policy<boost::math::policies::domain_error<boost::math::policies::ignore_error>>;
      const Real nan = std::numeric_limits<Real>::quiet_NaN();
      CHECK_TRUE((boost::math::isnan)(Real(ligamma(nan, Real(2), ignore_domain()))));
      CHECK_TRUE((boost::math::isnan)(Real(ligamma(Real(2), nan, ignore_domain()))));
      CHECK_TRUE((boost::math::isnan)(Real(ligamma_lower(nan, Real(2), ignore_domain()))));
      CHECK_TRUE((boost::math::isnan)(Real(ligamma_lower(Real(2), nan, ignore_domain()))));
#ifndef BOOST_NO_EXCEPTIONS
      CHECK_THROW(ligamma(nan, Real(2)), std::domain_error);
      CHECK_THROW(ligamma_lower(Real(2), nan), std::domain_error);
#endif
   }
#ifndef BOOST_NO_EXCEPTIONS
   CHECK_THROW(ligamma_lower(Real(2), Real(0)), std::overflow_error);
   CHECK_THROW(ligamma(Real(0), Real(1)), std::domain_error);
   CHECK_THROW(ligamma(Real(-1), Real(1)), std::domain_error);
   CHECK_THROW(ligamma_lower(Real(2), Real(-1)), std::domain_error);
#endif
}

int main()
{
   test_spots<float>(2, 2 * std::numeric_limits<float>::epsilon());
   test_spots<double>(4, 8 * std::numeric_limits<double>::epsilon());
#ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
   test_spots<long double>(8, 16 * std::numeric_limits<long double>::epsilon());
#endif
#ifndef BOOST_MATH_NO_REAL_CONCEPT_TESTS
   test_spots<boost::math::concepts::real_concept>(8, boost::math::concepts::real_concept(16 * std::numeric_limits<long double>::epsilon()));
#endif
   // At 50 digits ligamma inherits tgamma's errors of up to about 20 epsilon where it is log(tgamma):
   using boost::multiprecision::cpp_bin_float_50;
   test_spots<cpp_bin_float_50>(16, 64 * std::numeric_limits<cpp_bin_float_50>::epsilon());
   return boost::math::test::report_errors();
}
