//  Copyright Nick Thompson 2026.
//  Distributed under the Boost Software License, Version 1.0.
//  (See accompanying file LICENSE_1_0.txt or copy at
//  https://www.boost.org/LICENSE_1_0.txt)
//
// logcdf for the distributions whose cdf is a regularised incomplete gamma or beta function: these use
// lgamma_p/lgamma_q and libeta/libetac, so they stay finite and accurate deep in the tails, where
// log(cdf) underflows to -infinity. References: mpmath at 80 digits, printed to 40 for 128-bit long double, at
// the double values of the arguments, each tail computed directly and the larger one's log through log1p.

#include "math_unit_test.hpp"
#include <boost/math/distributions/beta.hpp>
#include <boost/math/distributions/binomial.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/extreme_value.hpp>
#include <boost/math/distributions/fisher_f.hpp>
#include <boost/math/distributions/gamma.hpp>
#include <boost/math/distributions/inverse_chi_squared.hpp>
#include <boost/math/distributions/inverse_gamma.hpp>
#include <boost/math/distributions/laplace.hpp>
#include <boost/math/distributions/logistic.hpp>
#include <boost/math/distributions/negative_binomial.hpp>
#include <boost/math/distributions/poisson.hpp>
#include <boost/math/distributions/students_t.hpp>
#include <boost/math/distributions/weibull.hpp>
#include <cmath>

using namespace boost::math;

// Some references (log(1 - tiny)) are below the range of T, where long double is double.
template <class T>
void compare(T got, long double expected, long double condition = 1)
{
   if (std::abs(expected) < (std::numeric_limits<T>::min)())
   {
      CHECK_LE(std::abs(got), (std::numeric_limits<T>::min)());
   }
   else if (std::numeric_limits<T>::digits > std::numeric_limits<double>::digits && std::abs(expected) < std::numeric_limits<T>::epsilon())
   {
      // log(1 - P) = -P for a tiny tail P, so the error is that of gamma_p/ibeta deep in the tail, which in
      // 80 and 128-bit long double reaches ~1000 ulps.
      CHECK_ULP_CLOSE(static_cast<T>(expected), got, static_cast<unsigned>((std::max)(2048.0L, 4 * condition)));
   }
   else
   {
      // ibeta's own error is the limit at large parameters (82 ulps for beta(200, 300) at 0.01).
      CHECK_ULP_CLOSE(static_cast<T>(expected), got, static_cast<unsigned>((std::max)(128.0L, 4 * condition)));
   }
}

// One tail is log1p(-exp(L)) of the other's log L. Where exp(L) is tiny, that multiplies the error in L
// by |L|, so the tolerance follows that condition number:
template <class Dist, class T>
void check(Dist const& dist, T x, long double log_p, long double log_q)
{
   compare(logcdf(dist, x), log_p, std::fabs(log_p) < 0.5L ? std::fabs(log_q) : 1);
   compare(logcdf(complement(dist, x)), log_q, std::fabs(log_q) < 0.5L ? std::fabs(log_p) : 1);
}

template <class T>
void test_logcdf()
{
   check(gamma_distribution<T>(3, 2), T(1e-5), -38.4099811558181071937949822254812246044L, -2.083325520848958844599724448666888704154e-17L);
   check(gamma_distribution<T>(3, 2), T(5), -0.7848527198374080791887526789161297924804L, -0.6091496281277140941072272253454719338665L);
   check(gamma_distribution<T>(3, 2), T(2000), -2.543060483631175388867961841076222370909e-429L, -986.8756366239270061386403350088440606221L);
   check(chi_squared_distribution<T>(4), T(1e-200), -923.113478739298109571247753439954693873L, -1.249999999999999955250655997706899390637e-401L);
   check(chi_squared_distribution<T>(4), T(3), -0.8160504531201044357014920530998701328173L, -0.5837092681258449348164727882319889285499L);
   check(chi_squared_distribution<T>(4), T(3000), -5.428224956160711803218770091166157414999e-649L, -1492.686113168366538045818027026630735541L);
   check(poisson_distribution<T>(3), T(0), -3.0L, -0.05106918094270158653872374052225925040403L);
   check(poisson_distribution<T>(3), T(2), -0.8599338365037292291676975035850500324873L, -0.5502424967772210739746239648958146543321L);
   check(poisson_distribution<T>(3), T(300), -2.241714473901493972281266342976153989923e-475L, -1092.920678208931809905521913588722779425L);
   check(poisson_distribution<T>(1000), T(10), -946.0168196296327941947927036687836193172L, -1.412913085111672570692003173433196756301e-411L);
   check(inverse_gamma_distribution<T>(3, 2), T(0.001), -1985.490342261642280711390591284337689917L, -5.158227394204309104543078555545342025769e-863L);
   check(inverse_gamma_distribution<T>(3, 2), T(1), -0.3905620875658996253992406667738123604744L, -1.129101649750928642851458680739580923656L);
   check(inverse_gamma_distribution<T>(3, 2), T(1e200), -1.333333333333333454400844483291894024106e-600L, -1381.263373723975629392555020442155780733L);
   check(inverse_chi_squared_distribution<T>(4, 1), T(0.001), -1992.398597665416224978057585027971434536L, -5.155648281795474993741236214897568167919e-866L);
   check(inverse_chi_squared_distribution<T>(4, 1), T(0.5), -2.390562087565899625399240666773812360474L, -0.0960464646835615600303621221018613777488L);
   check(inverse_chi_squared_distribution<T>(4, 1), T(1e100), -1.999999999999999936388435560960329330067e-200L, -459.8238714182491915259868410349345007929L);
   check(beta_distribution<T>(2, 3), T(1e-200), -919.242277728390218642183579717199461896L, -5.999999999999999785203148788993117075057e-400L);
   check(beta_distribution<T>(2, 3), T(0.3), -1.054691101610181603246491925280495985446L, -0.4281709456438023306957306177940229395512L);
   check(beta_distribution<T>(2, 3), T(0.999999), -3.999997000345067637136471505133501476859e-18L, -40.06023806268694597164937815153169172127L);
   check(beta_distribution<T>(200, 300), T(0.01), -591.3421539316795129251336275659566742894L, -1.525336387634633744977167522063499563125e-257L);
   check(beta_distribution<T>(200, 300), T(0.4), -0.6883017062147303609057176958307249733836L, -0.6980162478923903525229536515796548666059L);
   check(beta_distribution<T>(200, 300), T(0.99), -2.754801960637096324366520479106238786735e-457L, -1051.268041941874049727730852801512754922L);
   check(binomial_distribution<T>(1000, 0.5), T(10), -639.2090519107084677186820095331484419643L, -2.483338791489635329496411692609225536818e-278L);
   check(binomial_distribution<T>(1000, 0.5), T(500), -0.6682350621326347656354844687758574554535L, -0.7186958030550969098247646442306168424813L);
   check(binomial_distribution<T>(1000, 0.5), T(990), -2.503319603882796982500653673088143766104e-280L, -643.8062083457074886731345211669864578298L);
   check(negative_binomial_distribution<T>(5, 0.5), T(3), -1.012577951326306538013447410200529629507L, -0.4514272436728001412729246055446626679721L);
   check(negative_binomial_distribution<T>(5, 0.5), T(2000), -1.830934139019977188922228120403794231405e-592L, -1362.525548757380464465430200340333535489L);
   check(fisher_f_distribution<T>(5, 7), T(1e-100), -574.1019093019765557361639405185045065192L, -4.684990774515470876061643685681347995374e-250L);
   check(fisher_f_distribution<T>(5, 7), T(1), -0.6564178625060657254346602263095448269658L, -0.7312771498375994711240725311825173468129L);
   check(fisher_f_distribution<T>(5, 7), T(1e100), -2.519700478313008342266150540721178897428e-349L, -802.6780574182750595910799873265228880996L);
   check(students_t_distribution<T>(5), T(-1e70), -803.6545263121025131443400282712354083919L, -9.490167245562357338012859910823596491189e-350L);
   check(students_t_distribution<T>(5), T(-1.5), -2.333540917368008583533202342985570491282L, -0.101979393788242724150847655148827902832L);
   check(students_t_distribution<T>(5), T(1e-9), -0.6931471798007319300604456897198997608536L, -0.6931471813191586893504235090284194356903L);
   check(students_t_distribution<T>(5), T(1e70), -9.490167245562357338012859910823596491189e-350L, -803.6545263121025131443400282712354083919L);
   check(students_t_distribution<T>(1e12), T(-40), -804.6084413729540376061252065587134160101L, -3.655895883611449397336474100007682207859e-350L);
   check(students_t_distribution<T>(1e12), T(40), -3.655895883611449397336474100007682207859e-350L, -804.6084413729540376061252065587134160101L);
   // |x| > sqrt(max) in double, where df / x^2 may underflow:
   check(students_t_distribution<T>(1), T(-1e200), -461.6617484846585369474748405007362614333L, -3.183098861837906811720139504195578159356e-201L);
   check(students_t_distribution<T>(2.5), T(-1e160), -921.3634587421731073387873653405855378588L, -7.193397190831722297482216141085036842542e-401L);
   check(students_t_distribution<T>(7), T(-1e155), -2493.42162960194500325806622545904493395L, -1.320520676354676868206420818690757852191e-1083L);
   check(students_t_distribution<T>(0.5), T(-1e300), -346.5250100798783758198490653831970251013L, -3.207009754142228919211973798954086743678e-151L);
   compare(logcdf(students_t_distribution<T>(1e300), T(-1e200)), -1.151292546497022901892142285621235491416e+302L);
   check(students_t_distribution<T>(1e300), T(-40), -804.6084420137537881666068329186099362001L, -3.655893540915029703748985802688283665054e-350L);
}

// At an infinite x, the log of the tail that vanishes is an exact -infinity, which goes to the overflow
// policy, and the log of the other tail is zero:
template <template <class, class> class Distribution, class T>
void test_infinite_x()
{
   using Dist = Distribution<T, boost::math::policies::policy<>>;
   using Quiet = Distribution<T, boost::math::policies::policy<boost::math::policies::overflow_error<boost::math::policies::ignore_error>>>;
   const T inf = std::numeric_limits<T>::infinity();
   CHECK_EQUAL(boost::math::logcdf(Dist(), inf), T(0));
   CHECK_EQUAL(boost::math::logcdf(boost::math::complement(Dist(), -inf)), T(0));
   CHECK_EQUAL(boost::math::logcdf(Quiet(), -inf), -inf);
   CHECK_EQUAL(boost::math::logcdf(boost::math::complement(Quiet(), inf)), -inf);
   CHECK_THROW(boost::math::logcdf(Dist(), -inf), std::overflow_error);
   CHECK_THROW(boost::math::logcdf(boost::math::complement(Dist(), inf)), std::overflow_error);
}

int main()
{
   test_logcdf<double>();
   test_logcdf<long double>();
   test_infinite_x<boost::math::extreme_value_distribution, double>();
   test_infinite_x<boost::math::laplace_distribution, double>();
   test_infinite_x<boost::math::logistic_distribution, double>();
   // The Weibull pdf vanishes at 0 for shape > 1:
   CHECK_THROW(boost::math::logpdf(boost::math::weibull_distribution<double>(2), 0.0), std::overflow_error);
   return boost::math::test::report_errors();
}
