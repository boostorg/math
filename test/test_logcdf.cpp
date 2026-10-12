//  Copyright Nick Thompson 2026.
//  Distributed under the Boost Software License, Version 1.0.
//  (See accompanying file LICENSE_1_0.txt or copy at
//  https://www.boost.org/LICENSE_1_0.txt)
//
// logcdf for the distributions whose cdf is a regularised incomplete gamma or beta function: these use
// lgamma_p/lgamma_q and libeta/libetac, so they stay finite and accurate deep in the tails, where
// log(cdf) underflows to -infinity. References: mpmath at 80 digits at the double values of the arguments,
// each tail computed directly and the larger one's log through log1p.

#include "math_unit_test.hpp"
#include <boost/math/distributions/beta.hpp>
#include <boost/math/distributions/binomial.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/fisher_f.hpp>
#include <boost/math/distributions/gamma.hpp>
#include <boost/math/distributions/inverse_chi_squared.hpp>
#include <boost/math/distributions/inverse_gamma.hpp>
#include <boost/math/distributions/negative_binomial.hpp>
#include <boost/math/distributions/poisson.hpp>
#include <boost/math/distributions/students_t.hpp>
#include <cmath>

using namespace boost::math;

// Some references (log(1 - tiny)) are below the range of T, where long double is double.
template <class T>
void compare(T got, long double expected)
{
   if (std::abs(expected) < (std::numeric_limits<T>::min)())
   {
      CHECK_LE(std::abs(got), (std::numeric_limits<T>::min)());
   }
   else
   {
      // ibeta's own error is the limit at large parameters (82 ulps for beta(200, 300) at 0.01).
      CHECK_ULP_CLOSE(static_cast<T>(expected), got, 128);
   }
}

template <class Dist, class T>
void check(Dist const& dist, T x, long double log_p, long double log_q)
{
   compare(logcdf(dist, x), log_p);
   compare(logcdf(complement(dist, x)), log_q);
}

template <class T>
void test_logcdf()
{
   check(gamma_distribution<T>(3, 2), T(1e-5), -38.40998115581810719379498L, -2.083325520848958844599724e-17L);
   check(gamma_distribution<T>(3, 2), T(5), -0.7848527198374080791887527L, -0.6091496281277140941072272L);
   check(gamma_distribution<T>(3, 2), T(2000), -2.543060483631175388867962e-429L, -986.8756366239270061386403L);
   check(chi_squared_distribution<T>(4), T(1e-200), -923.1134787392981095712478L, -1.249999999999999955250656e-401L);
   check(chi_squared_distribution<T>(4), T(3), -0.8160504531201044357014921L, -0.5837092681258449348164728L);
   check(chi_squared_distribution<T>(4), T(3000), -5.42822495616071180321877e-649L, -1492.686113168366538045818L);
   check(poisson_distribution<T>(3), T(0), -3.0L, -0.05106918094270158653872374L);
   check(poisson_distribution<T>(3), T(2), -0.8599338365037292291676975L, -0.550242496777221073974624L);
   check(poisson_distribution<T>(3), T(300), -2.241714473901493972281266e-475L, -1092.920678208931809905522L);
   check(poisson_distribution<T>(1000), T(10), -946.0168196296327941947927L, -1.412913085111672570692003e-411L);
   check(inverse_gamma_distribution<T>(3, 2), T(0.001), -1985.490342261642280711391L, -5.158227394204309104543079e-863L);
   check(inverse_gamma_distribution<T>(3, 2), T(1), -0.3905620875658996253992407L, -1.129101649750928642851459L);
   check(inverse_gamma_distribution<T>(3, 2), T(1e200), -1.333333333333333454400844e-600L, -1381.263373723975629392555L);
   check(inverse_chi_squared_distribution<T>(4, 1), T(0.001), -1992.398597665416224978058L, -5.155648281795474993741236e-866L);
   check(inverse_chi_squared_distribution<T>(4, 1), T(0.5), -2.390562087565899625399241L, -0.09604646468356156003036212L);
   check(inverse_chi_squared_distribution<T>(4, 1), T(1e100), -1.999999999999999936388436e-200L, -459.8238714182491915259868L);
   check(beta_distribution<T>(2, 3), T(1e-200), -919.2422777283902186421836L, -5.999999999999999785203149e-400L);
   check(beta_distribution<T>(2, 3), T(0.3), -1.054691101610181603246492L, -0.4281709456438023306957306L);
   check(beta_distribution<T>(2, 3), T(0.999999), -3.999997000345067637136472e-18L, -40.06023806268694597164938L);
   check(beta_distribution<T>(200, 300), T(0.01), -591.3421539316795129251336L, -1.525336387634633744977168e-257L);
   check(beta_distribution<T>(200, 300), T(0.4), -0.6883017062147303609057177L, -0.6980162478923903525229537L);
   check(beta_distribution<T>(200, 300), T(0.99), -2.75480196063709632436652e-457L, -1051.268041941874049727731L);
   check(binomial_distribution<T>(1000, 0.5), T(10), -639.209051910708467718682L, -2.483338791489635329496412e-278L);
   check(binomial_distribution<T>(1000, 0.5), T(500), -0.6682350621326347656354845L, -0.7186958030550969098247646L);
   check(binomial_distribution<T>(1000, 0.5), T(990), -2.503319603882796982500654e-280L, -643.8062083457074886731345L);
   check(negative_binomial_distribution<T>(5, 0.5), T(3), -1.012577951326306538013447L, -0.4514272436728001412729246L);
   check(negative_binomial_distribution<T>(5, 0.5), T(2000), -1.830934139019977188922228e-592L, -1362.52554875738046446543L);
   check(fisher_f_distribution<T>(5, 7), T(1e-100), -574.1019093019765557361639L, -4.684990774515470876061644e-250L);
   check(fisher_f_distribution<T>(5, 7), T(1), -0.6564178625060657254346602L, -0.7312771498375994711240725L);
   check(fisher_f_distribution<T>(5, 7), T(1e100), -2.519700478313008342266151e-349L, -802.67805741827505959108L);
   check(students_t_distribution<T>(5), T(-1e70), -803.65452631210251314434L, -9.49016724556235733801286e-350L);
   check(students_t_distribution<T>(5), T(-1.5), -2.333540917368008583533202L, -0.1019793937882427241508477L);
   check(students_t_distribution<T>(5), T(1e-9), -0.6931471798007319300604457L, -0.6931471813191586893504235L);
   check(students_t_distribution<T>(5), T(1e70), -9.49016724556235733801286e-350L, -803.65452631210251314434L);
   check(students_t_distribution<T>(1e12), T(-40), -804.6084413729540376061252L, -3.655895883611449397336474e-350L);
   check(students_t_distribution<T>(1e12), T(40), -3.655895883611449397336474e-350L, -804.6084413729540376061252L);
   check(students_t_distribution<T>(1e300), T(-40), -804.6084420137537881666068L, -3.655893540915029703748986e-350L);
}

int main()
{
   test_logcdf<double>();
   test_logcdf<long double>();
   return boost::math::test::report_errors();
}
