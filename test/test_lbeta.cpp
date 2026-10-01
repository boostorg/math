// Copyright Nick Thompson, 2026
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "math_unit_test.hpp"
#include <cmath>
#include <stdexcept>
#include <boost/math/concepts/real_concept.hpp>
#include <boost/math/special_functions/beta.hpp>
#include <boost/math/special_functions/trunc.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>

using boost::math::lbeta;
using boost::multiprecision::cpp_bin_float_100;

// a = am * 2^ae and b = bm * 2^be are exact in every type whose exponent range holds them.
// The values of log(beta(a, b)) were computed as lgamma(a) + lgamma(b) - lgamma(a + b) with
// 800-digit cpp_bin_float, enough to absorb the cancellation when an argument is 2^1000.
struct lbeta_case
{
   double am;
   int ae;
   double bm;
   int be;
   const char* value;
};

const lbeta_case cases[] = {
      { 1, -60, 1, -50, "4.15898069195697740239298881123986100039915095416801455108115264409117033485238337627682905027475554197949082691e+01" },
      { 1, -100, 3, 0, "6.93147180559945309417232121458164735161921819183174078192200771589333072536381647225831806423195466923922439014e+01" },
      { 0.5, 0, 0.5, 0, "1.14472988584940017414342735135305871164729481291531157151362307147213776988482607978362327027548970770200981223e+00" },
      { 0.25, 0, 0.75, 0, "1.49130347612937282885204341208214699568504488009543919857396307621883458086967343758655493377369905147301055274e+00" },
      { 1.5, 0, 2.5, 0, "-1.62785883639038106352550113447964756065470572452570944496909696650143671799395278263983003771018504246599611185e+00" },
      { 2, 0, 0.75, 0, "-2.71933715483641758831669494532999161982574749963589623711364445601499668091899372237820945415928865613025799190e-01" },
      { 2.5, 0, 2.5, 0, "-2.60868808940210730038195226193165156023371556978372575559644266134412329068442796258380426388570901630403052585e+00" },
      { 2.5, 0, 7.25, 0, "-4.90533661883930395460834936695643873217655041103926652565285221548103484033393583597181887094577679099681463945e+00" },
      { 7.25, 0, 10, 0, "-1.15206093828576774238792925676023403474351969826024007014606103008266556359021233866230652434846476474713883000e+01" },
      { 33, 0, 70, 0, "-6.52309611058463901652530664244886314226939493764776246951104105272718684570163848327941904447963477233497622302e+01" },
      { 1000, 0, 2.5, 0, "-1.69865790781530187420413395240576683950617459492281475607293952992950240177344530711756692292901949091581441171e+01" },
      { 1, 20, 2.5, 0, "-3.43726779456625527055870837198099350327971371888199561621811864168868157894165670678814799899328577816876559438e+01" },
      { 1, 40, 0.125, 0, "-1.44631754524588046377436687345621730310702661961084719229448263705785310178427063985532264153346167268065642354e+00" },
      { 1, 60, 3, 0, "-1.24073345320230210388286634954978816325251619811257890584231209580161422603673162680077814295275683492847917281e+02" },
      { 7.25, 100, 1, -3, "-6.89254658305384294843871197074079636920399542538686737751765487854131887259626569180127928805023089821623285643e+00" },
      { 1, 20, 1, 20, "-1.45364066196521333105311408462318440517413066016099383479834829660518401041863896220590715872159106794242688144e+06" },
      { 1, 60, 1.5, 60, "-1.93982405936568553724082316071548785790406066925570806020585049192655840461166500901608060631298307929683783052e+18" },
      { 1, 100, 1, 90, "-9.81929078485523071471385110686110593997035850585013134716105089872660111623496248017922967392287894559854586350e+27" },
      { 2.5, 200, 2.5, 0, "-3.48579634239185123211942384088825609792055182534338388500569779067988144304556260885743443173501924555859593246e+02" },
      { 1, -1000, 1, 1000, "6.93147180559945309417232121458176568075500134360255254120680009493393621969694715605863326996418687542001481021e+02" },
      { 1, 1000, 1, 1000, "-1.48542634003375029367139526623878447128207418381951509618626907063551595543031141778517375655554942542595596811e+301" },
      { 1, -1000, 1, -1000, "6.93840327740505254726649353579634744643575634494615509374800689502887015591664410321469190323415106229543482502e+02" },
      { 1.5, 500, 7.25, 3, "-1.99483893641016901555392148959348216596399017908133451106724466122953422408704368366657864286454193673764307194e+04" }
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

// Where |log(beta)| < 1, it passes through zero and relative error is meaningless, so check the absolute error.
template <class Real>
void test_spots(std::size_t ulps, Real absolute_tolerance)
{
   using std::abs;
   using std::ldexp;
   const int max_exponent = boost::math::itrunc(boost::math::tools::log_max_value<Real>() / boost::math::constants::ln_two<Real>()) - 2;
   for (const lbeta_case& c : cases)
   {
      if ((abs(c.ae) > max_exponent) || (abs(c.be) > max_exponent))
         continue;
      Real a = ldexp(Real(c.am), c.ae);
      Real b = ldexp(Real(c.bm), c.be);
      Real expected = from_reference<Real>(cpp_bin_float_100(c.value));
      Real computed = lbeta(a, b);
      if (abs(expected) >= 1)
         CHECK_ULP_CLOSE(expected, computed, ulps);
      else
         CHECK_LE(Real(abs(expected - computed)), absolute_tolerance);
      CHECK_EQUAL(computed, Real(lbeta(b, a)));
   }
   CHECK_EQUAL(Real(lbeta(Real(1), Real(3))), Real(-log(Real(3))));
   CHECK_EQUAL(Real(lbeta(Real(0.25), Real(1))), Real(-log(Real(0.25))));
#ifndef BOOST_NO_EXCEPTIONS
   CHECK_THROW(lbeta(Real(0), Real(1)), std::domain_error);
   CHECK_THROW(lbeta(Real(1), Real(-2)), std::domain_error);
#endif
}

int main()
{
   test_spots<float>(2, 2 * std::numeric_limits<float>::epsilon());
   test_spots<double>(4, 2 * std::numeric_limits<double>::epsilon());
#ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
   test_spots<long double>(8, 4 * std::numeric_limits<long double>::epsilon());
#endif
#ifndef BOOST_MATH_NO_REAL_CONCEPT_TESTS
   test_spots<boost::math::concepts::real_concept>(8, boost::math::concepts::real_concept(4 * std::numeric_limits<long double>::epsilon()));
#endif
   using boost::multiprecision::cpp_bin_float_50;
   test_spots<cpp_bin_float_50>(8, 8 * std::numeric_limits<cpp_bin_float_50>::epsilon());
   // No Lanczos approximation at this precision, so this exercises the generic code path, whose
   // errors come from beta, lgamma and tgamma_delta_ratio at this precision:
   test_spots<cpp_bin_float_100>(64, 64 * std::numeric_limits<cpp_bin_float_100>::epsilon());
   return boost::math::test::report_errors();
}
