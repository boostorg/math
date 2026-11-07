//  (C) Copyright Nick Thompson 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/math/special_functions/fresnel.hpp>
#include <boost/math/constants/constants.hpp>
#include <cmath>
#include <complex>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include "math_unit_test.hpp"

#ifndef BOOST_MATH_STANDALONE
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/cpp_complex.hpp>
#endif

using boost::math::fresnel_c;
using boost::math::fresnel_s;
using boost::math::fresnel_cis;
using float_in_float = boost::math::policies::policy<boost::math::policies::promote_float<false> >;

struct reference
{
   double x;
   double c;
   double s;
};

// mpmath's fresnelc and fresnels at 80 digits, at the exact double arguments. The rows lie on both sides of each switch
// point: the series up to x = 1.8368 in double (1.2361 for float evaluated in float), the continued fraction, and the
// leading asymptotic terms from x = 2^26 (2^11.5).
static const reference references[] = {
   {                   1e-100, 1.000000000000000019992e-100, 5.235987755982989044803e-301 },
   {                    1e-20, 9.999999999999999451533e-21, 5.235987755982987869241e-61 },
   {                    1e-10, 1.000000000000000036432e-10, 5.235987755982989303047e-31 },
   {                    1e-05, 1.000000000000000081801e-5, 5.235987755982990015721e-16 },
   {                    0.001, 9.999999999997532807067e-4, 5.235987755982066251905e-10 },
   {                    0.125, 1.249924702994110995056e-1, 1.022609856621743053937e-3 },
   {                      0.5, 4.923442258714463928788e-1, 6.473243285999927761148e-2 },
   {                     0.75, 6.935259907871358974859e-1, 2.088771112333835702224e-1 },
   {                      1.0, 7.798934003768228294742e-1, 4.382591473903547660768e-1 },
   {                 1.234375, 6.918587752705619847781e-1, 6.483803784681175674913e-1 },
   {                1.2421875, 6.860530795326560933472e-1, 6.536076488212688474774e-1 },
   {                     1.25, 6.800907410754550836679e-1, 6.586555116366791305604e-1 },
   {                      1.5, 4.452611760398215350646e-1, 6.975049600820930130807e-1 },
   {                     1.75, 3.219350461495378595884e-1, 4.993846745906745515668e-1 },
   {                   1.8368, 3.506520996447624300799e-1, 4.183871562153951924773e-1 },
   {                   1.8369, 3.507075283352108602154e-1, 4.183039236496671791475e-1 },
   {                      2.0, 4.882534060753407545002e-1, 3.434156783636982421953e-1 },
   {                     2.25, 6.401242099465546405236e-1, 5.053022268236937035593e-1 },
   {                      2.5, 4.574130096417770452457e-1, 6.191817558195929361136e-1 },
   {                      3.0, 6.057207892976856295562e-1, 4.963129989673750360976e-1 },
   {                      3.5, 5.325724350280008453295e-1, 4.152480119724375238977e-1 },
   {                      4.0, 4.984260330381776155307e-1, 4.205157542469284244453e-1 },
   {                      5.0, 5.636311887040122311021e-1, 4.991913819171168867519e-1 },
   {                     7.25, 5.337665692478711423152e-1, 4.719447271538931635396e-1 },
   {                     10.0, 4.998986942055157236142e-1, 4.681699785848822404034e-1 },
   {                     15.5, 5.078336554052517314774e-1, 4.810167854210038870429e-1 },
   {                     37.5, 4.967534577323829425305e-1, 5.078428670985603923392e-1 },
   {                    100.0, 4.999998986788178975595e-1, 4.968169011478375532715e-1 },
   {                   1000.5, 5.001217509508100827182e-1, 4.997060669389542172314e-1 },
   {                   2895.5, 5.000420693865084798103e-1, 4.998984355056238997633e-1 },
   {                  2897.25, 5.000696983910289338272e-1, 5.000849277286293120142e-1 },
   {                 10000.25, 5.000031199043797985593e-1, 4.999683230782052267722e-1 },
   {               100000.125, 5.000000781170588317721e-1, 4.999968178638056951688e-1 },
   {               67108863.5, 5.000000018151390658039e-1, 4.999999956178666497429e-1 },
   {               67108864.5, 5.000000018151390387563e-1, 4.999999956178667150418e-1 },
   {              100000000.0, 5.0e-1, 4.999999968169011381621e-1 },
   {             1500000000.0, 5.0e-1, 4.999999997877934092108e-1 },
   {          1000000000000.0, 5.0e-1, 4.999999999996816901138e-1 },
   {       1000000000000000.5, 5.000000000000001218119e-1, 4.9999999999999970592e-1 },
   {                    1e+16, 5.0e-1, 4.99999999999999968169e-1 },
   {       9007199254740994.0, 5.0e-1, 4.999999999999999646605e-1 },
   {                   1e+300, 5.0e-1, 5.0e-1 },
};

template <class Real, class Policy>
void test_references(const Policy& pol)
{
   for (const reference& r : references)
   {
      const Real x = static_cast<Real>(r.x);
      // Only the rows whose argument and results the type can hold.
      if (static_cast<double>(x) != r.x || r.s < (std::numeric_limits<Real>::min)())
      {
         continue;
      }
      const Real c = fresnel_c(x, pol);
      const Real s = fresnel_s(x, pol);
      CHECK_ULP_CLOSE(static_cast<Real>(r.c), c, 2);
      CHECK_ULP_CLOSE(static_cast<Real>(r.s), s, 2);
      const std::complex<Real> z = fresnel_cis(x, pol);
      CHECK_EQUAL(z.real(), c);
      CHECK_EQUAL(z.imag(), s);
      CHECK_EQUAL(fresnel_c(Real(-x), pol), Real(-c));
      CHECK_EQUAL(fresnel_s(Real(-x), pol), Real(-s));
   }
}

template <class Real>
void test_special_values()
{
   static_assert(std::is_same<decltype(fresnel_cis(Real(1))), std::complex<Real> >::value, "fresnel_cis returns std::complex");
   const Real inf = std::numeric_limits<Real>::infinity();
   const Real nan = std::numeric_limits<Real>::quiet_NaN();
   CHECK_EQUAL(fresnel_c(Real(0)), Real(0));
   CHECK_EQUAL(fresnel_s(Real(0)), Real(0));
   CHECK_EQUAL(fresnel_c(inf), Real(0.5));
   CHECK_EQUAL(fresnel_s(inf), Real(0.5));
   CHECK_EQUAL(fresnel_c(-inf), Real(-0.5));
   CHECK_EQUAL(fresnel_s(-inf), Real(-0.5));
   CHECK_EQUAL(fresnel_s((std::numeric_limits<Real>::max)()), Real(0.5));
   CHECK_THROW(fresnel_c(nan), std::domain_error);
   CHECK_THROW(fresnel_cis(nan), std::domain_error);
   using ignore = boost::math::policies::policy<boost::math::policies::domain_error<boost::math::policies::ignore_error> >;
   CHECK_NAN(fresnel_c(nan, ignore()));
   CHECK_NAN(fresnel_s(nan, ignore()));
   CHECK_NAN(fresnel_cis(nan, ignore()).real());
   // Only the component asked for is narrowed: S(1e-20) = 5e-61 underflows float, C does not.
   using throw_on_underflow = boost::math::policies::policy<boost::math::policies::underflow_error<boost::math::policies::throw_on_error> >;
   CHECK_ULP_CLOSE(Real(1e-20), fresnel_c(Real(1e-20), throw_on_underflow()), 1);
}

// C' = cos(pi x^2 / 2) and S' = sin(pi x^2 / 2), by central differences, across the switch points.
template <class Real>
void test_derivatives()
{
   using std::cos;
   using std::sin;
   using std::cbrt;
   const Real h = cbrt(std::numeric_limits<Real>::epsilon());
   const Real pi = boost::math::constants::pi<Real>();
   for (Real x = Real(0.0625); x < 12; x += Real(0.3125))
   {
      const Real dc = (fresnel_c(Real(x + h)) - fresnel_c(Real(x - h))) / (2 * h);
      const Real ds = (fresnel_s(Real(x + h)) - fresnel_s(Real(x - h))) / (2 * h);
      // The differences lose a third of the digits, and the truncation error is about (pi x)^2 h^2 / 6.
      const Real tol = 2 * h * h * (1 + (pi * x) * (pi * x));
      CHECK_ABSOLUTE_ERROR(cos(pi * x * x / 2), dc, tol);
      CHECK_ABSOLUTE_ERROR(sin(pi * x * x / 2), ds, tol);
   }
}

#ifndef BOOST_MATH_STANDALONE
void test_multiprecision()
{
   using boost::multiprecision::cpp_bin_float_50;
   using boost::multiprecision::cpp_complex_50;
   struct mp_reference
   {
      const char* x;
      const char* c;
      const char* s;
   };
   // mpmath at 80 digits.
   static const mp_reference mp_references[] = {
      { "0.5", "4.92344225871446392878843665156681637766095145771501253294653e-1", "6.47324328599992776114805122306147676507259184935124927875889e-2" },
      { "1.25", "6.80090741075455083667922305104958983784241559388376099014171e-1", "6.58655511636679130560412368153016199567044695691829379208089e-1" },
      { "3.75", "4.93590462646024960474683108954835291386009734200941247437068e-1", "5.84533202150472704719459454924418335415065571573929306528568e-1" },
      { "4.5", "5.26025915053538741369606936920230966515796770826606306072205e-1", "4.34272975048703589471916132425111149908241084375481668703982e-1" },
      { "9.125", "4.68056738770710906150288013487707638412601164910222633017687e-1", "4.85986360751166388293159802143237262102860699440370784891478e-1" },
      { "1000.0625", "5.00318284001961536065373542684528038167853232628090415062788e-1", "5.0000195288938210967679264547887177771141615760004272894327e-1" },
   };
   static_assert(std::is_same<decltype(fresnel_cis(cpp_bin_float_50())), cpp_complex_50>::value, "fresnel_cis returns cpp_complex_50");
   for (const mp_reference& r : mp_references)
   {
      const cpp_bin_float_50 x(r.x);
      const cpp_bin_float_50 c = fresnel_c(x);
      const cpp_bin_float_50 s = fresnel_s(x);
      CHECK_ULP_CLOSE(cpp_bin_float_50(r.c), c, 2);
      CHECK_ULP_CLOSE(cpp_bin_float_50(r.s), s, 2);
      const cpp_complex_50 z = fresnel_cis(x);
      CHECK_EQUAL(z.real(), c);
      CHECK_EQUAL(z.imag(), s);
   }
}
#endif

int main()
{
   test_references<float>(boost::math::policies::policy<>());
   test_references<float>(float_in_float());
   test_references<double>(boost::math::policies::policy<>());
   test_special_values<float>();
   test_special_values<double>();
   test_derivatives<float>();
   test_derivatives<double>();
#ifndef BOOST_MATH_STANDALONE
   test_multiprecision();
#endif
   return boost::math::test::report_errors();
}
