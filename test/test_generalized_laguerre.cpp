//  Copyright Nick Thompson 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "math_unit_test.hpp"
#include <boost/math/special_functions/laguerre.hpp>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <type_traits>
#ifndef BOOST_MATH_STANDALONE
#include <boost/multiprecision/cpp_bin_float.hpp>
#endif

using boost::math::laguerre;

// 40-digit references, parsed by the type under test when it has a string constructor.
template <class Real>
Real parse(const char* s, const std::true_type&)
{
   return static_cast<Real>(std::strtod(s, nullptr));
}

template <class Real>
Real parse(const char* s, const std::false_type&)
{
   return Real(s);
}

// L_n^alpha(x) = sum_k (-1)^k C(n + alpha, n - k) x^k / k!, summed in 600-bit arithmetic at the double arguments.
template <class Real>
void test_references()
{
   struct row { unsigned n; double alpha; double x; const char* value; };
   const row rows[] = {
      { 0, 0.5, 2, "1.0000000000000000000000000000000000000000e+00" },
      { 1, -0.5, 0.25, "2.5000000000000000000000000000000000000000e-01" },
      { 3, 0.5, 1.2, "-8.3049999999999996025401571841939358485888e-01" },
      { 5, -0.75, 3.5, "-1.0129313151041666666666666666666666666667e+00" },
      { 8, 2.5, 0.125, "6.5788231699666865761317904033358134920635e+01" },
      { 10, -1.5, 4, "1.8836377098446800595238095238095238095238e-01" },
      { 12, 0.3, 30, "-1.9232443892011784385158674644505192876150e+05" },
      { 20, 7.25, 11, "-7.1416476650801495331821635594845989966468e+01" },
      { 30, -0.9, 0.01, "-7.9721964323995015941239219126087524923842e-03" },
      { 40, 15.5, 60, "-8.7492063742854258334316742125927228113629e+10" },
      { 50, 0.5, 100, "-2.7346062282742768260267194147297534490136e+20" },
      { 100, 0.5, 7, "3.0727157599619286066044931268521883933253e+00" },
      { 100, -0.5, 250, "-1.3917259244440939207872843963045265233663e+53" },
      { 200, 3.75, 21, "-1.4004797311279646567157372238235365117328e+05" },
      { 500, 1.5, 0, "8.4419979753762600656249253474057168928388e+03" },
      { 7, 1e-10, 2, "6.6349206359134920634853968615656955087901e-01" },
      { 6, 1000, 5, "1.3763199107250937430555555555555555555556e+15" },
      { 25, 0.1, 1e-5, "1.4531232486766974586801754772798603081638e+00" },
   };
   const Real tol = 4 * (std::max)(std::numeric_limits<Real>::epsilon(), Real(1e-38));
   for (const row& r : rows)
   {
      // Float uses only the rows whose arguments it represents exactly and whose value it can hold.
      if (std::is_same<Real, float>::value && (static_cast<double>(static_cast<float>(r.alpha)) != r.alpha || static_cast<double>(static_cast<float>(r.x)) != r.x
          || std::fabs(std::strtod(r.value, nullptr)) > (std::numeric_limits<float>::max)()))
      {
         continue;
      }
      CHECK_RELATIVE_ERROR(parse<Real>(r.value, std::is_floating_point<Real>()), laguerre(r.n, static_cast<Real>(r.alpha), static_cast<Real>(r.x)), tol);
   }
}

template <class Real>
void test_identities()
{
   // An integral alpha is the associated polynomial, and alpha = 0 the ordinary one.
   for (unsigned n : { 0u, 1u, 5u, 40u })
   {
      for (Real x : { Real(0), Real(0.5), Real(3), Real(40) })
      {
         CHECK_EQUAL(laguerre(n, Real(3), x), laguerre(n, 3u, x));
         CHECK_EQUAL(laguerre(n, Real(3), x), laguerre(n, 3, x));
         CHECK_ULP_CLOSE(laguerre(n, x), laguerre(n, Real(0), x), 2);
      }
   }
   // A negative integer order is a real alpha, not an unsigned m.
   CHECK_EQUAL(laguerre(5u, Real(-1), Real(2.5)), laguerre(5u, -1, Real(2.5)));
   // L_n^alpha(0) = C(n + alpha, n).
   CHECK_ULP_CLOSE(Real(35) / Real(128), laguerre(4u, Real(-0.5), Real(0)), 2);
   // Recurrence in n with alpha: L_n^alpha = L_n^{alpha+1} - L_{n-1}^{alpha+1}.
   for (Real x : { Real(0.25), Real(2), Real(9) })
   {
      CHECK_ULP_CLOSE(laguerre(7u, Real(0.75), x), laguerre(7u, Real(1.75), x) - laguerre(6u, Real(1.75), x), 64);
   }
   // Result types: a floating alpha promotes with x; an integer order does not.
   static_assert(std::is_same<decltype(laguerre(3u, 0.5f, 1.0f)), float>::value, "");
   static_assert(std::is_same<decltype(laguerre(3u, 0.5, 1.0f)), double>::value, "");
   static_assert(std::is_same<decltype(laguerre(3u, 2, 1.0f)), float>::value, "");
   static_assert(std::is_same<decltype(laguerre(3u, 0.5f, 1.0f, boost::math::policies::policy<>())), float>::value, "");
}

int main()
{
   test_references<float>();
   test_references<double>();
   test_identities<float>();
   test_identities<double>();
#ifndef BOOST_MATH_STANDALONE
   test_references<boost::multiprecision::cpp_bin_float_50>();
#endif
   return boost::math::test::report_errors();
}
