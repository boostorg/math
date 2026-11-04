// Copyright Nick Thompson, 2026
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "math_unit_test.hpp"
#include <array>
#include <cmath>
#include <limits>
#include <utility>
#include <vector>
#include <boost/math/filters/daubechies.hpp>
#include <boost/math/filters/symlets.hpp>
#if !defined(BOOST_MATH_STANDALONE)
#include <boost/multiprecision/cpp_bin_float.hpp>
#endif

using boost::math::filters::daubechies_scaling_filter;
using boost::math::filters::symlet_scaling_filter;

// sum h_k = sqrt(2), sum h_k h_{k+2m} = delta_m, and p vanishing moments: sum (-1)^k k^d h_k = 0, d < p.
template <class Real, unsigned p>
void test_filter(const Real& tol)
{
   using std::sqrt; using std::abs; using std::pow;
   const std::array<Real, 2 * p> h = symlet_scaling_filter<Real, p>();
   Real s = 0;
   for (const Real& v : h) s += v;
   CHECK_ABSOLUTE_ERROR(sqrt(Real(2)), s, tol);
   for (unsigned m = 0; m < p; ++m)
   {
      Real c = 0;
      for (std::size_t k = 0; k + 2 * m < h.size(); ++k) c += h[k] * h[k + 2 * m];
      CHECK_ABSOLUTE_ERROR(Real(m == 0 ? 1 : 0), c, tol);
   }
   // moments of the wavelet filter, relative to sum |k^d h_k|
   for (unsigned d = 0; d < p; ++d)
   {
      Real c = 0, scale = 0;
      for (std::size_t k = 0; k < h.size(); ++k)
      {
         const Real t = pow(Real(k), int(d)) * h[k];
         c += (k % 2 ? -t : t);
         scale += abs(t);
      }
      CHECK_ABSOLUTE_ERROR(Real(0), c / scale, tol);
   }
}

template <class Real>
void test_filters(const Real& tol)
{
   test_filter<Real, 2>(tol);
   test_filter<Real, 3>(tol);
   test_filter<Real, 4>(tol);
   test_filter<Real, 5>(tol);
   test_filter<Real, 6>(tol);
   test_filter<Real, 7>(tol);
   test_filter<Real, 8>(tol);
   test_filter<Real, 9>(tol);
   test_filter<Real, 10>(tol);
}

// For p = 2, 3 every root choice is the minimum-phase one or its mirror image, so the symlet is the Daubechies filter.
template <class Real, unsigned p>
void test_equals_daubechies()
{
   const auto h = symlet_scaling_filter<Real, p>();
   const auto d = daubechies_scaling_filter<Real, p>();
   for (std::size_t k = 0; k < h.size(); ++k)
      CHECK_ULP_CLOSE(d[k], h[k], 2);
}

// The filters are not the minimum-phase ones (p >= 4), and they agree with PyWavelets' sym4 and sym8
// (double precision values, rec_lo, from pywt.Wavelet('sym4').rec_lo).
void test_pywavelets()
{
   const auto s4 = symlet_scaling_filter<double, 4>();
   const auto d4 = daubechies_scaling_filter<double, 4>();
   const double pywt4[] = {0.0322231006040427, -0.012603967262037833, -0.09921954357684722, 0.29785779560527736,
                           0.8037387518059161, 0.49761866763201545, -0.02963552764599851, -0.07576571478927333};
   double dd = 0;
   for (std::size_t k = 0; k < 8; ++k)
   {
      CHECK_ABSOLUTE_ERROR(pywt4[k], s4[k], 1e-11);
      dd = (std::max)(dd, std::abs(s4[k] - d4[k]));
   }
   CHECK_LE(0.1, dd);
}

// The least-asymmetric choice: sym7 and sym10 are not PyWavelets' (whose first coefficients are
// 0.010268176708511255 and -0.0004593294210046588), but more symmetric root choices.
void test_convention()
{
   const auto s7 = symlet_scaling_filter<double, 7>();
   const auto s10 = symlet_scaling_filter<double, 10>();
   CHECK_ABSOLUTE_ERROR(0.012015419283549189, s7[0], 1e-15);
   CHECK_ABSOLUTE_ERROR(-0.00041011591580439834, s10[0], 1e-15);
}

int main()
{
   test_filters<float>(8 * std::numeric_limits<float>::epsilon());
   test_filters<double>(16 * std::numeric_limits<double>::epsilon());
   test_equals_daubechies<double, 2>();
   test_equals_daubechies<double, 3>();
   test_pywavelets();
   test_convention();
#if !defined(BOOST_MATH_STANDALONE)
   // the tables have 75 significant digits
   using mp = boost::multiprecision::cpp_bin_float_100;
   test_filters<mp>(mp(1e-73));
   test_equals_daubechies<boost::multiprecision::cpp_bin_float_50, 2>();
   test_equals_daubechies<boost::multiprecision::cpp_bin_float_50, 3>();
#endif
   return boost::math::test::report_errors();
}
