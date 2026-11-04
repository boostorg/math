// Copyright Nick Thompson, 2026
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "math_unit_test.hpp"
#include <cmath>
#include <limits>
#include <random>
#include <stdexcept>
#include <type_traits>
#include <vector>
#include <boost/math/filters/wavelet_denoise.hpp>
#if !defined(BOOST_MATH_STANDALONE)
#include <boost/multiprecision/cpp_bin_float.hpp>
#endif

using boost::math::filters::orthogonal_wavelet;
using boost::math::filters::threshold_rule;
using boost::math::filters::wavelet_denoise;
using boost::math::filters::denoise;
using boost::math::filters::wavelet_transform;
using boost::math::filters::inverse_wavelet_transform;
namespace wavelet = boost::math::filters::wavelet;
namespace bmd = boost::math::filters::detail;

// db2..db10 (index 0..8), sym2..sym10 (9..17)
template <class Real>
orthogonal_wavelet<Real> family(unsigned t)
{
   return t < 9 ? wavelet::daubechies<Real>(t + 2) : wavelet::symlet<Real>(t - 7);
}

// standard normal samples by Box-Muller, independent of the standard library's normal_distribution
template <class Real>
std::vector<Real> noise(std::size_t n, unsigned seed)
{
   std::mt19937_64 gen(seed);
   std::uniform_real_distribution<double> u(0, 1);
   std::vector<Real> z(n);
   for (Real& v : z)
   {
      const double a = u(gen), b = u(gen);
      v = Real(std::sqrt(-2 * std::log(1 - a)) * std::cos(2 * 3.141592653589793 * b));
   }
   return z;
}

template <class Real>
Real mse(const std::vector<Real>& a, const std::vector<Real>& f)
{
   Real s = 0;
   for (std::size_t i = 0; i < f.size(); ++i) s += (a[i] - f[i]) * (a[i] - f[i]);
   return s / Real(f.size());
}

// Round trip and energy at n samples. The transform of float is computed in double, so its errors are
// rounding of the output only.
template <class Real>
void test_round_trip(unsigned t, std::size_t n)
{
   using std::abs;
   const Real eps = std::numeric_limits<Real>::epsilon();
   const orthogonal_wavelet<Real> w = family<Real>(t);
   const std::vector<Real> x = noise<Real>(n, t);
   std::vector<Real> y = x;
   y = wavelet_transform(y, w);
   Real ex = 0, ey = 0;
   for (std::size_t i = 0; i < n; ++i) { ex += x[i] * x[i]; ey += y[i] * y[i]; }
   CHECK_LE(Real(abs(ey - ex) / ex), Real(256 * eps)); // the sums themselves
   y = inverse_wavelet_transform(y, w);
   Real e = 0, xm = 0;
   for (std::size_t i = 0; i < n; ++i) { e = (std::max)(e, Real(abs(y[i] - x[i]))); xm = (std::max)(xm, Real(abs(x[i]))); }
   CHECK_LE(Real(e / xm), Real(10 * eps));
}

// |Q Q^T - I| for the matrix Q of the transform at n samples.
template <class Real>
void test_matrix(unsigned t, std::size_t n)
{
   using std::abs;
   const orthogonal_wavelet<Real> w = family<Real>(t);
   std::vector<std::vector<Real>> Q(n);
   for (std::size_t j = 0; j < n; ++j)
   {
      Q[j].assign(n, Real(0));
      Q[j][j] = 1;
      Q[j] = wavelet_transform(Q[j], w);
   }
   Real o = 0;
   for (std::size_t i = 0; i < n; ++i)
      for (std::size_t k = i; k < n; ++k)
      {
         Real s = 0;
         for (std::size_t j = 0; j < n; ++j) s += Q[i][j] * Q[k][j];
         o = (std::max)(o, Real(abs(s - Real(i == k ? 1 : 0))));
      }
   CHECK_LE(o, 10 * std::numeric_limits<Real>::epsilon());
}

// Polynomials of degree < p have no detail coefficients, at every level and at both edges.
template <class Real>
void test_annihilation(unsigned t, std::size_t n)
{
   using std::abs; using std::pow;
   const orthogonal_wavelet<Real> w = family<Real>(t);
   const unsigned p = w.vanishing_moments();
   const std::size_t m = bmd::wavelet_implementation(w)->plan(n).back().nc;
   for (unsigned k = 0; k < p; ++k)
   {
      std::vector<Real> x(n);
      for (std::size_t i = 0; i < n; ++i) x[i] = pow(Real(i) / Real(n), int(k));
      x = wavelet_transform(x, w);
      Real d = 0;
      for (std::size_t i = m; i < n; ++i) d = (std::max)(d, Real(abs(x[i])));
      CHECK_LE(d, 1000 * std::numeric_limits<Real>::epsilon());
   }
}

// A partial transform: levels = 0 is the identity, k levels give the full transform's finest k levels of
// details (bit for bit) and invert exactly; more levels than n allows is a domain error.
template <class Real>
void test_levels(unsigned t, std::size_t n)
{
   using std::abs;
   const orthogonal_wavelet<Real> w = family<Real>(t);
   const std::vector<Real> x = noise<Real>(n, 7);
   const std::size_t L = bmd::wavelet_levels(w, n);
   const std::vector<Real> full = wavelet_transform(x, w);
   const std::vector<Real> all = wavelet_transform(x, w, L);
   CHECK_TRUE(all == full);
   std::size_t m = n; // n_k, the length of the part still to transform after k levels
   for (std::size_t k = 0; k <= L; ++k, m = (m + 1) / 2)
   {
      const std::vector<Real> c = wavelet_transform(x, w, k);
      if (k == 0) CHECK_TRUE(c == x);
      bool same = true;
      for (std::size_t i = m; i < n; ++i) same = same && c[i] == full[i];
      CHECK_TRUE(same);
      const std::vector<Real> y = inverse_wavelet_transform(c, w, k);
      Real e = 0;
      for (std::size_t i = 0; i < n; ++i) e = (std::max)(e, Real(abs(y[i] - x[i])));
      CHECK_LE(e, 64 * std::numeric_limits<Real>::epsilon());
   }
   CHECK_THROW(wavelet_transform(x, w, L + 1), std::domain_error);
   CHECK_THROW(inverse_wavelet_transform(x, w, L + 1), std::domain_error);
}

// A smooth signal with a jump; every rule ends far below the noise level.
template <class Real>
std::vector<Real> test_signal(std::size_t n)
{
   using std::sin;
   std::vector<Real> f(n);
   for (std::size_t i = 0; i < n; ++i)
   {
      const Real s = Real(i) / Real(n);
      f[i] = sin(4 * boost::math::constants::pi<Real>() * s) + 2 * s * s + (s > Real(0.6) ? Real(0.5) : Real(0));
   }
   return f;
}

template <class Real>
void test_rules(std::size_t n)
{
   const Real sigma = Real(1) / 10;
   const std::vector<Real> f = test_signal<Real>(n), z = noise<Real>(n, 99);
   std::vector<Real> y(n);
   for (std::size_t i = 0; i < n; ++i) y[i] = f[i] + sigma * z[i];
   // MSE / sigma^2 measured in double: 0.02-0.05, universal soft 0.06-0.14
   for (unsigned t : {2u, 15u})
   {
      const orthogonal_wavelet<Real> w = family<Real>(t);
      for (threshold_rule r : {threshold_rule::eb_posterior_median, threshold_rule::sure, threshold_rule::wiener})
      {
         CHECK_LE(mse(wavelet_denoise(y, w, r), f), sigma * sigma / 10);
         CHECK_LE(mse(wavelet_denoise(y, w, r, sigma), f), sigma * sigma / 10);
      }
      CHECK_LE(mse(wavelet_denoise(y, w, threshold_rule::universal_soft), f), sigma * sigma / 4);
   }
   // the default rule is the empirical Wiener filter
   const orthogonal_wavelet<Real> sym8 = wavelet::symlet<Real>(8);
   const std::vector<Real> a = wavelet_denoise(y, sym8), b = wavelet_denoise(y, sym8, threshold_rule::wiener);
   for (std::size_t i = 0; i < n; ++i) CHECK_EQUAL(a[i], b[i]);
   // sigma = 0: nothing to remove
   const std::vector<Real> c = wavelet_denoise(y, sym8, threshold_rule::sure, Real(0));
   for (std::size_t i = 0; i < n; ++i) CHECK_ABSOLUTE_ERROR(y[i], c[i], Real(64 * std::numeric_limits<Real>::epsilon()));
}

// The posterior median of the empirical-Bayes rule against R's EbayesThresh:
// ebayesthresh(x, prior = "laplace", a = 0.5, sdev = 1) on 32 values drawn after set.seed(42).
// SURE includes the zero threshold: with every coefficient large, keeping them all has the least risk.
void test_sure_zero()
{
   const std::vector<double> d(8, 10.0);
   CHECK_EQUAL(bmd::sure_threshold(d.data(), d.size(), 1.0), 0.0);
}

void test_ebayesthresh()
{
   namespace eb = boost::math::filters::detail::eb;
   const std::vector<double> x{-0.306638594078474602, -1.781308433980000494, -0.171917355759621360, 1.214674699172598960, 1.895193461264965240,
      -4.244562562695575636, -0.257269382768929633, -1.763163085194780155, 0.460097354831271399, -0.639994875960119280, 0.455450123241219362,
      0.704837337228819139, 7.085931561189633854, -0.608926375407211395, 0.504955123297970210, -1.717008679073342980, -0.784459008379496336,
      -0.850907594176518400, -2.414207649946632284, 0.036122606892255625, 0.205998600200253873, -0.361057298548666195, 0.758163235699517046,
      -0.726704827076575133, -1.368281044419294545, 0.432818025888717040, -0.811393176186671616, 1.444101261721252749, -0.431446202613345431,
      0.655647883402206588, 0.321925265203946542, -0.783838940880375423};
   const double w = eb::weight(x.data(), x.size(), 0.5);
   CHECK_RELATIVE_ERROR(0.20391319282278733, w, 4 * std::numeric_limits<double>::epsilon());
   for (std::size_t i = 0; i < x.size(); ++i)
   {
      const double mu = eb::posterior_median(x[i], w, 0.5);
      const double expected = i == 5 ? -3.7375413914174218 : i == 12 ? 6.5859315582154556 : 0.0;
      CHECK_ABSOLUTE_ERROR(expected, mu, 1e-14);
   }
   // and the 8-value case: w = 0.37378033273295874, one coefficient survives
   const std::vector<double> x8{1.370958447146668480, -1.679853238665574455, -0.170156934237292889, 0.632862604961040520,
                                0.404268323140999086, -10.731946199710588274, 1.511521997438938891, -0.094659038413097557};
   const double w8 = eb::weight(x8.data(), x8.size(), 0.5);
   CHECK_RELATIVE_ERROR(0.37378033273295874, w8, 4 * std::numeric_limits<double>::epsilon());
   CHECK_ABSOLUTE_ERROR(-10.231946199710588, eb::posterior_median(x8[5], w8, 0.5), 1e-14);
   CHECK_EQUAL(eb::posterior_median(x8[0], w8, 0.5), 0.0);
}

// Cross-validation: the selected wavelet denoiser beats the noise by far; denoise() returns one of its
// two candidates, wavelet_denoise(y) or trend_filter(y), and the trend filter for piecewise constant data.
template <class Real>
void test_cross_validation(std::size_t n)
{
   const Real sigma = Real(1) / 10;
   const std::vector<Real> z = noise<Real>(n, 99);
   std::vector<Real> f = test_signal<Real>(n), y(n);
   for (std::size_t i = 0; i < n; ++i) y[i] = f[i] + sigma * z[i];
   const std::vector<Real> wcv = wavelet_denoise(y);
   CHECK_LE(mse(wcv, f), sigma * sigma / 20);
   // the result is the fixed-wavelet denoiser at the chosen wavelet and rule
   bmd::wavelet_status status;
   using T = bmd::trend_filter_compute_type<Real>;
   const std::vector<T> yt(y.begin(), y.end());
   const bmd::wavelet_choice<T> choice = bmd::wavelet_select(yt, status);
   CHECK_TRUE(status == bmd::wavelet_status::ok);
   const std::vector<Real> fixed = wavelet_denoise(y, family<Real>(choice.wavelet), choice.rule);
   for (std::size_t i = 0; i < n; ++i) CHECK_EQUAL(wcv[i], fixed[i]);
   const std::vector<Real> d = denoise(y);
   CHECK_LE(mse(d, f), sigma * sigma / 20);
   {
      const std::vector<Real> tf = boost::math::filters::trend_filter(y);
      bool is_wavelet = true, is_tf = true;
      for (std::size_t i = 0; i < n; ++i) { is_wavelet = is_wavelet && d[i] == wcv[i]; is_tf = is_tf && d[i] == tf[i]; }
      CHECK_TRUE(is_wavelet || is_tf);
   }

   for (std::size_t i = 0; i < n; ++i)
   {
      const Real s = Real(i) / Real(n);
      f[i] = Real(s > Real(0.1)) - 2 * Real(s > Real(0.3)) + Real(1.5) * Real(s > Real(0.55)) + Real(0.7) * Real(s > Real(0.8));
      y[i] = f[i] + sigma * z[i];
   }
   const std::vector<Real> blocks = denoise(y), tf = boost::math::filters::trend_filter(y);
   for (std::size_t i = 0; i < n; ++i) CHECK_EQUAL(blocks[i], tf[i]);
   CHECK_LE(mse(blocks, f), sigma * sigma / 20);
   CHECK_LE(mse(blocks, f), mse(wavelet_denoise(y), f) / 2);
}

// The right-edge rows of odd lengths live in a small LRU cache: 40 odd lengths, then the first ones again
// (evicted by then), all match a fresh wavelet; the cache holds the all-even entry and at most 16 others.
void test_row_cache()
{
   const orthogonal_wavelet<double> w = wavelet::symlet<double>(6);
   const auto* W = bmd::wavelet_implementation(w);
   std::vector<std::size_t> lengths;
   for (std::size_t k = 0; k < 40; ++k) lengths.push_back(101 + 34 * k);
   for (std::size_t k = 0; k < 5; ++k) lengths.push_back(lengths[k]);
   std::vector<double> big(4096, 1.0);
   big = wavelet_transform(big, w); // the all-even entry
   for (std::size_t n : lengths)
   {
      const orthogonal_wavelet<double> fresh = bmd::make_orthogonal_wavelet<double>(bmd::wavelet_scaling_filter(w));
      std::vector<double> a = noise<double>(n, unsigned(n)), b = a;
      a = wavelet_transform(a, w);
      b = wavelet_transform(b, fresh);
      bool same = true;
      for (std::size_t i = 0; i < n; ++i) same = same && a[i] == b[i];
      CHECK_TRUE(same);
      CHECK_LE(W->right_cache_size(), std::size_t(17));
   }
   CHECK_EQUAL(W->right_cache_size(), std::size_t(17));
}

// denoise() drops a trend filter that fails (here: the interior point at a lambda near 1e102) and returns
// the wavelet result.
void test_denoise_fallback()
{
   std::vector<double> y = noise<double>(256, 1);
   for (std::size_t i = 0; i < y.size(); ++i) y[i] = 1e100 * (std::sin(double(i) / 20) + y[i]);
   CHECK_THROW(boost::math::filters::trend_filter(y), boost::math::evaluation_error);
   std::vector<double> d;
   CHECK_NO_THROW(d = denoise(y));
   const std::vector<double> w = wavelet_denoise(y);
   for (std::size_t i = 0; i < y.size(); ++i) CHECK_EQUAL(d[i], w[i]);
}

template <class Real>
void test_errors()
{
   const Real nan = std::numeric_limits<Real>::quiet_NaN();
   CHECK_THROW(wavelet::daubechies<Real>(1), std::domain_error);
   CHECK_THROW(wavelet::daubechies<Real>(11), std::domain_error);
   CHECK_THROW(wavelet::symlet<Real>(1), std::domain_error);
   CHECK_THROW(wavelet::symlet<Real>(11), std::domain_error);
   CHECK_EQUAL(orthogonal_wavelet<Real>().vanishing_moments(), 0u);
   CHECK_EQUAL(wavelet::symlet<Real>(7).vanishing_moments(), 7u);
   const orthogonal_wavelet<Real> db4 = wavelet::daubechies<Real>(4), sym8 = wavelet::symlet<Real>(8);
   // db4 needs n >= 4p - 1 = 15
   CHECK_THROW(wavelet_denoise(std::vector<Real>(14, Real(1)), db4), std::domain_error);
   CHECK_NO_THROW(wavelet_denoise(std::vector<Real>(15, Real(1)), db4));
   CHECK_NO_THROW(wavelet_denoise(std::vector<Real>(17, Real(1)), db4));
   CHECK_EQUAL(bmd::wavelet_levels(db4, 14), std::size_t(0));
   CHECK_EQUAL(bmd::wavelet_levels(db4, 15), std::size_t(1));
   CHECK_EQUAL(bmd::wavelet_levels(db4, 16), std::size_t(1));
   CHECK_EQUAL(bmd::wavelet_levels(db4, 48), std::size_t(2));
   CHECK_EQUAL(bmd::wavelet_levels(db4, 59), std::size_t(3)); // 59 -> 30 -> 15 -> 8
   // the documented count: L = bit_width((n - 1) / (4p - 2))
   for (unsigned t = 0; t < 18; ++t)
   {
      const orthogonal_wavelet<Real> w = family<Real>(t);
      bool same = true;
      for (std::size_t n = 1; n <= 3000; ++n)
      {
         std::size_t q = (n - 1) / (4 * w.vanishing_moments() - 2), L = 0;
         for (; q; q >>= 1) ++L;
         same = same && L == bmd::wavelet_levels(w, n);
      }
      CHECK_TRUE(same);
   }
   // the Wiener rule works whenever the wavelet does: the pilot falls back from sym8 to db4, sym4, db2, or the wavelet itself
   for (unsigned t : {0u, 2u, 9u, 15u})
   {
      const orthogonal_wavelet<Real> w = family<Real>(t);
      const std::size_t nmin = 4 * w.vanishing_moments() - 1;
      CHECK_NO_THROW(wavelet_denoise(noise<Real>(nmin, 3), w, threshold_rule::wiener));
      CHECK_THROW(wavelet_denoise(noise<Real>(nmin - 1, 3), w, threshold_rule::wiener), std::domain_error);
   }
   CHECK_THROW(wavelet_denoise(std::vector<Real>(64, Real(1)), db4, threshold_rule::sure, Real(-1)), std::domain_error);
   CHECK_THROW(wavelet_denoise(std::vector<Real>(64, Real(1)), db4, threshold_rule::sure, std::numeric_limits<Real>::infinity()), std::domain_error);
   CHECK_NO_THROW(wavelet_denoise(std::vector<Real>(64, Real(1)), db4, threshold_rule::sure, nan));
   CHECK_THROW(wavelet_denoise(std::vector<Real>(64, Real(1)), orthogonal_wavelet<Real>()), std::domain_error);
   CHECK_THROW(wavelet_denoise(std::vector<Real>(64, Real(1)), db4, threshold_rule(17)), std::domain_error);
   std::vector<Real> bad = noise<Real>(64, 2);
   bad[30] = nan;
   CHECK_THROW(wavelet_denoise(bad, db4), std::domain_error);
   CHECK_THROW(wavelet_denoise(bad), std::domain_error);
   CHECK_THROW(denoise(bad), std::domain_error);
   bad[30] = std::numeric_limits<Real>::infinity();
   CHECK_THROW(wavelet_denoise(bad, db4, threshold_rule::sure), std::domain_error);
   std::vector<Real> z(64, Real(1));
   CHECK_THROW(wavelet_transform(z, orthogonal_wavelet<Real>()), std::domain_error);
   // the halves must suit db2 (7 samples each)
   CHECK_THROW(wavelet_denoise(std::vector<Real>(13, Real(1))), std::domain_error);
   CHECK_THROW(denoise(std::vector<Real>(13, Real(1))), std::domain_error);
   CHECK_NO_THROW(wavelet_denoise(std::vector<Real>(14, Real(1))));
   CHECK_NO_THROW(denoise(noise<Real>(15, 1)));
   std::vector<Real> x(14, Real(1));
   CHECK_THROW(wavelet_transform(x, db4), std::domain_error);
}

int main()
{
   test_errors<double>();
   test_errors<float>();
   for (unsigned t = 0; t < 18; ++t)
   {
      const std::size_t p = t < 9 ? t + 2 : t - 7;
      for (std::size_t n : {4 * p - 1, 4 * p, 4 * p + 1, 8 * p - 1, 16 * p - 3})
      {
         test_matrix<double>(t, n);
         test_annihilation<double>(t, n);
      }
      for (std::size_t n : {4 * p - 1, 4 * p, std::size_t(1000), std::size_t(1001), std::size_t(2047), std::size_t(3001), std::size_t(1) << 16})
         test_round_trip<double>(t, n);
      for (std::size_t n : {std::size_t(1000), std::size_t(2047), std::size_t(3001), std::size_t(4096)})
         test_annihilation<double>(t, n);
      test_round_trip<float>(t, 4096);
      test_round_trip<float>(t, 3001);
   }
   for (std::size_t n : {2048, 2047, 1001})
   {
      test_rules<double>(n);
      test_cross_validation<double>(n);
   }
   for (unsigned t : {0u, 2u, 15u})
      for (std::size_t n : {std::size_t(1001), std::size_t(1024), std::size_t(4 * (t < 9 ? t + 2 : t - 7) - 1)})
         test_levels<double>(t, n);
   test_levels<float>(11, 777);
   test_rules<float>(2048);
   test_ebayesthresh();
   test_sure_zero();
   test_denoise_fallback();
   test_row_cache();
   test_cross_validation<float>(2001);
#if !defined(BOOST_MATH_STANDALONE)
   using mp = boost::multiprecision::cpp_bin_float_50;
   for (unsigned t : {0u, 11u, 8u})
   {
      test_round_trip<mp>(t, 509);
      test_annihilation<mp>(t, 509);
   }
   test_matrix<mp>(11, 31);
   // expression templates on (cpp_bin_float_50 has them off)
   using mp_et = boost::multiprecision::number<boost::multiprecision::cpp_bin_float<50>, boost::multiprecision::et_on>;
   test_round_trip<mp_et>(15, 301);
#endif
   return boost::math::test::report_errors();
}
