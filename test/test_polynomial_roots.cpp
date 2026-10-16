// Copyright Nick Thompson 2026.
// Distributed under the Boost Software License, Version 1.0.

#include "math_unit_test.hpp"
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/cpp_complex.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/complex_adaptor.hpp>
#include <boost/math/tools/polynomial_roots.hpp>
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <limits>
#include <random>
#include <stdexcept>
#include <vector>
#if defined(__cpp_lib_span) || (__cplusplus >= 202002L && __has_include(<span>))
#include <span>
#endif

namespace {

// |p(z)| / sum |p_i| |z|^i. For |z| > 1 the same ratio comes from the reversed polynomial in 1/z, so
// that neither sum overflows at high degree.
double backward_error(std::vector<double> const& p, std::complex<double> z)
{
   std::size_t n = p.size() - 1;
   bool reversed = std::abs(z) > 1;
   std::complex<double> x = reversed ? 1.0 / z : z;
   std::complex<double> value = reversed ? p[0] : p[n];
   double scale = std::abs(value);
   for(std::size_t i = 1; i <= n; ++i)
   {
      double c = reversed ? p[i] : p[n - i];
      value = value * x + c;
      scale = scale * std::abs(x) + std::abs(c);
   }
   return std::abs(value) / (std::max)(scale, (std::numeric_limits<double>::min)());
}

// The normwise backward error |p(z)| / (||p||_2 ||(1, z, ..., z^n)||_2), the quantity the companion QR keeps
// small; computed from the reversed polynomial when |z| > 1, as above.
double normwise_backward_error(std::vector<double> const& p, std::complex<double> z)
{
   std::size_t n = p.size() - 1;
   bool reversed = std::abs(z) > 1;
   std::complex<double> x = reversed ? 1.0 / z : z;
   std::complex<double> value = reversed ? p[0] : p[n];
   double powers = 1, norm = 0;
   for(std::size_t i = 1; i <= n; ++i)
   {
      value = value * x + (reversed ? p[i] : p[n - i]);
      powers = powers * std::norm(x) + 1;
   }
   for(double c : p)
      norm += c * c;
   return std::abs(value) / std::sqrt(norm * powers);
}

void check_residuals(std::vector<double> const& p,
                     std::vector<std::complex<double> > const& roots,
                     double tolerance = 1e-12)
{
   CHECK_EQUAL(roots.size() + 1, p.size());
   for(auto root : roots)
   {
      CHECK_TRUE(std::isfinite(root.real()));
      CHECK_TRUE(std::isfinite(root.imag()));
      CHECK_LE(backward_error(p, root), tolerance);
   }
}

std::vector<double> polynomial_from_real_roots(std::vector<double> const& roots)
{
   std::vector<double> p(1, 1);
   for(double root : roots)
   {
      std::vector<double> next(p.size() + 1);
      for(std::size_t i = 0; i < p.size(); ++i)
      {
         next[i] -= root * p[i];
         next[i + 1] += p[i];
      }
      p.swap(next);
   }
   return p;
}

// Both solvers with a multiprecision type, which returns the matching multiprecision complex type.
void test_multiprecision()
{
   using boost::multiprecision::cpp_bin_float_50;
   using boost::multiprecision::cpp_complex_50;
   using T = cpp_bin_float_50;
   T eps = std::numeric_limits<T>::epsilon();
   auto mp_backward_error = [](std::vector<T> const& p, cpp_complex_50 z)
   {
      std::size_t n = p.size() - 1;
      bool reversed = abs(z) > 1;
      cpp_complex_50 x = reversed ? cpp_complex_50(T(1) / z) : z;
      cpp_complex_50 value(reversed ? p[0] : p[n]);
      T scale = abs(value);
      for(std::size_t i = 1; i <= n; ++i)
      {
         T c = reversed ? p[i] : p[n - i];
         value = value * x + c;
         scale = scale * abs(x) + abs(c);
      }
      return T(abs(value) / scale);
   };

   // x^2 - 2: +-sqrt(2) to the last few bits.
   std::vector<T> two{-2, 0, 1};
   for(auto const& roots : {boost::math::tools::polynomial_roots(two), boost::math::tools::polynomial_roots_companion_matrix(two)})
   {
      CHECK_EQUAL(roots.size(), std::size_t(2));
      for(auto const& z : roots)
      {
         CHECK_LE(T(abs(abs(z.real()) - sqrt(T(2)))), 4 * eps * sqrt(T(2)));
         CHECK_LE(T(abs(z.imag())), 4 * eps);
      }
   }

   // A type with expression templates on, whose arithmetic yields expression objects rather than numbers.
   using D = boost::multiprecision::cpp_dec_float_50;
   std::vector<D> two_dec{-2, 0, 1};
   for(auto const& roots : {boost::math::tools::polynomial_roots(two_dec), boost::math::tools::polynomial_roots_companion_matrix(two_dec)})
   {
      CHECK_EQUAL(roots.size(), std::size_t(2));
      for(auto const& z : roots)
         CHECK_LE(D(abs(abs(z.real()) - sqrt(D(2)))), 16 * std::numeric_limits<D>::epsilon());
   }

   // Random coefficients and a repeated root: every root passes the same acceptance test as for double.
   std::mt19937_64 generator(0xbee5);
   std::normal_distribution<double> normal;
   std::vector<T> kac(31);
   for(auto& c : kac)
      c = normal(generator);
   std::vector<T> repeated{2, -3, 0, 1}; // (x - 1)^2 (x + 2)
   for(auto const& p : {kac, repeated})
   {
      auto roots = boost::math::tools::polynomial_roots(p);
      CHECK_EQUAL(roots.size() + 1, p.size());
      for(auto const& z : roots)
         CHECK_LE(mp_backward_error(p, z), T(16 * (p.size() - 1)) * eps);
   }

   // The companion QR is normwise backward stable: |p(z)| / (||p|| ||(1, z, ..., z^n)||).
   auto companion = boost::math::tools::polynomial_roots_companion_matrix(kac);
   CHECK_EQUAL(companion.size() + 1, kac.size());
   T norm = 0;
   for(auto const& c : kac)
      norm += c * c;
   for(auto const& z : companion)
   {
      bool reversed = abs(z) > 1;
      cpp_complex_50 x = reversed ? cpp_complex_50(T(1) / z) : z;
      cpp_complex_50 value(reversed ? kac[0] : kac.back());
      T powers = 1;
      for(std::size_t i = 1; i < kac.size(); ++i)
      {
         value = value * x + (reversed ? kac[i] : kac[kac.size() - 1 - i]);
         powers = powers * T(x.real() * x.real() + x.imag() * x.imag()) + 1;
      }
      CHECK_LE(T(abs(value) / sqrt(norm * powers)), T(64 * kac.size()) * eps);
   }
}

} // namespace

int main()
{
   using boost::math::tools::polynomial_roots;

   auto quadratic = polynomial_roots(std::vector<double>{1, 0, 1});
   check_residuals({1, 0, 1}, quadratic);
   CHECK_EQUAL(quadratic.size(), std::size_t(2));

   std::vector<double> with_zero_roots{0, 0, -2, 1};
   auto zero_roots = polynomial_roots(with_zero_roots);
   check_residuals(with_zero_roots, zero_roots);
   CHECK_EQUAL(std::count(zero_roots.begin(), zero_roots.end(), std::complex<double>(0, 0)), std::ptrdiff_t(2));

   std::mt19937_64 generator(0x5eed);
   std::uniform_real_distribution<double> distribution(-2, 2);
   for(unsigned degree = 2; degree <= 32; ++degree)
   {
      std::vector<double> known(degree);
      for(auto& root : known)
         root = distribution(generator);
      auto p = polynomial_from_real_roots(known);
      check_residuals(p, polynomial_roots(p));
   }

   // A sparse polynomial exercises the finite radial bounds needed when the
   // Newton polygon has a single edge.
   std::vector<double> roots_of_unity(257);
   roots_of_unity.front() = -1;
   roots_of_unity.back() = 1;
   check_residuals(roots_of_unity, polynomial_roots(roots_of_unity), 2e-8);

   // Random coefficients of high degree: an Aberth step can throw a root to |z| ~ 2.5, where Horner's
   // rule overflows at this degree and one non-finite root used to poison all the others.
   std::normal_distribution<double> normal;
   std::vector<double> kac(1001);
   for(auto& c : kac)
      c = normal(generator);
   check_residuals(kac, polynomial_roots(kac));

   // Wilkinson's polynomial of degree 15: every root found once, despite condition numbers ~1e10.
   std::vector<double> integers(15);
   for(std::size_t i = 0; i < integers.size(); ++i)
      integers[i] = double(i + 1);
   auto wilkinson = polynomial_roots(polynomial_from_real_roots(integers));
   CHECK_EQUAL(wilkinson.size(), integers.size());
   std::vector<int> found(integers.size() + 1, 0);
   for(auto z : wilkinson)
   {
      long k = std::lround(z.real());
      CHECK_LE(std::abs(z - std::complex<double>(double(k), 0)), 1e-4);
      if(k >= 1 && k <= long(integers.size()))
         ++found[k];
   }
   for(std::size_t k = 1; k <= integers.size(); ++k)
      CHECK_EQUAL(found[k], 1);

   // A close pair must not be merged into one root.
   auto close = polynomial_roots(polynomial_from_real_roots({1, 1 + 1e-5, -2}));
   std::vector<double> close_real;
   for(auto z : close)
      close_real.push_back(z.real());
   std::sort(close_real.begin(), close_real.end());
   CHECK_EQUAL(close.size(), std::size_t(3));
   CHECK_LE(std::abs(close_real[1] - 1), 1e-9);
   CHECK_LE(std::abs(close_real[2] - (1 + 1e-5)), 1e-9);

   // (z^2 + 1)^8: the piecewise pass leaves this cluster to the global fallback, whose roots must pass
   // the same backward-error test as the others.
   std::vector<double> cluster{1, 0, 8, 0, 28, 0, 56, 0, 70, 0, 56, 0, 28, 0, 8, 0, 1};
   check_residuals(cluster, polynomial_roots(cluster), 16 * 16 * std::numeric_limits<double>::epsilon());

   // Multiple roots of every kind, (z - c)^m, (z^2 + c)^m and (z^2 - z + c)^m: each root must pass the
   // acceptance test, which needs the fallback to stop on the residual and not on a noisy small step.
   for(int m = 2; m <= 10; ++m)
      for(double c : {0.25, 0.7, 1.0, 1.3, 3.0})
      {
         std::vector<double> factors[3] = {{-c, 1}, {c, 0, 1}, {c, -1, 1}};
         for(auto const& factor : factors)
         {
            std::vector<double> p(1, 1);
            for(int i = 0; i < m; ++i)
            {
               std::vector<double> next(p.size() + factor.size() - 1);
               for(std::size_t j = 0; j < p.size(); ++j)
                  for(std::size_t k = 0; k < factor.size(); ++k)
                     next[j + k] += p[j] * factor[k];
               p.swap(next);
            }
            check_residuals(p, polynomial_roots(p), 16 * double(p.size() - 1) * std::numeric_limits<double>::epsilon());
         }
      }

   // (z - 1)^4 + 2^-51 has no real roots, and its roots come in conjugate pairs; a cluster's inclusion
   // discs reach the real axis, so snapping on that alone made two of them real.
   std::vector<double> perturbed{1.0000000000000004, -4, 6, -4, 1};
   auto no_real = polynomial_roots(perturbed);
   check_residuals(perturbed, no_real);
   CHECK_EQUAL(std::count_if(no_real.begin(), no_real.end(), [](std::complex<double> z) { return z.imag() > 0; }), std::ptrdiff_t(2));
   CHECK_EQUAL(std::count_if(no_real.begin(), no_real.end(), [](std::complex<double> z) { return z.imag() < 0; }), std::ptrdiff_t(2));

   // Any random access container: (x - 1)(x - 2).
   std::array<double, 3> array_coefficients{2, -3, 1};
   for(auto const& roots : {polynomial_roots(array_coefficients), boost::math::tools::polynomial_roots_companion_matrix(array_coefficients)})
   {
      CHECK_EQUAL(roots.size(), std::size_t(2));
      for(auto z : roots)
         CHECK_LE((std::min)(std::abs(z - 1.0), std::abs(z - 2.0)), 16 * std::numeric_limits<double>::epsilon());
   }
#ifdef __cpp_lib_span
   std::span<const double> span_coefficients(array_coefficients);
   CHECK_EQUAL(polynomial_roots(span_coefficients).size(), std::size_t(2));
   CHECK_EQUAL(boost::math::tools::polynomial_roots_companion_matrix(span_coefficients).size(), std::size_t(2));
#endif

   CHECK_THROW(polynomial_roots(std::vector<double>{0, 0}), std::domain_error);
   CHECK_THROW(polynomial_roots(std::vector<double>{1, std::numeric_limits<double>::infinity()}), std::domain_error);

   // The companion QR: n roots, each with a small normwise backward error.
   using boost::math::tools::polynomial_roots_companion_matrix;
   auto check_companion = [](std::vector<double> const& p)
   {
      auto roots = polynomial_roots_companion_matrix(p);
      CHECK_EQUAL(roots.size() + 1, p.size());
      for(auto root : roots)
      {
         CHECK_TRUE(std::isfinite(root.real()));
         CHECK_TRUE(std::isfinite(root.imag()));
         CHECK_LE(normwise_backward_error(p, root), 64 * double(p.size()) * std::numeric_limits<double>::epsilon());
      }
      return roots;
   };
   check_companion({1, 0, 1});
   check_companion({-1, 1});
   auto companion_zero_roots = check_companion(with_zero_roots);
   CHECK_EQUAL(std::count(companion_zero_roots.begin(), companion_zero_roots.end(), std::complex<double>(0, 0)), std::ptrdiff_t(2));
   for(unsigned degree = 2; degree <= 32; ++degree)
   {
      std::vector<double> known(degree);
      for(auto& root : known)
         root = distribution(generator);
      check_companion(polynomial_from_real_roots(known));
   }
   // z^256 + 1: the trailing blocks are nilpotent, so only the exceptional shift makes progress.
   check_companion(roots_of_unity);
   check_companion(kac);
   check_companion(polynomial_from_real_roots(integers));
   // The scaling x = s w keeps the scaled coefficients finite when |a_0 / a_n|^(1/n) is far from 1.
   check_companion({1e300, 0, 0, 0, 1e-300});

   CHECK_THROW(polynomial_roots_companion_matrix(std::vector<double>{0, 0}), std::domain_error);
   CHECK_THROW(polynomial_roots_companion_matrix(std::vector<double>{1, std::numeric_limits<double>::quiet_NaN()}), std::domain_error);

   test_multiprecision();

   return boost::math::test::report_errors();
}
