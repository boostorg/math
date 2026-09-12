//  Copyright Nick Thompson 2026.
//  Distributed under the Boost Software License, Version 1.0.
//  (See accompanying file LICENSE_1_0.txt or copy at
//  https://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_TOOLS_POLYNOMIAL_ROOTS_HPP
#define BOOST_MATH_TOOLS_POLYNOMIAL_ROOTS_HPP

#include <boost/math/constants/constants.hpp>
#include <boost/math/tools/config.hpp>

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace boost { namespace math { namespace tools {
namespace polynomial_roots_detail {

template <class Real>
inline Real log2_abs(Real x)
{
   using std::abs;
   using std::log2;
   return x == 0 ? std::numeric_limits<Real>::infinity() : -log2(abs(x));
}

template <class Real>
struct hull_point
{
   std::size_t index;
   Real height;
};

template <class Real>
inline Real slope(hull_point<Real> const& a, hull_point<Real> const& b)
{
   return (b.height - a.height) / Real(b.index - a.index);
}

template <class Real>
std::vector<hull_point<Real> > lower_newton_hull(std::vector<Real> const& a)
{
   std::vector<hull_point<Real> > hull;
   hull.reserve(a.size());
   for(std::size_t j = 0; j < a.size(); ++j)
   {
      if(a[j] == 0)
         continue;
      hull_point<Real> p{j, log2_abs(a[j])};
      while(hull.size() >= 2 &&
            slope(hull[hull.size() - 2], hull.back()) >= slope(hull.back(), p))
         hull.pop_back();
      hull.push_back(p);
   }
   return hull;
}

template <class Real>
std::vector<Real> hull_heights(std::vector<hull_point<Real> > const& hull,
                               std::size_t degree)
{
   std::vector<Real> heights(degree + 1);
   std::size_t edge = 0;
   for(std::size_t j = 0; j <= degree; ++j)
   {
      while(edge + 1 < hull.size() && j > hull[edge + 1].index)
         ++edge;
      if(j == hull[edge].index || edge + 1 == hull.size())
         heights[j] = hull[edge].height;
      else
      {
         Real t = Real(j - hull[edge].index) /
                  Real(hull[edge + 1].index - hull[edge].index);
         heights[j] = hull[edge].height +
                      t * (hull[edge + 1].height - hull[edge].height);
      }
   }
   return heights;
}

template <class Real>
struct radial_piece
{
   Real log2_inner;
   Real log2_outer;
   std::size_t lower;
   std::size_t upper;
};

template <class Real>
std::vector<radial_piece<Real> > make_radial_pieces(std::vector<Real> const& a,
                                                    unsigned precision)
{
   std::size_t d = a.size() - 1;
   if(d == 1)
      return {{radial_piece<Real>{-std::numeric_limits<Real>::infinity(),
                                  std::numeric_limits<Real>::infinity(), 0, 1}}};
   auto hull = lower_newton_hull(a);
   auto heights = hull_heights(hull, d);
   Real m = static_cast<Real>(precision);
   std::vector<Real> lower_slope(d + 1), upper_slope(d + 1);
   upper_slope[0] = -std::numeric_limits<Real>::infinity();
   lower_slope[d] = std::numeric_limits<Real>::infinity();
   for(std::size_t j = 0; j <= d; ++j)
   {
      // The two supporting lines through (j, H(j)-m) touch the Newton
      // polygon on opposite sides of j.  Computing them from hull vertices
      // also handles absent coefficients without manufacturing finite
      // logarithms for them.
      Real shifted = heights[j] - m;
      if(j != 0)
      {
         Real high = -std::numeric_limits<Real>::infinity();
         for(auto const& p : hull)
            if(p.index < j)
               high = (std::max)(high, (shifted - p.height) /
                  Real(j - p.index));
         upper_slope[j] = high;
      }
      if(j != d)
      {
         Real low = std::numeric_limits<Real>::infinity();
         for(auto const& p : hull)
            if(p.index > j)
               low = (std::min)(low, (p.height - shifted) /
                  Real(p.index - j));
         lower_slope[j] = low;
      }
   }

   // Finite Fujiwara bounds avoid the degenerate start == end sweep for
   // sparse polynomials such as x^d-1.
   Real log_a0 = std::log2(std::abs(a[0]));
   Real log_ad = std::log2(std::abs(a[d]));
   Real start = std::numeric_limits<Real>::infinity();
   Real finish = -std::numeric_limits<Real>::infinity();
   for(std::size_t i = 1; i <= d; ++i)
   {
      if(a[i] != 0)
         start = (std::min)(start,
            (log_a0 - std::log2(std::abs(a[i]))) / Real(i));
      if(a[d - i] != 0)
         finish = (std::max)(finish,
            (std::log2(std::abs(a[d - i])) - log_ad) / Real(i));
   }
   start -= Real(1);
   finish += Real(1);
   std::vector<radial_piece<Real> > pieces;
   std::size_t first_upper = 1;
   while(first_upper <= d && start > upper_slope[first_upper])
      ++first_upper;
   first_upper = (std::min)(first_upper, d);
   pieces.push_back({-std::numeric_limits<Real>::infinity(), start, 0, first_upper});
   Real s = start;
   std::size_t previous_lower = 0;
   std::size_t previous_upper = first_upper;
   std::size_t guard = 0;
   while(s < finish && guard++ < 66 * d + 16)
   {
      std::size_t l = previous_lower;
      while(l < d && !(lower_slope[l] > s))
         ++l;
      std::size_t u = (std::max)(l, previous_upper);
      while(u < d && (upper_slope[u] - s) < m / Real(u - l + 1))
         ++u;
      Real next = s + m / Real(u - l + 1);
      if(!(next > s))
         next = std::nextafter(s, std::numeric_limits<Real>::infinity());
      pieces.push_back({s, next, l, u});
      previous_lower = l;
      previous_upper = u;
      s = next;
   }
   std::size_t last_lower = previous_lower;
   while(last_lower < d && !(lower_slope[last_lower] > s))
      ++last_lower;
   pieces.push_back({s, std::numeric_limits<Real>::infinity(), last_lower, d});
   return pieces;
}

template <class Real>
void fft(std::vector<std::complex<Real> >& a, bool inverse)
{
   std::size_t n = a.size();
   for(std::size_t i = 1, j = 0; i < n; ++i)
   {
      std::size_t bit = n >> 1;
      for(; j & bit; bit >>= 1)
         j ^= bit;
      j ^= bit;
      if(i < j)
         std::swap(a[i], a[j]);
   }
   Real pi = boost::math::constants::pi<Real>();
   for(std::size_t len = 2; len <= n; len <<= 1)
   {
      Real angle = (inverse ? Real(2) : Real(-2)) * pi / Real(len);
      std::complex<Real> step(std::cos(angle), std::sin(angle));
      for(std::size_t i = 0; i < n; i += len)
      {
         std::complex<Real> w(1, 0);
         for(std::size_t j = 0; j < len / 2; ++j)
         {
            auto u = a[i + j];
            auto v = a[i + j + len / 2] * w;
            a[i + j] = u + v;
            a[i + j + len / 2] = u - v;
            w *= step;
         }
      }
   }
   if(inverse)
      for(auto& z : a)
         z /= Real(n);
}

// Positive-sign DFT of arbitrary length, using Bluestein's reduction.
template <class Real>
std::vector<std::complex<Real> > positive_dft(std::vector<std::complex<Real> > const& x)
{
   std::size_t n = x.size();
   if(n <= 1)
      return x;
   std::size_t size = 1;
   while(size < 2 * n - 1)
      size <<= 1;
   std::vector<std::complex<Real> > a(size), b(size);
   Real pi = boost::math::constants::pi<Real>();
   for(std::size_t j = 0; j < n; ++j)
   {
      Real angle = pi * Real((j * j) % (2 * n)) / Real(n);
      std::complex<Real> positive(std::cos(angle), std::sin(angle));
      std::complex<Real> negative = std::conj(positive);
      a[j] = x[j] * positive;
      b[j] = negative;
      if(j != 0)
         b[size - j] = negative;
   }
   fft(a, false);
   fft(b, false);
   for(std::size_t j = 0; j < size; ++j)
      a[j] *= b[j];
   fft(a, true);
   a.resize(n);
   for(std::size_t j = 0; j < n; ++j)
   {
      Real angle = pi * Real((j * j) % (2 * n)) / Real(n);
      a[j] *= std::complex<Real>(std::cos(angle), std::sin(angle));
   }
   return a;
}

template <class Real>
std::pair<std::complex<Real>, std::complex<Real> >
evaluate_with_derivative(std::vector<std::complex<Real> > const& a,
                         std::complex<Real> z)
{
   std::complex<Real> f = a.back();
   std::complex<Real> df(0, 0);
   for(std::size_t i = a.size() - 1; i-- > 0;)
   {
      df = df * z + f;
      f = f * z + a[i];
   }
   return {f, df};
}

template <class Real>
std::vector<std::complex<Real> > aberth(std::vector<std::complex<Real> > a)
{
   while(a.size() > 1 && std::abs(a.back()) == 0)
      a.pop_back();
   std::size_t n = a.size() - 1;
   if(n == 0)
      return {};
   if(n == 1)
      return {-a[0] / a[1]};
   std::complex<Real> leading = a.back();
   Real log_leading = std::log(std::abs(leading));
   Real log_radius = -std::numeric_limits<Real>::infinity();
   for(std::size_t i = 0; i < n; ++i)
      if(a[i] != std::complex<Real>(0, 0))
         log_radius = (std::max)(log_radius,
            (std::log(std::abs(a[i])) - log_leading) / Real(n - i));
   // Fujiwara's bound is far tighter than 1+max|a_i/a_n| for a
   // translated high-degree polynomial, and computing it logarithmically
   // prevents the initial circle itself from overflowing.
   Real radius = Real(2) * std::exp(log_radius);
   if(!(radius > 0) || !std::isfinite(radius))
      radius = Real(1);
   Real pi = boost::math::constants::pi<Real>();
   std::vector<std::complex<Real> > roots(n);
   std::vector<hull_point<Real> > newton;
   for(std::size_t i = 0; i <= n; ++i)
      if(a[i] != std::complex<Real>(0, 0))
      {
         hull_point<Real> p{i, -std::log2(std::abs(a[i]))};
         while(newton.size() >= 2 &&
               slope(newton[newton.size() - 2], newton.back()) >=
               slope(newton.back(), p))
            newton.pop_back();
         newton.push_back(p);
      }
   std::size_t root_index = 0;
   for(std::size_t edge = 1; edge < newton.size(); ++edge)
   {
      std::size_t count = newton[edge].index - newton[edge - 1].index;
      Real edge_radius = std::exp2(slope(newton[edge - 1], newton[edge]));
      if(!(edge_radius > 0) || !std::isfinite(edge_radius))
         edge_radius = radius;
      for(std::size_t k = 0; k < count; ++k, ++root_index)
      {
         Real scramble = std::fmod(Real(root_index + 1) * Real(0.6180339887498948482), Real(1));
         Real angle = Real(2) * pi *
            (Real(k) + Real(0.25) + Real(edge) / Real(newton.size()) +
             Real(0.2) * scramble) / Real(count);
         Real radial_jitter = Real(0.95) + Real(0.1) * scramble;
         roots[root_index] = radial_jitter * edge_radius *
            std::complex<Real>(std::cos(angle), std::sin(angle));
      }
   }
   Real eps = std::numeric_limits<Real>::epsilon();
   for(unsigned iteration = 0; iteration < 20; ++iteration)
   {
      bool converged = true;
      for(std::size_t i = 0; i < n; ++i)
      {
         auto fd = evaluate_with_derivative(a, roots[i]);
         if(fd.second == std::complex<Real>(0, 0))
            continue;
         std::complex<Real> newton = fd.first / fd.second;
         std::complex<Real> repulsion(0, 0);
         for(std::size_t j = 0; j < n; ++j)
            if(i != j && roots[i] != roots[j])
               repulsion += Real(1) / (roots[i] - roots[j]);
         std::complex<Real> denominator = Real(1) - newton * repulsion;
         std::complex<Real> correction = denominator == std::complex<Real>(0, 0)
            ? newton : newton / denominator;
         roots[i] -= correction;
         if(std::abs(correction) > Real(64) * eps * (Real(1) + std::abs(roots[i])))
            converged = false;
      }
      if(converged)
         break;
   }
   return roots;
}

template <class Real>
std::complex<Real> root_of_unity(std::size_t k, std::size_t n)
{
   Real angle = Real(2) * boost::math::constants::pi<Real>() * Real(k) / Real(n);
   return {std::cos(angle), std::sin(angle)};
}

template <class Real>
Real scaled_residual(std::vector<Real> const& a, std::complex<Real> z)
{
   std::complex<Real> f = a.back();
   Real scale = std::abs(a.back());
   for(std::size_t i = a.size() - 1; i-- > 0;)
   {
      f = f * z + a[i];
      scale = scale * std::abs(z) + std::abs(a[i]);
   }
   return std::abs(f) / ((std::max)(scale, std::numeric_limits<Real>::min()));
}

} // namespace polynomial_roots_detail

BOOST_MATH_EXPORT template <class Real>
std::vector<std::complex<Real> > polynomial_roots(std::vector<Real> coefficients,
                                                  unsigned precision = std::numeric_limits<Real>::digits)
{
   static_assert(std::is_floating_point<Real>::value,
      "polynomial_roots requires a built-in floating-point coefficient type");
   while(!coefficients.empty() && coefficients.back() == 0)
      coefficients.pop_back();
   if(precision == 0)
      throw std::domain_error("polynomial_roots: precision must be positive");
   if(coefficients.empty())
      throw std::domain_error("polynomial_roots: the zero polynomial has no finite root set");
   for(auto coefficient : coefficients)
      if(!std::isfinite(coefficient))
         throw std::domain_error("polynomial_roots: coefficients must be finite");
   if(coefficients.size() == 1)
      return {};
   std::size_t zero_roots = 0;
   while(zero_roots + 1 < coefficients.size() && coefficients[zero_roots] == 0)
      ++zero_roots;
   if(zero_roots)
      coefficients.erase(coefficients.begin(), coefficients.begin() + zero_roots);

   using complex_type = std::complex<Real>;
   using namespace polynomial_roots_detail;
   std::size_t degree = coefficients.size() - 1;
   auto pieces = make_radial_pieces(coefficients, precision);
   std::vector<complex_type> candidates;
   std::size_t taylor_degree = (std::min)(degree, std::size_t(4) * precision);
   Real log_min = std::log2((std::numeric_limits<Real>::min)());
   Real log_max = std::log2((std::numeric_limits<Real>::max)());

   for(auto const& piece : pieces)
   {
      if(piece.upper <= piece.lower || !std::isfinite(piece.log2_inner) ||
         !std::isfinite(piece.log2_outer))
         continue;
      Real inner = std::exp2(piece.log2_inner);
      Real outer = std::exp2(piece.log2_outer);
      Real gamma = (inner + outer) / Real(2);
      Real rho = Real(0.75) * (outer - inner);
      if(!(gamma > 0) || !(rho > 0) || !std::isfinite(gamma + rho))
         continue;
      std::size_t delta = piece.upper - piece.lower;
      std::size_t sectors = static_cast<std::size_t>(std::ceil(
         Real(2) * boost::math::constants::pi<Real>() * gamma / rho));
      sectors = (std::max)(std::size_t(1), sectors);
      std::size_t local_degree = (std::min)(delta, taylor_degree);

      Real log_gamma = std::log2(gamma);
      Real common_log = -std::numeric_limits<Real>::infinity();
      for(std::size_t j = piece.lower; j <= piece.upper; ++j)
         if(coefficients[j] != 0)
            common_log = (std::max)(common_log,
               std::log2(std::abs(coefficients[j])) + Real(j) * log_gamma);

      std::vector<std::vector<complex_type> > local(
         sectors, std::vector<complex_type>(local_degree + 1));
      std::vector<std::vector<complex_type> > residues(
         local_degree + 1, std::vector<complex_type>(sectors));
      for(std::size_t exponent = piece.lower; exponent <= piece.upper; ++exponent)
      {
         if(coefficients[exponent] == 0)
            continue;
         Real log_term = std::log2(std::abs(coefficients[exponent]))
                       + Real(exponent) * log_gamma - common_log;
         if(log_term < log_min)
            continue;
         Real term = std::copysign(std::exp2((std::min)(log_term, log_max)),
                                   coefficients[exponent]);
         std::size_t limit = (std::min)(exponent, local_degree);
         for(std::size_t derivative = 0; derivative <= limit; ++derivative)
         {
            residues[derivative][exponent % sectors] += term;
            if(derivative != limit)
               term *= (Real(exponent - derivative) / Real(derivative + 1)) *
                       (rho / gamma);
         }
      }
      for(std::size_t derivative = 0; derivative <= local_degree; ++derivative)
      {
         auto values = positive_dft(residues[derivative]);
         for(std::size_t k = 0; k < sectors; ++k)
            local[k][derivative] = values[k];
      }
      for(std::size_t k = 0; k < sectors; ++k)
      {
         Real coefficient_norm = 0;
         for(auto const& coefficient : local[k])
            coefficient_norm += std::abs(coefficient);
         Real discarded = 0;
         std::size_t effective_size = local[k].size();
         Real discard_limit = Real(8) * std::numeric_limits<Real>::epsilon() * coefficient_norm;
         while(effective_size > 1 &&
               discarded + std::abs(local[k][effective_size - 1]) <= discard_limit)
         {
            discarded += std::abs(local[k][effective_size - 1]);
            --effective_size;
         }
         local[k].resize(effective_size);
         auto local_roots = aberth(local[k]);
         complex_type phase = root_of_unity<Real>(k, sectors);
         for(auto t : local_roots)
         {
            if(std::abs(t) > Real(1.05))
               continue;
            complex_type local_value = local[k].back();
            Real local_scale = std::abs(local[k].back());
            for(std::size_t q = local[k].size() - 1; q-- > 0;)
            {
               local_value = local_value * t + local[k][q];
               local_scale = local_scale * std::abs(t) + std::abs(local[k][q]);
            }
            if(!std::isfinite(local_value.real()) || !std::isfinite(local_value.imag()) ||
               std::abs(local_value) > std::sqrt(std::numeric_limits<Real>::epsilon()) *
                  (std::max)(local_scale, std::numeric_limits<Real>::min()))
               continue;
            complex_type z = phase * (gamma + rho * t);
            Real modulus = std::abs(z);
            if(modulus >= inner * (Real(1) - Real(32) * std::numeric_limits<Real>::epsilon()) &&
               modulus <= outer * (Real(1) + Real(32) * std::numeric_limits<Real>::epsilon()))
               candidates.push_back(z);
         }
      }
   }

   // Coalesce the deliberately overlapping rings and sectors before Horner
   // evaluation against the full input.  Otherwise validating O(d) local
   // representatives with an O(d) Horner evaluation would reintroduce a
   // quadratic step.
   std::sort(candidates.begin(), candidates.end(), [](complex_type a, complex_type b) {
      return std::arg(a) < std::arg(b) || (std::arg(a) == std::arg(b) && std::abs(a) < std::abs(b));
   });
   std::vector<complex_type> approximate;
   Real approximate_merge = Real(64) * std::sqrt(std::numeric_limits<Real>::epsilon());
   for(auto z : candidates)
   {
      if(!std::isfinite(z.real()) || !std::isfinite(z.imag()))
         continue;
      if(approximate.empty() ||
         std::abs(z - approximate.back()) > approximate_merge *
            (Real(1) + (std::max)(std::abs(z), std::abs(approximate.back()))))
         approximate.push_back(z);
   }
   candidates.swap(approximate);

   // Polish candidates against the original polynomial and keep the best
   // representative from overlapping sectors.
   std::vector<complex_type> complex_coefficients(coefficients.begin(), coefficients.end());
   for(auto& z : candidates)
      for(unsigned i = 0; i < 8; ++i)
      {
         auto fd = evaluate_with_derivative(complex_coefficients, z);
         if(fd.second == complex_type(0, 0))
            break;
         complex_type step = fd.first / fd.second;
         z -= step;
         if(std::abs(step) <= Real(16) * std::numeric_limits<Real>::epsilon() *
                              (Real(1) + std::abs(z)))
            break;
      }
   std::sort(candidates.begin(), candidates.end(), [](complex_type a, complex_type b) {
      return std::arg(a) < std::arg(b) || (std::arg(a) == std::arg(b) && std::abs(a) < std::abs(b));
   });
   std::vector<complex_type> roots;
   Real merge = std::sqrt(std::numeric_limits<Real>::epsilon());
   for(auto z : candidates)
   {
      // A root of a truncated local Taylor polynomial need not be a root of
      // the input polynomial.  Algorithm 3 validates local factors against
      // f and f'; scaled backward error is the floating-point analogue of
      // that acceptance test.
      Real residual = scaled_residual(coefficients, z);
      if(!std::isfinite(z.real()) || !std::isfinite(z.imag()) ||
         !std::isfinite(residual) ||
         residual > std::sqrt(std::numeric_limits<Real>::epsilon()))
         continue;
      auto duplicate = std::find_if(roots.begin(), roots.end(), [=](complex_type w) {
         return std::abs(z - w) <= merge * (Real(1) + (std::max)(std::abs(z), std::abs(w)));
      });
      if(duplicate == roots.end())
         roots.push_back(z);
      else if(scaled_residual(coefficients, z) < scaled_residual(coefficients, *duplicate))
         *duplicate = z;
   }

   // The piecewise algorithm intentionally resolves only roots meaningful at
   // the requested precision. Fall back to a global factorization when an
   // ill-conditioned cluster leaves fewer than degree candidates.
   if(roots.size() != degree)
      roots = aberth(complex_coefficients);
   roots.insert(roots.end(), zero_roots, complex_type(0, 0));
   return roots;
}

}}} // namespaces

#endif
