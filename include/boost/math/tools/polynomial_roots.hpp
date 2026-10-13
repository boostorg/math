//  Copyright Nick Thompson 2026.
//  Distributed under the Boost Software License, Version 1.0.
//  (See accompanying file LICENSE_1_0.txt or copy at
//  https://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_TOOLS_POLYNOMIAL_ROOTS_HPP
#define BOOST_MATH_TOOLS_POLYNOMIAL_ROOTS_HPP

#include <boost/math/tools/config.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/tools/fft.hpp>
#include <boost/math/tools/precision.hpp>
#include <boost/math/tools/detail/complex_from_real.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/math/special_functions/next.hpp>
#include <boost/math/special_functions/sign.hpp>

#ifndef BOOST_MATH_BUILD_MODULE
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>
#endif

BOOST_MATH_NAMESPACE_BEGIN namespace tools {
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
   using std::abs;
   using std::log2;
   Real log_a0 = log2(abs(a[0]));
   Real log_ad = log2(abs(a[d]));
   Real start = std::numeric_limits<Real>::infinity();
   Real finish = -std::numeric_limits<Real>::infinity();
   for(std::size_t i = 1; i <= d; ++i)
   {
      if(a[i] != 0)
         start = (std::min)(start,
            (log_a0 - log2(abs(a[i]))) / Real(i));
      if(a[d - i] != 0)
         finish = (std::max)(finish,
            (log2(abs(a[d - i])) - log_ad) / Real(i));
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
         next = BOOST_MATH_NAMESPACE::float_next(s);
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

// The Newton correction p(z)/p'(z). For |z| > 1 Horner would overflow at high degree, so it uses the
// reversed polynomial q(w) = w^n p(1/w), w = 1/z: p/p' = z / (n - w q'(w) / q(w)). Returns zero where
// p'(z) vanishes and p(z) does not, so callers take no step there.
template <class Coefficient, class Complex>
Complex newton_correction(std::vector<Coefficient> const& a, Complex z)
{
   using std::abs;
   using complex_type = Complex;
   using Real = detail::real_from_complex_t<Complex>;
   std::size_t n = a.size() - 1;
   if(abs(z) <= 1)
   {
      complex_type f(a[n]), df(0, 0);
      for(std::size_t i = n; i-- > 0;)
      {
         df = df * z + f;
         f = f * z + complex_type(a[i]);
      }
      if(f == complex_type(0, 0))
         return f;
      return df == complex_type(0, 0) ? complex_type(0, 0) : f / df;
   }
   complex_type w = Real(1) / z;
   complex_type q(a[0]), dq(0, 0);
   for(std::size_t i = 1; i <= n; ++i)
   {
      dq = dq * w + q;
      q = q * w + complex_type(a[i]);
   }
   if(q == complex_type(0, 0))
      return q;
   complex_type denominator = Real(n) - w * dq / q;
   return denominator == complex_type(0, 0) ? complex_type(0, 0) : z / denominator;
}

// All roots of the polynomial a by Aberth's method. Only the roots inside |z| <= interest are wanted
// (the default wants all): a root still moving well outside that disc after a sweep is a root of a
// truncated Taylor polynomial that the caller will discard, so it is frozen rather than iterated.
// |p(z)| / sum |a_i| |z|^i, the backward error of z as a root. For |z| > 1 the same ratio is formed from
// the reversed polynomial in w = 1/z, so that neither sum overflows.
template <class Coefficient, class Complex>
detail::real_from_complex_t<Complex> scaled_residual(std::vector<Coefficient> const& a, Complex z)
{
   using std::abs;
   using Real = detail::real_from_complex_t<Complex>;
   std::size_t n = a.size() - 1;
   bool reversed = abs(z) > 1;
   Complex x = reversed ? Complex(Real(1) / z) : z;
   Real ax = abs(x);
   Complex f(reversed ? a[0] : a[n]);
   Real scale = abs(f);
   for(std::size_t i = 1; i <= n; ++i)
   {
      Coefficient c = reversed ? a[i] : a[n - i];
      f = f * x + c;
      scale = scale * ax + abs(c);
   }
   return abs(f) / ((std::max)(scale, tools::min_value<Real>()));
}

template <class Complex>
std::vector<Complex> aberth(std::vector<Complex> a,
                           detail::real_from_complex_t<Complex> interest = std::numeric_limits<detail::real_from_complex_t<Complex> >::infinity())
{
   using std::abs;
   using std::log;
   using std::log2;
   using std::exp;
   using std::exp2;
   using std::fmod;
   using std::cos;
   using std::sin;
   using Real = detail::real_from_complex_t<Complex>;
   while(a.size() > 1 && abs(a.back()) == 0)
      a.pop_back();
   std::size_t n = a.size() - 1;
   if(n == 0)
      return {};
   if(n == 1)
      return {-a[0] / a[1]};
   Complex leading = a.back();
   Real log_leading = log(abs(leading));
   Real log_radius = -std::numeric_limits<Real>::infinity();
   for(std::size_t i = 0; i < n; ++i)
      if(a[i] != Complex(0, 0))
         log_radius = (std::max)(log_radius,
            (log(abs(a[i])) - log_leading) / Real(n - i));
   // Fujiwara's bound is far tighter than 1+max|a_i/a_n| for a
   // translated high-degree polynomial, and computing it logarithmically
   // prevents the initial circle itself from overflowing.
   Real radius = Real(2) * exp(log_radius);
   if(!(radius > 0) || !(BOOST_MATH_NAMESPACE::isfinite)(radius))
      radius = Real(1);
   Real pi = BOOST_MATH_NAMESPACE::constants::pi<Real>();
   std::vector<Complex > roots(n);
   std::vector<hull_point<Real> > polygon;
   for(std::size_t i = 0; i <= n; ++i)
      if(a[i] != Complex(0, 0))
      {
         hull_point<Real> p{i, -log2(abs(a[i]))};
         while(polygon.size() >= 2 &&
               slope(polygon[polygon.size() - 2], polygon.back()) >=
               slope(polygon.back(), p))
            polygon.pop_back();
         polygon.push_back(p);
      }
   std::size_t root_index = 0;
   for(std::size_t edge = 1; edge < polygon.size(); ++edge)
   {
      std::size_t count = polygon[edge].index - polygon[edge - 1].index;
      Real edge_radius = exp2(slope(polygon[edge - 1], polygon[edge]));
      if(!(edge_radius > 0) || !(BOOST_MATH_NAMESPACE::isfinite)(edge_radius))
         edge_radius = radius;
      for(std::size_t k = 0; k < count; ++k, ++root_index)
      {
         Real scramble = fmod(Real(root_index + 1) * Real(0.6180339887498948482), Real(1));
         Real angle = Real(2) * pi *
            (Real(k) + Real(0.25) + Real(edge) / Real(polygon.size()) +
             Real(0.2) * scramble) / Real(count);
         Real radial_jitter = Real(0.95) + Real(0.1) * scramble;
         roots[root_index] = radial_jitter * edge_radius *
            Complex(cos(angle), sin(angle));
      }
   }
   Real eps = tools::epsilon<Real>();
   // A root whose last correction was negligible is frozen: it still repels the others but is no longer
   // updated. Each sweep then costs O(n * active) instead of O(n^2).
   const bool global = !(BOOST_MATH_NAMESPACE::isfinite)(interest);
   const Real acceptable = Real(16) * Real(n) * eps;
   std::vector<Complex > last_acceptable(global ? n : 0);
   std::vector<char> has_acceptable(global ? n : 0, 0);
   std::vector<char> active(n, 1);
   std::size_t remaining = n;
   // Simple roots converge in a few sweeps; a cluster of multiplicity m only linearly, at a rate of about
   // 1 - 1/m. The global fallback (no disc of interest) is what such clusters reach, so it gets a far
   // larger budget, which costs little because converged roots are frozen.
   unsigned max_iterations = 20 + static_cast<unsigned>(log2(Real(n)));
   if(!(BOOST_MATH_NAMESPACE::isfinite)(interest))
      max_iterations = 400 + static_cast<unsigned>(log2(Real(n)));
   for(unsigned iteration = 0; iteration < max_iterations && remaining; ++iteration)
   {
      for(std::size_t i = 0; i < n; ++i)
      {
         if(!active[i])
            continue;
         Complex newton = newton_correction(a, roots[i]);
         // sum 1/(z_i - z_j) in real arithmetic: 1/d = conj(d)/|d|^2, without the library's slow
         // general-purpose complex division.
         Real repulsion_re = 0, repulsion_im = 0;
         Real xi = roots[i].real(), yi = roots[i].imag();
         for(std::size_t j = 0; j < n; ++j)
         {
            Real dx = xi - roots[j].real(), dy = yi - roots[j].imag();
            Real norm = dx * dx + dy * dy;
            if(j != i && norm > 0)
            {
               repulsion_re += dx / norm;
               repulsion_im -= dy / norm;
            }
         }
         Complex repulsion(repulsion_re, repulsion_im);
         Complex denominator = Real(1) - newton * repulsion;
         Complex correction = denominator == Complex(0, 0)
            ? newton : newton / denominator;
         // One non-finite root would poison every other root through the repulsion sum.
         if(!(BOOST_MATH_NAMESPACE::isfinite)(correction.real()) || !(BOOST_MATH_NAMESPACE::isfinite)(correction.imag()))
            correction = newton;
         if(!(BOOST_MATH_NAMESPACE::isfinite)(correction.real()) || !(BOOST_MATH_NAMESPACE::isfinite)(correction.imag()))
            continue;
         roots[i] -= correction;
         // A root stops once its step is negligible. In the global fallback, which is where multiple roots
         // end up, a small step can be noise long before z is acceptable, so there it must also pass the
         // backward-error test the caller accepts roots by; the last iterate that passed is kept in case the
         // noise carries the root out again.
         bool converged = abs(correction) <= Real(4) * eps * abs(roots[i]);
         if(global)
         {
            bool acceptable_now = scaled_residual(a, roots[i]) <= acceptable;
            if(acceptable_now)
            {
               last_acceptable[i] = roots[i];
               has_acceptable[i] = 1;
            }
            converged = converged && acceptable_now;
         }
         if(converged || (iteration >= 2 && abs(roots[i]) > Real(2) * interest))
         {
            active[i] = 0;
            --remaining;
         }
      }
      // Done once every root still moving lies well outside the disc of interest.
      if(iteration >= 2 && (BOOST_MATH_NAMESPACE::isfinite)(interest))
      {
         bool wanted_moving = false;
         for(std::size_t i = 0; i < n && !wanted_moving; ++i)
            wanted_moving = active[i] && abs(roots[i]) <= Real(1.25) * interest;
         if(!wanted_moving)
            break;
      }
   }
   if(global)
      for(std::size_t i = 0; i < n; ++i)
         if(has_acceptable[i] && !(scaled_residual(a, roots[i]) <= acceptable))
            roots[i] = last_acceptable[i];
   return roots;
}

template <class Real>
detail::complex_from_real_t<Real> root_of_unity(std::size_t k, std::size_t n)
{
   using std::cos;
   using std::sin;
   Real angle = Real(2) * BOOST_MATH_NAMESPACE::constants::pi<Real>() * Real(k) / Real(n);
   return detail::complex_from_real_t<Real>(cos(angle), sin(angle));
}


// The radius n max(|p(z)|, eps sum |a_i| |z|^i) / |p'(z)| of a disc about z that contains a root of p:
// the Newton inclusion radius, floored at the rounding error of evaluating p, which is eps cond(z) |z|.
// Copies of one root found in overlapping sectors have overlapping discs; distinct roots have disjoint
// ones unless they are closer than the data can resolve. For |z| > 1 it uses the reversed polynomial
// q(w) = w^n p(1/w), w = 1/z, where p'(z) = z^(n-1) (n q - w q').
template <class Real, class Complex>
Real inclusion_radius(std::vector<Real> const& a, Complex z)
{
   using std::abs;
   using complex_type = Complex;
   std::size_t n = a.size() - 1;
   Real eps = tools::epsilon<Real>();
   bool reversed = abs(z) > 1;
   complex_type x = reversed ? complex_type(Real(1) / z) : z;
   Real ax = abs(x);
   complex_type f(reversed ? a[0] : a[n]), df(0, 0);
   Real scale = abs(f);
   for(std::size_t i = 1; i <= n; ++i)
   {
      Real c = reversed ? a[i] : a[n - i];
      df = df * x + f;
      f = f * x + c;
      scale = scale * ax + abs(c);
   }
   // In the reversed case f, df and scale are q(w), q'(w) and |z|^-n sum |a_i||z|^i, and |p'(z)| is
   // |z|^(n-1) |n q - w q'|, so the radius is |z| max(|q|, eps scale) / |n q - w q'|.
   Real derivative = reversed ? abs(Real(n) * f - x * df) : abs(df);
   Real numerator = (std::max)(Real(abs(f)), Real(eps * scale)) * (reversed ? Real(abs(z)) : Real(1));
   if(!(derivative > 0))
      return std::numeric_limits<Real>::infinity();
   return Real(n) * numerator / derivative;
}

} // namespace polynomial_roots_detail

namespace polynomial_roots_detail {

// The piecewise method, with local Taylor expansions to the full precision of Real.
template <class Real>
std::vector<detail::complex_from_real_t<Real> > piecewise(std::vector<Real> coefficients)
{
   const unsigned precision = static_cast<unsigned>(tools::digits<Real>());
   static_assert(!std::is_integral<Real>::value, "polynomial_roots requires a floating-point coefficient type");
   using std::abs;
   using std::log2;
   using std::exp2;
   using std::ceil;
   using std::sqrt;
   using std::conj;
   while(!coefficients.empty() && coefficients.back() == 0)
      coefficients.pop_back();
   if(coefficients.empty())
      throw std::domain_error("polynomial_roots: the zero polynomial has no finite root set");
   for(auto coefficient : coefficients)
      if(!(BOOST_MATH_NAMESPACE::isfinite)(coefficient))
         throw std::domain_error("polynomial_roots: coefficients must be finite");
   if(coefficients.size() == 1)
      return {};
   std::size_t zero_roots = 0;
   while(zero_roots + 1 < coefficients.size() && coefficients[zero_roots] == 0)
      ++zero_roots;
   if(zero_roots)
      coefficients.erase(coefficients.begin(), coefficients.begin() + zero_roots);

   using complex_type = detail::complex_from_real_t<Real>;
   using namespace polynomial_roots_detail;
   std::size_t degree = coefficients.size() - 1;
   auto pieces = make_radial_pieces(coefficients, precision);
   // A local root mapped back to z, with a merge radius (a few Newton steps on the local Taylor polynomial,
   // floored at its rounding and truncation error) and its local backward error.
   struct candidate { complex_type z; Real radius; Real residual; };
   std::vector<candidate> candidates;
   std::size_t taylor_degree = (std::min)(degree, std::size_t(4) * precision);
   Real log_min = log2(tools::min_value<Real>());
   Real log_max = log2(tools::max_value<Real>());

   for(auto const& piece : pieces)
   {
      if(piece.upper <= piece.lower || !(BOOST_MATH_NAMESPACE::isfinite)(piece.log2_inner) ||
         !(BOOST_MATH_NAMESPACE::isfinite)(piece.log2_outer))
         continue;
      Real inner = exp2(piece.log2_inner);
      Real outer = exp2(piece.log2_outer);
      Real gamma = (inner + outer) / Real(2);
      Real rho = Real(0.75) * (outer - inner);
      if(!(gamma > 0) || !(rho > 0) || !(BOOST_MATH_NAMESPACE::isfinite)(gamma + rho))
         continue;
      std::size_t delta = piece.upper - piece.lower;
      std::size_t sectors = static_cast<std::size_t>(ceil(
         Real(2) * BOOST_MATH_NAMESPACE::constants::pi<Real>() * gamma / rho));
      sectors = (std::max)(std::size_t(1), sectors);
      std::size_t local_degree = (std::min)(delta, taylor_degree);

      Real log_gamma = log2(gamma);
      Real common_log = -std::numeric_limits<Real>::infinity();
      for(std::size_t j = piece.lower; j <= piece.upper; ++j)
         if(coefficients[j] != 0)
            common_log = (std::max)(common_log,
               log2(abs(coefficients[j])) + Real(j) * log_gamma);

      std::vector<std::vector<complex_type> > local(
         sectors, std::vector<complex_type>(local_degree + 1));
      std::vector<std::vector<complex_type> > residues(
         local_degree + 1, std::vector<complex_type>(sectors));
      for(std::size_t exponent = piece.lower; exponent <= piece.upper; ++exponent)
      {
         if(coefficients[exponent] == 0)
            continue;
         Real log_term = log2(abs(coefficients[exponent]))
                       + Real(exponent) * log_gamma - common_log;
         if(log_term < log_min)
            continue;
         Real term = BOOST_MATH_NAMESPACE::copysign(Real(exp2((std::min)(log_term, log_max))),
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
      // Positive-sign DFT of each derivative's residues, as conj(F(conj(x))).
      fft<complex_type> dft(sectors);
      for(std::size_t derivative = 0; derivative <= local_degree; ++derivative)
      {
         auto& values = residues[derivative];
         for(auto& v : values)
            v = conj(v);
         dft.transform(values);
         for(std::size_t k = 0; k < sectors; ++k)
            local[k][derivative] = conj(values[k]);
      }
      for(std::size_t k = 0; k < sectors; ++k)
      {
         Real coefficient_norm = 0;
         for(auto const& coefficient : local[k])
            coefficient_norm += abs(coefficient);
         Real discarded = 0;
         std::size_t effective_size = local[k].size();
         Real discard_limit = Real(8) * tools::epsilon<Real>() * coefficient_norm;
         while(effective_size > 1 &&
               discarded + abs(local[k][effective_size - 1]) <= discard_limit)
         {
            discarded += abs(local[k][effective_size - 1]);
            --effective_size;
         }
         local[k].resize(effective_size);
         auto local_roots = aberth(local[k], Real(1.05));
         complex_type phase = root_of_unity<Real>(k, sectors);
         for(auto t : local_roots)
         {
            if(abs(t) > Real(1.05))
               continue;
            complex_type local_value = local[k].back(), local_derivative(0, 0);
            Real local_scale = abs(local[k].back());
            for(std::size_t q = local[k].size() - 1; q-- > 0;)
            {
               local_derivative = local_derivative * t + local_value;
               local_value = local_value * t + local[k][q];
               local_scale = local_scale * abs(t) + abs(local[k][q]);
            }
            local_scale = (std::max)(local_scale, tools::min_value<Real>());
            if(!(BOOST_MATH_NAMESPACE::isfinite)(local_value.real()) || !(BOOST_MATH_NAMESPACE::isfinite)(local_value.imag()) ||
               abs(local_value) > sqrt(tools::epsilon<Real>()) * local_scale)
               continue;
            complex_type z = phase * (gamma + rho * t);
            Real modulus = abs(z);
            if(modulus >= inner * (Real(1) - Real(32) * tools::epsilon<Real>()) &&
               modulus <= outer * (Real(1) + Real(32) * tools::epsilon<Real>()))
            {
               Real error = abs(local_value) + Real(8) * tools::epsilon<Real>() * local_scale;
               // A few Newton steps' worth: copies of one root agree to about the local accuracy, and the
               // guaranteed radius (degree times the step) would also merge neighbouring roots. Copies
               // this misses are merged after polishing, with the guaranteed radius on the full polynomial.
               Real radius = abs(local_derivative) > 0
                  ? Real(4) * error / abs(local_derivative) * rho
                  : std::numeric_limits<Real>::infinity();
               candidates.push_back({z, radius, abs(local_value) / local_scale});
            }
         }
      }
   }

   // Overlapping rings and sectors find most roots many times. Merge two copies when each lies within the
   // other's radius (sorted by real part, so there is no branch cut to split a pair), keeping the one with
   // the smaller local backward error, so that each root is polished against the full polynomial only
   // once. Using the smaller radius keeps a poor candidate from absorbing a good one for a different root.
   Real eps = tools::epsilon<Real>();
   auto by_real = [](complex_type a, complex_type b) { return a.real() < b.real() || (a.real() == b.real() && a.imag() < b.imag()); };
   candidates.erase(std::remove_if(candidates.begin(), candidates.end(), [](candidate const& c) {
      return !(BOOST_MATH_NAMESPACE::isfinite)(c.z.real()) || !(BOOST_MATH_NAMESPACE::isfinite)(c.z.imag()) || !(c.radius < std::numeric_limits<Real>::infinity()); }),
      candidates.end());
   std::sort(candidates.begin(), candidates.end(), [&](candidate const& a, candidate const& b) { return by_real(a.z, b.z); });
   std::vector<complex_type> representatives;
   {
      std::vector<candidate> kept;
      for(auto const& c : candidates)
      {
         bool copy = false;
         for(std::size_t k = kept.size(); k-- > 0 && c.z.real() - kept[k].z.real() <= c.radius;)
            if(abs(c.z - kept[k].z) <= (std::min)(c.radius, kept[k].radius))
            {
               if(c.residual < kept[k].residual)
                  kept[k] = c;
               copy = true;
               break;
            }
         if(!copy)
            kept.push_back(c);
      }
      for(auto const& c : kept)
         representatives.push_back(c.z);
   }
   std::vector<complex_type>& polished = representatives;

   // Polish against the input polynomial, while Newton's steps keep shrinking: near a close neighbour
   // a rough local estimate can need more than a handful of iterations. A candidate on which Newton is
   // only converging linearly after a few steps is not near a simple root (it is a spurious root of a
   // truncated Taylor polynomial); it is left for the residual test below to reject.
   auto polish = [&](complex_type& z)
   {
      Real previous = std::numeric_limits<Real>::infinity();
      for(unsigned i = 0; i < 64; ++i)
      {
         complex_type step = newton_correction(coefficients, z);
         Real size = abs(step);
         if(!(size < previous))
            break;
         z -= step;
         if(size <= Real(4) * eps * abs(z) || (i >= 8 && size > previous / 2))
            break;
         previous = size;
      }
   };
   for(auto& z : polished)
      polish(z);

   // A root of a truncated local Taylor polynomial need not be a root of the input polynomial: Algorithm 3
   // validates local factors against f and f', and scaled backward error is the floating-point analogue.
   // A polished root has a backward error of a few eps even in a cluster; a larger one means Newton has
   // not converged, and such a candidate can sit a long way from the root it was headed for.
   // Copies of one root are then merged when their Newton inclusion discs overlap, which separates
   // distinct roots as finely as their conditioning allows.
   struct accepted_root { complex_type z; Real radius; Real residual; };
   std::vector<accepted_root> found;
   for(auto z : polished)
   {
      Real residual = scaled_residual(coefficients, z);
      if((BOOST_MATH_NAMESPACE::isfinite)(z.real()) && (BOOST_MATH_NAMESPACE::isfinite)(z.imag()) && (BOOST_MATH_NAMESPACE::isfinite)(residual) &&
         residual <= Real(16) * Real(degree) * eps)
         found.push_back({z, inclusion_radius(coefficients, z), residual});
   }
   std::sort(found.begin(), found.end(), [&](accepted_root const& a, accepted_root const& b) { return by_real(a.z, b.z); });
   std::vector<accepted_root> unique;
   Real largest_radius = 0;
   for(auto const& c : found)
   {
      bool copy = false;
      for(std::size_t k = unique.size(); k-- > 0 &&
          c.z.real() - unique[k].z.real() <= c.radius + largest_radius;)
         if(abs(c.z - unique[k].z) <= c.radius + unique[k].radius)
         {
            if(c.residual < unique[k].residual)
               unique[k] = c;
            copy = true;
            break;
         }
      if(!copy)
      {
         unique.push_back(c);
         largest_radius = (std::max)(largest_radius, c.radius);
      }
   }
   std::vector<complex_type> roots;
   roots.reserve(degree + zero_roots);
   for(auto const& c : unique)
      roots.push_back(c.z);

   // The piecewise algorithm resolves only roots meaningful at the precision of Real. Fall back to a
   // global factorization when an ill-conditioned cluster leaves a different number of roots than the degree.
   if(roots.size() != degree)
   {
      roots = aberth(std::vector<complex_type>(coefficients.begin(), coefficients.end()));
      // Polish each root, but keep the result only if it is still acceptable or no worse: near a multiple
      // root Newton's steps are noise and can move an acceptable root to an unacceptable place. (Comparing
      // residuals alone would discard good polishing of simple roots, whose residuals are both rounding.)
      const Real acceptable = Real(16) * Real(degree) * eps;
      for(auto& z : roots)
      {
         complex_type candidate = z;
         polish(candidate);
         Real after = scaled_residual(coefficients, candidate);
         if(after <= acceptable || after <= scaled_residual(coefficients, z))
            z = candidate;
      }
   }

   // The coefficients are real, so the roots are closed under conjugation: an inclusion disc that reaches
   // the real axis and overlaps no other root's disc can't hold a non-real root without its conjugate, so
   // that root is real. In a cluster the discs overlap and nothing is snapped: a complex pair there would
   // otherwise become two reals.
   std::vector<Real> radii(roots.size());
   for(std::size_t k = 0; k < roots.size(); ++k)
      radii[k] = inclusion_radius(coefficients, roots[k]);
   std::vector<complex_type> snapped = roots;
   for(std::size_t k = 0; k < roots.size(); ++k)
   {
      complex_type z = roots[k];
      if(z.imag() == 0 || !(abs(z.imag()) <= radii[k]))
         continue;
      // An imaginary part below the rounding of the real part is noise; a conjugate pair that close
      // snaps together.
      if(abs(z.imag()) <= eps * abs(z.real()))
      {
         snapped[k] = complex_type(z.real(), Real(0));
         continue;
      }
      bool isolated = true;
      for(std::size_t j = 0; j < roots.size() && isolated; ++j)
         if(j != k && abs(roots[j] - z) <= radii[j] + radii[k])
            isolated = false;
      if(isolated && scaled_residual(coefficients, complex_type(z.real(), Real(0))) <= scaled_residual(coefficients, z))
         snapped[k] = complex_type(z.real(), Real(0));
   }
   roots.swap(snapped);
   roots.insert(roots.end(), zero_roots, complex_type(0, 0));
   return roots;
}

} // namespace polynomial_roots_detail


// Roots of a polynomial as the eigenvalues of its companion matrix, by the companion QR algorithm of
// Aurentz, Mach, Vandebril and Watkins ("Fast and backward stable computation of roots of polynomials",
// SIAM J. Matrix Anal. Appl. 36 (2015); Part II with Robol, SIMAX 39 (2018)): O(n^2) time, O(n) memory,
// backward stable. This is the complex single-shift variant.
//
// The monic companion matrix is A = QR, Q = Q_1 ... Q_{n-1} a product of core transformations (2x2 unitary
// matrices acting on rows i, i+1) and R upper triangular, unitary plus rank one, stored through its
// (n+1)x(n+1) extension Rbar = C^* (B + alpha e_1 y^T) as the core sequences C = C_1 ... C_n and
// B = B_1 ... B_n; y is implied by the zero last row of Rbar and is never formed. All cores are in SU(2),
// [[a, -conj(b)], [b, conj(a)]].
namespace companion_detail {


// The std overloads for the built-in types; argument-dependent lookup finds the multiprecision ones.
using std::abs;
using std::conj;
using std::cos;
using std::exp;
using std::hypot;
using std::log;
using std::sin;
using std::sqrt;

template <class Real>
struct core
{
   detail::complex_from_real_t<Real> a, b;
};

template <class Complex>
inline core<detail::real_from_complex_t<Complex> > normalized(Complex a, Complex b)
{
   using Real = detail::real_from_complex_t<Complex>;
   Real r = hypot(abs(a), abs(b));
   return r == 0 ? core<Real>{Complex(1, 0), Complex(0, 0)} : core<Real>{a / r, b / r};
}

// The core U with U^* (x1, x2) = (|x|, 0).
template <class Complex>
inline core<detail::real_from_complex_t<Complex> > eliminating(Complex x1, Complex x2)
{
   return normalized(x1, x2);
}

template <class Real>
inline core<Real> adjoint(core<Real> const& c)
{
   return {conj(c.a), -c.b};
}

// c1 c2, both acting on the same rows.
template <class Real>
inline core<Real> fuse(core<Real> const& c1, core<Real> const& c2)
{
   using complex_type = detail::complex_from_real_t<Real>;
   return normalized(complex_type(c1.a * c2.a - conj(c1.b) * c2.b), complex_type(c1.b * c2.a + conj(c1.a) * c2.b));
}

// F_i G_{i+1} H_i = F'_{i+1} G'_i H'_{i+1}.
template <class Real>
inline void turnover(core<Real> const& F, core<Real> const& G, core<Real> const& H,
                     core<Real>& Fp, core<Real>& Gp, core<Real>& Hp)
{
   using C = detail::complex_from_real_t<Real>;
   // M = F (rows 1,2) * G (rows 2,3) * H (rows 1,2), 3x3.
   C f[2][2] = {{F.a, -conj(F.b)}, {F.b, conj(F.a)}};
   C g[2][2] = {{G.a, -conj(G.b)}, {G.b, conj(G.a)}};
   C h[2][2] = {{H.a, -conj(H.b)}, {H.b, conj(H.a)}};
   // FG = F * G
   C fg[3][3] = {{f[0][0], f[0][1] * g[0][0], f[0][1] * g[0][1]},
                 {f[1][0], f[1][1] * g[0][0], f[1][1] * g[0][1]},
                 {C(0), g[1][0], g[1][1]}};
   C m[3][3];
   for(int i = 0; i < 3; ++i)
   {
      m[i][0] = fg[i][0] * h[0][0] + fg[i][1] * h[1][0];
      m[i][1] = fg[i][0] * h[0][1] + fg[i][1] * h[1][1];
      m[i][2] = fg[i][2];
   }
   // F' on rows (2,3) zeroes m[2][0]; G' on rows (1,2) zeroes m[1][0]; what remains is H' on rows (2,3).
   Fp = eliminating(m[1][0], m[2][0]);
   C fa = conj(Fp.a), fb = conj(Fp.b);
   C m1[3][3];
   for(int j = 0; j < 3; ++j)
   {
      m1[0][j] = m[0][j];
      m1[1][j] = fa * m[1][j] + fb * m[2][j];
      m1[2][j] = -Fp.b * m[1][j] + Fp.a * m[2][j];
   }
   Gp = eliminating(m1[0][0], m1[1][0]);
   // Column 2 of G'^* m1, rows 2 and 3, is H' e_1.
   C h11 = -Gp.b * m1[0][1] + Gp.a * m1[1][1];
   C h21 = m1[2][1];
   Hp = normalized(h11, h21);
}

template <class Real>
class companion
{
public:
   using complex_type = detail::complex_from_real_t<Real>;

   // Coefficients of a monic polynomial, a_0 ... a_{n-1} (a_n = 1 implied), n >= 2.
   explicit companion(std::vector<complex_type> const& a) : n_(a.size())
   {
      std::size_t n = n_;
      Q_.assign(n - 1, core<Real>{complex_type(0, 0), complex_type(1, 0)});
      // R = I except its last column r: r_k = -a_k (k < n), r_n = -a_0 (-1)^(n-1). With Rbar = Y + z e_n^T
      // (Y swaps columns n and n+1, with a sign to stay in SU(2)), z = (r, -1); C zeroes z from the bottom.
      std::vector<complex_type> z(n + 1);
      for(std::size_t k = 1; k < n; ++k)
         z[k - 1] = -a[k];
      z[n - 1] = (n % 2 == 1) ? -a[0] : a[0];
      z[n] = complex_type(-1, 0);
      C_.resize(n);
      for(std::size_t m = n; m-- > 0;)
      {
         complex_type vm = z[m], vn = z[m + 1];
         Real rho = hypot(abs(vm), abs(vn));
         if(rho == 0)
            C_[m] = {complex_type(1, 0), complex_type(0, 0)};
         else
            C_[m] = {conj(vm) / rho, -vn / rho};
         z[m] = rho;
         z[m + 1] = 0;
      }
      B_ = C_;
      core<Real> cn = C_[n - 1];
      B_[n - 1] = {-conj(cn.b), conj(cn.a)};
   }

   std::vector<complex_type> eigenvalues()
   {
      std::size_t n = n_;
      Real eps = tools::epsilon<Real>();
      std::size_t hi = n - 1;
      unsigned its = 0, restarts = 0;
      while(hi > 0)
      {
         for(std::size_t i = hi; i-- > 0;)
            if(abs(Q_[i].b) < eps)
            {
               Real m = abs(Q_[i].a);
               Q_[i] = {m > 0 ? Q_[i].a / m : complex_type(1, 0), complex_type(0, 0)};
            }
         std::size_t lo = hi;
         while(lo > 0 && Q_[lo - 1].b != complex_type(0, 0))
            --lo;
         if(lo == hi)
         {
            --hi;
            its = 0;
            continue;
         }
         if(++its > 30 * (hi - lo + 1))
            throw std::runtime_error("polynomial_roots_companion_matrix: the companion QR iteration did not converge");
         complex_type mu = shift(lo, hi, its, restarts);
         step(lo, hi, mu);
      }
      std::vector<complex_type> lambda(n);
      for(std::size_t k = 0; k < n; ++k)
      {
         complex_type ak = k + 1 < n ? Q_[k].a : complex_type(1, 0);
         complex_type akm1 = k > 0 ? Q_[k - 1].a : complex_type(1, 0);
         complex_type rkk[3];
         r_column(k, 1, rkk);
         lambda[k] = ak * conj(akm1) * rkk[0];
      }
      return lambda;
   }

private:
   // Entries R_{k,k}, R_{k-1,k}, ... (depth of them) of column k, by back substitution in
   // C Rbar e_k = B e_k + alpha y_k e_1, whose rows below the first involve only the cores.
   void r_column(std::size_t k, std::size_t depth, complex_type* out) const
   {
      complex_type beta[4];
      complex_type w = B_[k].a;
      beta[0] = B_[k].b;                           // row k+1
      std::size_t d = 1;
      for(std::size_t m = k; m-- > 0 && d < depth;)
      {
         beta[d++] = conj(B_[m].a) * w;      // row m+1
         w = -conj(B_[m].b) * w;
      }
      complex_type u(0, 0);
      for(std::size_t j = 0; j < depth && j <= k; ++j)
      {
         std::size_t m = k - j;
         complex_type v = (beta[j] - conj(C_[m].a) * u) / C_[m].b;
         out[j] = v;
         u = C_[m].a * v - conj(C_[m].b) * u;
      }
   }

   // R G_i = G'_i R', i <= n - 2.
   core<Real> pass_through(std::size_t i, core<Real> const& G)
   {
      core<Real> Gt, Gp, c1, c0;
      turnover(B_[i], B_[i + 1], G, Gt, B_[i], B_[i + 1]);
      turnover(adjoint(C_[i + 1]), adjoint(C_[i]), Gt, Gp, c1, c0);
      C_[i + 1] = adjoint(c1);
      C_[i] = adjoint(c0);
      return Gp;
   }

   complex_type shift(std::size_t lo, std::size_t hi, unsigned its, unsigned& restarts) const
   {
      complex_type col_hi[3] = {}, col_h1[3] = {};
      r_column(hi, 3, col_hi);
      r_column(hi - 1, 2, col_h1);
      // R entries by (row, column) relative to hi: col_hi[j] = R_{hi-j, hi}, col_h1[j] = R_{hi-1-j, hi-1}.
      complex_type dq = hi + 1 < n_ ? Q_[hi].a : complex_type(1, 0);
      core<Real> q1 = Q_[hi - 1];
      auto rows = [&](complex_type x_hm2, complex_type x_hm1, complex_type x_hi, complex_type& y_hm1, complex_type& y_hi)
      {
         x_hi *= dq;
         y_hm1 = q1.a * x_hm1 - conj(q1.b) * x_hi;
         y_hi = q1.b * x_hm1 + conj(q1.a) * x_hi;
         if(hi >= 2)
         {
            core<Real> q2 = Q_[hi - 2];
            y_hm1 = hi - 2 >= lo ? complex_type(q2.b * x_hm2 + conj(q2.a) * y_hm1) : complex_type(conj(q2.a) * y_hm1);
         }
      };
      complex_type p, q, s, t;
      rows(hi >= 2 ? col_h1[1] : complex_type(0), col_h1[0], complex_type(0), p, q);
      rows(hi >= 2 ? col_hi[2] : complex_type(0), col_hi[1], col_hi[0], s, t);
      if(its % 10 == 0)
      {
         // Exceptional shift, nonzero even when the trailing block is nilpotent (as for x^n - 1).
         ++restarts;
         Real theta = Real(2.4) * Real(restarts);
         Real scale = (std::max)((std::max)(abs(q), abs(t)), Real(1e-3));
         return t + Real(0.75) * scale * complex_type(cos(theta), sin(theta));
      }
      // Wilkinson: the eigenvalue of [[p, s], [q, t]] nearer t.
      complex_type half_trace = (p + t) / Real(2);
      complex_type disc = sqrt(half_trace * half_trace - (p * t - s * q));
      complex_type l1 = half_trace + disc, l2 = half_trace - disc;
      return abs(l1 - t) < abs(l2 - t) ? l1 : l2;
   }

   void step(std::size_t lo, std::size_t hi, complex_type mu)
   {
      complex_type rlo[3];
      r_column(lo, 1, rlo);
      complex_type dlo = lo > 0 ? conj(Q_[lo - 1].a) : complex_type(1, 0);
      complex_type x1 = dlo * Q_[lo].a * rlo[0] - mu, x2 = Q_[lo].b * rlo[0];
      core<Real> U = eliminating(x1, x2);
      core<Real> Ustar = adjoint(U);
      // U^* passes the deflated core Q_{lo-1} = diag(a, conj(a)), which multiplies its b by conj(a).
      if(lo > 0)
         Ustar.b *= conj(Q_[lo - 1].a);
      Q_[lo] = fuse(Ustar, Q_[lo]);
      core<Real> X = U;
      for(std::size_t i = lo; i < hi; ++i)
      {
         core<Real> Xp = pass_through(i, X);
         if(i + 1 == hi)
         {
            // Xp passes the deflated core Q_hi = diag(a, conj(a)), which multiplies its b by a.
            if(hi + 1 < n_)
               Xp.b *= Q_[hi].a;
            Q_[i] = fuse(Q_[i], Xp);
         }
         else
         {
            core<Real> q0, q1;
            turnover(Q_[i], Q_[i + 1], Xp, X, q0, q1);
            Q_[i] = q0;
            Q_[i + 1] = q1;
         }
      }
   }

   std::size_t n_;
   std::vector<core<Real> > Q_, C_, B_;
};

} // namespace companion_detail

namespace polynomial_roots_detail {

// All roots of a_0 + a_1 x + ... + a_n x^n, as eigenvalues of the companion matrix.
template <class Real>
std::vector<detail::complex_from_real_t<Real> > companion_qr(std::vector<Real> coefficients)
{
   using std::abs;
   using std::exp;
   using std::log;
   using complex_type = detail::complex_from_real_t<Real>;
   while(!coefficients.empty() && coefficients.back() == 0)
      coefficients.pop_back();
   if(coefficients.empty())
      throw std::domain_error("polynomial_roots_companion_matrix: the zero polynomial has no finite root set");
   for(Real c : coefficients)
      if(!(BOOST_MATH_NAMESPACE::isfinite)(c))
         throw std::domain_error("polynomial_roots_companion_matrix: coefficients must be finite");
   std::size_t zero_roots = 0;
   while(zero_roots + 1 < coefficients.size() && coefficients[zero_roots] == 0)
      ++zero_roots;
   std::vector<complex_type> roots;
   std::size_t n = coefficients.size() - 1 - zero_roots;
   if(n == 1)
      roots.push_back(complex_type(-coefficients[zero_roots] / coefficients.back(), Real(0)));
   else if(n >= 2)
   {
      // The companion QR is normwise backward stable, so scale the variable, x = s w, to make the
      // constant and leading coefficients equal in magnitude: s = |a_0 / a_n|^(1/n), formed in logs so
      // that neither it nor the scaled coefficients a_k s^(k-n) / a_n overflow.
      Real const* a = coefficients.data() + zero_roots;
      Real log_s = (log(abs(a[0])) - log(abs(a[n]))) / Real(n);
      std::vector<complex_type> monic(n);
      for(std::size_t k = 0; k < n; ++k)
         monic[k] = a[k] == 0 ? complex_type(0, 0)
            : complex_type(Real((a[k] < 0) != (a[n] < 0) ? -1 : 1), Real(0)) *
              exp(log(abs(a[k])) - log(abs(a[n])) + (Real(k) - Real(n)) * log_s);
      roots = companion_detail::companion<Real>(monic).eigenvalues();
      Real s = exp(log_s);
      for(auto& z : roots)
         z *= s;
   }
   roots.insert(roots.end(), zero_roots, complex_type(0, 0));
   return roots;
}


} // namespace polynomial_roots_detail

// All roots of a_0 + a_1 x + ... + a_n x^n, coefficients in ascending order, by the piecewise method.
// RandomAccessContainer is std::vector, std::array, std::span, ...; the roots come back as
// std::vector<std::complex<Real>> for built-in Real, or of the matching Boost.Multiprecision complex type.
BOOST_MATH_EXPORT template <class RandomAccessContainer>
auto polynomial_roots(RandomAccessContainer const& coefficients)
{
   using Real = typename std::remove_cv<typename RandomAccessContainer::value_type>::type;
   return polynomial_roots_detail::piecewise(std::vector<Real>(coefficients.begin(), coefficients.end()));
}

// The same roots as the eigenvalues of the companion matrix, by the companion QR algorithm of Aurentz, Mach,
// Vandebril and Watkins: O(n^2), and normwise backward stable.
BOOST_MATH_EXPORT template <class RandomAccessContainer>
auto polynomial_roots_companion_matrix(RandomAccessContainer const& coefficients)
{
   using Real = typename std::remove_cv<typename RandomAccessContainer::value_type>::type;
   return polynomial_roots_detail::companion_qr(std::vector<Real>(coefficients.begin(), coefficients.end()));
}

} BOOST_MATH_NAMESPACE_END

#endif
