// Copyright 2026 Nick Thompson
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// Embind bindings for the WebAssembly build: interpolators, and the Daubechies scaling functions and wavelets.

#include "bindings.hpp"

#include <boost/math/interpolators/barycentric_rational.hpp>
#include <boost/math/interpolators/cardinal_cubic_b_spline.hpp>
#include <boost/math/interpolators/catmull_rom.hpp>
#include <boost/math/interpolators/makima.hpp>
#include <boost/math/interpolators/pchip.hpp>
#include <boost/math/quadrature/wavelet_transforms.hpp>
#include <boost/math/special_functions/daubechies_scaling.hpp>
#include <boost/math/special_functions/daubechies_wavelet.hpp>

#include <array>
#include <functional>
#include <limits>
#include <memory>
#include <utility>

namespace
{
using vec = std::vector<double>;

bool finite(const vec& v)
{
   return std::all_of(v.begin(), v.end(), [](double x) { return std::isfinite(x); });
}

// The interpolators take (x, y, t): the data, and the abscissas t at which to evaluate. Some report bad data through
// the error policy, which would give NaN everywhere and no reason, so the data are checked here.
std::pair<vec, vec> data(const val& x, const val& y, std::size_t least)
{
   std::pair<vec, vec> d(numbers(x), numbers(y));
   if (d.first.size() != d.second.size())
      range_error("x and y must have the same length.");
   if (d.first.size() < least)
      range_error("Need at least " + std::to_string(least) + " points.");
   if (!finite(d.first) || !finite(d.second))
      range_error("The data must be finite.");
   if (std::adjacent_find(d.first.begin(), d.first.end(), std::greater_equal<double>()) != d.first.end())
      range_error("The abscissas must be strictly increasing.");
   return d;
}

template <class Interpolator>
val evaluate(const Interpolator& f, const val& t)
{
   vec v = numbers(t);
   for (double& u : v)
      u = f(u);
   return to_js(v);
}

// Points are rows of 2 or 3 numbers: catmull_rom wants points of fixed size.
using point = std::array<double, 3>;
std::vector<point> points(const val& p, std::size_t& dimension)
{
   std::vector<point> m;
   dimension = 0;
   for (std::size_t i = 0, n = p["length"].as<std::size_t>(); i < n; ++i)
   {
      const vec v = numbers(p[i]);
      dimension = i ? dimension : v.size();
      if (v.size() != dimension || dimension < 2 || dimension > 3 || !finite(v))
         range_error("The points must be arrays of 2 or of 3 finite numbers.");
      m.push_back({v[0], v[1], dimension == 3 ? v[2] : 0});
   }
   return m;
}

bm::catmull_rom<point> curve(const val& p, const val& o, std::size_t& dimension)
{
   options opt(o);
   const bool closed = opt.flag("closed");
   const double alpha = opt.real("alpha", 0.5);
   opt.done();
   return bm::catmull_rom<point>(points(p, dimension), closed, alpha);
}

// The curve at the parameters s in [0, max_parameter], as rows flattened: s.length * dimension numbers. A parameter
// outside that range, or NaN, throws.
val catmull_rom(const val& p, const val& s, const val& o)
{
   std::size_t d;
   const bm::catmull_rom<point> c = curve(p, o, d);
   vec out;
   for (double u : numbers(s))
   {
      const point q = c(u);
      out.insert(out.end(), q.begin(), q.begin() + d);
   }
   return to_js(out);
}

// The parameter at each point, then max_parameter: for a closed curve, where it returns to the first point.
val catmull_rom_parameters(const val& p, const val& o)
{
   std::size_t d;
   const std::size_t n = p["length"].as<std::size_t>();
   const bm::catmull_rom<point> c = curve(p, o, d);
   vec s(n);
   for (std::size_t i = 0; i < n; ++i)
      s[i] = c.parameter_at_point(i);
   if (s.back() != c.max_parameter())
      s.push_back(c.max_parameter());
   return to_js(s);
}

// The Daubechies scaling function and wavelet with p = 2..10 vanishing moments, zero outside their supports
// [0, 2p - 1] and [1 - p, p]. They interpolate a grid of 2^grid_refinements points per unit, an option. The error
// falls with the refinement at a rate set by the smoothness, which grows with p: at 10, 14 and 16 refinements it is
// 1e-2, 2e-3 and 9e-4 for p = 2; 6e-4, 3e-5 and 4e-6 for p = 3; 2e-5, 2e-7 and 2e-8 for p = 4; and 7e-8, 1e-10 and 5e-12
// for p = 6. The default is 14 for p <= 5, 12 for p = 6 and 7, and 10 above, which is 4e-9 or better from p = 5 on, in
// at most 10 MB; the C++ default of up to 21 takes hundreds, so a grid of more than 2^23 numbers is refused. The grid
// of each function and p is built on first use and kept until another refinement is asked for; whoever uses one holds
// it, so that it outlives its replacement.
int refinements(options& opt, int p)
{
   const int r = int(opt.whole("grid_refinements", p <= 5 ? 14 : p <= 7 ? 12 : 10, 1, 21));
   // The function and one to three derivatives, at each point.
   if (std::ldexp((2 * p - 1) * (p < 6 ? 2 : p < 10 ? 3 : 4), r) > 0x1p23)
      range_error("grid_refinements: the grid would hold more than 2^23 numbers.");
   return r;
}

template <template <class, int> class W, int p>
std::shared_ptr<const W<double, p>> grid(int r)
{
   static std::shared_ptr<const W<double, p>> kept;
   static int at = -1;
   if (at != r)
   {
      kept = std::make_shared<const W<double, p>>(r);
      at = r;
   }
   return kept;
}

template <template <class, int> class W, int p>
std::function<double(double)> daubechies_at(int r)
{
   return [w = grid<W, p>(r)](double x) {
      const auto support = w->support();
      return x >= support.first && x <= support.second ? (*w)(x) : x == x ? 0 : x;
   };
}

template <template <class, int> class W, std::size_t... I>
std::function<double(double)> daubechies(double p, const val& o, std::index_sequence<I...>)
{
   static std::function<double(double)> (*const f[])(int) = {&daubechies_at<W, int(I) + 2>...};
   const int moments = int(whole(p, "p", 2, 10));
   options opt(o);
   const int r = refinements(opt, moments);
   opt.done();
   return f[moments - 2](r);
}

// name(p, x) and name_array(p, xs), each with an optional {grid_refinements}.
template <template <class, int> class W>
void bind_daubechies(const char* name, const char* array_name)
{
   static const auto f = [](double p, const val& o) { return daubechies<W>(p, o, std::make_index_sequence<9>()); };
   entry(name, +[](double p, double x) { return f(p, val::undefined())(x); }, +[](double p, double x, const val& o) { return f(p, o)(x); });
   entry(array_name, +[](double p, const val& x) { return evaluate(f(p, val::undefined()), x); },
         +[](double p, const val& x, const val& o) { return evaluate(f(p, o), x); });
}

// W(s, t) = sqrt|s| int f(s u + t) psi_p(u) du for each scale s and each time t, s.length * t.length numbers row major,
// where f interpolates the samples y_i at i dt linearly and is zero before the first and after the last.
template <int p>
vec wavelet_transform(const vec& y, double dt, const vec& s, const vec& t, double tolerance, int max_refinements, int r)
{
   const std::size_t n = y.size();
   auto f = [&y, dt, n](double u) {
      const double i = u / dt;
      if (!(i >= 0 && i <= double(n - 1)))
         return 0.0;
      const std::size_t k = (std::min)(std::size_t(i), n - 2);
      return y[k] + (i - double(k)) * (y[k + 1] - y[k]);
   };
   const auto psi = grid<bm::daubechies_wavelet, p>(r);
   const bm::quadrature::daubechies_wavelet_transform<decltype(f), double, p> W(f, *psi, tolerance, max_refinements);
   vec out;
   out.reserve(s.size() * t.size());
   for (double scale : s)
      for (double time : t)
         out.push_back(W(scale, time));
   return out;
}

template <std::size_t... I>
val wavelet_transform(const val& y, double dt, double p, const val& s, const val& t, const val& o, std::index_sequence<I...>)
{
   static vec (*const f[])(const vec&, double, const vec&, const vec&, double, int, int) = {&wavelet_transform<int(I) + 2>...};
   options opt(o);
   const double tolerance = opt.real("tolerance", 100 * std::numeric_limits<double>::epsilon());
   const int moments = int(whole(p, "p", 2, 10)), max_refinements = int(opt.whole("max_refinements", 12, 1, 20)), r = refinements(opt, moments);
   opt.done();
   const vec v = numbers(y), scales = numbers(s), times = numbers(t);
   check_size(double(scales.size()) * double(times.size()), "The transform");
   if (v.size() < 2 || !finite(v) || !(dt > 0 && std::isfinite(dt)) || !finite(scales) || !finite(times) || !(tolerance > 0))
      range_error("Need at least 2 finite samples, dt > 0, finite scales and times, and a tolerance > 0.");
   return to_js(f[moments - 2](v, dt, scales, times, tolerance, max_refinements, r));
}

} // namespace

EMSCRIPTEN_BINDINGS(boost_math_interpolation_wavelets)
{
   entry("makima", +[](const val& x, const val& y, const val& t) {
      auto d = data(x, y, 4);
      return evaluate(bm::interpolators::makima<vec>(std::move(d.first), std::move(d.second)), t);
   });
   entry("pchip", +[](const val& x, const val& y, const val& t) {
      auto d = data(x, y, 4);
      return evaluate(bm::interpolators::pchip<vec>(std::move(d.first), std::move(d.second)), t);
   });
   // The approximation order, 3 by default, is at most n - 1.
   entry("barycentric_rational", +[](const val& x, const val& y, const val& t) {
      auto d = data(x, y, 4);
      return evaluate(bm::interpolators::barycentric_rational<double>(std::move(d.first), std::move(d.second)), t);
   });
   entry("barycentric_rational", +[](const val& x, const val& y, const val& t, double order) {
      auto d = data(x, y, 2);
      if (!(order == std::trunc(order) && order >= 0 && order < double(d.first.size())))
         range_error("The approximation order must be an integer in [0, n).");
      return evaluate(bm::interpolators::barycentric_rational<double>(std::move(d.first), std::move(d.second), std::size_t(order)), t);
   });
   // The samples y_i are at t0 + i h.
   entry("cardinal_cubic_b_spline", +[](const val& y, double t0, double h, const val& t) {
      const vec v = numbers(y);
      if (v.size() < 5 || !finite(v) || !std::isfinite(t0) || !(h > 0 && std::isfinite(h)))
         range_error("Need at least 5 finite samples, a finite t0 and h > 0.");
      return evaluate(bm::interpolators::cardinal_cubic_b_spline<double>(v.data(), v.size(), t0, h), t);
   });
   entry("catmull_rom", +[](const val& p, const val& s) { return catmull_rom(p, s, val::undefined()); });
   entry("catmull_rom", &catmull_rom);
   entry("catmull_rom_parameters", +[](const val& p) { return catmull_rom_parameters(p, val::undefined()); });
   entry("catmull_rom_parameters", &catmull_rom_parameters);

   bind_daubechies<bm::daubechies_scaling>("daubechies_scaling", "daubechies_scaling_array");
   bind_daubechies<bm::daubechies_wavelet>("daubechies_wavelet", "daubechies_wavelet_array");
   entry("daubechies_wavelet_transform", +[](const val& y, double dt, double p, const val& s, const val& t) {
      return wavelet_transform(y, dt, p, s, t, val::undefined(), std::make_index_sequence<9>());
   });
   entry("daubechies_wavelet_transform", +[](const val& y, double dt, double p, const val& s, const val& t, const val& o) {
      return wavelet_transform(y, dt, p, s, t, o, std::make_index_sequence<9>());
   });
}
