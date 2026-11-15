// Copyright 2026 Nick Thompson
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// Embind bindings for the WebAssembly build: quadrature and root finding. The integrand or function is a JavaScript
// callback, which is synchronous, and the last argument of most functions is an optional object of options.

#include "bindings.hpp"

#include <emscripten/em_js.h>

#include <boost/math/quadrature/exp_sinh.hpp>
#include <boost/math/quadrature/gauss.hpp>
#include <boost/math/quadrature/gauss_hermite.hpp>
#include <boost/math/quadrature/gauss_jacobi.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/math/quadrature/gauss_laguerre.hpp>
#include <boost/math/quadrature/leja.hpp>
#include <boost/math/quadrature/ooura_fourier_integrals.hpp>
#include <boost/math/quadrature/sinh_sinh.hpp>
#include <boost/math/quadrature/sparse_grid.hpp>
#include <boost/math/quadrature/tanh_sinh.hpp>
#include <boost/math/tools/multivariate_roots.hpp>
#include <boost/math/tools/polynomial_roots.hpp>
#include <boost/math/tools/roots.hpp>
#include <boost/math/tools/toms748_solve.hpp>

#include <array>
#include <cstdint>
#include <memory>
#include <tuple>
#include <type_traits>
#include <utility>

using emscripten::EM_VAL;
namespace bq = boost::math::quadrature;
namespace bt = boost::math::tools;

namespace
{
const double root_epsilon = bt::root_epsilon<double>();

// x => f(x) for a callback f: a direct call, a third of the time of f(x).as<double>(). Callbacks are synchronous: a
// promise, or anything else with a then, is refused.
EM_JS_DEPS(boost_math_call_real, "$Emval");
EM_JS(double, call_real, (EM_VAL f, double x), {
   const y = Emval.toValue(f)(x);
   if (typeof y === 'object' && y !== null && typeof y.then === 'function')
      throw new TypeError('The callback returned a promise: callbacks must be synchronous.');
   return y;
});

auto real_function(const val& f)
{
   if (f.typeOf().as<std::string>() != "function")
      range_error("The function must be a callback.");
   return [h = f.as_handle()](double x) { return call_real(h, x); };
}

// {result, error, L1}, without what the rule does not report.
val integral(double result, const double* error, const double* L1)
{
   val o = val::object();
   o.set("result", result);
   if (error)
      o.set("error", *error);
   if (L1)
      o.set("L1", *L1);
   return o;
}

// g(N) for the rule with n points, N an integral_constant, where the number of points is a template parameter.
template <unsigned... Ns, class G>
val with_points(double n, G g)
{
   val r = val::undefined();
   if (!((n == Ns && (r = g(std::integral_constant<unsigned, Ns>()), true)) || ...))
      range_error("points must be one of" + ((" " + std::to_string(Ns)) + ...) + ".");
   return r;
}

// tanh_sinh over [a, b], exp_sinh over a half-infinite interval, sinh_sinh over the real line.
template <class Rule, class... Bounds>
val double_exponential(const val& f, const val& o, std::size_t refinements, Bounds... bounds)
{
   options opt(o);
   const double tol = opt.real("tolerance", root_epsilon);
   const Rule rule(std::size_t(opt.whole("max_refinements", double(refinements), 0, 30)));
   opt.done();
   double error = NAN, L1 = NAN;
   std::size_t levels = 0;
   val r = integral(rule.integrate(real_function(f), bounds..., tol, &error, &L1, &levels), &error, &L1);
   r.set("levels", double(levels));
   return r;
}

val gauss_kronrod(const val& f, double a, double b, const val& o)
{
   options opt(o);
   const double points = opt.real("points", 15), tol = opt.real("tolerance", root_epsilon);
   const unsigned depth = unsigned(opt.whole("max_depth", 15, 0, 60));
   opt.done();
   return with_points<15, 21, 31, 41, 51, 61>(points, [&](auto N) {
      double error = NAN, L1 = NAN;
      const double r = bq::gauss_kronrod<double, decltype(N)::value>::integrate(real_function(f), a, b, depth, tol, &error, &L1);
      return integral(r, &error, &L1);
   });
}

val gauss(const val& f, double a, double b, const val& o)
{
   options opt(o);
   const double points = opt.real("points", 7);
   opt.done();
   return with_points<7, 10, 15, 20, 25, 30>(points, [&](auto N) {
      double L1 = NAN;
      const double r = bq::gauss<double, decltype(N)::value>::integrate(real_function(f), a, b, &L1);
      return integral(r, nullptr, &L1);
   });
}

// A Leja rule computes its points as they are first used, so the rules are kept: one per weight, and for a weight
// with exponents, the one last used. An integration holds its rule, which an integrand may replace here by
// integrating with other exponents.
template <class Rule, class... Exponents>
std::shared_ptr<const Rule> kept(Exponents... exponents)
{
   static std::tuple<Exponents...> last;
   static std::shared_ptr<const Rule> rule;
   if (!rule || last != std::make_tuple(exponents...))
   {
      rule = std::make_shared<const Rule>(exponents...);
      last = std::make_tuple(exponents...);
   }
   return rule;
}

// The exponent of a weight (1 - x)^alpha, x^alpha exp(-x), ...
double exponent(options& opt, const char* name)
{
   const double a = opt.real(name, 0);
   if (!(a > -1 && std::isfinite(a)))
      range_error(std::string(name) + " must be finite and greater than -1.");
   return a;
}

// The integral of f(x) w(x) for a nested Leja rule with weight w; bounds are [a, b] for leja_legendre and leja_jacobi.
template <class Rule, class... Bounds>
val leja(const val& f, options& opt, std::shared_ptr<const Rule> rule, Bounds... bounds)
{
   const double tol = opt.real("tolerance", root_epsilon);
   opt.done();
   double error = NAN, L1 = NAN;
   return integral(rule->integrate(real_function(f), bounds..., tol, &error, &L1), &error, &L1);
}

// A Gauss rule for a weight; averaged: true adds Laurie's anti-Gauss rule, points + 1 more evaluations, for an error.
template <class Rule>
val weighted_gauss(const val& f, const Rule& rule, bool averaged)
{
   double error = NAN, L1 = NAN;
   const double r = averaged ? rule.integrate_averaged(real_function(f), &error, &L1) : rule.integrate(real_function(f), &L1);
   return integral(r, averaged ? &error : nullptr, &L1);
}

val gauss_jacobi(const val& f, const val& o)
{
   options opt(o);
   const unsigned n = unsigned(opt.whole("points", 8, 1, 4096));
   const double alpha = exponent(opt, "alpha"), beta = exponent(opt, "beta");
   const bool averaged = opt.flag("averaged");
   opt.done();
   return weighted_gauss(f, bq::gauss_jacobi<double>(n, alpha, beta), averaged);
}

template <bool Laguerre>
val gauss_on_line(const val& f, const val& o)
{
   options opt(o);
   const double points = opt.real("points", 8), alpha = Laguerre ? exponent(opt, "alpha") : 0;
   const bool averaged = opt.flag("averaged");
   opt.done();
   return with_points<2, 4, 8, 16, 32, 64, 128>(points, [&](auto N) {
      if constexpr (Laguerre)
         return weighted_gauss(f, bq::gauss_laguerre<double, decltype(N)::value>(alpha), averaged);
      else
         return weighted_gauss(f, bq::gauss_hermite<double, decltype(N)::value>(), averaged);
   });
}

// The integral of f(x) sin(omega x) or f(x) cos(omega x) over [0, infinity): {result, relative_error}.
template <class Rule>
val ooura(const val& f, double omega, const val& o)
{
   options opt(o);
   const double tol = opt.real("tolerance", root_epsilon);
   Rule rule(tol, std::size_t(opt.whole("levels", 8, 1, 20)));
   opt.done();
   const std::pair<double, double> r = rule.integrate(real_function(f), omega);
   val result = val::object();
   result.set("result", r.first);
   result.set("relative_error", r.second);
   return result;
}

// The integral of f over the box [a, b], by a dimension-adaptive sparse grid of Leja points; f takes a point.
val sparse_grid(const val& f, const val& a, const val& b, const val& o)
{
   options opt(o);
   const std::vector<double> lower = numbers(a), upper = numbers(b);
   const std::size_t calls = std::size_t(opt.whole("max_calls", 100000, 1, 1 << 24));
   const double tol = opt.real("tolerance", root_epsilon);
   opt.done();
   if (lower.empty() || lower.size() != upper.size())
      range_error("The box needs as many lower as upper bounds, and at least one.");
   double error = NAN, L1 = NAN;
   const auto rule = kept<bq::leja_legendre<double>>();
   const double r = bq::sparse_grid([&f](const std::vector<double>& x) { return sync(f(to_js(x))).as<double>(); },
                                    *rule, lower, upper, calls, tol, &error, &L1);
   return integral(r, &error, &L1);
}

// f(x) returns [f, f'], and [f, f', f''] for the methods of third order.
template <std::size_t N>
auto derivatives(const val& f)
{
   return [&f](double x) {
      const std::vector<double> y = numbers(sync(f(x)));
      if (y.size() < N)
         range_error("f(x) must return an array of " + std::to_string(N) + " numbers.");
      return N == 2 ? std::make_tuple(y[0], y[1], 0.0) : std::make_tuple(y[0], y[1], y[2]);
   };
}

// The options of every root finder: binary digits wanted, and the most evaluations.
struct root_options
{
   int digits;
   std::uintmax_t iterations;
   explicit root_options(const val& o)
   {
      options opt(o);
      digits = int(opt.whole("digits", 53, 1, 53));
      iterations = std::uintmax_t(opt.whole("max_iterations", 200, 1, 1e9));
      opt.done();
   }
};

// 0: Newton-Raphson, 1: Halley, 2: Schroder. {root, iterations}.
template <int Method>
val iterate(const val& f, double guess, double lower, double upper, const val& o)
{
   root_options opt(o);
   const double x = Method == 0 ? bt::newton_raphson_iterate([g = derivatives<2>(f)](double x) { auto [y, d, _] = g(x); return std::make_pair(y, d); },
                                                             guess, lower, upper, opt.digits, opt.iterations)
                    : Method == 1 ? bt::halley_iterate(derivatives<3>(f), guess, lower, upper, opt.digits, opt.iterations)
                                  : bt::schroder_iterate(derivatives<3>(f), guess, lower, upper, opt.digits, opt.iterations);
   val r = val::object();
   r.set("root", x);
   r.set("iterations", double(opt.iterations));
   return r;
}

// {lower, upper, iterations}: a bracket of the root, NaN unless f changes sign over [a, b].
template <bool Bisect>
val bracket(const val& f, double a, double b, const val& o)
{
   root_options opt(o);
   const bt::eps_tolerance<double> tol(opt.digits);
   const std::pair<double, double> x = Bisect ? bt::bisect(real_function(f), a, b, tol, opt.iterations)
                                              : bt::toms748_solve(real_function(f), a, b, tol, opt.iterations);
   val r = val::object();
   r.set("lower", x.first);
   r.set("upper", x.second);
   r.set("iterations", double(opt.iterations));
   return r;
}

// Coefficients a_0, ..., a_n of a_0 + a_1 z + ... + a_n z^n in; its roots out, as (re, im) pairs.
template <bool Companion>
val polynomial_roots(const val& coefficients)
{
   const std::vector<double> a = numbers(coefficients);
   return to_js(Companion ? bt::polynomial_roots_companion_matrix(a) : bt::polynomial_roots(a));
}

// A batch of multivariate_newton or multivariate_halley solves, too many for a callback each: p(z) = 0 as the system
// (Re p, Im p)(x, y) = 0, from the centres of the n by n cells of [-extent, extent]^2, the top row first. Per start,
// two bytes: the index of the root reached, or 255 for none, and the iterations. The roots are in order of argument.
val polynomial_basins(const val& coefficients, bool halley, double n, double extent, const val& o)
{
   options opt(o);
   const int digits = int(opt.whole("digits", 53, 1, 53));
   const std::uintmax_t most = std::uintmax_t(opt.whole("max_iterations", 100, 1, 255));
   opt.done();
   const std::vector<double> a = numbers(coefficients);
   if (!(n == std::trunc(n) && n >= 1 && n <= 4096 && extent > 0 && std::isfinite(extent)))
      range_error("The grid must have 1 to 4096 points a side and a positive extent.");
   if (a.size() < 2 || a.size() > 255 || a.back() == 0)
      range_error("The polynomial must have degree 1 to 254 and a nonzero leading coefficient.");
   std::vector<cplx> roots = bt::polynomial_roots(a);
   std::sort(roots.begin(), roots.end(), [](const cplx& u, const cplx& v) { return std::arg(u) < std::arg(v); });
   using point = std::array<double, 2>;
   const auto system = [&a](const point& x) {
      const cplx z(x[0], x[1]);
      cplx p = 0, d1 = 0, d2 = 0;
      for (std::size_t i = a.size(); i-- > 0;)
      {
         d2 = d2 * z + d1;
         d1 = d1 * z + p;
         p = p * z + a[i];
      }
      const auto second = [d2](const point& v) {
         const cplx w = 2.0 * d2 * cplx(v[0], v[1]) * cplx(v[0], v[1]);
         return point{w.real(), w.imag()};
      };
      return std::make_tuple(point{p.real(), p.imag()}, std::array<double, 4>{d1.real(), -d1.imag(), d1.imag(), d1.real()}, second);
   };
   const std::size_t side = std::size_t(n);
   std::vector<std::uint8_t> basin(2 * side * side);
   for (std::size_t j = 0, k = 0; j < side; ++j)
      for (std::size_t i = 0; i < side; ++i, k += 2)
      {
         const point start{-extent + (i + 0.5) * 2 * extent / n, extent - (j + 0.5) * 2 * extent / n};
         std::uintmax_t it = most;
         const point x = halley ? bt::multivariate_halley(system, start, digits, it) : bt::multivariate_newton(system, start, digits, it);
         std::size_t nearest = 255;
         double least = 1e-6;
         for (std::size_t r = 0; r < roots.size(); ++r)
            if (const double d = std::abs(cplx(x[0], x[1]) - roots[r]); d <= least)
               least = d, nearest = r;
         basin[k] = std::uint8_t(nearest);
         basin[k + 1] = std::uint8_t(it);
      }
   val r = val::object();
   r.set("basin", val(typed_memory_view(basin.size(), basin.data())).call<val>("slice"));
   r.set("roots", to_js(roots));
   return r;
}
// f(x) for the multivariate solvers returns [F, J] or [F, J, d2]: F an array of n, J n arrays of n or n*n numbers
// row major, and d2(v) = F''(x)[v, v] an array of n. A wrong shape throws a RangeError.
struct js_system
{
   val f;
   std::size_t n;

   static std::vector<double> numbers(const val& a, std::size_t count, const char* what)
   {
      std::vector<double> v;
      if (a.isArray() && a["length"].as<std::size_t>() > 0 && a[0].isArray())
         v = ::numbers(a.call<val>("flat"));
      else if (!a.isUndefined() && !a.isNull())
         v = ::numbers(a);
      if (v.size() != count)
         range_error(std::string(what) + " has " + std::to_string(v.size()) + " numbers, expected " + std::to_string(count) + ".");
      return v;
   }

   auto operator()(const std::vector<double>& x) const
   {
      val r = sync(f(to_js(x)));
      if (!r.isArray() || r["length"].as<std::size_t>() < 2)
         range_error("f(x) must return [F, J] or [F, J, d2].");
      val d2 = r["length"].as<std::size_t>() > 2 ? r[2] : val::undefined();
      auto second = [d2, n = n](const std::vector<double>& v) {
         if (d2.typeOf().as<std::string>() != "function")
            range_error("multivariate_halley: f(x) must return [F, J, d2], with d2 a function.");
         return numbers(sync(d2(to_js(v))), n, "d2(v)");
      };
      return std::make_tuple(numbers(r[0], n, "F(x)"), numbers(r[1], n * n, "J(x)"), second);
   }
};

// {root, iterations}: the root is all NaN if the solver failed, and iterations == max_iterations if it did not converge.
template <bool Halley>
val multivariate_root(const val& f, const val& x0, double digits, double max_iterations)
{
   std::vector<double> x = ::numbers(x0);
   std::uintmax_t it = count(max_iterations, "max_iterations", 1);
   const int d = int(whole(digits, "digits", 1, 53));
   if (x.empty())
      range_error("The starting point is empty.");
   const js_system sys{f, x.size()};
   x = Halley ? bm::tools::multivariate_halley(sys, std::move(x), d, it) : bm::tools::multivariate_newton(sys, std::move(x), d, it);
   val o = val::object();
   o.set("root", to_js(x));
   o.set("iterations", double(it));
   return o;
}

using V = const val&;
val tanh_sinh(V f, double a, double b, V o) { return double_exponential<bq::tanh_sinh<double>>(f, o, 15, a, b); }
val exp_sinh(V f, double a, double b, V o) { return double_exponential<bq::exp_sinh<double>>(f, o, 9, a, b); }
val sinh_sinh(V f, V o) { return double_exponential<bq::sinh_sinh<double>>(f, o, 9); }
val leja_legendre(V f, double a, double b, V o)
{
   options opt(o);
   return leja(f, opt, kept<bq::leja_legendre<double>>(), a, b);
}
// The integral of f(x) |b - x|^alpha |x - a|^beta.
val leja_jacobi(V f, double a, double b, V o)
{
   options opt(o);
   const double alpha = exponent(opt, "alpha"), beta = exponent(opt, "beta");
   return leja(f, opt, kept<bq::leja_jacobi<double>>(alpha, beta), a, b);
}
// The integral of f(x) exp(-x^2) over the real line, and of f(x) x^alpha exp(-x) over [0, infinity).
val leja_hermite(V f, V o)
{
   options opt(o);
   return leja(f, opt, kept<bq::leja_hermite<double>>());
}
val leja_laguerre(V f, V o)
{
   options opt(o);
   const double alpha = exponent(opt, "alpha");
   return leja(f, opt, kept<bq::leja_laguerre<double>>(alpha));
}
} // namespace

// name(args...) and name(args..., options).
#define WITH_OPTIONS(name, f, ...) entry(name, &f, +[](__VA_ARGS__) { return f(ARGS, val::undefined()); })

EMSCRIPTEN_BINDINGS(boost_math_quadrature_roots)
{
#define ARGS f, a, b
   WITH_OPTIONS("tanh_sinh", tanh_sinh, V f, double a, double b);
   WITH_OPTIONS("exp_sinh", exp_sinh, V f, double a, double b);
   WITH_OPTIONS("gauss_kronrod", gauss_kronrod, V f, double a, double b);
   WITH_OPTIONS("gauss", gauss, V f, double a, double b);
   WITH_OPTIONS("leja_legendre", leja_legendre, V f, double a, double b);
   WITH_OPTIONS("leja_jacobi", leja_jacobi, V f, double a, double b);
   WITH_OPTIONS("toms748_solve", bracket<false>, V f, double a, double b);
   WITH_OPTIONS("bisect", bracket<true>, V f, double a, double b);
   WITH_OPTIONS("sparse_grid", sparse_grid, V f, V a, V b);
#undef ARGS
#define ARGS f
   WITH_OPTIONS("sinh_sinh", sinh_sinh, V f);
   WITH_OPTIONS("leja_hermite", leja_hermite, V f);
   WITH_OPTIONS("leja_laguerre", leja_laguerre, V f);
   // Gauss rules with a given number of points, for the weights of the Leja rules; gauss_jacobi is over [-1, 1].
   WITH_OPTIONS("gauss_hermite", gauss_on_line<false>, V f);
   WITH_OPTIONS("gauss_laguerre", gauss_on_line<true>, V f);
   WITH_OPTIONS("gauss_jacobi", gauss_jacobi, V f);
#undef ARGS
#define ARGS f, omega
   WITH_OPTIONS("ooura_fourier_sin", ooura<bq::ooura_fourier_sin<double>>, V f, double omega);
   WITH_OPTIONS("ooura_fourier_cos", ooura<bq::ooura_fourier_cos<double>>, V f, double omega);
#undef ARGS
#define ARGS f, guess, lower, upper
   WITH_OPTIONS("newton_raphson_iterate", iterate<0>, V f, double guess, double lower, double upper);
   WITH_OPTIONS("halley_iterate", iterate<1>, V f, double guess, double lower, double upper);
   WITH_OPTIONS("schroder_iterate", iterate<2>, V f, double guess, double lower, double upper);
#undef ARGS
#define ARGS a, halley, n, extent
   WITH_OPTIONS("polynomial_basins", polynomial_basins, V a, bool halley, double n, double extent);
#undef ARGS
   entry("polynomial_roots", &polynomial_roots<false>);
   entry("polynomial_roots_companion_matrix", &polynomial_roots<true>);
   entry("multivariate_newton", +[](V f, V x0) { return multivariate_root<false>(f, x0, 53, 1000); },
         +[](V f, V x0, double d) { return multivariate_root<false>(f, x0, d, 1000); }, &multivariate_root<false>);
   entry("multivariate_halley", +[](V f, V x0) { return multivariate_root<true>(f, x0, 53, 1000); },
         +[](V f, V x0, double d) { return multivariate_root<true>(f, x0, d, 1000); }, &multivariate_root<true>);
}
