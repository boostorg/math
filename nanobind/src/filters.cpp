// Copyright 2026 Nick Thompson
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// boost.math.filters: filters::trend_filter, wavelet_denoise, wavelet_transform and denoise on 1-D float64 data
// (anything np.ascontiguousarray understands). The work releases the GIL; errors follow boost.math.errstate,
// as for a scalar argument, and a failed result is NaN.

#include "bindings.hpp"

#include <boost/math/filters/trend_filter.hpp>
#include <boost/math/filters/wavelet_denoise.hpp>

#include <nanobind/stl/optional.h>
#include <nanobind/stl/string.h>

#include <cmath>
#include <functional>
#include <iterator>
#include <limits>
#include <optional>
#include <string>
#include <vector>

namespace bmpy {

namespace {

namespace bf = boost::math::filters;

// An orthogonal_wavelet<double> with the factory call that made it, for __repr__.
struct py_wavelet {
   bf::orthogonal_wavelet<double> w;
   std::string family;
   unsigned p;
};

// Donoho and Johnstone's (1994) test signals at t_i = i/n, scaled to standard deviation 7 (their convention:
// adding unit-variance noise gives a signal-to-noise ratio of 7).
std::vector<double> test_signal(std::string const& name, std::size_t n) {
   static double const tj[] = {0.1, 0.13, 0.15, 0.23, 0.25, 0.40, 0.44, 0.65, 0.76, 0.78, 0.81};
   static double const hb[] = {4, -5, 3, -4, 5, -4.2, 2.1, 4.3, -3.1, 2.1, -4.2};
   static double const hu[] = {4, 5, 3, 4, 5, 4.2, 2.1, 4.3, 3.1, 5.1, 4.2};
   static double const wu[] = {0.005, 0.005, 0.006, 0.01, 0.01, 0.03, 0.01, 0.01, 0.005, 0.008, 0.005};
   const double pi = 3.14159265358979323846;
   auto sgn = [](double x) { return double((x > 0) - (x < 0)); };
   std::function<double(double)> f;
   if (name == "blocks")
      f = [&](double t) { double v = 0; for (int j = 0; j < 11; ++j) v += hb[j] * (1 + sgn(t - tj[j])) / 2; return v; };
   else if (name == "bumps")
      f = [&](double t) { double v = 0; for (int j = 0; j < 11; ++j) v += hu[j] / std::pow(1 + std::abs(t - tj[j]) / wu[j], 4); return v; };
   else if (name == "heavisine")
      f = [&](double t) { return 4 * std::sin(4 * pi * t) - sgn(t - 0.3) - sgn(0.72 - t); };
   else if (name == "doppler")
      f = [&](double t) { return std::sqrt(t * (1 - t)) * std::sin(2 * pi * 1.05 / (t + 0.05)); };
   else
      throw nb::value_error(("unknown signal '" + name + "'; expected 'blocks', 'bumps', 'heavisine' or 'doppler'").c_str());
   if (n < 2)
      throw nb::value_error(("signals." + name + ": n must be at least 2").c_str());
   std::vector<double> v(n);
   double mean = 0, var = 0;
   for (std::size_t i = 0; i < n; ++i) { v[i] = f(double(i) / double(n)); mean += v[i]; }
   mean /= double(n);
   for (double x : v) var += (x - mean) * (x - mean);
   const double scale = var > 0 ? 7 / std::sqrt(var / double(n)) : 1;
   for (double& x : v) x *= scale;
   return v;
}

std::vector<double> to_vector(nb::handle x, char const* name) {
   nb::object a = numpy().attr("ascontiguousarray")(x, nb::arg("dtype") = numpy().attr("float64"));
   if (nb::cast<int>(a.attr("ndim")) != 1)
      throw nb::value_error((std::string(name) + ": expected a 1-D array").c_str());
   auto d = nb::cast<nb::ndarray<double const, nb::ndim<1>, nb::c_contig, nb::device::cpu>>(a);
   return std::vector<double>(d.data(), d.data() + d.shape(0));
}

// f() -> std::vector<double> without the GIL; Boost.Math errors are recorded and reported per errstate,
// and an exception that escapes the policy gives a NaN result of length nan_size.
template <class F>
auto run(char const* name, std::size_t nan_size, F f) {
   error_actions actions = current_error_actions;
   error_record errors;
   std::vector<double> out;
   {
      recording rec(errors);
      nb::gil_scoped_release release;
      try {
         out = f();
      } catch (std::exception const& e) {
         record_exception(e);
         out.assign(nan_size, std::numeric_limits<double>::quiet_NaN());
      }
   }
   if (errors.any())
      report(name, errors, 1, true, false, actions);
   return vector_to_array(std::move(out));
}

constexpr char const* rule_names[] = {"wiener", "universal_soft", "sure", "eb_posterior_median"};

static_assert(std::size(rule_names) == std::size_t(bf::threshold_rule::eb_posterior_median) + 1, "rule_names must name every threshold_rule");

bf::threshold_rule parse_rule(std::string const& r) {
   for (std::size_t i = 0; i < std::size(rule_names); ++i)
      if (r == rule_names[i]) return static_cast<bf::threshold_rule>(i);
   std::string msg = "wavelet_denoise: unknown rule '" + r + "'; expected one of";
   for (char const* n : rule_names) msg += std::string(" '") + n + "'";
   throw nb::value_error(msg.c_str());
}

py_wavelet make_wavelet(int p, bool sym) {
   const std::string family = sym ? "symlet" : "daubechies";
   if (p < 2 || p > 10)
      throw nb::value_error((family + ": p must be in [2, 10], but got " + std::to_string(p)).c_str());
   nb::gil_scoped_release release;
   return py_wavelet{sym ? bf::wavelet::symlet<double>(unsigned(p)) : bf::wavelet::daubechies<double>(unsigned(p)), family, unsigned(p)};
}

}  // namespace

void bind_filters(nb::module_& m) {
   m.def("trend_filter", [](nb::handle x, std::optional<int> k, std::optional<double> lambda) {
      if (k.has_value() != lambda.has_value())
         throw nb::value_error("trend_filter: give both k and lambda_, or neither for cross-validation");
      std::vector<double> y = to_vector(x, "trend_filter");
      return run("trend_filter", y.size(), [&] { return k ? bf::trend_filter(y, *k, *lambda) : bf::trend_filter(y); });
   }, nb::arg("x"), nb::arg("k") = nb::none(), nb::arg("lambda_") = nb::none(),
      nb::sig("def trend_filter(x, /, k=None, lambda_=None)"),
      "The trend filter of x: argmin_z ½ Σ (xᵢ − zᵢ)² + λ Σ |(Dᵏ⁺¹ z)ᵢ|, Dᵏ⁺¹ the (k + 1)-th difference.\n\n"
      "The fit is piecewise polynomial of degree k with knots chosen by the data: k = 0 is piecewise constant\n"
      "(total-variation denoising), k = 1 piecewise linear, k = 2 piecewise quadratic. λ = lambda_ ≥ 0; λ = 0\n"
      "returns x, and a large enough λ the least-squares polynomial of degree k. Give both k and lambda_, or\n"
      "neither: then k ∈ {0, 1, 2} and λ are chosen by cross-validation on the even- and odd-indexed samples\n"
      "(at least 4). Cross-validation assumes the noise is independent from sample to sample; with correlated\n"
      "noise it chooses too small a λ. If the noise estimate is 0 (noise-free data), x is returned unchanged.\n"
      "Returns a new float64 array."
      BMPY_URL("trend_filter"));

   nb::class_<py_wavelet>(m, "OrthogonalWavelet",
      "An orthonormal wavelet basis on an interval: a Daubechies or least-asymmetric (symlet) filter in the interior and\n"
      "orthogonal boundary rows at both ends, so polynomials of degree < p have no detail coefficients. Make one with\n"
      "daubechies(p) or symlet(p); copies share the boundary rows, which are generated on first use."
      BMPY_URL("wavelet_denoise"))
      .def_ro("vanishing_moments", &py_wavelet::p, "p, the number of vanishing moments.")
      .def("__repr__", [](py_wavelet const& w) { return w.family + "(" + std::to_string(w.p) + ")"; });
   m.def("daubechies", [](int p) { return make_wavelet(p, false); }, nb::arg("p"), nb::sig("def daubechies(p, /) -> OrthogonalWavelet"),
      "Daubechies' minimum-phase wavelet with p ∈ [2, 10] vanishing moments, on an interval." BMPY_URL("wavelet_denoise"));
   m.def("symlet", [](int p) { return make_wavelet(p, true); }, nb::arg("p"), nb::sig("def symlet(p, /) -> OrthogonalWavelet"),
      "The least-asymmetric Daubechies wavelet (symlet) with p ∈ [2, 10] vanishing moments, on an interval. sym7 and sym10\n"
      "are more symmetric than PyWavelets' (a different root choice); the others agree." BMPY_URL("wavelet_denoise"));

   m.def("wavelet_denoise", [](nb::handle x, std::optional<py_wavelet> const& w, std::optional<std::string> const& rule,
                               std::optional<double> sigma) {
      std::vector<double> y = to_vector(x, "wavelet_denoise");
      if (!w && (rule || sigma))
         throw nb::value_error("wavelet_denoise: rule and sigma need a wavelet; give none of the three for cross-validation");
      const bf::threshold_rule r = rule ? parse_rule(*rule) : bf::threshold_rule::wiener;
      const double s = sigma ? *sigma : std::numeric_limits<double>::quiet_NaN();
      const std::size_t n = y.size();
      return run("wavelet_denoise", n, [&] { return w ? bf::wavelet_denoise(y, w->w, r, s) : bf::wavelet_denoise(y); });
   }, nb::arg("x"), nb::arg("wavelet") = nb::none(), nb::arg("rule") = nb::none(), nb::arg("sigma") = nb::none(),
      nb::sig("def wavelet_denoise(x, /, wavelet=None, rule=None, sigma=None)"),
      "Denoises samples of a function on an interval by shrinking their wavelet coefficients in the basis `wavelet`.\n\n"
      "rule: 'wiener' (the default with a wavelet; empirical Wiener weights from a pilot: the first of sym8, db4, sym4,\n"
      "db2 that differs from the wavelet and fits), 'universal_soft', 'sure', or 'eb_posterior_median' (empirical-Bayes\n"
      "posterior median, Laplace prior of scale 1/2, as R's EbayesThresh).\n"
      "sigma: the noise standard deviation; None estimates it from the finest details.\n"
      "len(x) ≥ 4p − 1 (any length). Without a wavelet (and without rule and sigma), the wavelet (db2..db10,\n"
      "sym4..sym10) and the rule are chosen by cross-validation on the even- and odd-indexed samples (len(x) ≥ 14).\n"
      "Returns a new float64 array."
      BMPY_URL("wavelet_denoise"));

   for (bool forward : {true, false}) {
      char const* name = forward ? "wavelet_transform" : "inverse_wavelet_transform";
      m.def(name, [forward, name](nb::handle x, py_wavelet const& w, std::optional<std::size_t> levels) {
         std::vector<double> y = to_vector(x, name);
         const std::size_t n = y.size(), L = levels ? *levels : std::size_t(-1);
         return run(name, n, [&] { return forward ? bf::wavelet_transform(y, w.w, L) : bf::inverse_wavelet_transform(y, w.w, L); });
      }, nb::arg(forward ? "x" : "c"), nb::arg("wavelet"), nb::arg("levels") = nb::none(),
         nb::sig(forward ? "def wavelet_transform(x, wavelet, /, levels=None)" : "def inverse_wavelet_transform(c, wavelet, /, levels=None)"),
         forward
         ? "The coefficients of x in the orthonormal basis `wavelet`, as one float64 array of the same length, coarsest first:\n"
           "[a_L | d_L | d_{L-1} | ... | d_1], d_1 the finest details. Level l takes n_{l-1} values to ceil(n_{l-1}/2) scaling\n"
           "and floor(n_{l-1}/2) wavelet coefficients (n_0 = len(x)), so d_l is c[n_l:n_{l-1}].\n"
           "levels: L, at most (and by default, None) as many as len(x) allows, ((n - 1) // (4*p - 2)).bit_length() for\n"
           "n = len(x); 0 returns x. len(x) ≥ 4p − 1."
         : "The inverse of wavelet_transform(x, wavelet, levels): x from its coefficients c, with the same levels.",
         BMPY_URL("wavelet_denoise"));
   }

   m.def("denoise", [](nb::handle x) {
      std::vector<double> y = to_vector(x, "denoise");
      const std::size_t n = y.size();
      return run("denoise", n, [&] { return bf::denoise(y); });
   }, nb::arg("x"), nb::sig("def denoise(x, /)"),
      "Whichever of wavelet_denoise(x) and trend_filter(x) (both cross-validated) has the lower even/odd\n"
      "cross-validation score. Returns a new float64 array."
      BMPY_URL("wavelet_denoise"));

   // boost.math.filters.signals: the test signals, as a nested module (submodule() makes top-level ones)
   std::string const full = nb::cast<std::string>(m.attr("__name__")) + ".signals";
   nb::module_ signals = nb::steal<nb::module_>(PyModule_New(full.c_str()));
   signals.doc() = "Donoho and Johnstone's (1994) test signals for denoising, sampled at t_i = i/n, i = 0..n-1, and scaled to\n"
                   "standard deviation 7, so that adding unit-variance noise gives their signal-to-noise ratio of 7:\n\n"
                   "    f = signals.doppler()\n"
                   "    y = f + np.random.default_rng(1).normal(size=f.size)\n"
                   "    x = filters.denoise(y)  # np.mean((x - f)**2) is the MSE in units of the noise variance";
   m.attr("signals") = signals;
   nb::module_::import_("sys").attr("modules")[full.c_str()] = signals;
   struct { char const* name; char const* doc; } const list[] = {
      {"blocks", "Piecewise constant, with 11 jumps."},
      {"bumps", "Eleven narrow peaks of different widths on a zero baseline."},
      {"heavisine", "A sinusoid with jumps at t = 0.3 and t = 0.72."},
      {"doppler", "sqrt(t(1 - t)) sin(2.1 pi / (t + 0.05)): its frequency rises steeply towards t = 0."},
   };
   for (auto const& s : list) {
      std::string const name = s.name;
      signals.def(s.name, [name](std::size_t n) { return vector_to_array(test_signal(name, n)); }, nb::arg("n") = 2048,
         nb::sig(("def " + name + "(n: int = 2048) -> numpy.ndarray[dtype=float64]").c_str()),
         (std::string(s.doc) + " n ≥ 2 samples, scaled to standard deviation 7; returns a new float64 array.").c_str());
   }
}

}  // namespace bmpy
