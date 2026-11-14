// Copyright 2026 Nick Thompson
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// Embind bindings for the WebAssembly build: statistics, and the trend and wavelet filters.

#include "bindings.hpp"

#include <boost/math/filters/trend_filter.hpp>
#include <boost/math/filters/wavelet_denoise.hpp>
#include <boost/math/statistics/anderson_darling.hpp>
#include <boost/math/statistics/linear_regression.hpp>
#include <boost/math/statistics/ljung_box.hpp>
#include <boost/math/statistics/runs_test.hpp>
#include <boost/math/statistics/signal_statistics.hpp>
#include <boost/math/statistics/t_test.hpp>
#include <boost/math/statistics/z_test.hpp>

#include <cstdint>
#include <utility>

namespace
{
namespace st = boost::math::statistics;
namespace bf = bm::filters;
using vec = std::vector<double>;

// The filters report a bad argument (a length, lambda, sigma, or a sample that is not finite) through the error
// policy, which sets errno and fills the result with NaN: throw instead, as for any other bad argument.
template <class F>
val filtered(const char* name, F f)
{
   errno = 0;
   const vec y = f();
   if (errno == EDOM)
      range_error(std::string(name) + ": the samples must be finite and long enough for the filter, and its parameters in range.");
   return to_js(y);
}

// A wavelet is named by its family, "daubechies" or "symlet", and p = 2..10 vanishing moments.
bf::orthogonal_wavelet<double> wavelet(const val& family, double p)
{
   const std::string name = text(family, "The family");
   if (name != "daubechies" && name != "symlet")
      range_error("The family must be daubechies or symlet.");
   const unsigned moments = unsigned(whole(p, "p", 2, 10));
   return name == "symlet" ? bf::wavelet::symlet<double>(moments) : bf::wavelet::daubechies<double>(moments);
}

const char* const rule_names[] = {"wiener", "universal_soft", "sure", "eb_posterior_median"};

bf::threshold_rule rule(const val& name)
{
   const std::string s = text(name, "The rule");
   for (unsigned i = 0; i < 4; ++i)
      if (s == rule_names[i])
         return static_cast<bf::threshold_rule>(i);
   range_error("The rule must be wiener, universal_soft, sure or eb_posterior_median.");
}

// sigma, the standard deviation of the noise, is estimated from the data if it is NaN.
val wavelet_denoise(const val& x, const val& family, double p, const val& rule_name, double sigma)
{
   const vec v = numbers(x);
   const auto w = wavelet(family, p);
   const bf::threshold_rule r = rule(rule_name);
   return filtered("wavelet_denoise", [&] { return bf::wavelet_denoise(v, w, r, sigma); });
}

// levels = -1 means as many as the length allows.
template <bool Forward>
val wavelet_transform(const val& x, const val& family, double p, double levels)
{
   const vec v = numbers(x);
   const auto w = wavelet(family, p);
   const std::size_t most = bf::detail::wavelet_levels(w, v.size());
   const std::size_t L = levels == -1 ? std::size_t(-1) : count(levels, "levels", 0, double(most));
   return filtered("wavelet_transform", [&] { return Forward ? bf::wavelet_transform(v, w, L) : bf::inverse_wavelet_transform(v, w, L); });
}

vec finite_samples(const val& x)
{
   vec y = numbers(x);
   if (!std::all_of(y.begin(), y.end(), [](double v) { return std::isfinite(v); }))
      range_error("The samples must be finite.");
   return y;
}

// What trend_filter(x) chooses by cross-validation, and its score: {k, lambda, score}. Data without noise, which
// trend_filter(x) returns unchanged, give k = -1, lambda = 0 and no score.
val trend_filter_select(const val& x)
{
   const vec y = finite_samples(x);
   if (y.size() < 4)
      range_error("Need at least 4 samples.");
   const double sigma = bf::detail::trend_filter_sigma(y);
   double score = NAN;
   // A fit that fails scores as the error policy's NaN does, rather than throwing.
   const std::pair<int, double> kl = sigma == 0 ? std::pair<int, double>(-1, 0.0)
      : bf::detail::trend_filter_select(y, sigma, &score, bm::policies::policy<bm::policies::evaluation_error<bm::policies::ignore_error>>());
   val o = val::object();
   o.set("k", kl.first);
   o.set("lambda", kl.second);
   o.set("score", score);
   return o;
}

// What wavelet_denoise(x) chooses by cross-validation, and its score: {family, p, rule, score}.
val wavelet_denoise_select(const val& x)
{
   const vec y = finite_samples(x);
   bf::detail::wavelet_status status;
   const bf::detail::wavelet_choice<double> c = bf::detail::wavelet_select(y, status);
   if (status != bf::detail::wavelet_status::ok)
      range_error("Need at least 14 samples.");
   const bool symlet = c.wavelet >= bf::detail::wavelet_index(true, 2);
   val o = val::object();
   o.set("family", symlet ? "symlet" : "daubechies");
   o.set("p", c.wavelet - bf::detail::wavelet_index(symlet, 2) + 2);
   o.set("rule", rule_names[std::size_t(c.rule)]);
   o.set("score", c.score);
   return o;
}

// The statistics assert or throw on too few values, so check first.
std::vector<double> sample(const val& x, std::size_t at_least)
{
   std::vector<double> v = numbers(x);
   if (v.size() < at_least)
      range_error("The array needs at least " + std::to_string(at_least) + (at_least == 1 ? " value." : " values."));
   return v;
}

std::vector<double> same_size(const val& y, const std::vector<double>& x)
{
   std::vector<double> v = sample(y, x.size());
   if (v.size() != x.size())
      range_error("The arrays must have the same length.");
   return v;
}

// Sorting and selection need an order, which NaN breaks.
bool has_nan(const std::vector<double>& v)
{
   return std::any_of(v.begin(), v.end(), [](double x) { return std::isnan(x); });
}

val test_result(std::pair<double, double> r)
{
   val o = val::object();
   o.set("statistic", r.first);
   o.set("p_value", r.second);
   return o;
}

// The statistics of one array, which needs at least n values.
#define S1(f, n) entry(#f, +[](const val& x) { std::vector<double> v = sample(x, n); return double(st::f(v)); })
#define S1_ORDER(f, n) entry(#f, +[](const val& x) { std::vector<double> v = sample(x, n); return has_nan(v) ? NAN : double(st::f(v)); })
} // namespace

EMSCRIPTEN_BINDINGS(boost_math_statistics_filters)
{
   // k is 0, 1 or 2, and lambda >= 0; without them, both are chosen by cross-validation.
   entry("trend_filter", +[](const val& x, double k, double lambda) {
      const vec v = numbers(x);
      const int order = int(whole(k, "k", 0, 2));
      return filtered("trend_filter", [&] { return bf::trend_filter(v, order, lambda); }); },
      +[](const val& x) { const vec v = numbers(x); return filtered("trend_filter", [&] { return bf::trend_filter(v); }); });
   // The wavelet and the rule are chosen by cross-validation unless given.
   entry("wavelet_denoise", +[](const val& x) { const vec v = numbers(x); return filtered("wavelet_denoise", [&] { return bf::wavelet_denoise(v); }); },
      +[](const val& x, const val& f, double p) { return wavelet_denoise(x, f, p, val("wiener"), NAN); },
      +[](const val& x, const val& f, double p, const val& r) { return wavelet_denoise(x, f, p, r, NAN); }, &wavelet_denoise);
   entry("wavelet_transform", +[](const val& x, const val& f, double p) { return wavelet_transform<true>(x, f, p, -1); }, &wavelet_transform<true>);
   entry("inverse_wavelet_transform", +[](const val& c, const val& f, double p) { return wavelet_transform<false>(c, f, p, -1); },
      &wavelet_transform<false>);
   entry("denoise", +[](const val& x) { const vec v = numbers(x); return filtered("denoise", [&] { return bf::denoise(v); }); });
   entry("trend_filter_select", &trend_filter_select);
   entry("wavelet_denoise_select", &wavelet_denoise_select);
   // The number of levels wavelet_transform(x, family, p) uses for n samples: 0 if n is too short.
   entry("wavelet_levels", +[](double n, const val& family, double p) {
      return double(bf::detail::wavelet_levels(wavelet(family, p), count(n, "n"))); });

   S1(mean, 1); S1(variance, 1); S1(sample_variance, 2); S1(skewness, 1); S1(kurtosis, 1); S1(excess_kurtosis, 1);
   S1_ORDER(median, 1); S1_ORDER(median_absolute_deviation, 1); S1_ORDER(interquartile_range, 3);
   S1_ORDER(gini_coefficient, 1); S1_ORDER(sample_gini_coefficient, 2); S1(hoyer_sparsity, 2);
   entry("mean_and_sample_variance", +[](const val& x) {
      std::vector<double> v = sample(x, 2);
      auto [m, s2] = st::mean_and_sample_variance(v);
      val o = val::object();
      o.set("mean", m);
      o.set("sample_variance", s2);
      return o;
   });
   entry("covariance", +[](const val& x, const val& y) { std::vector<double> u = sample(x, 1); return st::covariance(u, same_size(y, u)); });
   entry("correlation_coefficient", +[](const val& x, const val& y) {
      std::vector<double> u = sample(x, 1);
      return st::correlation_coefficient(u, same_size(y, u));
   });
   entry("simple_ordinary_least_squares", +[](const val& x, const val& y) {
      std::vector<double> u = sample(x, 2), v = same_size(y, u);
      if (st::variance(u) <= 0)
         range_error("x must not be constant.");
      auto [c0, c1, r2] = st::simple_ordinary_least_squares_with_R_squared(u, v);
      val o = val::object();
      o.set("intercept", c0);
      o.set("slope", c1);
      o.set("r_squared", r2);
      return o;
   });
   entry("one_sample_t_test", +[](const val& x, double mu) { return test_result(st::one_sample_t_test(sample(x, 2), mu)); });
   entry("one_sample_z_test", +[](const val& x, double mu) { return test_result(st::one_sample_z_test(sample(x, 2), mu)); });
   entry("two_sample_t_test", +[](const val& x, const val& y) { return test_result(st::two_sample_t_test(sample(x, 2), sample(y, 2))); });
   entry("two_sample_z_test", +[](const val& x, const val& y) { return test_result(st::two_sample_z_test(sample(x, 2), sample(y, 2))); });
   entry("paired_samples_t_test", +[](const val& x, const val& y) {
      std::vector<double> u = sample(x, 2);
      return test_result(st::paired_samples_t_test(u, same_size(y, u)));
   });
   entry("runs_above_and_below_threshold", +[](const val& x, double t) {
      std::vector<double> v = sample(x, 2);
      return test_result(has_nan(v) ? std::pair<double, double>(NAN, NAN) : st::runs_above_and_below_threshold(v, t));
   });
   entry("runs_above_and_below_median", +[](const val& x) {
      std::vector<double> v = sample(x, 2);
      return test_result(has_nan(v) ? std::pair<double, double>(NAN, NAN) : st::runs_above_and_below_median(v));
   });
   // Sorts a copy of the data, and estimates the mean and standard deviation from it.
   entry("anderson_darling_normality_statistic", +[](const val& x) -> double {
      std::vector<double> v = sample(x, 2);
      if (has_nan(v))
         return NAN;
      std::sort(v.begin(), v.end());
      return st::anderson_darling_normality_statistic(v);
   });
   // lags = -1 takes ceil(log(n)) lags.
   entry("ljung_box", +[](const val& x, double lags, double fit_dof) {
      std::vector<double> v = sample(x, 2);
      const double n = double(v.size()), l = lags == -1 ? std::ceil(std::log(n)) : lags;
      if (!(l == std::trunc(l) && l >= 1 && l < n))
         range_error("The lags must be -1 or an integer in [1, n).");
      if (!(fit_dof == std::trunc(fit_dof) && std::abs(fit_dof) < 1e15))
         range_error("fit_dof must be an integer.");
      return test_result(st::ljung_box(v, int64_t(l), int64_t(fit_dof)));
   });
}
