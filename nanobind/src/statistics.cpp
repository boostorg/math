// Copyright 2026 Nick Thompson
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// boost.math.statistics: the sequential functions of boost/math/statistics/ on 1-D float64 data. Anything
// np.ascontiguousarray understands is accepted; a contiguous float64 array is used without a copy. Functions that
// sort or partition their argument (median, ...) work on a private copy, so the caller's array is never modified.
// Several values come back as a NamedTuple. Errors follow boost.math.errstate, as for a scalar argument.

#include "bindings.hpp"

#include <boost/math/statistics/anderson_darling.hpp>
#include <boost/math/statistics/bivariate_statistics.hpp>
#include <boost/math/statistics/chatterjee_correlation.hpp>
#include <boost/math/statistics/linear_regression.hpp>
#include <boost/math/statistics/ljung_box.hpp>
#include <boost/math/statistics/runs_test.hpp>
#include <boost/math/statistics/signal_statistics.hpp>
#include <boost/math/statistics/t_test.hpp>
#include <boost/math/statistics/univariate_statistics.hpp>
#include <boost/math/statistics/z_test.hpp>

#include <algorithm>
#include <cmath>
#include <span>
#include <string>
#include <tuple>
#include <vector>

namespace bmpy {

namespace {

namespace bs = boost::math::statistics;
using span = std::span<double const>;
using vec = std::vector<double>;
constexpr double nan = std::numeric_limits<double>::quiet_NaN();

// x as a 1-D float64 view (a converted copy unless it already is one), with at least min_size elements.
// Boost.Math asserts on too-short input, so the length is checked here.
struct data {
   nb::ndarray<double const, nb::ndim<1>, nb::c_contig, nb::device::cpu> a;
   span s() const { return {a.data(), a.shape(0)}; }
   vec copy() const { return {a.data(), a.data() + a.shape(0)}; }
};

data as_data(char const* name, nb::handle x, std::size_t min_size) {
   nb::module_& np = numpy();
   nb::object a = np.attr("ascontiguousarray")(x, nb::arg("dtype") = np.attr("float64"));
   if (nb::cast<int>(a.attr("ndim")) != 1)
      throw nb::value_error((std::string(name) + ": expected a 1-D array").c_str());
   data d{nb::cast<decltype(data::a)>(a)};
   if (d.a.shape(0) < min_size)
      throw nb::value_error((std::string(name) + ": expected at least " + std::to_string(min_size) + " elements").c_str());
   return d;
}

std::pair<data, data> as_pair(char const* name, nb::handle x, nb::handle y, std::size_t min_size, bool same_length = true) {
   std::pair<data, data> d{as_data(name, x, min_size), as_data(name, y, min_size)};
   if (same_length && d.first.a.shape(0) != d.second.a.shape(0))
      throw nb::value_error((std::string(name) + ": the two samples must have the same length").c_str());
   return d;
}

// f() without the GIL, with the errors handled as for a scalar argument. A float comes back as a Python float, a
// pair or tuple as an instance of the NamedTuple `type`, NaN-filled if the evaluation failed.
template <class F>
nb::object evaluate(char const* name, nb::handle type, F const& f) {
   using R = decltype(f());
   error_actions actions = current_error_actions;
   error_record errors;
   R r;
   {
      recording rec(errors);
      nb::gil_scoped_release release;
      try {
         r = f();
      } catch (std::exception const& e) {
         record_exception(e);
         if constexpr (std::is_floating_point_v<R>)
            r = nan;
         else
            std::apply([](auto&... e) { ((e = nan), ...); }, r);
      }
   }
   if (errors.any())
      report(name, errors, 1, true, false, actions);
   if constexpr (std::is_floating_point_v<R>)
      return nb::float_(r);
   else
      return std::apply([&](auto... e) { return type(e...); }, r);
}

// name(x, /) and name(x, y, /), calling f(data) or f(data, data).
template <class F>
void def_x(nb::module_& m, char const* name, std::size_t min_size, nb::handle type, F f, char const* doc) {
   m.def(name, [=](nb::handle x) {
      data d = as_data(name, x, min_size);
      return evaluate(name, type, [&] { return f(d); });
   }, nb::arg("x"), nb::sig((std::string("def ") + name + "(x, /)").c_str()), doc);
}

template <class F>
void def_xy(nb::module_& m, char const* name, std::size_t min_size, nb::handle type, F f, char const* doc, bool same_length = true) {
   m.def(name, [=](nb::handle x, nb::handle y) {
      auto [u, v] = as_pair(name, x, y, min_size, same_length);
      return evaluate(name, type, [&] { return f(u, v); });
   }, nb::arg("x"), nb::arg("y"), nb::sig((std::string("def ") + name + "(x, y, /)").c_str()), doc);
}

}  // namespace

void bind_statistics(nb::module_& m) {
   nb::handle real;
   nb::object test = named_tuple(m, "TestResult", "statistic pvalue");
   nb::object mean_var = named_tuple(m, "MeanVariance", "mean sample_variance");
   nb::object moments = named_tuple(m, "Moments", "mean m2 m3 m4");
   nb::object means_cov = named_tuple(m, "MeansCovariance", "mean_x mean_y covariance");
   nb::object fit = named_tuple(m, "LinearFit", "intercept slope");
   nb::object fit_r2 = named_tuple(m, "LinearFitR2", "intercept slope r_squared");

   // univariate_statistics.hpp
#define BMPY_UNIVARIATE BMPY_URL("univariate_statistics")
   def_x(m, "mean", 1, real, [](data const& d) { return bs::mean(d.s()); }, "Arithmetic mean." BMPY_UNIVARIATE);
   def_x(m, "variance", 1, real, [](data const& d) { return bs::variance(d.s()); },
         "Population variance, dividing by n." BMPY_UNIVARIATE);
   def_x(m, "sample_variance", 2, real, [](data const& d) { return bs::sample_variance(d.s()); },
         "Sample variance, dividing by n - 1." BMPY_UNIVARIATE);
   def_x(m, "mean_and_sample_variance", 2, mean_var, [](data const& d) { return bs::mean_and_sample_variance(d.s()); },
         "MeanVariance(mean, sample_variance), in one pass." BMPY_UNIVARIATE);
   def_x(m, "first_four_moments", 1, moments, [](data const& d) { return bs::first_four_moments(d.s()); },
         "Moments(mean, m2, m3, m4): the mean and the central moments sum((x - mean)^k)/n." BMPY_UNIVARIATE);
   def_x(m, "skewness", 1, real, [](data const& d) { return bs::skewness(d.s()); }, "Population skewness." BMPY_UNIVARIATE);
   def_x(m, "kurtosis", 1, real, [](data const& d) { return bs::kurtosis(d.s()); },
         "Population kurtosis (3 for a normal distribution)." BMPY_UNIVARIATE);
   def_x(m, "excess_kurtosis", 1, real, [](data const& d) { return bs::excess_kurtosis(d.s()); },
         "kurtosis - 3." BMPY_UNIVARIATE);
   def_x(m, "median", 1, real, [](data const& d) { vec w = d.copy(); return bs::median(w); }, "Median." BMPY_UNIVARIATE);
   def_x(m, "interquartile_range", 3, real, [](data const& d) { vec w = d.copy(); return bs::interquartile_range(w); },
         "Interquartile range Q3 - Q1." BMPY_UNIVARIATE);
   def_x(m, "gini_coefficient", 1, real, [](data const& d) { vec w = d.copy(); return bs::gini_coefficient(w); },
         "Gini coefficient of non-negative data." BMPY_UNIVARIATE);
   def_x(m, "sample_gini_coefficient", 2, real, [](data const& d) { vec w = d.copy(); return bs::sample_gini_coefficient(w); },
         "n/(n - 1) times gini_coefficient." BMPY_UNIVARIATE);
   m.def("median_absolute_deviation", [](nb::handle x, double center) {
      data d = as_data("median_absolute_deviation", x, 1);
      return evaluate("median_absolute_deviation", nb::handle(), [&] {
         vec w = d.copy();
         return bs::median_absolute_deviation(w, center);
      });
   }, nb::arg("x"), nb::arg("center") = nan, nb::sig("def median_absolute_deviation(x, /, center=nan)"),
         "median(|x - center|); center NaN means the median of x." BMPY_UNIVARIATE);
   m.def("mode", [](nb::handle x) {
      vec w = as_data("mode", x, 1).copy();
      if (std::any_of(w.begin(), w.end(), [](double v) { return v != v; }))  // NaN != NaN makes Boost loop forever.
         throw nb::value_error("mode: x contains NaN");
      nb::list out;
      for (double v : bs::mode(w))
         out.append(v);
      return out;
   }, nb::arg("x"), nb::sig("def mode(x, /)"), "The most frequent values of x, as a sorted list." BMPY_UNIVARIATE);
#undef BMPY_UNIVARIATE

   // bivariate_statistics.hpp
   def_xy(m, "covariance", 1, real, [](data const& u, data const& v) { return bs::covariance(u.s(), v.s()); },
          "Population covariance." BMPY_URL("bivariate_statistics"));
   def_xy(m, "means_and_covariance", 1, means_cov,
          [](data const& u, data const& v) { return bs::means_and_covariance(u.s(), v.s()); },
          "MeansCovariance(mean_x, mean_y, covariance), in one pass." BMPY_URL("bivariate_statistics"));
   def_xy(m, "correlation_coefficient", 1, real,
          [](data const& u, data const& v) { return bs::correlation_coefficient(u.s(), v.s()); },
          "Pearson correlation coefficient." BMPY_URL("bivariate_statistics"));

   // signal_statistics.hpp
#define BMPY_SIGNAL BMPY_URL("signal_statistics")
   def_x(m, "absolute_gini_coefficient", 2, real, [](data const& d) { vec w = d.copy(); return bs::absolute_gini_coefficient(w); },
         "Gini coefficient of |x|, a measure of sparsity in [0, 1]." BMPY_SIGNAL);
   def_x(m, "sample_absolute_gini_coefficient", 2, real,
         [](data const& d) { vec w = d.copy(); return bs::sample_absolute_gini_coefficient(w); },
         "n/(n - 1) times absolute_gini_coefficient." BMPY_SIGNAL);
   def_x(m, "hoyer_sparsity", 2, real, [](data const& d) { return bs::hoyer_sparsity(d.s().begin(), d.s().end()); },
         "Hoyer sparsity, in [0, 1]." BMPY_SIGNAL);
   def_xy(m, "oracle_snr", 1, real, [](data const& u, data const& v) { return bs::oracle_snr(u.s(), v.s()); },
          "Signal-to-noise ratio of the noisy signal y against the known signal x." BMPY_SIGNAL);
   def_xy(m, "oracle_snr_db", 1, real, [](data const& u, data const& v) { return bs::oracle_snr_db(u.s(), v.s()); },
          "oracle_snr in decibels." BMPY_SIGNAL);
   def_xy(m, "mean_invariant_oracle_snr", 1, real,
          [](data const& u, data const& v) { return bs::mean_invariant_oracle_snr(u.s(), v.s()); },
          "oracle_snr with the mean of the signal x removed." BMPY_SIGNAL);
   def_xy(m, "mean_invariant_oracle_snr_db", 1, real,
          [](data const& u, data const& v) { return bs::mean_invariant_oracle_snr_db(u.s(), v.s()); },
          "mean_invariant_oracle_snr in decibels." BMPY_SIGNAL);
   for (bool db : {false, true}) {
      char const* name = db ? "m2m4_snr_estimator_db" : "m2m4_snr_estimator";
      m.def(name, [name, db](nb::handle x, double ks, double kn) {
         if (!(ks > 0 && kn > 0))
            throw nb::value_error((std::string(name) + ": the kurtoses must be positive").c_str());
         data d = as_data(name, x, 1);
         return evaluate(name, nb::handle(), [&] {
            span s = d.s();
            return db ? bs::m2m4_snr_estimator_db(s.begin(), s.end(), ks, kn) : bs::m2m4_snr_estimator(s.begin(), s.end(), ks, kn);
         });
      }, nb::arg("x"), nb::arg("signal_kurtosis") = 1.0, nb::arg("noise_kurtosis") = 3.0,
            nb::sig((std::string("def ") + name + "(x, /, signal_kurtosis=1.0, noise_kurtosis=3.0)").c_str()),
            db ? "m2m4_snr_estimator in decibels." BMPY_SIGNAL
               : "Blind SNR estimate of a noisy signal x from its second and fourth moments." BMPY_SIGNAL);
   }
#undef BMPY_SIGNAL

   // anderson_darling.hpp, ljung_box.hpp, runs_test.hpp
   m.def("anderson_darling_normality_statistic", [](nb::handle x, double mu, double sd) {
      // sample_variance asserts n >= 2.
      data d = as_data("anderson_darling_normality_statistic", x, std::isnan(sd) ? 2 : 1);
      return evaluate("anderson_darling_normality_statistic", nb::handle(), [&] {
         vec w = d.copy();
         std::sort(w.begin(), w.end());
         return bs::anderson_darling_normality_statistic(w, mu, sd);
      });
   }, nb::arg("x"), nb::arg("mu") = nan, nb::arg("sd") = nan,
         nb::sig("def anderson_darling_normality_statistic(x, /, mu=nan, sd=nan)"),
         "Anderson-Darling statistic of x against N(mu, sd); NaN means the sample mean or standard deviation."
         BMPY_URL("anderson_darling"));
   m.def("ljung_box", [test](nb::handle x, std::int64_t lags, std::int64_t fit_dof) {
      data d = as_data("ljung_box", x, 2);
      return evaluate("ljung_box", test, [&] { return bs::ljung_box(d.s(), lags, fit_dof); });
   }, nb::arg("x"), nb::arg("lags") = -1, nb::arg("fit_dof") = 0, nb::sig("def ljung_box(x, /, lags=-1, fit_dof=0)"),
         "Ljung-Box test for autocorrelation: TestResult(statistic, pvalue); lags = -1 means ceil(log n)." BMPY_URL("ljung_box"));
   m.def("runs_above_and_below_threshold", [test](nb::handle x, double threshold) {
      data d = as_data("runs_above_and_below_threshold", x, 2);
      return evaluate("runs_above_and_below_threshold", test,
                      [&] { return bs::runs_above_and_below_threshold(d.s(), threshold); });
   }, nb::arg("x"), nb::arg("threshold"), nb::sig("def runs_above_and_below_threshold(x, threshold, /)"),
         "Wald-Wolfowitz runs test of x about threshold: TestResult(statistic, pvalue)." BMPY_URL("runs_test"));
   def_x(m, "runs_above_and_below_median", 2, test, [](data const& d) { return bs::runs_above_and_below_median(d.copy()); },
         "runs_above_and_below_threshold about the median of x." BMPY_URL("runs_test"));

   // t_test.hpp, z_test.hpp
   for (bool z : {false, true}) {
      char const* name = z ? "one_sample_z_test" : "one_sample_t_test";
      m.def(name, [test, name, z](nb::handle x, double assumed_mean) {
         data d = as_data(name, x, 2);
         return evaluate(name, test, [&] {
            return z ? bs::one_sample_z_test(d.s(), assumed_mean) : bs::one_sample_t_test(d.s(), assumed_mean);
         });
      }, nb::arg("x"), nb::arg("assumed_mean"), nb::sig((std::string("def ") + name + "(x, assumed_mean, /)").c_str()),
            z ? "Two-sided one-sample z-test of mean(x) = assumed_mean: TestResult(statistic, pvalue)." BMPY_URL("z_test")
              : "Two-sided one-sample t-test of mean(x) = assumed_mean: TestResult(statistic, pvalue)." BMPY_URL("t_test"));
   }
   def_xy(m, "two_sample_t_test", 2, test, [](data const& u, data const& v) { return bs::two_sample_t_test(u.s(), v.s()); },
          "Two-sided two-sample t-test of equal means, Welch's when one standard deviation exceeds twice the other: "
          "TestResult(statistic, pvalue)." BMPY_URL("t_test"), false);
   def_xy(m, "paired_samples_t_test", 2, test, [](data const& u, data const& v) { return bs::paired_samples_t_test(u.s(), v.s()); },
          "Two-sided t-test of mean(x - y) = 0: TestResult(statistic, pvalue)." BMPY_URL("t_test"));
   def_xy(m, "two_sample_z_test", 2, test, [](data const& u, data const& v) { return bs::two_sample_z_test(u.s(), v.s()); },
          "Two-sided two-sample z-test of equal means: TestResult(statistic, pvalue)." BMPY_URL("z_test"), false);

   // linear_regression.hpp
   def_xy(m, "simple_ordinary_least_squares", 2, fit,
          [](data const& u, data const& v) { return bs::simple_ordinary_least_squares(u.s(), v.s()); },
          "Least-squares line y ~ intercept + slope x: LinearFit(intercept, slope)." BMPY_URL("linear_regression"));
   def_xy(m, "simple_ordinary_least_squares_with_R_squared", 2, fit_r2,
          [](data const& u, data const& v) { return bs::simple_ordinary_least_squares_with_R_squared(u.s(), v.s()); },
          "LinearFitR2(intercept, slope, r_squared)." BMPY_URL("linear_regression"));

   // chatterjee_correlation.hpp
   m.def("chatterjee_correlation", [](nb::handle x, nb::handle y, std::size_t M) {
      auto [u, v] = as_pair("chatterjee_correlation", x, y, 2);
      if (!std::is_sorted(u.s().begin(), u.s().end()))
         throw nb::value_error("chatterjee_correlation: x must be sorted in increasing order");
      if (M < 1 || M > u.a.shape(0))
         throw nb::value_error("chatterjee_correlation: M must be in [1, n]");
      return evaluate("chatterjee_correlation", nb::handle(), [&] {
         return M == 1 ? bs::chatterjee_correlation(u.s(), v.s()) : bs::chatterjee_correlation_mnn(u.s(), v.s(), M);
      });
   }, nb::arg("x"), nb::arg("y"), nb::arg("M") = 1, nb::sig("def chatterjee_correlation(x, y, /, M=1)"),
         "Chatterjee's rank correlation of y on x, x sorted increasing; M > 1 gives Lin and Han's M-nearest-neighbour "
         "version." BMPY_URL("chatterjee_correlation"));
}

}  // namespace bmpy
