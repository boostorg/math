// Copyright 2026 Nick Thompson
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// Embind bindings for the WebAssembly build: special functions, distributions and FFTs.

#include "bindings.hpp"

#include <boost/math/distributions.hpp>
#include <boost/math/special_functions.hpp>
#include <boost/math/special_functions/spheroidal.hpp>
#include <boost/math/special_functions/struve.hpp>
#include <boost/math/tools/fft.hpp>

#include <limits>
#include <type_traits>

// f(a..., x) for each x of an array, returned as a Float64Array: one call from JavaScript instead of one per x.
template <class F, class... A>
val map_last(F f, const val& xs, A... a)
{
   std::vector<double> x = numbers(xs);
   for (double& v : x)
      v = f(a..., v);
   return to_js(x);
}

// An integer order or degree n of type I: NaN unless n is an integer that I can hold.
template <class I, class F>
double with_order(double n, F f)
{
   return n == std::trunc(n) && n >= std::numeric_limits<I>::min() && n <= std::numeric_limits<I>::max() ? double(f(I(n))) : NAN;
}

// Each function f also has f_array, whose last argument is an array: erf_array(xs), gamma_p_array(a, xs).
#define F1(f) function(#f, +[](double x) { return bm::f(x); }); \
   function(#f "_array", +[](const val& x) { return map_last([](double v) { return double(bm::f(v)); }, x); })
#define F2(f) function(#f, +[](double a, double x) { return bm::f(a, x); }); \
   function(#f "_array", +[](double a, const val& x) { return map_last([](double a, double v) { return double(bm::f(a, v)); }, x, a); })
#define F3(f) function(#f, +[](double a, double b, double x) { return bm::f(a, b, x); }); \
   function(#f "_array", +[](double a, double b, const val& x) { \
      return map_last([](double a, double b, double v) { return double(bm::f(a, b, v)); }, x, a, b); })
#define FN(f, I) function(#f, +[](double n, double x) { return with_order<I>(n, [x](I k) { return bm::f(k, x); }); }); \
   function(#f "_array", +[](double n, const val& x) { \
      return map_last([](double n, double v) { return with_order<I>(n, [v](I k) { return bm::f(k, v); }); }, x, n); })

// Struve functions f(v, x), with f_prime the derivative in x (the C++ derivative pointer).
#define STRUVE(f) F2(f); \
   function(#f "_prime", +[](double v, double x) { double d = NAN; bm::f(v, x, &d); return d; }); \
   function(#f "_prime_array", +[](double v, const val& x) { \
      return map_last([](double v, double x) { double d = NAN; bm::f(v, x, &d); return d; }, x, v); })

// Spheroidal wave functions f(n, m, c, x), with f_prime the derivative in x (the C++ derivative pointer).
#define SPHEROIDAL(f) \
   function(#f, +[](double n, double m, double c, double x) { \
      return with_order<unsigned>(n, [&](unsigned k) { return with_order<unsigned>(m, [&](unsigned j) { return bm::f(k, j, c, x); }); }); }); \
   function(#f "_prime", +[](double n, double m, double c, double x) { \
      return with_order<unsigned>(n, [&](unsigned k) { return with_order<unsigned>(m, [&](unsigned j) { double d = NAN; bm::f(k, j, c, x, &d); return d; }); }); }); \
   function(#f "_array", +[](double n, double m, double c, const val& x) { \
      return map_last([](double n, double m, double c, double v) { \
         return with_order<unsigned>(n, [&](unsigned k) { return with_order<unsigned>(m, [&](unsigned j) { return bm::f(k, j, c, v); }); }); }, x, n, m, c); }); \
   function(#f "_prime_array", +[](double n, double m, double c, const val& x) { \
      return map_last([](double n, double m, double c, double v) { \
         return with_order<unsigned>(n, [&](unsigned k) { return with_order<unsigned>(m, [&](unsigned j) { double d = NAN; bm::f(k, j, c, v, &d); return d; }); }); }, x, n, m, c); })
#define SPHEROIDAL_EIGENVALUE(f) \
   function(#f, +[](double n, double m, double c) { \
      return with_order<unsigned>(n, [&](unsigned k) { return with_order<unsigned>(m, [&](unsigned j) { return bm::f(k, j, c); }); }); }); \
   function(#f "_array", +[](double n, double m, const val& c) { \
      return map_last([](double n, double m, double v) { \
         return with_order<unsigned>(n, [&](unsigned k) { return with_order<unsigned>(m, [&](unsigned j) { return bm::f(k, j, v); }); }); }, c, n, m); })

// Only the constructor checks the parameters, so remember whether it set errno, and give NaN for everything
// about an invalid distribution: mean(geometric_distribution(2.2)) would otherwise be -0.55.
template <class D>
struct checked
{
   D d;
   bool valid;
};

template <class D, class... Params>
checked<D> make(Params... params)
{
   errno = 0;
   D d(params...);
   return {d, errno == 0};
}

// Not every distribution has an entropy.
template <class D>
auto entropy_of(const D& d, int) -> decltype(entropy(d)) { return entropy(d); }
template <class D>
double entropy_of(const D&, long) { return NAN; }

template <class D, class... Params>
void distribution(const char* name)
{
   using C = checked<D>;
   // pdf(x) and pdf_array(xs), and so on.
#define AT(f, e) function(f, +[](const C& c, double x) { return c.valid ? double(e) : NAN; }) \
      .function(f "_array", +[](const C& c, const val& xs) { return map_last([&c](double x) { return c.valid ? double(e) : NAN; }, xs); })
   class_<C>(name)
      .constructor(&make<D, Params...>)
      .AT("pdf", pdf(c.d, x))
      .AT("cdf", cdf(c.d, x))
      .AT("ccdf", cdf(complement(c.d, x)))
      .AT("quantile", quantile(c.d, x))
#undef AT
      .function("mean", +[](const C& c) { return c.valid ? mean(c.d) : NAN; })
      .function("variance", +[](const C& c) { return c.valid ? variance(c.d) : NAN; })
      .function("entropy", +[](const C& c) { return c.valid ? entropy_of(c.d, 0) : NAN; })
      .function("support", +[](const C& c) {
         auto s = support(c.d);
         const double values[] = {c.valid ? s.first : NAN, c.valid ? s.second : NAN};
         return val::array(values, values + 2);
      });
}

std::size_t length(double n)
{
   if (!(n >= 1 && n <= (1 << 24) && n == std::trunc(n)))
      range_error("The length must be an integer in [1, 2^24].");
   return static_cast<std::size_t>(n);
}

// Complex data are interleaved (re, im) Float64Arrays.
struct fft
{
   explicit fft(double n) : f(length(n)), z(f.size()) {}
   val run(const val& x, bool inverse)
   {
      std::vector<double> v = numbers(x);
      if (v.size() != 2 * z.size())
         range_error("The input must hold 2n numbers.");
      for (std::size_t k = 0; k < z.size(); ++k)
         z[k] = cplx(v[2 * k], v[2 * k + 1]);
      inverse ? f.inverse_transform(z) : f.transform(z);
      return to_js(z);
   }
   bm::tools::fft<cplx> f;
   std::vector<cplx> z;
};

// n real values in, X[0], ..., X[n/2] out.
struct real_fft
{
   explicit real_fft(double n) : f(length(n)), z(f.size() / 2 + 1) {}
   val transform(const val& x)
   {
      std::vector<double> v = numbers(x);
      if (v.size() != f.size())
         range_error("The input must hold n numbers.");
      f.transform(v, z);
      return to_js(z);
   }
   bm::tools::real_fft<cplx> f;
   std::vector<cplx> z;
};

// Transforms each lane along one axis of a C-order array with the given shape, whose elements are In (double or
// cplx). make(n) gives the transform of one lane of length n, a function from std::vector<In> to std::vector<Out>.
template <class In, class Make>
val along_axis(const val& x, const val& shape, double axis, Make make)
{
   std::vector<double> v = numbers(x);
   std::vector<double> s = numbers(shape);
   const double d = double(s.size());
   if (!(axis == std::trunc(axis) && axis >= -d && axis < d))
      range_error("The axis must be an integer in [-ndim, ndim).");
   const std::size_t a = std::size_t(axis < 0 ? axis + d : axis);
   double before = 1, after = 1;
   for (std::size_t i = 0; i < s.size(); ++i)
   {
      if (!(s[i] >= 1 && s[i] == std::trunc(s[i])))
         range_error("The shape must hold positive integers.");
      if (i != a)
         (i < a ? before : after) *= s[i];
   }
   if (before * s[a] * after * (sizeof(In) / sizeof(double)) != double(v.size()))
      range_error(sizeof(In) == sizeof(double) ? "The data must hold product(shape) numbers."
                                               : "The data must hold 2 product(shape) numbers.");
   const std::size_t outer = std::size_t(before), n = length(s[a]), inner = std::size_t(after);
   auto t = make(n);
   const In* in = reinterpret_cast<const In*>(v.data());
   std::vector<In> lane(n);
   std::vector<std::decay_t<decltype(t(lane)[0])>> out;
   for (std::size_t o = 0; o < outer; ++o)
      for (std::size_t i = 0; i < inner; ++i)
      {
         for (std::size_t j = 0; j < n; ++j)
            lane[j] = in[(o * n + j) * inner + i];
         const auto& y = t(lane);
         out.resize(outer * y.size() * inner);
         for (std::size_t k = 0; k < y.size(); ++k)
            out[(o * y.size() + k) * inner + i] = y[k];
      }
   return to_js(out);
}

val fft_axis(const val& x, const val& shape, double axis, bool inverse)
{
   return along_axis<cplx>(x, shape, axis, [inverse](std::size_t n) {
      return [f = bm::tools::fft<cplx>(n), inverse](std::vector<cplx>& z) mutable -> std::vector<cplx>& {
         inverse ? f.inverse_transform(z) : f.transform(z);
         return z;
      };
   });
}

val real_fft_axis(const val& x, const val& shape, double axis)
{
   return along_axis<double>(x, shape, axis, [](std::size_t n) {
      return [f = bm::tools::real_fft<cplx>(n), z = std::vector<cplx>(n / 2 + 1)](std::vector<double>& v) mutable -> std::vector<cplx>& {
         f.transform(v, z);
         return z;
      };
   });
}

// P[k] = c |X[k]|^2 / n^2 for k = 0, ..., n/2, with c = 1 at k = 0 and the Nyquist k = n/2 and c = 2 otherwise:
// sum P = mean(x^2), and A cos(2 pi k j / n) gives P[k] = A^2 / 2.
val power_spectrum_axis(const val& x, const val& shape, double axis)
{
   return along_axis<double>(x, shape, axis, [](std::size_t n) {
      return [f = bm::tools::real_fft<cplx>(n), z = std::vector<cplx>(n / 2 + 1), p = std::vector<double>(n / 2 + 1), n](
                std::vector<double>& v) mutable -> std::vector<double>& {
         f.transform(v, z);
         for (std::size_t k = 0; k < z.size(); ++k)
            p[k] = (k == 0 || 2 * k == n ? 1 : 2) * std::norm(z[k]) / (double(n) * double(n));
         return p;
      };
   });
}

EMSCRIPTEN_BINDINGS(boost_math)
{
   F1(tgamma); F2(tgamma); F1(lgamma); F1(digamma); F1(trigamma); FN(polygamma, int); F1(tgamma1pm1);
   F2(tgamma_lower); F2(tgamma_ratio); F2(tgamma_delta_ratio);
   F2(gamma_p); F2(gamma_q); F2(gamma_p_inv); F2(gamma_q_inv); F2(gamma_p_derivative);
   F2(beta); F3(beta); F3(betac); F3(ibeta); F3(ibetac); F3(ibeta_inv); F3(ibetac_inv); F3(ibeta_derivative);
   F1(erf); F1(erfc); F1(erf_inv); F1(erfc_inv);
   F1(expint); FN(expint, unsigned); F1(zeta);
   F2(cyl_bessel_j); F2(cyl_neumann); F2(cyl_bessel_i); F2(cyl_bessel_k); FN(sph_bessel, unsigned); FN(sph_neumann, unsigned);
   F1(airy_ai); F1(airy_bi); F1(airy_ai_prime); F1(airy_bi_prime);
   F1(ellint_1); F2(ellint_1); F1(ellint_2); F2(ellint_2); F2(ellint_3); F3(ellint_3);
   F1(lambert_w0); F1(lambert_wm1); F2(owens_t);
   FN(legendre_p, int); FN(laguerre, unsigned); FN(hermite, unsigned); FN(chebyshev_t, unsigned);
   F1(log1p); F1(expm1); F1(cbrt); F1(sqrt1pm1); F2(powm1); F1(sinc_pi); F1(sin_pi); F1(cos_pi);
   F2(hypergeometric_0F1); F3(hypergeometric_1F1); F3(hypergeometric_2F0);
   F2(jacobi_sn); F2(jacobi_cn); F2(jacobi_dn); F2(jacobi_zeta); F2(heuman_lambda); F1(ellint_d); F2(ellint_d);
   F2(jacobi_theta1); F2(jacobi_theta2); F2(jacobi_theta3); F2(jacobi_theta4);
   SPHEROIDAL_EIGENVALUE(prolate_spheroidal_eigenvalue); SPHEROIDAL_EIGENVALUE(oblate_spheroidal_eigenvalue);
   SPHEROIDAL(prolate_spheroidal_angular); SPHEROIDAL(oblate_spheroidal_angular);
   SPHEROIDAL(prolate_spheroidal_radial_1); SPHEROIDAL(prolate_spheroidal_radial_2);
   SPHEROIDAL(oblate_spheroidal_radial_1); SPHEROIDAL(oblate_spheroidal_radial_2);
   STRUVE(struve_h); STRUVE(struve_l);
   // Fresnel C(x) and S(x); fresnel_cis gives both, [C, S], for the cost of one, and fresnel_cis_array interleaves them.
   F1(fresnel_c); F1(fresnel_s);
   function("fresnel_cis", +[](double x) { std::vector<cplx> z{bm::fresnel_cis(x)}; return to_js(z); });
   function("fresnel_cis_array", +[](const val& xs) {
      std::vector<double> x = numbers(xs);
      std::vector<cplx> z(x.size());
      for (std::size_t k = 0; k < x.size(); ++k)
         z[k] = bm::fresnel_cis(x[k]);
      return to_js(z); });
   F3(ellint_rf); F3(ellint_rd); F2(ellint_rc); F3(ellint_rg);
   F2(cyl_bessel_j_prime); F2(cyl_neumann_prime); F2(cyl_bessel_i_prime); F2(cyl_bessel_k_prime);
   // Zeros are numbered m = 1, 2, ...; the order n of gegenbauer and spherical_harmonic_* comes first.
   function("cyl_bessel_j_zero", +[](double v, double m) { return with_order<int>(m, [v](int k) { return bm::cyl_bessel_j_zero(v, k); }); });
   function("cyl_neumann_zero", +[](double v, double m) { return with_order<int>(m, [v](int k) { return bm::cyl_neumann_zero(v, k); }); });
   function("airy_ai_zero", +[](double m) { return with_order<int>(m, [](int k) { return bm::airy_ai_zero<double>(k); }); });
   function("airy_bi_zero", +[](double m) { return with_order<int>(m, [](int k) { return bm::airy_bi_zero<double>(k); }); });
   function("gegenbauer", +[](double n, double l, double x) { return with_order<unsigned>(n, [l, x](unsigned k) { return bm::gegenbauer(k, l, x); }); });
   function("spherical_harmonic_r", +[](double n, double m, double theta, double phi) {
      return with_order<unsigned>(n, [&](unsigned k) { return with_order<int>(m, [&](int j) { return bm::spherical_harmonic_r(k, j, theta, phi); }); }); });
   function("spherical_harmonic_i", +[](double n, double m, double theta, double phi) {
      return with_order<unsigned>(n, [&](unsigned k) { return with_order<int>(m, [&](int j) { return bm::spherical_harmonic_i(k, j, theta, phi); }); }); });

   distribution<bm::normal_distribution<double>, double, double>("normal_distribution");
   distribution<bm::students_t_distribution<double>, double>("students_t_distribution");
   distribution<bm::chi_squared_distribution<double>, double>("chi_squared_distribution");
   distribution<bm::beta_distribution<double>, double, double>("beta_distribution");
   distribution<bm::gamma_distribution<double>, double, double>("gamma_distribution");
   distribution<bm::exponential_distribution<double>, double>("exponential_distribution");
   distribution<bm::lognormal_distribution<double>, double, double>("lognormal_distribution");
   distribution<bm::weibull_distribution<double>, double, double>("weibull_distribution");
   distribution<bm::cauchy_distribution<double>, double, double>("cauchy_distribution");
   distribution<bm::fisher_f_distribution<double>, double, double>("fisher_f_distribution");
   distribution<bm::uniform_distribution<double>, double, double>("uniform_distribution");
   distribution<bm::laplace_distribution<double>, double, double>("laplace_distribution");
   distribution<bm::logistic_distribution<double>, double, double>("logistic_distribution");
   distribution<bm::non_central_chi_squared_distribution<double>, double, double>("non_central_chi_squared_distribution");
   distribution<bm::non_central_t_distribution<double>, double, double>("non_central_t_distribution");
   distribution<bm::binomial_distribution<double>, double, double>("binomial_distribution");
   distribution<bm::poisson_distribution<double>, double>("poisson_distribution");
   distribution<bm::negative_binomial_distribution<double>, double, double>("negative_binomial_distribution");
   distribution<bm::geometric_distribution<double>, double>("geometric_distribution");
   distribution<bm::arcsine_distribution<double>, double, double>("arcsine_distribution");
   distribution<bm::bernoulli_distribution<double>, double>("bernoulli_distribution");
   distribution<bm::extreme_value_distribution<double>, double, double>("extreme_value_distribution");
   distribution<bm::generalized_extreme_value_distribution<double>, double, double, double>("generalized_extreme_value_distribution");
   distribution<bm::holtsmark_distribution<double>, double, double>("holtsmark_distribution");
   distribution<bm::hypergeometric_distribution<double>, unsigned, unsigned, unsigned>("hypergeometric_distribution");
   distribution<bm::inverse_chi_squared_distribution<double>, double, double>("inverse_chi_squared_distribution");
   distribution<bm::inverse_gamma_distribution<double>, double, double>("inverse_gamma_distribution");
   distribution<bm::inverse_gaussian_distribution<double>, double, double>("inverse_gaussian_distribution");
   distribution<bm::kolmogorov_smirnov_distribution<double>, double>("kolmogorov_smirnov_distribution");
   distribution<bm::landau_distribution<double>, double, double>("landau_distribution");
   distribution<bm::mapairy_distribution<double>, double, double>("mapairy_distribution");
   distribution<bm::non_central_beta_distribution<double>, double, double, double>("non_central_beta_distribution");
   distribution<bm::non_central_f_distribution<double>, double, double, double>("non_central_f_distribution");
   distribution<bm::pareto_distribution<double>, double, double>("pareto_distribution");
   distribution<bm::rayleigh_distribution<double>, double>("rayleigh_distribution");
   distribution<bm::saspoint5_distribution<double>, double, double>("saspoint5_distribution");
   distribution<bm::skew_normal_distribution<double>, double, double, double>("skew_normal_distribution");
   distribution<bm::triangular_distribution<double>, double, double, double>("triangular_distribution");
   distribution<bm::von_mises_distribution<double>, double, double>("von_mises_distribution");

   class_<fft>("fft")
      .constructor<double>()
      .function("transform", +[](fft& f, const val& x) { return f.run(x, false); })
      .function("inverse_transform", +[](fft& f, const val& x) { return f.run(x, true); });
   class_<real_fft>("real_fft").constructor<double>().function("transform", &real_fft::transform);
   for (const char* name : {"new fft", "fft.transform", "fft.inverse_transform", "new real_fft", "real_fft.transform"})
      throws(name);
   entry("fft_axis", +[](const val& x, const val& shape, double axis) { return fft_axis(x, shape, axis, false); });
   entry("inverse_fft_axis", +[](const val& x, const val& shape, double axis) { return fft_axis(x, shape, axis, true); });
   entry("real_fft_axis", &real_fft_axis);
   entry("power_spectrum_axis", &power_spectrum_axis);
   entry("power_spectrum", +[](const val& x) { return power_spectrum_axis(x, val::array(std::vector<double>{x["length"].as<double>()}), 0); });
}
