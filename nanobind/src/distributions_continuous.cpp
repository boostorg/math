// Copyright 2026 Nick Thompson
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/math/distributions/arcsine.hpp>
#include <boost/math/distributions/beta.hpp>
#include <boost/math/distributions/cauchy.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/exponential.hpp>
#include <boost/math/distributions/extreme_value.hpp>
#include <boost/math/distributions/fisher_f.hpp>
#include <boost/math/distributions/gamma.hpp>
#include <boost/math/distributions/generalized_extreme_value.hpp>
#include <boost/math/distributions/inverse_chi_squared.hpp>
#include <boost/math/distributions/inverse_gamma.hpp>
#include <boost/math/distributions/inverse_gaussian.hpp>
#include <boost/math/distributions/laplace.hpp>
#include <boost/math/distributions/logistic.hpp>
#include <boost/math/distributions/lognormal.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/pareto.hpp>
#include <boost/math/distributions/rayleigh.hpp>
#include <boost/math/distributions/skew_normal.hpp>
#include <boost/math/distributions/students_t.hpp>
#include <boost/math/distributions/triangular.hpp>
#include <boost/math/distributions/uniform.hpp>
#include <boost/math/distributions/weibull.hpp>

#include "bindings.hpp"
#include "distributions.hpp"

namespace bm = boost::math;
namespace nb = nanobind;
using nb::arg;
using bmpy::add_dist_methods;
using pol = bmpy::dist_policy;

void bmpy::bind_distributions_continuous(nb::module_& m) {
   {
      using D = bm::arcsine_distribution<double, pol>;
      nb::class_<D> c(m, "arcsine_distribution",
                     "Arcsine distribution on [x_min, x_max]: pdf 1 / (pi sqrt((x - x_min)(x_max - x)))." BMPY_DIST_URL("arcine_dist"));
      c.def(bmpy::init<double, double>(), arg("x_min") = 0.0, arg("x_max") = 1.0);
      add_dist_methods(c);
      bmpy::add_dist_params(c, {BMPY_PARAM("x_min", x_min), BMPY_PARAM("x_max", x_max)});
   }
   {
      using D = bm::beta_distribution<double, pol>;
      nb::class_<D> c(m, "beta_distribution",
                     "Beta distribution on [0, 1]: pdf x^(alpha - 1) (1 - x)^(beta - 1) / B(alpha, beta)." BMPY_DIST_URL("beta_dist"));
      c.def(bmpy::init<double, double>(), arg("alpha") = 1.0, arg("beta") = 1.0);
      add_dist_methods(c);
      bmpy::add_dist_params(c, {BMPY_PARAM("alpha", alpha), BMPY_PARAM("beta", beta)});
   }
   {
      using D = bm::cauchy_distribution<double, pol>;
      nb::class_<D> c(m, "cauchy_distribution",
                     "Cauchy-Lorentz distribution: pdf 1 / (pi scale (1 + ((x - location) / scale)^2))." BMPY_DIST_URL("cauchy_dist"));
      c.def(bmpy::init<double, double>(), arg("location") = 0.0, arg("scale") = 1.0);
      add_dist_methods(c);
      bmpy::add_dist_params(c, {BMPY_PARAM("location", location), BMPY_PARAM("scale", scale)});
   }
   {
      using D = bm::chi_squared_distribution<double, pol>;
      nb::class_<D> c(m, "chi_squared_distribution",
                     "Chi-squared distribution with df degrees of freedom: the sum of df squared standard normals." BMPY_DIST_URL("chi_squared_dist"));
      c.def(bmpy::init<double>(), arg("df") = 1.0);
      add_dist_methods(c);
      bmpy::add_dist_params(c, {BMPY_PARAM("df", degrees_of_freedom)});
      bmpy::def_static_scalar(c, "find_degrees_of_freedom", "find_degrees_of_freedom(difference_from_variance, alpha, beta, variance, /)",
         +[](double difference, double alpha, double beta, double variance) {
            return D::find_degrees_of_freedom(difference, alpha, beta, variance);
         },
         "The degrees of freedom a chi-squared test for the variance needs to detect a difference_from_variance from the "
         "nominal variance, with risks alpha (type I) and beta (type II).");
   }
   {
      using D = bm::exponential_distribution<double, pol>;
      nb::class_<D> c(m, "exponential_distribution",
                     "Exponential distribution with rate lambda_: pdf lambda_ exp(-lambda_ x), x >= 0." BMPY_DIST_URL("exp_dist"));
      c.def(bmpy::init<double>(), arg("lambda_") = 1.0);
      add_dist_methods(c);
      bmpy::add_dist_params(c, {BMPY_PARAM("lambda_", lambda)});
   }
   {
      using D = bm::extreme_value_distribution<double, pol>;
      nb::class_<D> c(m, "extreme_value_distribution",
                     "Extreme value (Gumbel) distribution: cdf exp(-exp(-(x - location) / scale))." BMPY_DIST_URL("extreme_dist"));
      c.def(bmpy::init<double, double>(), arg("location") = 0.0, arg("scale") = 1.0);
      add_dist_methods(c);
      bmpy::add_dist_params(c, {BMPY_PARAM("location", location), BMPY_PARAM("scale", scale)});
   }
   {
      using D = bm::generalized_extreme_value_distribution<double, pol>;
      nb::class_<D> c(m, "generalized_extreme_value_distribution",
                     "Generalized extreme value distribution: cdf exp(-(1 + shape (x - location) / scale)^(-1 / shape))." BMPY_DIST_URL("gev_dist"));
      c.def(bmpy::init<double, double, double>(), arg("location") = 0.0, arg("scale") = 1.0, arg("shape") = 0.0);
      add_dist_methods(c);
      bmpy::add_dist_params(c, {BMPY_PARAM("location", location), BMPY_PARAM("scale", scale), BMPY_PARAM("shape", shape)});
   }
   {
      using D = bm::fisher_f_distribution<double, pol>;
      nb::class_<D> c(m, "fisher_f_distribution",
                     "Fisher F distribution with df1 and df2 degrees of freedom: the ratio of chi-squared variates, each over its df." BMPY_DIST_URL("f_dist"));
      c.def(bmpy::init<double, double>(), arg("df1"), arg("df2"));
      add_dist_methods(c);
      bmpy::add_dist_params(c, {BMPY_PARAM("df1", degrees_of_freedom1), BMPY_PARAM("df2", degrees_of_freedom2)});
   }
   {
      using D = bm::gamma_distribution<double, pol>;
      nb::class_<D> c(m, "gamma_distribution",
                     "Gamma distribution: pdf x^(shape - 1) exp(-x / scale) / (Gamma(shape) scale^shape), x > 0." BMPY_DIST_URL("gamma_dist"));
      c.def(bmpy::init<double, double>(), arg("shape"), arg("scale") = 1.0);
      add_dist_methods(c);
      bmpy::add_dist_params(c, {BMPY_PARAM("shape", shape), BMPY_PARAM("scale", scale)});
   }
   {
      using D = bm::inverse_chi_squared_distribution<double, pol>;
      nb::class_<D> c(m, "inverse_chi_squared_distribution",
                     "Inverse chi-squared distribution with df degrees of freedom and scale; scale defaults to 1 / df." BMPY_DIST_URL("inverse_chi_squared_dist"));
      c.def(bmpy::init<double>(), arg("df") = 1.0);
      c.def(bmpy::init<double, double>(), arg("df"), arg("scale"));
      add_dist_methods(c);
      bmpy::add_dist_params(c, {BMPY_PARAM("df", degrees_of_freedom), BMPY_PARAM("scale", scale)});
   }
   {
      using D = bm::inverse_gamma_distribution<double, pol>;
      nb::class_<D> c(m, "inverse_gamma_distribution",
                     "Inverse gamma distribution: pdf scale^shape x^(-shape - 1) exp(-scale / x) / Gamma(shape), x > 0." BMPY_DIST_URL("inverse_gamma_dist"));
      c.def(bmpy::init<double, double>(), arg("shape") = 1.0, arg("scale") = 1.0);
      add_dist_methods(c);
      bmpy::add_dist_params(c, {BMPY_PARAM("shape", shape), BMPY_PARAM("scale", scale)});
   }
   {
      using D = bm::inverse_gaussian_distribution<double, pol>;
      nb::class_<D> c(m, "inverse_gaussian_distribution",
                     "Inverse Gaussian (Wald) distribution: pdf sqrt(scale / (2 pi x^3)) exp(-scale (x - mean)^2 / (2 mean^2 x))." BMPY_DIST_URL("inverse_gaussian_dist"));
      c.def(bmpy::init<double, double>(), arg("mean") = 1.0, arg("scale") = 1.0);
      add_dist_methods(c);
      bmpy::add_dist_params(c, {BMPY_PARAM("mean", mean), BMPY_PARAM("scale", scale)});
   }
   {
      using D = bm::laplace_distribution<double, pol>;
      nb::class_<D> c(m, "laplace_distribution",
                     "Laplace (double exponential) distribution: pdf exp(-|x - location| / scale) / (2 scale)." BMPY_DIST_URL("laplace_dist"));
      c.def(bmpy::init<double, double>(), arg("location") = 0.0, arg("scale") = 1.0);
      add_dist_methods(c);
      bmpy::add_dist_params(c, {BMPY_PARAM("location", location), BMPY_PARAM("scale", scale)});
   }
   {
      using D = bm::logistic_distribution<double, pol>;
      nb::class_<D> c(m, "logistic_distribution",
                     "Logistic distribution: cdf 1 / (1 + exp(-(x - location) / scale))." BMPY_DIST_URL("logistic_dist"));
      c.def(bmpy::init<double, double>(), arg("location") = 0.0, arg("scale") = 1.0);
      add_dist_methods(c);
      bmpy::add_dist_params(c, {BMPY_PARAM("location", location), BMPY_PARAM("scale", scale)});
   }
   {
      using D = bm::lognormal_distribution<double, pol>;
      nb::class_<D> c(m, "lognormal_distribution",
                     "Log-normal distribution: log(X) is normal with mean location and standard deviation scale." BMPY_DIST_URL("lognormal_dist"));
      c.def(bmpy::init<double, double>(), arg("location") = 0.0, arg("scale") = 1.0);
      add_dist_methods(c);
      bmpy::add_dist_params(c, {BMPY_PARAM("location", location), BMPY_PARAM("scale", scale)});
   }
   {
      using D = bm::normal_distribution<double, pol>;
      nb::class_<D> c(m, "normal_distribution",
                     "Normal (Gaussian) distribution: pdf exp(-(x - mean)^2 / (2 sd^2)) / (sd sqrt(2 pi))." BMPY_DIST_URL("normal_dist"));
      c.def(bmpy::init<double, double>(), arg("mean") = 0.0, arg("sd") = 1.0);
      add_dist_methods(c);
      bmpy::add_dist_params(c, {BMPY_PARAM("mean", mean), BMPY_PARAM("sd", standard_deviation)});
   }
   {
      using D = bm::pareto_distribution<double, pol>;
      nb::class_<D> c(m, "pareto_distribution",
                     "Pareto distribution: pdf shape scale^shape / x^(shape + 1), x >= scale." BMPY_DIST_URL("pareto"));
      c.def(bmpy::init<double, double>(), arg("scale") = 1.0, arg("shape") = 1.0);
      add_dist_methods(c);
      bmpy::add_dist_params(c, {BMPY_PARAM("scale", scale), BMPY_PARAM("shape", shape)});
   }
   {
      using D = bm::rayleigh_distribution<double, pol>;
      nb::class_<D> c(m, "rayleigh_distribution",
                     "Rayleigh distribution: pdf x exp(-x^2 / (2 sigma^2)) / sigma^2, x >= 0." BMPY_DIST_URL("rayleigh"));
      c.def(bmpy::init<double>(), arg("sigma") = 1.0);
      add_dist_methods(c);
      bmpy::add_dist_params(c, {BMPY_PARAM("sigma", sigma)});
   }
   {
      using D = bm::skew_normal_distribution<double, pol>;
      nb::class_<D> c(m, "skew_normal_distribution",
                     "Skew normal distribution: pdf 2 phi(z) Phi(shape z) / scale, z = (x - location) / scale." BMPY_DIST_URL("skew_normal_dist"));
      c.def(bmpy::init<double, double, double>(), arg("location") = 0.0, arg("scale") = 1.0, arg("shape") = 0.0);
      add_dist_methods(c);
      bmpy::add_dist_params(c, {BMPY_PARAM("location", location), BMPY_PARAM("scale", scale), BMPY_PARAM("shape", shape)});
   }
   {
      using D = bm::students_t_distribution<double, pol>;
      nb::class_<D> c(m, "students_t_distribution",
                     "Student's t distribution with df degrees of freedom: pdf proportional to (1 + x^2 / df)^(-(df + 1) / 2)." BMPY_DIST_URL("students_t_dist"));
      c.def(bmpy::init<double>(), arg("df") = 1.0);
      add_dist_methods(c);
      bmpy::add_dist_params(c, {BMPY_PARAM("df", degrees_of_freedom)});
      bmpy::def_static_scalar(c, "find_degrees_of_freedom", "find_degrees_of_freedom(difference_from_mean, alpha, beta, sd, /)",
         +[](double difference, double alpha, double beta, double sd) { return D::find_degrees_of_freedom(difference, alpha, beta, sd); },
         "The degrees of freedom a t test needs to detect a mean differing by difference_from_mean, with risks alpha (type I) "
         "and beta (type II), for standard deviation sd.");
      bmpy::def_static_scalar(c, "find_degrees_of_freedom", "find_degrees_of_freedom(t, p, /)",
         +[](double t, double p) { return D::find_degrees_of_freedom(t, p); }, "The degrees of freedom df with cdf(t; df) = p.");
   }
   {
      using D = bm::triangular_distribution<double, pol>;
      nb::class_<D> c(m, "triangular_distribution",
                     "Triangular distribution on [lower, upper] with its peak at mode." BMPY_DIST_URL("triangular_dist"));
      c.def(bmpy::init<double, double, double>(), arg("lower") = -1.0, arg("mode") = 0.0, arg("upper") = 1.0);
      add_dist_methods(c);
      bmpy::add_dist_params(c, {BMPY_PARAM("lower", lower), BMPY_PARAM("mode", mode), BMPY_PARAM("upper", upper)});
   }
   {
      using D = bm::uniform_distribution<double, pol>;
      nb::class_<D> c(m, "uniform_distribution",
                     "Uniform distribution on [lower, upper]: pdf 1 / (upper - lower)." BMPY_DIST_URL("uniform_dist"));
      c.def(bmpy::init<double, double>(), arg("lower") = 0.0, arg("upper") = 1.0);
      add_dist_methods(c);
      bmpy::add_dist_params(c, {BMPY_PARAM("lower", lower), BMPY_PARAM("upper", upper)});
   }
   {
      using D = bm::weibull_distribution<double, pol>;
      nb::class_<D> c(m, "weibull_distribution",
                     "Weibull distribution: cdf 1 - exp(-(x / scale)^shape), x >= 0." BMPY_DIST_URL("weibull_dist"));
      c.def(bmpy::init<double, double>(), arg("shape"), arg("scale") = 1.0);
      add_dist_methods(c);
      bmpy::add_dist_params(c, {BMPY_PARAM("shape", shape), BMPY_PARAM("scale", scale)});
   }
}
