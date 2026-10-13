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
      nb::class_<D> c(m, "arcsine_distribution");
      c.def(nb::init<double, double>(), arg("x_min") = 0.0, arg("x_max") = 1.0);
      add_dist_methods(c);
      bmpy::add_dist_params(c, {BMPY_PARAM("x_min", x_min), BMPY_PARAM("x_max", x_max)});
   }
   {
      using D = bm::beta_distribution<double, pol>;
      nb::class_<D> c(m, "beta_distribution");
      c.def(nb::init<double, double>(), arg("alpha") = 1.0, arg("beta") = 1.0);
      add_dist_methods(c);
      bmpy::add_dist_params(c, {BMPY_PARAM("alpha", alpha), BMPY_PARAM("beta", beta)});
   }
   {
      using D = bm::cauchy_distribution<double, pol>;
      nb::class_<D> c(m, "cauchy_distribution");
      c.def(nb::init<double, double>(), arg("location") = 0.0, arg("scale") = 1.0);
      add_dist_methods(c);
      bmpy::add_dist_params(c, {BMPY_PARAM("location", location), BMPY_PARAM("scale", scale)});
   }
   {
      using D = bm::chi_squared_distribution<double, pol>;
      nb::class_<D> c(m, "chi_squared_distribution");
      c.def(nb::init<double>(), arg("df") = 1.0);
      add_dist_methods(c);
      bmpy::add_dist_params(c, {BMPY_PARAM("df", degrees_of_freedom)});
   }
   {
      using D = bm::exponential_distribution<double, pol>;
      nb::class_<D> c(m, "exponential_distribution");
      c.def(nb::init<double>(), arg("lambda_") = 1.0);
      add_dist_methods(c);
      bmpy::add_dist_params(c, {BMPY_PARAM("lambda_", lambda)});
   }
   {
      using D = bm::extreme_value_distribution<double, pol>;
      nb::class_<D> c(m, "extreme_value_distribution");
      c.def(nb::init<double, double>(), arg("location") = 0.0, arg("scale") = 1.0);
      add_dist_methods(c);
      bmpy::add_dist_params(c, {BMPY_PARAM("location", location), BMPY_PARAM("scale", scale)});
   }
   {
      using D = bm::generalized_extreme_value_distribution<double, pol>;
      nb::class_<D> c(m, "generalized_extreme_value_distribution");
      c.def(nb::init<double, double, double>(), arg("location") = 0.0, arg("scale") = 1.0, arg("shape") = 0.0);
      add_dist_methods(c);
      bmpy::add_dist_params(c, {BMPY_PARAM("location", location), BMPY_PARAM("scale", scale), BMPY_PARAM("shape", shape)});
   }
   {
      using D = bm::fisher_f_distribution<double, pol>;
      nb::class_<D> c(m, "fisher_f_distribution");
      c.def(nb::init<double, double>(), arg("df1"), arg("df2"));
      add_dist_methods(c);
      bmpy::add_dist_params(c, {BMPY_PARAM("df1", degrees_of_freedom1), BMPY_PARAM("df2", degrees_of_freedom2)});
   }
   {
      using D = bm::gamma_distribution<double, pol>;
      nb::class_<D> c(m, "gamma_distribution");
      c.def(nb::init<double, double>(), arg("shape"), arg("scale") = 1.0);
      add_dist_methods(c);
      bmpy::add_dist_params(c, {BMPY_PARAM("shape", shape), BMPY_PARAM("scale", scale)});
   }
   {
      using D = bm::inverse_chi_squared_distribution<double, pol>;
      nb::class_<D> c(m, "inverse_chi_squared_distribution");
      c.def(nb::init<double>(), arg("df") = 1.0);
      c.def(nb::init<double, double>(), arg("df"), arg("scale"));
      add_dist_methods(c);
      bmpy::add_dist_params(c, {BMPY_PARAM("df", degrees_of_freedom), BMPY_PARAM("scale", scale)});
   }
   {
      using D = bm::inverse_gamma_distribution<double, pol>;
      nb::class_<D> c(m, "inverse_gamma_distribution");
      c.def(nb::init<double, double>(), arg("shape") = 1.0, arg("scale") = 1.0);
      add_dist_methods(c);
      bmpy::add_dist_params(c, {BMPY_PARAM("shape", shape), BMPY_PARAM("scale", scale)});
   }
   {
      using D = bm::inverse_gaussian_distribution<double, pol>;
      nb::class_<D> c(m, "inverse_gaussian_distribution");
      c.def(nb::init<double, double>(), arg("mean") = 1.0, arg("scale") = 1.0);
      add_dist_methods(c);
      bmpy::add_dist_params(c, {BMPY_PARAM("mean", mean), BMPY_PARAM("scale", scale)});
   }
   {
      using D = bm::laplace_distribution<double, pol>;
      nb::class_<D> c(m, "laplace_distribution");
      c.def(nb::init<double, double>(), arg("location") = 0.0, arg("scale") = 1.0);
      add_dist_methods(c);
      bmpy::add_dist_params(c, {BMPY_PARAM("location", location), BMPY_PARAM("scale", scale)});
   }
   {
      using D = bm::logistic_distribution<double, pol>;
      nb::class_<D> c(m, "logistic_distribution");
      c.def(nb::init<double, double>(), arg("location") = 0.0, arg("scale") = 1.0);
      add_dist_methods(c);
      bmpy::add_dist_params(c, {BMPY_PARAM("location", location), BMPY_PARAM("scale", scale)});
   }
   {
      using D = bm::lognormal_distribution<double, pol>;
      nb::class_<D> c(m, "lognormal_distribution");
      c.def(nb::init<double, double>(), arg("location") = 0.0, arg("scale") = 1.0);
      add_dist_methods(c);
      bmpy::add_dist_params(c, {BMPY_PARAM("location", location), BMPY_PARAM("scale", scale)});
   }
   {
      using D = bm::normal_distribution<double, pol>;
      nb::class_<D> c(m, "normal_distribution");
      c.def(nb::init<double, double>(), arg("mean") = 0.0, arg("sd") = 1.0);
      add_dist_methods(c);
      bmpy::add_dist_params(c, {BMPY_PARAM("mean", mean), BMPY_PARAM("sd", standard_deviation)});
   }
   {
      using D = bm::pareto_distribution<double, pol>;
      nb::class_<D> c(m, "pareto_distribution");
      c.def(nb::init<double, double>(), arg("scale") = 1.0, arg("shape") = 1.0);
      add_dist_methods(c);
      bmpy::add_dist_params(c, {BMPY_PARAM("scale", scale), BMPY_PARAM("shape", shape)});
   }
   {
      using D = bm::rayleigh_distribution<double, pol>;
      nb::class_<D> c(m, "rayleigh_distribution");
      c.def(nb::init<double>(), arg("sigma") = 1.0);
      add_dist_methods(c);
      bmpy::add_dist_params(c, {BMPY_PARAM("sigma", sigma)});
   }
   {
      using D = bm::skew_normal_distribution<double, pol>;
      nb::class_<D> c(m, "skew_normal_distribution");
      c.def(nb::init<double, double, double>(), arg("location") = 0.0, arg("scale") = 1.0, arg("shape") = 0.0);
      add_dist_methods(c);
      bmpy::add_dist_params(c, {BMPY_PARAM("location", location), BMPY_PARAM("scale", scale), BMPY_PARAM("shape", shape)});
   }
   {
      using D = bm::students_t_distribution<double, pol>;
      nb::class_<D> c(m, "students_t_distribution");
      c.def(nb::init<double>(), arg("df") = 1.0);
      add_dist_methods(c);
      bmpy::add_dist_params(c, {BMPY_PARAM("df", degrees_of_freedom)});
   }
   {
      using D = bm::triangular_distribution<double, pol>;
      nb::class_<D> c(m, "triangular_distribution");
      c.def(nb::init<double, double, double>(), arg("lower") = -1.0, arg("mode") = 0.0, arg("upper") = 1.0);
      add_dist_methods(c);
      bmpy::add_dist_params(c, {BMPY_PARAM("lower", lower), BMPY_PARAM("mode", mode), BMPY_PARAM("upper", upper)});
   }
   {
      using D = bm::uniform_distribution<double, pol>;
      nb::class_<D> c(m, "uniform_distribution");
      c.def(nb::init<double, double>(), arg("lower") = 0.0, arg("upper") = 1.0);
      add_dist_methods(c);
      bmpy::add_dist_params(c, {BMPY_PARAM("lower", lower), BMPY_PARAM("upper", upper)});
   }
   {
      using D = bm::weibull_distribution<double, pol>;
      nb::class_<D> c(m, "weibull_distribution");
      c.def(nb::init<double, double>(), arg("shape"), arg("scale") = 1.0);
      add_dist_methods(c);
      bmpy::add_dist_params(c, {BMPY_PARAM("shape", shape), BMPY_PARAM("scale", scale)});
   }
}
