// Copyright 2026 Nick Thompson
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/math/distributions/non_central_beta.hpp>
#include <boost/math/distributions/non_central_chi_squared.hpp>
#include <boost/math/distributions/non_central_f.hpp>
#include <boost/math/distributions/non_central_t.hpp>

#include "bindings.hpp"
#include "distributions.hpp"

namespace bm = boost::math;
namespace nb = nanobind;
using nb::arg;
using bmpy::add_dist_methods;
using pol = bmpy::dist_policy;

void bmpy::bind_distributions_noncentral(nb::module_& m) {
   {
      using D = bm::non_central_beta_distribution<double, pol>;
      nb::class_<D> c(m, "non_central_beta_distribution");
      c.def(nb::init<double, double, double>(), arg("alpha"), arg("beta"), arg("lambda_"));
      add_dist_methods(c);
      bmpy::add_dist_params(c, {BMPY_PARAM("alpha", alpha), BMPY_PARAM("beta", beta), BMPY_PARAM("lambda_", non_centrality)});
   }
   {
      using D = bm::non_central_chi_squared_distribution<double, pol>;
      nb::class_<D> c(m, "non_central_chi_squared_distribution");
      c.def(nb::init<double, double>(), arg("df"), arg("lambda_"));
      add_dist_methods(c);
      bmpy::add_dist_params(c, {BMPY_PARAM("df", degrees_of_freedom), BMPY_PARAM("lambda_", non_centrality)});
   }
   {
      using D = bm::non_central_f_distribution<double, pol>;
      nb::class_<D> c(m, "non_central_f_distribution");
      c.def(nb::init<double, double, double>(), arg("df1"), arg("df2"), arg("lambda_"));
      add_dist_methods(c);
      bmpy::add_dist_params(c, {BMPY_PARAM("df1", degrees_of_freedom1), BMPY_PARAM("df2", degrees_of_freedom2), BMPY_PARAM("lambda_", non_centrality)});
   }
   {
      using D = bm::non_central_t_distribution<double, pol>;
      nb::class_<D> c(m, "non_central_t_distribution");
      c.def(nb::init<double, double>(), arg("df"), arg("delta"));
      add_dist_methods(c);
      bmpy::add_dist_params(c, {BMPY_PARAM("df", degrees_of_freedom), BMPY_PARAM("delta", non_centrality)});
   }
}
