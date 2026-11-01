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
      nb::class_<D> c(m, "non_central_beta_distribution",
                     "Noncentral beta distribution with shapes alpha, beta and non-centrality lambda_." BMPY_DIST_URL("nc_beta_dist"));
      c.def(bmpy::init<double, double, double>(), arg("alpha"), arg("beta"), arg("lambda_"));
      add_dist_methods(c);
      bmpy::add_dist_params(c, {BMPY_PARAM("alpha", alpha), BMPY_PARAM("beta", beta), BMPY_PARAM("lambda_", non_centrality)});
   }
   {
      using D = bm::non_central_chi_squared_distribution<double, pol>;
      nb::class_<D> c(m, "non_central_chi_squared_distribution",
                     "Noncentral chi-squared distribution with df degrees of freedom and non-centrality lambda_: the sum of df squared normals of means mu_i, lambda_ = sum mu_i^2." BMPY_DIST_URL("nc_chi_squared_dist"));
      c.def(bmpy::init<double, double>(), arg("df"), arg("lambda_"));
      add_dist_methods(c);
      bmpy::add_dist_params(c, {BMPY_PARAM("df", degrees_of_freedom), BMPY_PARAM("lambda_", non_centrality)});
      bmpy::def_static_scalar(c, "find_degrees_of_freedom", "find_degrees_of_freedom(lambda_, x, p, /)",
         +[](double lambda, double x, double p) { return D::find_degrees_of_freedom(lambda, x, p); },
         "The degrees of freedom df with cdf(x; df, lambda_) = p.");
      bmpy::def_static_scalar(c, "find_non_centrality", "find_non_centrality(df, x, p, /)",
         +[](double df, double x, double p) { return D::find_non_centrality(df, x, p); },
         "The non-centrality lambda_ with cdf(x; df, lambda_) = p.");
   }
   {
      using D = bm::non_central_f_distribution<double, pol>;
      nb::class_<D> c(m, "non_central_f_distribution",
                     "Noncentral F distribution with df1 and df2 degrees of freedom and non-centrality lambda_ of the numerator." BMPY_DIST_URL("nc_f_dist"));
      c.def(bmpy::init<double, double, double>(), arg("df1"), arg("df2"), arg("lambda_"));
      add_dist_methods(c);
      bmpy::add_dist_params(c, {BMPY_PARAM("df1", degrees_of_freedom1), BMPY_PARAM("df2", degrees_of_freedom2), BMPY_PARAM("lambda_", non_centrality)});
      bmpy::def_static_scalar(c, "find_non_centrality", "find_non_centrality(x, df1, df2, p, /)",
         +[](double x, double df1, double df2, double p) { return D::find_non_centrality(x, df1, df2, p); },
         "The non-centrality lambda_ with cdf(x; df1, df2, lambda_) = p.");
   }
   {
      using D = bm::non_central_t_distribution<double, pol>;
      nb::class_<D> c(m, "non_central_t_distribution",
                     "Noncentral t distribution with df degrees of freedom and non-centrality delta: (Z + delta) / sqrt(V / df)." BMPY_DIST_URL("nc_t_dist"));
      c.def(bmpy::init<double, double>(), arg("df"), arg("delta"));
      add_dist_methods(c);
      bmpy::add_dist_params(c, {BMPY_PARAM("df", degrees_of_freedom), BMPY_PARAM("delta", non_centrality)});
      bmpy::def_static_scalar(c, "find_degrees_of_freedom", "find_degrees_of_freedom(delta, x, p, /)",
         +[](double delta, double x, double p) { return D::find_degrees_of_freedom(delta, x, p); },
         "The degrees of freedom df with cdf(x; df, delta) = p.");
      bmpy::def_static_scalar(c, "find_non_centrality", "find_non_centrality(df, x, p, /)",
         +[](double df, double x, double p) { return D::find_non_centrality(df, x, p); },
         "The non-centrality delta with cdf(x; df, delta) = p.");
   }
}
