// Copyright 2026 Nick Thompson
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// Discrete distributions use the default discrete_quantile policy.

#include <boost/math/distributions/bernoulli.hpp>
#include <boost/math/distributions/binomial.hpp>
#include <boost/math/distributions/geometric.hpp>
#include <boost/math/distributions/hypergeometric.hpp>
#include <boost/math/distributions/negative_binomial.hpp>
#include <boost/math/distributions/poisson.hpp>

#include <cstdint>

#include "bindings.hpp"
#include "distributions.hpp"

namespace bm = boost::math;
namespace nb = nanobind;
using nb::arg;
using bmpy::add_dist_methods;
using pol = bmpy::dist_policy;

void bmpy::bind_distributions_discrete(nb::module_& m) {
   {
      using D = bm::bernoulli_distribution<double, pol>;
      nb::class_<D> c(m, "bernoulli_distribution");
      c.def(nb::init<double>(), arg("p") = 0.5);
      add_dist_methods(c);
   }
   {
      using D = bm::binomial_distribution<double, pol>;
      nb::class_<D> c(m, "binomial_distribution");
      c.def(nb::init<double, double>(), arg("n") = 1.0, arg("p") = 0.5);
      add_dist_methods(c);
   }
   {
      using D = bm::geometric_distribution<double, pol>;
      nb::class_<D> c(m, "geometric_distribution");
      c.def(nb::init<double>(), arg("p"));
      add_dist_methods(c);
   }
   {
      using D = bm::hypergeometric_distribution<double, pol>;
      nb::class_<D> c(m, "hypergeometric_distribution");
      c.def(nb::init<std::uint64_t, std::uint64_t, std::uint64_t>(), arg("r"), arg("n"), arg("N"));
      add_dist_methods(c);
   }
   {
      using D = bm::negative_binomial_distribution<double, pol>;
      nb::class_<D> c(m, "negative_binomial_distribution");
      c.def(nb::init<double, double>(), arg("r"), arg("p"));
      add_dist_methods(c);
   }
   {
      using D = bm::poisson_distribution<double, pol>;
      nb::class_<D> c(m, "poisson_distribution");
      c.def(nb::init<double>(), arg("mean") = 1.0);
      add_dist_methods(c);
   }
}
