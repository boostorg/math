// Copyright 2026 Nick Thompson
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// Stable-law and other distributions with large tables or unusual constructors.

#include <boost/math/distributions/holtsmark.hpp>
#include <boost/math/distributions/landau.hpp>
#include <boost/math/distributions/mapairy.hpp>
#include <boost/math/distributions/saspoint5.hpp>
#include <boost/math/distributions/von_mises.hpp>
#include <boost/math/distributions/kolmogorov_smirnov.hpp>
#include <boost/math/distributions/hyperexponential.hpp>

#include <nanobind/stl/vector.h>

#include <vector>

#include "bindings.hpp"
#include "distributions.hpp"

namespace bm = boost::math;
namespace nb = nanobind;
using nb::arg;
using bmpy::add_dist_methods;
using pol = bmpy::dist_policy;

void bmpy::bind_distributions_heavy(nb::module_& m) {
   {
      using D = bm::holtsmark_distribution<double, pol>;
      nb::class_<D> c(m, "holtsmark_distribution");
      c.def(nb::init<double, double>(), arg("location") = 0.0, arg("scale") = 1.0);
      add_dist_methods(c);
   }
   {
      using D = bm::landau_distribution<double, pol>;
      nb::class_<D> c(m, "landau_distribution");
      c.def(nb::init<double, double>(), arg("location") = 0.0, arg("scale") = 1.0);
      add_dist_methods(c);
   }
   {
      using D = bm::mapairy_distribution<double, pol>;
      nb::class_<D> c(m, "mapairy_distribution");
      c.def(nb::init<double, double>(), arg("location") = 0.0, arg("scale") = 1.0);
      add_dist_methods(c);
   }
   {
      using D = bm::saspoint5_distribution<double, pol>;
      nb::class_<D> c(m, "saspoint5_distribution");
      c.def(nb::init<double, double>(), arg("location") = 0.0, arg("scale") = 1.0);
      add_dist_methods(c);
   }
   {
      using D = bm::von_mises_distribution<double, pol>;
      nb::class_<D> c(m, "von_mises_distribution");
      c.def(nb::init<double, double>(), arg("mean") = 0.0, arg("concentration") = 1.0);
      add_dist_methods(c);
   }
   {
      using D = bm::kolmogorov_smirnov_distribution<double, pol>;
      nb::class_<D> c(m, "kolmogorov_smirnov_distribution");
      c.def(nb::init<double>(), arg("n"));
      add_dist_methods(c);
   }
   {
      using D = bm::hyperexponential_distribution<double, pol>;
      nb::class_<D> c(m, "hyperexponential_distribution");
      c.def(nb::init<const std::vector<double>&>(), arg("rates"));
      c.def(
          "__init__",
          [](D* self, const std::vector<double>& probs, const std::vector<double>& rates) { new (self) D(probs, rates); },
          arg("probabilities"), arg("rates"));
      add_dist_methods(c);
   }
}
