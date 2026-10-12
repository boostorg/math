// Copyright 2026 Nick Thompson
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// Binds a Boost.Math distribution (double, default policy) as a Python class
// with vectorised pdf/cdf/quantile-style methods and scalar moment methods.

#ifndef BOOST_MATH_NANOBIND_DISTRIBUTIONS_HPP
#define BOOST_MATH_NANOBIND_DISTRIBUTIONS_HPP

#include <boost/math/policies/policy.hpp>
#include <boost/math/distributions/complement.hpp>
#include <boost/math/distributions/detail/derived_accessors.hpp>

#include <nanobind/stl/pair.h>

#include <utility>

#include "vectorize.hpp"

namespace bmpy {

// The default policy, except that moments the distribution lacks (the mean of a
// Cauchy, ...) raise ValueError at run time instead of failing to compile.
using dist_policy = boost::math::policies::policy<boost::math::policies::assert_undefined<false>>;

template <class D, class G>
void def_dist_vec(nb::class_<D>& cls, const char* name, const char* argname, G g) {
   cls.def(name, [g](const D& d, double x) { return g(d, x); }, nb::arg(argname));
   cls.def(
       name,
       [g, name, argname](const D& d, nb::object x) {
          auto f = [&d, &g]<class T>(double v) { return g(d, v); };
          return call_array<false, Real>(name, f, {x}, {argname});
       },
       nb::arg(argname));
}

template <class D>
nb::class_<D>& add_dist_methods(nb::class_<D>& cls) {
   namespace bm = boost::math;
   def_dist_vec(cls, "pdf", "x", [](const D& d, double x) { return bm::pdf(d, x); });
   def_dist_vec(cls, "logpdf", "x", [](const D& d, double x) { return bm::logpdf(d, x); });
   def_dist_vec(cls, "cdf", "x", [](const D& d, double x) { return bm::cdf(d, x); });
   def_dist_vec(cls, "logcdf", "x", [](const D& d, double x) { return bm::logcdf(d, x); });
   def_dist_vec(cls, "cdf_complement", "x", [](const D& d, double x) { return bm::cdf(bm::complement(d, x)); });
   def_dist_vec(cls, "logcdf_complement", "x", [](const D& d, double x) { return bm::logcdf(bm::complement(d, x)); });
   def_dist_vec(cls, "quantile", "p", [](const D& d, double p) { return bm::quantile(d, p); });
   def_dist_vec(cls, "quantile_complement", "q", [](const D& d, double q) { return bm::quantile(bm::complement(d, q)); });
   def_dist_vec(cls, "hazard", "x", [](const D& d, double x) { return bm::hazard(d, x); });
   def_dist_vec(cls, "chf", "x", [](const D& d, double x) { return bm::chf(d, x); });

   // Unqualified calls (found by ADL) so that accessors only some distributions
   // declare, such as entropy, need no declaration here.
#define BMPY_DIST_ACCESSOR(fn)                       \
   if constexpr (requires(const D& d) { fn(d); })    \
      cls.def(#fn, [](const D& d) { return fn(d); });
   BMPY_DIST_ACCESSOR(mean)
   BMPY_DIST_ACCESSOR(variance)
   BMPY_DIST_ACCESSOR(standard_deviation)
   BMPY_DIST_ACCESSOR(skewness)
   BMPY_DIST_ACCESSOR(kurtosis)
   BMPY_DIST_ACCESSOR(kurtosis_excess)
   BMPY_DIST_ACCESSOR(mode)
   BMPY_DIST_ACCESSOR(median)
   BMPY_DIST_ACCESSOR(entropy)
   BMPY_DIST_ACCESSOR(support)
   BMPY_DIST_ACCESSOR(range)
#undef BMPY_DIST_ACCESSOR
   return cls;
}

}  // namespace bmpy

#endif
