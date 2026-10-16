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

#include <initializer_list>
#include <string>
#include <utility>
#include <vector>

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

// A parameter: the constructor's argument name, and its value read back from the distribution.
template <class D>
using dist_param = std::pair<const char*, nb::object (*)(const D&)>;

#define BMPY_PARAM(name, accessor) \
   bmpy::dist_param<D>{name, +[](const D& d) { return nb::cast(d.accessor()); }}

// Exposes the parameters as read-only properties, under the constructor's argument names, and a repr that
// constructs the same distribution: normal_distribution(mean=0.0, sd=1.0). A parameter named like a method
// (the mean of a normal, the mode of a triangular) stays the method, which returns the same value.
template <class D>
nb::class_<D>& add_dist_params(nb::class_<D>& cls, std::initializer_list<dist_param<D>> params) {
   std::vector<dist_param<D>> ps(params);
   for (auto const& [name, get] : ps) {
      if (!nb::hasattr(cls, name)) {
         cls.def_prop_ro(name, get);
      }
   }
   std::string type = nb::cast<std::string>(cls.attr("__name__"));
   cls.def("__repr__", [type, ps](const D& d) {
      std::string s = type + "(";
      for (std::size_t i = 0; i < ps.size(); ++i) {
         s += (i ? ", " : "") + std::string(ps[i].first) + "=" + nb::cast<std::string>(nb::repr(ps[i].second(d)));
      }
      return s + ")";
   });
   return cls;
}

}  // namespace bmpy

#endif
