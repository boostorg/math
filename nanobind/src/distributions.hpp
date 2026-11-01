// Copyright 2026 Nick Thompson
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// Binds a Boost.Math distribution (double, default policy) as a Python class
// with vectorised pdf/cdf/quantile-style methods, scalar moment methods, rvs and interval.

#ifndef BOOST_MATH_NANOBIND_DISTRIBUTIONS_HPP
#define BOOST_MATH_NANOBIND_DISTRIBUTIONS_HPP

#include <boost/math/policies/policy.hpp>
#include <boost/math/distributions/complement.hpp>
#include <boost/math/distributions/detail/derived_accessors.hpp>

#include <nanobind/stl/pair.h>

#include <cstdint>
#include <initializer_list>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "vectorize.hpp"

namespace bmpy {

// The default policy, except that moments the distribution lacks (the mean of a
// Cauchy, ...) raise ValueError at run time instead of failing to compile.
using dist_policy = boost::math::policies::policy<boost::math::policies::assert_undefined<false>>;

// A distribution parameter of arithmetic type T. Anything that does not convert to T is let through, so that
// the constructor can say what is wrong with it instead of nanobind's generic "incompatible arguments".
template <class T>
struct param {
   T value{};
   PyObject* bad = nullptr;
};

template <class A>
decltype(auto) param_value(std::string const& type, A&& a) {
   if constexpr (requires { a.bad; }) {
      if (a.bad) {
         nb::handle h(a.bad);
         if (PyFloat_Check(a.bad) || PyLong_Check(a.bad))
            throw nb::value_error((type + ": invalid parameter " + nb::cast<std::string>(nb::repr(h))).c_str());
         throw nb::type_error((type + ": parameters must be scalars, not " + nb::cast<std::string>(nb::type_name(h.type()))).c_str());
      }
      return a.value;
   } else
      return std::forward<A>(a);
}

// nb::init<Args...>, with param<A> for the arithmetic Args.
template <class... Args>
struct init : nb::def_visitor<init<Args...>> {
   template <class C, class... Extra>
   static void execute(C& cls, const Extra&... extra) {
      using D = typename C::Type;
      std::string type = nb::cast<std::string>(cls.attr("__name__"));
      cls.def(
          "__init__",
          [type](D* self, std::conditional_t<std::is_arithmetic_v<Args>, param<Args>, Args>... args) {
             new (self) D(param_value(type, args)...);
          },
          extra...);
   }
};

// The docs page of a distribution, for its docstring.
#define BMPY_DIST_URL(page) "\n\nhttps://www.boost.org/doc/libs/release/libs/math/doc/html/math_toolkit/dist_ref/dists/" page ".html"

// The variate of d at u, a draw of numpy's Generator.random, k 2^-53 for 0 <= k < 2^53: the quantile of the
// midpoint of u's bin, exactly (2k + 1) 2^-54, which is in (0, 1). The upper half goes through the complement,
// whose argument, 1 minus the midpoint, is then exact too. A discrete quantile rounds up by default, to the
// least k with cdf(k) >= p.
template <bool Discrete, class D>
auto variate(D const& d, double u) {
   namespace bm = boost::math;
   constexpr double half_bin = 0x1p-54;
   double x = u < 0.5 ? bm::quantile(d, u + half_bin) : bm::quantile(bm::complement(d, (1 - u) - half_bin));
   if constexpr (Discrete) {
      if (!(x < 0x1p63))
         throw std::overflow_error("variate too large for a 64-bit integer");
      return static_cast<std::int64_t>(x);
   } else
      return x;
}

// A static method of the class whose arguments are positional only: sig is its Python signature, without the
// "@staticmethod\ndef ", and f is called with the arguments under the errors.hpp handling. The docstring ends
// with the last paragraph (the docs URL) of the class's.
template <class D, class... A>
void def_static_scalar(nb::class_<D>& cls, const char* name, const char* sig, double (*f)(A...), const char* doc) {
   nb::object class_doc = cls.attr("__doc__");
   std::string url = class_doc.is_none() ? "" : nb::cast<std::string>(class_doc);
   url = url.rfind("\n\n") == std::string::npos ? "" : url.substr(url.rfind("\n\n"));
   cls.def_static(name, [name, f](A... a) { return call_scalar(name, [&] { return f(a...); }); },
                  nb::sig((std::string("@staticmethod\ndef ") + sig).c_str()), (doc + url).c_str());
}

template <class D, bool Discrete = false>
nb::class_<D>& add_dist_methods(nb::class_<D>& cls) {
   namespace bm = boost::math;
   def_method_vec(cls, "pdf", "x", [](const D& d, double x) { return bm::pdf(d, x); });
   def_method_vec(cls, "logpdf", "x", [](const D& d, double x) { return bm::logpdf(d, x); });
   def_method_vec(cls, "cdf", "x", [](const D& d, double x) { return bm::cdf(d, x); });
   def_method_vec(cls, "logcdf", "x", [](const D& d, double x) { return bm::logcdf(d, x); });
   def_method_vec(cls, "cdf_complement", "x", [](const D& d, double x) { return bm::cdf(bm::complement(d, x)); });
   def_method_vec(cls, "logcdf_complement", "x", [](const D& d, double x) { return bm::logcdf(bm::complement(d, x)); });
   def_method_vec(cls, "quantile", "p", [](const D& d, double p) { return bm::quantile(d, p); });
   def_method_vec(cls, "quantile_complement", "q", [](const D& d, double q) { return bm::quantile(bm::complement(d, q)); });
   def_method_vec(cls, "hazard", "x", [](const D& d, double x) { return bm::hazard(d, x); });
   def_method_vec(cls, "chf", "x", [](const D& d, double x) { return bm::chf(d, x); });
   cls.def(
       "rvs",
       [](const D& d, nb::object size, nb::object rng) {
          if (rng.is_none())
             rng = numpy().attr("random").attr("default_rng")();
          auto f = [&d]<class T>(double u) { return variate<Discrete>(d, u); };
          return call_array<false, Real>("rvs", f, {rng.attr("random")(size)}, {"u"});
       },
       nb::arg("size") = nb::none(), nb::kw_only(), nb::arg("rng") = nb::none(), nb::sig("def rvs(self, size=None, *, rng=None)"),
       "Random variates by inversion: quantile(u) for uniform u from rng.random(size) (a numpy.random.Generator, default "
       "numpy.random.default_rng()), the upper half through quantile_complement(1 - u) for accuracy in the upper tail. "
       "size None gives a float (an int for a discrete distribution), else an array of that shape (of int64)."
       "\n\nhttps://www.boost.org/doc/libs/release/libs/math/doc/html/math_toolkit/stat_tut/variates.html");
   cls.def(
       "interval",
       [](const D& d, double confidence) {
          if (!(confidence >= 0 && confidence <= 1))
             throw nb::value_error("interval: confidence must be in [0, 1]");
          double a = (1 - confidence) / 2;
          return std::make_pair(call_scalar("interval", [&] { return bm::quantile(d, a); }),
                                call_scalar("interval", [&] { return bm::quantile(bm::complement(d, a)); }));
       },
       nb::arg(), nb::sig("def interval(self, confidence, /)"),
       "The equal-tailed interval (lower, upper) = (quantile(a), quantile_complement(a)), a = (1 - confidence) / 2, holding "
       "the given fraction of the probability; a discrete distribution rounds both ends up."
       "\n\nhttps://www.boost.org/doc/libs/release/libs/math/doc/html/math_toolkit/dist_ref/nmp.html");

   // Unqualified calls (found by ADL) so that accessors only some distributions
   // declare, such as entropy, need no declaration here.
#define BMPY_DIST_ACCESSOR(fn)                       \
   if constexpr (requires(const D& d) { fn(d); })    \
      cls.def(#fn, [](const D& d) { return call_scalar(#fn, [&] { return fn(d); }); });
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
   // Equal, and hashed alike, when the type and the parameters are; list parameters hash as tuples.
   auto key = [ps](const D& d) {
      nb::list k;
      for (auto const& p : ps) {
         nb::object v = p.second(d);
         k.append(nb::isinstance<nb::list>(v) ? nb::tuple(v) : v);
      }
      return nb::tuple(k);
   };
   cls.def("__eq__", [key](const D& a, nb::handle b) -> nb::object {
      if (!nb::isinstance<D>(b))
         return nb::borrow(Py_NotImplemented);
      return nb::bool_(key(a).equal(key(nb::cast<const D&>(b))));
   }, nb::sig("def __eq__(self, other: object, /) -> bool"));
   cls.def("__hash__", [key](const D& d) { return nb::hash(key(d)); });
   return cls;
}

}  // namespace bmpy

namespace nanobind::detail {

template <class T>
struct type_caster<bmpy::param<T>> {
   NB_TYPE_CASTER(bmpy::param<T>, make_caster<T>::Name)
   template <class Flags>
   bool from_python(handle src, Flags flags, cleanup_list* cleanup) noexcept {
      make_caster<T> c;
      value.bad = nullptr;
      if (c.from_python(src, flags, cleanup)) {
         value.value = c.value;
         return true;
      }
      // Rejected outright only in the first, no-conversion pass over the overloads.
      if (!(flags & static_cast<Flags>(cast_flags::convert)))
         return false;
      value.bad = src.ptr();
      return true;
   }
   static handle from_cpp(bmpy::param<T> const& p, rv_policy policy, cleanup_list* cleanup) noexcept {
      return make_caster<T>::from_cpp(p.value, policy, cleanup);
   }
};

}  // namespace nanobind::detail

#endif
