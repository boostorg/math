// Copyright 2026 Nick Thompson
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// daubechies_scaling and daubechies_wavelet for p chosen at run time, through a virtual interface over the
// instantiations p = 2, ..., 19. (p = 1, the Haar functions, is left out: daubechies_scaling<double, 1>
// dereferences a null interpolator.)

#include <boost/math/special_functions/daubechies_scaling.hpp>
#include <boost/math/special_functions/daubechies_wavelet.hpp>

#include <nanobind/stl/pair.h>

#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>

#include "bindings.hpp"
#include "macros.hpp"

namespace {

namespace nb = nanobind;

struct daubechies_base {
   virtual ~daubechies_base() = default;
   virtual double value(double x) const = 0;
   virtual double prime(double x) const = 0;
   virtual double double_prime(double x) const = 0;
   virtual std::pair<double, double> support() const = 0;
   virtual std::int64_t bytes() const = 0;
};

template <class D, int p>
struct daubechies_impl final : daubechies_base {
   D d;
   explicit daubechies_impl(int grid_refinements) : d(grid_refinements) {}
   double value(double x) const override { return d(x); }
   double prime(double x) const override {
      if constexpr (p > 2)
         return d.prime(x);
      else
         throw std::domain_error("prime needs p >= 3");
   }
   double double_prime(double x) const override {
      if constexpr (p >= 6)
         return d.double_prime(x);
      else
         throw std::domain_error("double_prime needs p >= 6");
   }
   std::pair<double, double> support() const override { return d.support(); }
   std::int64_t bytes() const override { return d.bytes(); }
};

template <bool Wavelet>
struct py_daubechies {
   int p, grid_refinements;
   std::unique_ptr<daubechies_base> impl;
   py_daubechies(int p_, int grid_refinements_) : p(p_), grid_refinements(grid_refinements_) {
      // Beyond 23 the grid outgrows memory (and 1 << grid_refinements, int).
      if (grid_refinements != -1 && (grid_refinements < 1 || grid_refinements > 23))
         throw nb::value_error("grid_refinements must be -1 (the default) or in [1, 23]");
      nb::gil_scoped_release release;
      impl = bmpy::with_int<2, 19>("p", p, [&](auto P) -> std::unique_ptr<daubechies_base> {
         using D = std::conditional_t<Wavelet, bm::daubechies_wavelet<double, P()>, bm::daubechies_scaling<double, P()>>;
         return std::make_unique<daubechies_impl<D, P()>>(grid_refinements);
      });
   }
};

template <bool Wavelet>
void bind_one(nb::module_& m, char const* name, char const* doc) {
   using C = py_daubechies<Wavelet>;
   nb::class_<C> c(m, name, doc);
   c.def(nb::init<int, int>(), nb::arg("p"), nb::kw_only(), nb::arg("grid_refinements") = -1,
         nb::sig("def __init__(self, p, /, *, grid_refinements=-1) -> None"))
      .def_ro("p", &C::p)
      .def_prop_ro("support", [](C const& d) { return d.impl->support(); }, "The interval outside which the function is 0.")
      .def_prop_ro("bytes", [](C const& d) { return d.impl->bytes(); }, "The memory used by the interpolation table.")
      .def("__repr__", [name](C const& d) {
         std::string r = std::string(name) + "(" + std::to_string(d.p);
         return r + (d.grid_refinements == -1 ? ")" : ", grid_refinements=" + std::to_string(d.grid_refinements) + ")");
      });
   bmpy::def_method_vec(c, "__call__", "x", [](C const& d, double x) { return d.impl->value(x); }, "The function at x.");
   bmpy::def_method_vec(c, "prime", "x", [](C const& d, double x) { return d.impl->prime(x); }, "The derivative at x; p >= 3.");
   bmpy::def_method_vec(c, "double_prime", "x", [](C const& d, double x) { return d.impl->double_prime(x); },
                        "The second derivative at x; p >= 6.");
}

}  // namespace

void bmpy::bind_daubechies(nb::module_& m) {
   bind_one<false>(m, "daubechies_scaling",
      "The Daubechies scaling function with p vanishing moments, p in [2, 19], interpolated from a dyadic grid of spacing "
      "2^-grid_refinements (-1: Boost's default, which may take hundreds of MB)." BMPY_URL("daubechies"));
   bind_one<true>(m, "daubechies_wavelet",
      "The Daubechies wavelet with p vanishing moments, p in [2, 19]; grid_refinements as for daubechies_scaling."
      BMPY_URL("daubechies"));
}
