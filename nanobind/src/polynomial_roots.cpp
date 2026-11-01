// Copyright 2026 Nick Thompson
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// tools::polynomial_roots<double>, coefficients in increasing degree (a_0 first). Also the global Aberth
// iteration it falls back to, for comparison: when the two agree bit for bit, the fallback was taken.

#include "bindings.hpp"

#include <boost/math/tools/polynomial_roots.hpp>
#include <nanobind/stl/vector.h>

#include <complex>
#include <vector>

namespace bmpy {

namespace {

using C = std::complex<double>;

nb::ndarray<nb::numpy, C, nb::ndim<1>> to_array(std::vector<C>&& roots)
{
   auto* data = new std::vector<C>(std::move(roots));
   nb::capsule owner(data, [](void* p) noexcept { delete static_cast<std::vector<C>*>(p); });
   return nb::ndarray<nb::numpy, C, nb::ndim<1>>(data->data(), {data->size()}, owner);
}

}  // namespace

void bind_polynomial_roots(nb::module_& m)
{
   namespace bm = boost::math;
   m.def("polynomial_roots", [](std::vector<double> a) {
      std::vector<C> roots;
      {
         nb::gil_scoped_release release;
         roots = bm::tools::polynomial_roots(a);
      }
      return to_array(std::move(roots));
   }, nb::arg("coefficients"),
      "The n complex roots of a_0 + a_1 x + ... + a_n x^n, given the coefficients in increasing order (the reverse of "
      "numpy.roots)." BMPY_URL("polynomial_roots"));
   m.def("polynomial_roots_aberth", [](std::vector<double> a) {
      if (a.empty())
         throw nb::value_error("polynomial_roots_aberth: no coefficients");
      std::vector<C> roots;
      {
         nb::gil_scoped_release release;
         roots = bm::tools::polynomial_roots_detail::aberth(std::vector<C>(a.begin(), a.end()));
      }
      return to_array(std::move(roots));
   }, nb::arg("coefficients"),
      "The roots by the Aberth-Ehrlich iteration alone, the fallback of polynomial_roots, for comparison; coefficients "
      "in increasing order." BMPY_URL("polynomial_roots"));
   m.def("polynomial_roots_companion_matrix", [](std::vector<double> a) {
      std::vector<C> roots;
      {
         nb::gil_scoped_release release;
         roots = bm::tools::polynomial_roots_companion_matrix(a);
      }
      return to_array(std::move(roots));
   }, nb::arg("coefficients"),
      "The roots as the eigenvalues of the companion matrix, by the QR algorithm of Aurentz, Mach, Vandebril and "
      "Watkins; coefficients in increasing order." BMPY_URL("polynomial_roots"));
}

}  // namespace bmpy
