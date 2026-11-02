// Copyright 2026 Nick Thompson
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// boost.math.filters: filters::trend_filter on 1-D float64 data (anything np.ascontiguousarray understands).
// The solve releases the GIL; errors follow boost.math.errstate, as for a scalar argument, and a failed fit is NaN.

#include "bindings.hpp"

#include <boost/math/filters/trend_filter.hpp>

#include <nanobind/stl/optional.h>

#include <algorithm>
#include <limits>
#include <optional>
#include <vector>

namespace bmpy {

void bind_filters(nb::module_& m) {
   m.def("trend_filter", [](nb::handle x, std::optional<int> k, std::optional<double> lambda) {
      if (k.has_value() != lambda.has_value())
         throw nb::value_error("trend_filter: give both k and lambda_, or neither for cross-validation");
      nb::object a = numpy().attr("ascontiguousarray")(x, nb::arg("dtype") = numpy().attr("float64"));
      if (nb::cast<int>(a.attr("ndim")) != 1)
         throw nb::value_error("trend_filter: expected a 1-D array");
      auto d = nb::cast<nb::ndarray<double const, nb::ndim<1>, nb::c_contig, nb::device::cpu>>(a);
      std::vector<double> y(d.data(), d.data() + d.shape(0));
      error_actions actions = current_error_actions;
      error_record errors;
      {
         recording rec(errors);
         nb::gil_scoped_release release;
         try {
            y = k ? boost::math::filters::trend_filter(y, *k, *lambda) : boost::math::filters::trend_filter(y);
         } catch (std::exception const& e) {
            record_exception(e);
            std::fill(y.begin(), y.end(), std::numeric_limits<double>::quiet_NaN());
         }
      }
      if (errors.any())
         report("trend_filter", errors, 1, true, false, actions);
      return vector_to_array(std::move(y));
   }, nb::arg("x"), nb::arg("k") = nb::none(), nb::arg("lambda_") = nb::none(),
      nb::sig("def trend_filter(x, /, k=None, lambda_=None)"),
      "The trend filter of x: argmin_z ½ Σ (xᵢ − zᵢ)² + λ Σ |(Dᵏ⁺¹ z)ᵢ|, Dᵏ⁺¹ the (k + 1)-th difference.\n\n"
      "The fit is piecewise polynomial of degree k with knots chosen by the data: k = 0 is piecewise constant\n"
      "(total-variation denoising), k = 1 piecewise linear, k = 2 piecewise quadratic. λ = lambda_ ≥ 0; λ = 0\n"
      "returns x, and a large enough λ the least-squares polynomial of degree k. Give both k and lambda_, or\n"
      "neither: then k ∈ {0, 1, 2} and λ are chosen by cross-validation on the even- and odd-indexed samples\n"
      "(at least 4). Cross-validation assumes the noise is independent from sample to sample; with correlated\n"
      "noise it chooses too small a λ. If the noise estimate is 0 (noise-free data), x is returned unchanged.\n"
      "Returns a new float64 array."
      BMPY_URL("trend_filter"));
}

}  // namespace bmpy
