// Copyright 2026 Nick Thompson
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/math/special_functions/hypergeometric_1F0.hpp>
#include <boost/math/special_functions/hypergeometric_0F1.hpp>
#include <boost/math/special_functions/hypergeometric_2F0.hpp>
#include <boost/math/special_functions/hypergeometric_1F1.hpp>
#include <boost/math/special_functions/hypergeometric_pFq.hpp>

#include <nanobind/stl/vector.h>

#include "bindings.hpp"
#include "macros.hpp"

void bmpy::bind_hypergeometric(nb::module_& m) {
   BMPY_R2(hypergeometric_1F0, "a", "z");
   BMPY_R2(hypergeometric_0F1, "b", "z");
   BMPY_R3(hypergeometric_2F0, "a1", "a2", "z");
   BMPY_R3(hypergeometric_1F1, "a", "b", "z");
   BMPY_R3(hypergeometric_1F1_regularized, "a", "b", "z");
   BMPY_R3(log_hypergeometric_1F1, "a", "b", "z");

   // pFq: the parameter lists are fixed per call, z is vectorised; double only.
   m.def(
       "hypergeometric_pFq",
       [](std::vector<double> a, std::vector<double> b, nb::object z) {
          auto f = [a = std::move(a), b = std::move(b)]<class T>(double x) {
             return bm::hypergeometric_pFq(a, b, x);
          };
          return bmpy::call_array<false, bmpy::Real>("hypergeometric_pFq", f, {z}, {"z"});
       },
       nb::arg("a"), nb::arg("b"), nb::arg("z"), "boost::math::hypergeometric_pFq(a, b, z); vectorised over z, double only");
}
