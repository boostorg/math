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
   BMPY_R2(hypergeometric_1F0, "a", "z", "₁F₀(a; ; z) = (1 − z)^(−a)." BMPY_URL("hypergeometric/hypergeometric_1f0"));
   BMPY_R2(hypergeometric_0F1, "b", "z", "₀F₁(; b; z) = Σ zᵏ/((b)ₖ k!); as scipy.special.hyp0f1(b, z)."
      BMPY_URL("hypergeometric/hypergeometric_0f1"));
   BMPY_R3(hypergeometric_2F0, "a1", "a2", "z",
      "₂F₀(a1, a2; ; z) = Σ (a1)ₖ (a2)ₖ zᵏ/k!; only the polynomial case (a1 or a2 a non-positive integer) converges."
      BMPY_URL("hypergeometric/hypergeometric_2f0"));
   BMPY_R3(hypergeometric_1F1, "a", "b", "z",
      "Kummer's function ₁F₁(a; b; z) = M(a, b, z) = Σ (a)ₖ zᵏ/((b)ₖ k!); as scipy.special.hyp1f1(a, b, z)."
      BMPY_URL("hypergeometric/hypergeometric_1f1"));
   BMPY_R3(hypergeometric_1F1_regularized, "a", "b", "z", "₁F₁(a; b; z)/Γ(b), finite at b = 0, −1, −2, ..."
      BMPY_URL("hypergeometric/hypergeometric_1f1"));
   BMPY_R3(log_hypergeometric_1F1, "a", "b", "z", "log |₁F₁(a; b; z)|, finite where ₁F₁ overflows; the sign is not returned."
      BMPY_URL("hypergeometric/hypergeometric_1f1"));

   // pFq: the parameter lists are fixed per call, z is vectorised; double only.
   m.def(
       "hypergeometric_pFq",
       [](std::vector<double> a, std::vector<double> b, nb::handle z, int threads) {
          auto f = [a = std::move(a), b = std::move(b)]<class T>(double x) {
             return bm::hypergeometric_pFq(a, b, x);
          };
          return bmpy::call_array<false, bmpy::Real>("hypergeometric_pFq", f, {z}, {"z"}, threads);
       },
       nb::arg("a"), nb::arg("b"), nb::arg("z"), nb::kw_only(), nb::arg("threads") = 1,
       "The generalised hypergeometric function ₚF_q(a; b; z), for sequences a and b; elementwise in z, double only."
       BMPY_URL("hypergeometric/hypergeometric_pfq"));
}
