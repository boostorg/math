// Copyright 2026 Nick Thompson
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/math/special_functions/beta.hpp>

#include "bindings.hpp"
#include "macros.hpp"

void bmpy::bind_beta(nb::module_& m) {
   BMPY_R2(beta, "a", "b", "The beta function B(a, b) = Γ(a) Γ(b)/Γ(a + b)." BMPY_URL("sf_beta/beta_function"));
   BMPY_R3(beta, "a", "b", "x", "The incomplete beta function B_x(a, b), not normalised (ibeta is).");
   BMPY_R3(betac, "a", "b", "x", "B(a, b) − B_x(a, b)." BMPY_URL("sf_beta/ibeta_function"));
   BMPY_R2(lbeta, "a", "b", "log B(a, b), finite where B(a, b) under- or overflows." BMPY_URL("sf_beta/beta_function"));
   BMPY_R3(ibeta, "a", "b", "x",
      "The regularised incomplete beta function I_x(a, b) = B_x(a, b)/B(a, b); as scipy.special.betainc(a, b, x)."
      BMPY_URL("sf_beta/ibeta_function"));
   BMPY_R3(ibetac, "a", "b", "x", "1 − I_x(a, b), computed directly; as scipy.special.betaincc(a, b, x)."
      BMPY_URL("sf_beta/ibeta_function"));
   BMPY_R3(libeta, "a", "b", "x", "log I_x(a, b), finite where I_x(a, b) underflows." BMPY_URL("sf_beta/ibeta_function"));
   BMPY_R3(libetac, "a", "b", "x", "log(1 − I_x(a, b)), finite where 1 − I_x(a, b) underflows."
      BMPY_URL("sf_beta/ibeta_function"));
   BMPY_R3(ibeta_inv, "a", "b", "p", "The x with I_x(a, b) = p; as scipy.special.betaincinv(a, b, p)."
      BMPY_URL("sf_beta/ibeta_inv_function"));
   BMPY_R3(ibetac_inv, "a", "b", "q", "The x with 1 − I_x(a, b) = q; as scipy.special.betainccinv(a, b, q)."
      BMPY_URL("sf_beta/ibeta_inv_function"));
   BMPY_R3(ibeta_inva, "b", "x", "p", "The a with I_x(a, b) = p; note the argument order (b, x, p)."
      BMPY_URL("sf_beta/ibeta_inv_function"));
   BMPY_R3(ibeta_invb, "a", "x", "p", "The b with I_x(a, b) = p; note the argument order (a, x, p)."
      BMPY_URL("sf_beta/ibeta_inv_function"));
   BMPY_R3(ibetac_inva, "b", "x", "q", "The a with 1 − I_x(a, b) = q; note the argument order (b, x, q)."
      BMPY_URL("sf_beta/ibeta_inv_function"));
   BMPY_R3(ibetac_invb, "a", "x", "q", "The b with 1 − I_x(a, b) = q; note the argument order (a, x, q)."
      BMPY_URL("sf_beta/ibeta_inv_function"));
   BMPY_R3(ibeta_derivative, "a", "b", "x",
      "∂I_x(a, b)/∂x = x^(a−1) (1 − x)^(b−1)/B(a, b), the density of the beta distribution." BMPY_URL("sf_beta/beta_derivative"));
}
