// Copyright 2026 Nick Thompson
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/math/special_functions/beta.hpp>

#include "bindings.hpp"
#include "macros.hpp"

void bmpy::bind_beta(nb::module_& m) {
   BMPY_R2(beta, "a", "b");
   BMPY_R3(beta, "a", "b", "x");  // non-normalised incomplete
   BMPY_R3(betac, "a", "b", "x");
   BMPY_R2(lbeta, "a", "b");
   BMPY_R3(ibeta, "a", "b", "x");
   BMPY_R3(ibetac, "a", "b", "x");
   BMPY_R3(libeta, "a", "b", "x");
   BMPY_R3(libetac, "a", "b", "x");
   BMPY_R3(ibeta_inv, "a", "b", "p");
   BMPY_R3(ibetac_inv, "a", "b", "q");
   BMPY_R3(ibeta_inva, "b", "x", "p");
   BMPY_R3(ibeta_invb, "a", "x", "p");
   BMPY_R3(ibetac_inva, "b", "x", "q");
   BMPY_R3(ibetac_invb, "a", "x", "q");
   BMPY_R3(ibeta_derivative, "a", "b", "x");
}
