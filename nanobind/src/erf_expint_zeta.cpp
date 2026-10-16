// Copyright 2026 Nick Thompson
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/math/special_functions/erf.hpp>
#include <boost/math/special_functions/expint.hpp>
#include <boost/math/special_functions/zeta.hpp>
#include <boost/math/special_functions/owens_t.hpp>
#include <boost/math/special_functions/lambert_w.hpp>

#include "bindings.hpp"
#include "macros.hpp"

void bmpy::bind_erf_expint_zeta(nb::module_& m) {
   BMPY_R1(erf, "z");
   BMPY_R1(erfc, "z");
   BMPY_R1(erf_inv, "p");
   BMPY_R1(erfc_inv, "q");
   BMPY_R1(expint, "z");                  // Ei
   BMPY_I1R1(expint, unsigned, "n", "z");  // E_n
   BMPY_R1(zeta, "s");
   BMPY_R2(owens_t, "h", "a");
   BMPY_R1(lambert_w0, "z");
   BMPY_R1(lambert_wm1, "z");
   BMPY_R1(lambert_w0_prime, "z");
   BMPY_R1(lambert_wm1_prime, "z");
}
