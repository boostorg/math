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
   BMPY_R1(erf, "z", "The error function erf(z) = (2/√π) ∫₀^z e^(−t²) dt." BMPY_URL("sf_erf/error_function"));
   BMPY_R1(erfc, "z", "The complementary error function erfc(z) = 1 − erf(z), computed directly."
      BMPY_URL("sf_erf/error_function"));
   BMPY_R1(erf_inv, "p", "The z with erf(z) = p, for p in [−1, 1]." BMPY_URL("sf_erf/error_inv"));
   BMPY_R1(erfc_inv, "q", "The z with erfc(z) = q, for q in [0, 2]." BMPY_URL("sf_erf/error_inv"));
   BMPY_R1(expint, "z", "The exponential integral Ei(z); as scipy.special.expi(z)." BMPY_URL("expint/expint_i"));
   BMPY_I1R1(expint, unsigned, "n", "z", "The exponential integral Eₙ(z) = ∫₁^∞ e^(−zt)/tⁿ dt; as scipy.special.expn(n, z)."
      BMPY_URL("expint/expint_n"));
   BMPY_R1(zeta, "s", "The Riemann zeta function ζ(s) = Σ_{k≥1} k^(−s), continued to s ≠ 1." BMPY_URL("zetas/zeta"));
   BMPY_R2(owens_t, "h", "a", "Owen's T function T(h, a) = (1/2π) ∫₀^a e^(−h²(1+x²)/2)/(1 + x²) dx." BMPY_URL("owens_t"));
   BMPY_R1(lambert_w0, "z", "The principal branch W₀ of the Lambert W function, W(z) e^W(z) = z, for z >= −1/e; W₀ >= −1."
      BMPY_URL("lambert_w"));
   BMPY_R1(lambert_wm1, "z", "The branch W₋₁ of the Lambert W function, for −1/e <= z < 0; W₋₁ <= −1." BMPY_URL("lambert_w"));
   BMPY_R1(lambert_w0_prime, "z", "W₀′(z) = W₀(z)/(z (1 + W₀(z)))." BMPY_URL("lambert_w"));
   BMPY_R1(lambert_wm1_prime, "z", "W₋₁′(z) = W₋₁(z)/(z (1 + W₋₁(z)))." BMPY_URL("lambert_w"));
}
