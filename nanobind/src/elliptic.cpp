// Copyright 2026 Nick Thompson
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/ellint_2.hpp>
#include <boost/math/special_functions/ellint_3.hpp>
#include <boost/math/special_functions/ellint_d.hpp>
#include <boost/math/special_functions/ellint_rc.hpp>
#include <boost/math/special_functions/ellint_rd.hpp>
#include <boost/math/special_functions/ellint_rf.hpp>
#include <boost/math/special_functions/ellint_rg.hpp>
#include <boost/math/special_functions/ellint_rj.hpp>
#include <boost/math/special_functions/jacobi_zeta.hpp>
#include <boost/math/special_functions/heuman_lambda.hpp>
#include <boost/math/special_functions/jacobi_elliptic.hpp>

#include "bindings.hpp"
#include "macros.hpp"

void bmpy::bind_elliptic(nb::module_& m) {
   BMPY_R1(ellint_1, "k", "The complete elliptic integral of the first kind K(k) = F(π/2, k).\n\n"
      "Takes the modulus k, not the parameter m = k² of scipy.special.ellipk(m) and mpmath." BMPY_URL("ellint/ellint_1"));
   BMPY_R2(ellint_1, "k", "phi",
      "The incomplete integral F(φ, k) = ∫₀^φ dθ/√(1 − k² sin²θ); scipy.special.ellipkinc(phi, k**2).");
   BMPY_R1(ellint_2, "k", "The complete elliptic integral of the second kind E(k) = E(π/2, k).\n\n"
      "Takes the modulus k, not the parameter m = k² of scipy.special.ellipe(m) and mpmath." BMPY_URL("ellint/ellint_2"));
   BMPY_R2(ellint_2, "k", "phi",
      "The incomplete integral E(φ, k) = ∫₀^φ √(1 − k² sin²θ) dθ; scipy.special.ellipeinc(phi, k**2).");
   BMPY_R2(ellint_3, "k", "n", "The complete elliptic integral of the third kind Π(n, k) = Π(n, π/2, k).\n\n"
      "Takes the modulus k (not m = k²) first; mpmath.ellippi(n, k**2) is the same integral." BMPY_URL("ellint/ellint_3"));
   BMPY_R3(ellint_3, "k", "n", "phi", "The incomplete integral Π(n, φ, k).");
   BMPY_R1(ellint_d, "k", "The complete elliptic integral D(k) = (K(k) − E(k))/k², modulus k." BMPY_URL("ellint/ellint_d"));
   BMPY_R2(ellint_d, "k", "phi", "The incomplete integral D(φ, k) = ∫₀^φ sin²θ/√(1 − k² sin²θ) dθ.");
   BMPY_R2(ellint_rc, "x", "y", "Carlson's degenerate integral R_C(x, y) = R_F(x, y, y)." BMPY_URL("ellint/ellint_carlson"));
   BMPY_R3(ellint_rd, "x", "y", "z", "Carlson's integral R_D(x, y, z) = R_J(x, y, z, z)." BMPY_URL("ellint/ellint_carlson"));
   BMPY_R3(ellint_rf, "x", "y", "z",
      "Carlson's symmetric integral of the first kind R_F(x, y, z)."
      BMPY_URL("ellint/ellint_carlson"));
   BMPY_R3(ellint_rg, "x", "y", "z", "Carlson's symmetric integral of the second kind R_G(x, y, z)."
      BMPY_URL("ellint/ellint_carlson"));
   BMPY_R4(ellint_rj, "x", "y", "z", "p",
      "Carlson's symmetric integral of the third kind R_J(x, y, z, p)."
      BMPY_URL("ellint/ellint_carlson"));
   BMPY_R2(jacobi_zeta, "k", "phi", "Jacobi's zeta function Z(φ, k) = E(φ, k) − E(k) F(φ, k)/K(k), modulus k."
      BMPY_URL("ellint/jacobi_zeta"));
   BMPY_R2(heuman_lambda, "k", "phi", "Heuman's lambda function Λ₀(φ, k), modulus k." BMPY_URL("ellint/heuman_lambda"));

   // A NaN modulus makes Boost's descending Landen recursion unbounded (a stack overflow), so it is NaN here.
#define BMPY_JACOBI(fn, doc)                                                                                          \
   bmpy::def_vec<bmpy::Real, bmpy::Real>(m, #fn, []<class T>(T k, T u) {                                             \
      return k != k ? std::numeric_limits<T>::quiet_NaN() : bm::fn(k, u); }, {"k", "u"}, doc)
   BMPY_JACOBI(jacobi_sn,
      "The Jacobi elliptic function sn(u, k) = sin(am(u, k)).\n\nNote the order (k, u) and the modulus k: scipy.special.ellipj(u, k**2)."
      BMPY_URL("jacobi/jacobi_sn"));
   BMPY_JACOBI(jacobi_cn, "The Jacobi elliptic function cn(u, k) = cos(am(u, k)), modulus k first."
      BMPY_URL("jacobi/jacobi_cn"));
   BMPY_JACOBI(jacobi_dn, "The Jacobi elliptic function dn(u, k) = √(1 − k² sn²(u, k)), modulus k first."
      BMPY_URL("jacobi/jacobi_dn"));
   BMPY_JACOBI(jacobi_cd, "The Jacobi elliptic function cd = cn/dn, modulus k first." BMPY_URL("jacobi/jacobi_cd"));
   BMPY_JACOBI(jacobi_dc, "The Jacobi elliptic function dc = dn/cn, modulus k first." BMPY_URL("jacobi/jacobi_dc"));
   BMPY_JACOBI(jacobi_ns, "The Jacobi elliptic function ns = 1/sn, modulus k first." BMPY_URL("jacobi/jacobi_ns"));
   BMPY_JACOBI(jacobi_sd, "The Jacobi elliptic function sd = sn/dn, modulus k first." BMPY_URL("jacobi/jacobi_sd"));
   BMPY_JACOBI(jacobi_ds, "The Jacobi elliptic function ds = dn/sn, modulus k first." BMPY_URL("jacobi/jacobi_ds"));
   BMPY_JACOBI(jacobi_nc, "The Jacobi elliptic function nc = 1/cn, modulus k first." BMPY_URL("jacobi/jacobi_nc"));
   BMPY_JACOBI(jacobi_nd, "The Jacobi elliptic function nd = 1/dn, modulus k first." BMPY_URL("jacobi/jacobi_nd"));
   BMPY_JACOBI(jacobi_sc, "The Jacobi elliptic function sc = sn/cn, modulus k first." BMPY_URL("jacobi/jacobi_sc"));
   BMPY_JACOBI(jacobi_cs, "The Jacobi elliptic function cs = cn/sn, modulus k first." BMPY_URL("jacobi/jacobi_cs"));
#undef BMPY_JACOBI
}
