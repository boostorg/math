// Copyright 2026 Nick Thompson
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/bessel_prime.hpp>
#include <boost/math/special_functions/airy.hpp>
#include <boost/math/special_functions/hankel.hpp>

#include "bindings.hpp"
#include "macros.hpp"

void bmpy::bind_bessel(nb::module_& m) {
   BMPY_R2(cyl_bessel_j, "v", "x", "The Bessel function of the first kind J_v(x); as scipy.special.jv."
      BMPY_URL("bessel/bessel_first"));
   BMPY_R2(cyl_neumann, "v", "x", "The Bessel function of the second kind Y_v(x); as scipy.special.yv."
      BMPY_URL("bessel/bessel_first"));
   BMPY_R2(cyl_bessel_i, "v", "x", "The modified Bessel function of the first kind I_v(x); as scipy.special.iv."
      BMPY_URL("bessel/mbessel"));
   BMPY_R2(cyl_bessel_k, "v", "x", "The modified Bessel function of the second kind K_v(x); as scipy.special.kv."
      BMPY_URL("bessel/mbessel"));
   BMPY_R2(cyl_bessel_j_prime, "v", "x", "J_v′(x), the derivative in x." BMPY_URL("bessel/bessel_derivatives"));
   BMPY_R2(cyl_neumann_prime, "v", "x", "Y_v′(x), the derivative in x." BMPY_URL("bessel/bessel_derivatives"));
   BMPY_R2(cyl_bessel_i_prime, "v", "x", "I_v′(x), the derivative in x." BMPY_URL("bessel/bessel_derivatives"));
   BMPY_R2(cyl_bessel_k_prime, "v", "x", "K_v′(x), the derivative in x." BMPY_URL("bessel/bessel_derivatives"));
   BMPY_I1R1(sph_bessel, unsigned, "n", "x",
      "The spherical Bessel function of the first kind jₙ(x) = √(π/2x) J_(n+½)(x); as scipy.special.spherical_jn."
      BMPY_URL("bessel/sph_bessel"));
   BMPY_I1R1(sph_neumann, unsigned, "n", "x",
      "The spherical Bessel function of the second kind yₙ(x) = √(π/2x) Y_(n+½)(x); as scipy.special.spherical_yn."
      BMPY_URL("bessel/sph_bessel"));
   BMPY_I1R1(sph_bessel_prime, unsigned, "n", "x", "jₙ′(x), the derivative in x." BMPY_URL("bessel/bessel_derivatives"));
   BMPY_I1R1(sph_neumann_prime, unsigned, "n", "x", "yₙ′(x), the derivative in x." BMPY_URL("bessel/bessel_derivatives"));
   BMPY_R1I1(cyl_bessel_j_zero, int, "v", "m", "The mth positive zero of J_v, m = 1, 2, ... (1-based)."
      BMPY_URL("bessel/bessel_root"));
   BMPY_R1I1(cyl_neumann_zero, int, "v", "m", "The mth positive zero of Y_v, m = 1, 2, ... (1-based)."
      BMPY_URL("bessel/bessel_root"));

   // Complex results.
   BMPY_R2(cyl_hankel_1, "v", "x", "The Hankel function H_v⁽¹⁾(x) = J_v(x) + i Y_v(x), complex." BMPY_URL("hankel/cyl_hankel"));
   BMPY_R2(cyl_hankel_2, "v", "x", "The Hankel function H_v⁽²⁾(x) = J_v(x) − i Y_v(x), complex." BMPY_URL("hankel/cyl_hankel"));
   BMPY_R2(sph_hankel_1, "v", "x", "The spherical Hankel function h_v⁽¹⁾(x) = √(π/2x) H_(v+½)⁽¹⁾(x), complex."
      BMPY_URL("hankel/sph_hankel"));
   BMPY_R2(sph_hankel_2, "v", "x", "The spherical Hankel function h_v⁽²⁾(x) = √(π/2x) H_(v+½)⁽²⁾(x), complex."
      BMPY_URL("hankel/sph_hankel"));

   BMPY_R1(airy_ai, "x", "The Airy function Ai(x), the solution of y″ = x y that decays as x → ∞." BMPY_URL("airy/ai"));
   BMPY_R1(airy_bi, "x", "The Airy function Bi(x), the solution of y″ = x y that grows as x → ∞." BMPY_URL("airy/bi"));
   BMPY_R1(airy_ai_prime, "x", "Ai′(x)." BMPY_URL("airy/aip"));
   BMPY_R1(airy_bi_prime, "x", "Bi′(x)." BMPY_URL("airy/bip"));
   BMPY_T_I1(airy_ai_zero, int, "m", "The mth zero of Ai, m = 1, 2, ... (1-based); the zeros are negative."
      BMPY_URL("airy/airy_root"));
   BMPY_T_I1(airy_bi_zero, int, "m", "The mth zero of Bi, m = 1, 2, ... (1-based); the zeros are negative."
      BMPY_URL("airy/airy_root"));
}
