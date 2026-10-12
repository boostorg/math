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
   BMPY_R2(cyl_bessel_j, "v", "x");
   BMPY_R2(cyl_neumann, "v", "x");
   BMPY_R2(cyl_bessel_i, "v", "x");
   BMPY_R2(cyl_bessel_k, "v", "x");
   BMPY_R2(cyl_bessel_j_prime, "v", "x");
   BMPY_R2(cyl_neumann_prime, "v", "x");
   BMPY_R2(cyl_bessel_i_prime, "v", "x");
   BMPY_R2(cyl_bessel_k_prime, "v", "x");
   BMPY_I1R1(sph_bessel, unsigned, "n", "x");
   BMPY_I1R1(sph_neumann, unsigned, "n", "x");
   BMPY_I1R1(sph_bessel_prime, unsigned, "n", "x");
   BMPY_I1R1(sph_neumann_prime, unsigned, "n", "x");
   BMPY_R1I1(cyl_bessel_j_zero, int, "v", "m");
   BMPY_R1I1(cyl_neumann_zero, int, "v", "m");

   // Complex results.
   BMPY_R2(cyl_hankel_1, "v", "x");
   BMPY_R2(cyl_hankel_2, "v", "x");
   BMPY_R2(sph_hankel_1, "v", "x");
   BMPY_R2(sph_hankel_2, "v", "x");

   BMPY_R1(airy_ai, "x");
   BMPY_R1(airy_bi, "x");
   BMPY_R1(airy_ai_prime, "x");
   BMPY_R1(airy_bi_prime, "x");
   BMPY_T_I1(airy_ai_zero, int, "m");
   BMPY_T_I1(airy_bi_zero, int, "m");
}
