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
   BMPY_R1(ellint_1, "k");
   BMPY_R2(ellint_1, "k", "phi");
   BMPY_R1(ellint_2, "k");
   BMPY_R2(ellint_2, "k", "phi");
   BMPY_R2(ellint_3, "k", "n");
   BMPY_R3(ellint_3, "k", "n", "phi");
   BMPY_R1(ellint_d, "k");
   BMPY_R2(ellint_d, "k", "phi");
   BMPY_R2(ellint_rc, "x", "y");
   BMPY_R3(ellint_rd, "x", "y", "z");
   BMPY_R3(ellint_rf, "x", "y", "z");
   BMPY_R3(ellint_rg, "x", "y", "z");
   BMPY_R4(ellint_rj, "x", "y", "z", "p");
   BMPY_R2(jacobi_zeta, "k", "phi");
   BMPY_R2(heuman_lambda, "k", "phi");

   BMPY_R2(jacobi_sn, "k", "u");
   BMPY_R2(jacobi_cn, "k", "u");
   BMPY_R2(jacobi_dn, "k", "u");
   BMPY_R2(jacobi_cd, "k", "u");
   BMPY_R2(jacobi_dc, "k", "u");
   BMPY_R2(jacobi_ns, "k", "u");
   BMPY_R2(jacobi_sd, "k", "u");
   BMPY_R2(jacobi_ds, "k", "u");
   BMPY_R2(jacobi_nc, "k", "u");
   BMPY_R2(jacobi_nd, "k", "u");
   BMPY_R2(jacobi_sc, "k", "u");
   BMPY_R2(jacobi_cs, "k", "u");
}
