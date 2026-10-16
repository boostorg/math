// Copyright 2026 Nick Thompson
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/ellint_2.hpp>
#include <boost/math/special_functions/ellint_3.hpp>
#include <boost/math/special_functions/ellint_d.hpp>
#include <boost/math/special_functions/ellint_rc.hpp>
#include <boost/math/special_functions/jacobi_zeta.hpp>
#include <boost/math/special_functions/heuman_lambda.hpp>
#include <boost/math/special_functions/jacobi_elliptic.hpp>
#include <boost/math/special_functions/jacobi_theta.hpp>

#include "ulps.hpp"

namespace bm = boost::math;

void bmpy::ulps::register_elliptic(registry& r) {
   BMPY_ULPS1(ellint_1);
   BMPY_ULPS2(ellint_1);
   BMPY_ULPS1(ellint_2);
   BMPY_ULPS2(ellint_2);
   BMPY_ULPS2(ellint_3);
   BMPY_ULPS1(ellint_d);
   BMPY_ULPS2(ellint_d);
   BMPY_ULPS2(ellint_rc);
   BMPY_ULPS2(jacobi_zeta);
   BMPY_ULPS2(heuman_lambda);

   BMPY_ULPS2(jacobi_sn);
   BMPY_ULPS2(jacobi_cn);
   BMPY_ULPS2(jacobi_dn);
   BMPY_ULPS2(jacobi_cd);
   BMPY_ULPS2(jacobi_dc);
   BMPY_ULPS2(jacobi_ns);
   BMPY_ULPS2(jacobi_sd);
   BMPY_ULPS2(jacobi_ds);
   BMPY_ULPS2(jacobi_nc);
   BMPY_ULPS2(jacobi_nd);
   BMPY_ULPS2(jacobi_sc);
   BMPY_ULPS2(jacobi_cs);

   BMPY_ULPS2(jacobi_theta1);
   BMPY_ULPS2(jacobi_theta2);
   BMPY_ULPS2(jacobi_theta3);
   BMPY_ULPS2(jacobi_theta4);
   BMPY_ULPS2(jacobi_theta3m1);
   BMPY_ULPS2(jacobi_theta4m1);
   BMPY_ULPS2(jacobi_theta1tau);
   BMPY_ULPS2(jacobi_theta2tau);
   BMPY_ULPS2(jacobi_theta3tau);
   BMPY_ULPS2(jacobi_theta4tau);
   BMPY_ULPS2(jacobi_theta3m1tau);
   BMPY_ULPS2(jacobi_theta4m1tau);
}
