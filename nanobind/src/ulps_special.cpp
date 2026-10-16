// Copyright 2026 Nick Thompson
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/math/special_functions/erf.hpp>
#include <boost/math/special_functions/expint.hpp>
#include <boost/math/special_functions/zeta.hpp>
#include <boost/math/special_functions/owens_t.hpp>
#include <boost/math/special_functions/lambert_w.hpp>
#include <boost/math/special_functions/log1p.hpp>
#include <boost/math/special_functions/expm1.hpp>
#include <boost/math/special_functions/cbrt.hpp>
#include <boost/math/special_functions/sqrt1pm1.hpp>
#include <boost/math/special_functions/powm1.hpp>
#include <boost/math/special_functions/pow1p.hpp>
#include <boost/math/special_functions/sin_pi.hpp>
#include <boost/math/special_functions/cos_pi.hpp>
#include <boost/math/special_functions/sinc.hpp>
#include <boost/math/special_functions/sinhc.hpp>
#include <boost/math/special_functions/asinh.hpp>
#include <boost/math/special_functions/acosh.hpp>
#include <boost/math/special_functions/atanh.hpp>
#include <boost/math/special_functions/hypot.hpp>
#include <boost/math/special_functions/rsqrt.hpp>
#include <boost/math/special_functions/logit.hpp>
#include <boost/math/special_functions/logistic_sigmoid.hpp>
#include <boost/math/special_functions/logaddexp.hpp>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/bessel_prime.hpp>
#include <boost/math/special_functions/airy.hpp>
#include <boost/math/special_functions/hypergeometric_1F0.hpp>
#include <boost/math/special_functions/hypergeometric_0F1.hpp>

#include "ulps.hpp"

namespace bm = boost::math;

void bmpy::ulps::register_special(registry& r) {
   BMPY_ULPS1(erf);
   BMPY_ULPS1(erfc);
   BMPY_ULPS1(erf_inv);
   BMPY_ULPS1(erfc_inv);
   BMPY_ULPS1(expint);
   BMPY_ULPS1(zeta);
   BMPY_ULPS1(lambert_w0);
   BMPY_ULPS1(lambert_wm1);
   BMPY_ULPS1(lambert_w0_prime);
   BMPY_ULPS1(lambert_wm1_prime);
   BMPY_ULPS2(owens_t);

   BMPY_ULPS1(log1p);
   BMPY_ULPS1(log1pmx);
   BMPY_ULPS1(expm1);
   BMPY_ULPS1(cbrt);
   BMPY_ULPS1(sqrt1pm1);
   BMPY_ULPS1(sin_pi);
   BMPY_ULPS1(cos_pi);
   BMPY_ULPS1(sinc_pi);
   BMPY_ULPS1(sinhc_pi);
   BMPY_ULPS1(asinh);
   BMPY_ULPS1(acosh);
   BMPY_ULPS1(atanh);
   BMPY_ULPS1(rsqrt);
   BMPY_ULPS1(logit);
   BMPY_ULPS1(logistic_sigmoid);
   BMPY_ULPS2(powm1);
   BMPY_ULPS2(pow1p);
   BMPY_ULPS2(hypot);
   BMPY_ULPS2(logaddexp);

   BMPY_ULPS2(cyl_bessel_j);
   BMPY_ULPS2(cyl_neumann);
   BMPY_ULPS2(cyl_bessel_i);
   BMPY_ULPS2(cyl_bessel_k);
   BMPY_ULPS2(cyl_bessel_j_prime);
   BMPY_ULPS2(cyl_neumann_prime);
   BMPY_ULPS2(cyl_bessel_i_prime);
   BMPY_ULPS2(cyl_bessel_k_prime);
   BMPY_ULPS1(airy_ai);
   BMPY_ULPS1(airy_bi);
   BMPY_ULPS1(airy_ai_prime);
   BMPY_ULPS1(airy_bi_prime);

   BMPY_ULPS2(hypergeometric_1F0);
   BMPY_ULPS2(hypergeometric_0F1);
}
