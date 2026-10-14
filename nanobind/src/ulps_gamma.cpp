// Copyright 2026 Nick Thompson
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/trigamma.hpp>
#include <boost/math/special_functions/beta.hpp>

#include "ulps.hpp"

namespace bm = boost::math;

void bmpy::ulps::register_gamma(registry& r) {
   BMPY_ULPS1(tgamma);
   BMPY_ULPS1(tgamma1pm1);
   BMPY_ULPS1(lgamma);
   BMPY_ULPS1(digamma);
   BMPY_ULPS1(trigamma);

   BMPY_ULPS2(tgamma);  // upper incomplete
   BMPY_ULPS2(tgamma_ratio);
   BMPY_ULPS2(tgamma_delta_ratio);
   BMPY_ULPS2(tgamma_lower);
   BMPY_ULPS2(gamma_p);
   BMPY_ULPS2(gamma_q);
   BMPY_ULPS2(gamma_p_inv);
   BMPY_ULPS2(gamma_q_inv);
   BMPY_ULPS2(gamma_p_inva);
   BMPY_ULPS2(gamma_q_inva);
   BMPY_ULPS2(gamma_p_derivative);
   BMPY_ULPS2(lgamma_p);
   BMPY_ULPS2(lgamma_q);
   BMPY_ULPS2(ligamma);
   BMPY_ULPS2(ligamma_lower);
   BMPY_ULPS2(beta);
   BMPY_ULPS2(lbeta);
}
