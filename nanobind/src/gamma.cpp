// Copyright 2026 Nick Thompson
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/trigamma.hpp>
#include <boost/math/special_functions/polygamma.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/special_functions/binomial.hpp>
#include <boost/math/special_functions/bernoulli.hpp>

#include "bindings.hpp"
#include "macros.hpp"

void bmpy::bind_gamma(nb::module_& m) {
   BMPY_R1(tgamma, "z");
   BMPY_R2(tgamma, "a", "z");  // upper incomplete
   BMPY_R1(tgamma1pm1, "z");
   BMPY_R1(lgamma, "z");
   BMPY_R1(digamma, "z");
   BMPY_R1(trigamma, "z");
   BMPY_I1R1(polygamma, int, "n", "z");
   BMPY_R2(tgamma_ratio, "a", "b");
   BMPY_R2(tgamma_delta_ratio, "a", "delta");
   BMPY_R2(tgamma_lower, "a", "z");
   BMPY_R2(gamma_p, "a", "z");
   BMPY_R2(gamma_q, "a", "z");
   BMPY_R2(gamma_p_inv, "a", "p");
   BMPY_R2(gamma_q_inv, "a", "q");
   BMPY_R2(gamma_p_inva, "z", "p");
   BMPY_R2(gamma_q_inva, "z", "q");
   BMPY_R2(gamma_p_derivative, "a", "z");
   BMPY_R2(lgamma_p, "a", "z");
   BMPY_R2(lgamma_q, "a", "z");
   BMPY_R2(ligamma, "a", "z");
   BMPY_R2(ligamma_lower, "a", "z");

   BMPY_T_I1(factorial, unsigned, "n");
   BMPY_T_I1(unchecked_factorial, unsigned, "n");
   BMPY_T_I1(double_factorial, unsigned, "n");
   BMPY_R1I1(rising_factorial, int, "x", "n");
   BMPY_R1I1(falling_factorial, unsigned, "x", "n");
   bmpy::def_vec<bmpy::Int<unsigned>, bmpy::Int<unsigned>>(
       m, "binomial_coefficient", []<class T>(unsigned n, unsigned k) { return bm::binomial_coefficient<double>(n, k); },
       {"n", "k"}, "boost::math::binomial_coefficient<double>");
   BMPY_T_I1(bernoulli_b2n, int, "n");
   BMPY_T_I1(tangent_t2n, int, "n");
}
