// Copyright 2026 Nick Thompson
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/math/constants/constants.hpp>

#include "bindings.hpp"
#include "macros.hpp"

// Every constant of constants.hpp, in its order there (pow23_four_minus_pi is commented out there).
#define BMPY_CONSTANTS(X)                                                                                               \
   X(half) X(third) X(twothirds) X(two_thirds) X(sixth) X(three_quarters) X(root_two) X(root_three) X(half_root_two)    \
   X(ln_two) X(ln_ln_two) X(root_ln_four) X(one_div_root_two) X(pi) X(half_pi) X(third_pi) X(sixth_pi) X(two_pi)      \
   X(two_thirds_pi) X(three_quarters_pi) X(four_thirds_pi) X(one_div_two_pi) X(one_div_root_two_pi) X(log_pi)          \
   X(root_pi) X(root_half_pi) X(root_two_pi) X(log_root_two_pi) X(one_div_root_pi) X(root_one_div_pi)                  \
   X(pi_minus_three) X(four_minus_pi) X(pi_pow_e) X(pi_sqr) X(pi_sqr_div_six) X(pi_cubed)       \
   X(cbrt_pi) X(one_div_cbrt_pi) X(log2_e) X(e) X(exp_minus_half) X(exp_minus_one) X(e_pow_pi) X(root_e) X(log10_e)    \
   X(one_div_log10_e) X(ln_ten) X(degree) X(radian) X(sin_one) X(cos_one) X(sinh_one) X(cosh_one) X(phi) X(ln_phi)     \
   X(one_div_ln_phi) X(euler) X(one_div_euler) X(euler_sqr) X(zeta_two) X(zeta_three) X(catalan) X(glaisher)          \
   X(khinchin) X(extreme_value_skewness) X(rayleigh_skewness) X(rayleigh_kurtosis) X(rayleigh_kurtosis_excess)          \
   X(two_div_pi) X(root_two_div_pi) X(quarter_pi) X(one_div_pi) X(two_div_root_pi) X(first_feigenbaum) X(plastic)      \
   X(gauss) X(dottie) X(reciprocal_fibonacci) X(laplace_limit)

void bmpy::bind_constants(nb::module_& m) {
#define BMPY_CONSTANT(name) m.attr(#name) = bm::constants::name<double>();
   BMPY_CONSTANTS(BMPY_CONSTANT)
#undef BMPY_CONSTANT
}
