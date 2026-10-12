// Copyright 2026 Nick Thompson
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_MATH_NANOBIND_BINDINGS_HPP
#define BOOST_MATH_NANOBIND_BINDINGS_HPP

#include "vectorize.hpp"

namespace bmpy {

void bind_gamma(nb::module_& m);
void bind_beta(nb::module_& m);
void bind_erf_expint_zeta(nb::module_& m);
void bind_elementary(nb::module_& m);
void bind_bessel(nb::module_& m);
void bind_elliptic(nb::module_& m);
void bind_jacobi_theta(nb::module_& m);
void bind_polynomials(nb::module_& m);
void bind_hypergeometric(nb::module_& m);
void bind_distributions_continuous(nb::module_& m);
void bind_distributions_heavy(nb::module_& m);
void bind_distributions_discrete(nb::module_& m);
void bind_distributions_noncentral(nb::module_& m);
void bind_fft(nb::module_& m);

}  // namespace bmpy

#endif
