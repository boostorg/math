// Copyright 2026 Nick Thompson
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include "bindings.hpp"

NB_MODULE(_ext, m) {
   m.doc() = "Development-only nanobind bindings for Boost.Math, using the default policy.";
   bmpy::bind_gamma(m);
   bmpy::bind_beta(m);
   bmpy::bind_erf_expint_zeta(m);
   bmpy::bind_elementary(m);
   bmpy::bind_bessel(m);
   bmpy::bind_elliptic(m);
   bmpy::bind_jacobi_theta(m);
   bmpy::bind_polynomials(m);
   bmpy::bind_hypergeometric(m);
   bmpy::bind_distributions_continuous(m);
   bmpy::bind_distributions_heavy(m);
   bmpy::bind_distributions_discrete(m);
   bmpy::bind_distributions_noncentral(m);
   bmpy::bind_fft(m);
}
