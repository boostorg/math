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
   BMPY_R1(tgamma, "z", "The gamma function Γ(z)." BMPY_URL("sf_gamma/tgamma"));
   BMPY_R2(tgamma, "a", "z", "The upper incomplete gamma function Γ(a, z), not normalised (gamma_q is).");
   BMPY_R1(tgamma1pm1, "z", "Γ(1 + z) − 1, accurate for small z." BMPY_URL("sf_gamma/tgamma"));
   BMPY_R1(lgamma, "z", "log |Γ(z)|; as scipy.special.gammaln." BMPY_URL("sf_gamma/lgamma"));
   BMPY_R1(digamma, "z", "The digamma function ψ(z) = Γ′(z)/Γ(z)." BMPY_URL("sf_gamma/digamma"));
   BMPY_R1(trigamma, "z", "The trigamma function ψ′(z) = d²/dz² log Γ(z)." BMPY_URL("sf_gamma/trigamma"));
   BMPY_I1R1(polygamma, int, "n", "z", "The polygamma function ψ⁽ⁿ⁾(z) = dⁿ⁺¹/dzⁿ⁺¹ log Γ(z); n = 0 is digamma."
      BMPY_URL("sf_gamma/polygamma"));
   BMPY_R2(tgamma_ratio, "a", "b", "Γ(a)/Γ(b), without intermediate overflow." BMPY_URL("sf_gamma/gamma_ratios"));
   BMPY_R2(tgamma_delta_ratio, "a", "delta", "Γ(a)/Γ(a + delta), accurate when delta is small relative to a."
      BMPY_URL("sf_gamma/gamma_ratios"));
   BMPY_R2(tgamma_lower, "a", "z",
      "The lower incomplete gamma function γ(a, z), not normalised (gamma_p is)."
      BMPY_URL("sf_gamma/igamma"));
   BMPY_R2(gamma_p, "a", "z",
      "The regularised lower incomplete gamma function P(a, z) = γ(a, z)/Γ(a); as scipy.special.gammainc(a, z)."
      BMPY_URL("sf_gamma/igamma"));
   BMPY_R2(gamma_q, "a", "z",
      "The regularised upper incomplete gamma function Q(a, z) = Γ(a, z)/Γ(a) = 1 − P(a, z); as scipy.special.gammaincc(a, z)."
      BMPY_URL("sf_gamma/igamma"));
   BMPY_R2(gamma_p_inv, "a", "p", "The z with P(a, z) = p; as scipy.special.gammaincinv(a, p)." BMPY_URL("sf_gamma/igamma_inv"));
   BMPY_R2(gamma_q_inv, "a", "q", "The z with Q(a, z) = q; as scipy.special.gammainccinv(a, q)." BMPY_URL("sf_gamma/igamma_inv"));
   BMPY_R2(gamma_p_inva, "z", "p", "The a with P(a, z) = p; note the argument order (z, p)." BMPY_URL("sf_gamma/igamma_inv"));
   BMPY_R2(gamma_q_inva, "z", "q", "The a with Q(a, z) = q; note the argument order (z, q)." BMPY_URL("sf_gamma/igamma_inv"));
   BMPY_R2(gamma_p_derivative, "a", "z",
      "∂P(a, z)/∂z = z^(a−1) e^(−z)/Γ(a), the density of the gamma distribution." BMPY_URL("sf_gamma/gamma_derivatives"));
   BMPY_R2(lgamma_p, "a", "z", "log P(a, z), finite where P(a, z) underflows." BMPY_URL("sf_gamma/igamma"));
   BMPY_R2(lgamma_q, "a", "z", "log Q(a, z), finite where Q(a, z) underflows." BMPY_URL("sf_gamma/igamma"));
   BMPY_R2(ligamma, "a", "z", "log Γ(a, z), finite where Γ(a, z) underflows." BMPY_URL("sf_gamma/igamma"));
   BMPY_R2(ligamma_lower, "a", "z", "log γ(a, z), finite where γ(a, z) underflows." BMPY_URL("sf_gamma/igamma"));

   BMPY_T_I1(factorial, unsigned, "n", "n! as a float; n > 170 overflows." BMPY_URL("factorials/sf_factorial"));
   BMPY_T_I1(double_factorial, unsigned, "n", "The double factorial n!!, as a float." BMPY_URL("factorials/sf_double_factorial"));
   BMPY_R1I1(rising_factorial, int, "x", "n",
      "The rising factorial (Pochhammer symbol) Γ(x + n)/Γ(x); as scipy.special.poch(x, n)."
      BMPY_URL("factorials/sf_rising_factorial"));
   BMPY_R1I1(falling_factorial, unsigned, "x", "n", "The falling factorial x (x − 1) ... (x − n + 1)."
      BMPY_URL("factorials/sf_falling_factorial"));
   bmpy::def_vec<bmpy::Int<unsigned>, bmpy::Int<unsigned>>(
       m, "binomial_coefficient", []<class T>(unsigned n, unsigned k) { return bm::binomial_coefficient<double>(n, k); },
       {"n", "k"}, "The binomial coefficient n!/(k! (n − k)!), as a float; k > n is a domain error."
      BMPY_URL("factorials/sf_binomial"));
   BMPY_T_I1(bernoulli_b2n, int, "n", "The Bernoulli number B₂ₙ: note the index is 2n, so bernoulli_b2n(1) = B₂ = 1/6."
      BMPY_URL("number_series/bernoulli_numbers"));
   BMPY_T_I1(tangent_t2n, int, "n", "The tangent number T₂ₙ; tangent_t2n(1) = 1." BMPY_URL("number_series/tangent_numbers"));
}
