// Copyright 2026 Nick Thompson
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/math/special_functions/legendre.hpp>
#include <boost/math/special_functions/laguerre.hpp>
#include <boost/math/special_functions/hermite.hpp>
#include <boost/math/special_functions/chebyshev.hpp>
#include <boost/math/special_functions/gegenbauer.hpp>
#include <boost/math/special_functions/jacobi.hpp>
#include <boost/math/special_functions/spherical_harmonic.hpp>

#include "bindings.hpp"
#include "macros.hpp"

void bmpy::bind_polynomials(nb::module_& m) {
   using bmpy::Int;
   using bmpy::Real;
   BMPY_I1R1(legendre_p, int, "l", "x", "The Legendre polynomial P_l(x); P_(−l−1) = P_l." BMPY_URL("sf_poly/legendre"));
   bmpy::def_vec<Int<int>, Int<int>, Real>(
       m, "legendre_p", []<class T>(int l, int mm, T x) { return bm::legendre_p(l, mm, x); }, {"l", "m", "x"},
       "The associated Legendre function P_l^m(x), with the Condon-Shortley phase (−1)^m.\n\n"
       "As scipy.special.lpmv(m, l, x): note the argument order (l, m, x)." BMPY_URL("sf_poly/legendre"));
   BMPY_I1R1(legendre_p_prime, int, "l", "x", "P_l′(x), the derivative of the Legendre polynomial." BMPY_URL("sf_poly/legendre"));
   BMPY_I1R1(legendre_q, unsigned, "l", "x", "The Legendre function of the second kind Q_l(x), |x| < 1."
      BMPY_URL("sf_poly/legendre"));
   BMPY_I1R1(laguerre, unsigned, "n", "x", "The Laguerre polynomial Lₙ(x), orthogonal for e^(−x) on [0, ∞)."
      BMPY_URL("sf_poly/laguerre"));
   bmpy::def_vec<Int<unsigned>, Int<unsigned>, Real>(
       m, "laguerre", []<class T>(unsigned n, unsigned mm, T x) { return bm::laguerre(n, mm, x); }, {"n", "m", "x"},
       "The associated Laguerre polynomial Lₙ^m(x); as scipy.special.eval_genlaguerre(n, m, x)." BMPY_URL("sf_poly/laguerre"));
   BMPY_I1R1(hermite, unsigned, "n", "x",
      "The (physicists') Hermite polynomial Hₙ(x), orthogonal for e^(−x²); as scipy.special.eval_hermite."
      BMPY_URL("sf_poly/hermite"));
   BMPY_I1R1(chebyshev_t, unsigned, "n", "x", "The Chebyshev polynomial of the first kind Tₙ(x), Tₙ(cos θ) = cos(nθ)."
      BMPY_URL("sf_poly/chebyshev"));
   BMPY_I1R1(chebyshev_u, unsigned, "n", "x",
      "The Chebyshev polynomial of the second kind Uₙ(x), Uₙ(cos θ) = sin((n+1)θ)/sin θ." BMPY_URL("sf_poly/chebyshev"));
   BMPY_I1R1(chebyshev_t_prime, unsigned, "n", "x", "Tₙ′(x), the derivative of the Chebyshev polynomial."
      BMPY_URL("sf_poly/chebyshev"));

   bmpy::def_vec<Int<unsigned>, Real, Real>(
       m, "gegenbauer", []<class T>(unsigned n, T lambda, T x) { return bm::gegenbauer(n, lambda, x); },
       {"n", "lambda_", "x"}, "The Gegenbauer (ultraspherical) polynomial Cₙ^(λ)(x), orthogonal for (1 − x²)^(λ−½); λ > −½."
       BMPY_URL("sf_poly/gegenbauer"));
   bmpy::def_vec<Int<unsigned>, Real, Real>(
       m, "gegenbauer_prime", []<class T>(unsigned n, T lambda, T x) { return bm::gegenbauer_prime(n, lambda, x); },
       {"n", "lambda_", "x"}, "d/dx Cₙ^(λ)(x)." BMPY_URL("sf_poly/gegenbauer"));
   bmpy::def_vec<Int<unsigned>, Real, Real, Int<unsigned>>(
       m, "gegenbauer_derivative",
       []<class T>(unsigned n, T lambda, T x, unsigned k) { return bm::gegenbauer_derivative(n, lambda, x, k); },
       {"n", "lambda_", "x", "k"}, "dᵏ/dxᵏ Cₙ^(λ)(x), the kth derivative." BMPY_URL("sf_poly/gegenbauer"));
   bmpy::def_vec<Int<unsigned>, Real, Real, Real>(
       m, "jacobi", []<class T>(unsigned n, T a, T b, T x) { return bm::jacobi(n, a, b, x); },
       {"n", "alpha", "beta", "x"}, "The Jacobi polynomial Pₙ^(α,β)(x), orthogonal for (1 − x)^α (1 + x)^β on [−1, 1].\n\n"
       "As scipy.special.eval_jacobi(n, alpha, beta, x)." BMPY_URL("sf_poly/jacobi"));
   bmpy::def_vec<Int<unsigned>, Real, Real, Real>(
       m, "jacobi_prime", []<class T>(unsigned n, T a, T b, T x) { return bm::jacobi_prime(n, a, b, x); },
       {"n", "alpha", "beta", "x"}, "d/dx Pₙ^(α,β)(x)." BMPY_URL("sf_poly/jacobi"));
   bmpy::def_vec<Int<unsigned>, Real, Real, Real, Int<unsigned>>(
       m, "jacobi_derivative",
       []<class T>(unsigned n, T a, T b, T x, unsigned k) { return bm::jacobi_derivative(n, a, b, x, k); },
       {"n", "alpha", "beta", "x", "k"}, "dᵏ/dxᵏ Pₙ^(α,β)(x), the kth derivative." BMPY_URL("sf_poly/jacobi"));

   bmpy::def_vec<Int<unsigned>, Int<int>, Real, Real>(
       m, "spherical_harmonic",
       []<class T>(unsigned n, int mm, T theta, T phi) { return bm::spherical_harmonic(n, mm, theta, phi); },
       {"n", "m", "theta", "phi"}, "The spherical harmonic Yₙ^m(θ, φ) ∝ P_n^m(cos θ) e^(imφ), complex.\n\n"
       "theta is the polar angle and phi the azimuth: the reverse of scipy.special.sph_harm(m, n, azimuth, polar)."
       BMPY_URL("sf_poly/sph_harm"));
   bmpy::def_vec<Int<unsigned>, Int<int>, Real, Real>(
       m, "spherical_harmonic_r",
       []<class T>(unsigned n, int mm, T theta, T phi) { return bm::spherical_harmonic_r(n, mm, theta, phi); },
       {"n", "m", "theta", "phi"}, "The real part of spherical_harmonic(n, m, theta, phi); theta polar, phi azimuthal."
      BMPY_URL("sf_poly/sph_harm"));
   bmpy::def_vec<Int<unsigned>, Int<int>, Real, Real>(
       m, "spherical_harmonic_i",
       []<class T>(unsigned n, int mm, T theta, T phi) { return bm::spherical_harmonic_i(n, mm, theta, phi); },
       {"n", "m", "theta", "phi"}, "The imaginary part of spherical_harmonic(n, m, theta, phi); theta polar, phi azimuthal."
      BMPY_URL("sf_poly/sph_harm"));
}
