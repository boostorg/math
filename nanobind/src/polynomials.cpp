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
   BMPY_I1R1(legendre_p, int, "l", "x");
   bmpy::def_vec<Int<int>, Int<int>, Real>(
       m, "legendre_p", []<class T>(int l, int mm, T x) { return bm::legendre_p(l, mm, x); }, {"l", "m", "x"},
       "boost::math::legendre_p (associated)");
   BMPY_I1R1(legendre_p_prime, int, "l", "x");
   BMPY_I1R1(legendre_q, unsigned, "l", "x");
   BMPY_I1R1(laguerre, unsigned, "n", "x");
   bmpy::def_vec<Int<unsigned>, Int<unsigned>, Real>(
       m, "laguerre", []<class T>(unsigned n, unsigned mm, T x) { return bm::laguerre(n, mm, x); }, {"n", "m", "x"},
       "boost::math::laguerre (associated)");
   BMPY_I1R1(hermite, unsigned, "n", "x");
   BMPY_I1R1(chebyshev_t, unsigned, "n", "x");
   BMPY_I1R1(chebyshev_u, unsigned, "n", "x");
   BMPY_I1R1(chebyshev_t_prime, unsigned, "n", "x");

   bmpy::def_vec<Int<unsigned>, Real, Real>(
       m, "gegenbauer", []<class T>(unsigned n, T lambda, T x) { return bm::gegenbauer(n, lambda, x); },
       {"n", "lambda_", "x"}, "boost::math::gegenbauer");
   bmpy::def_vec<Int<unsigned>, Real, Real>(
       m, "gegenbauer_prime", []<class T>(unsigned n, T lambda, T x) { return bm::gegenbauer_prime(n, lambda, x); },
       {"n", "lambda_", "x"}, "boost::math::gegenbauer_prime");
   bmpy::def_vec<Int<unsigned>, Real, Real, Int<unsigned>>(
       m, "gegenbauer_derivative",
       []<class T>(unsigned n, T lambda, T x, unsigned k) { return bm::gegenbauer_derivative(n, lambda, x, k); },
       {"n", "lambda_", "x", "k"}, "boost::math::gegenbauer_derivative");
   bmpy::def_vec<Int<unsigned>, Real, Real, Real>(
       m, "jacobi", []<class T>(unsigned n, T a, T b, T x) { return bm::jacobi(n, a, b, x); },
       {"n", "alpha", "beta", "x"}, "boost::math::jacobi");
   bmpy::def_vec<Int<unsigned>, Real, Real, Real>(
       m, "jacobi_prime", []<class T>(unsigned n, T a, T b, T x) { return bm::jacobi_prime(n, a, b, x); },
       {"n", "alpha", "beta", "x"}, "boost::math::jacobi_prime");
   bmpy::def_vec<Int<unsigned>, Real, Real, Real, Int<unsigned>>(
       m, "jacobi_derivative",
       []<class T>(unsigned n, T a, T b, T x, unsigned k) { return bm::jacobi_derivative(n, a, b, x, k); },
       {"n", "alpha", "beta", "x", "k"}, "boost::math::jacobi_derivative");

   bmpy::def_vec<Int<unsigned>, Int<int>, Real, Real>(
       m, "spherical_harmonic",
       []<class T>(unsigned n, int mm, T theta, T phi) { return bm::spherical_harmonic(n, mm, theta, phi); },
       {"n", "m", "theta", "phi"}, "boost::math::spherical_harmonic (complex)");
   bmpy::def_vec<Int<unsigned>, Int<int>, Real, Real>(
       m, "spherical_harmonic_r",
       []<class T>(unsigned n, int mm, T theta, T phi) { return bm::spherical_harmonic_r(n, mm, theta, phi); },
       {"n", "m", "theta", "phi"}, "boost::math::spherical_harmonic_r");
   bmpy::def_vec<Int<unsigned>, Int<int>, Real, Real>(
       m, "spherical_harmonic_i",
       []<class T>(unsigned n, int mm, T theta, T phi) { return bm::spherical_harmonic_i(n, mm, theta, phi); },
       {"n", "m", "theta", "phi"}, "boost::math::spherical_harmonic_i");
}
