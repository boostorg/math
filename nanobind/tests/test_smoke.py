#  Copyright 2026 Nick Thompson
#  Distributed under the Boost Software License, Version 1.0.
#  https://www.boost.org/LICENSE_1_0.txt
"""Smoke tests: each binding once against mpmath (loose 1e-12), plus the conventions."""

import inspect
import math
import warnings

import mpmath as mp
import numpy as np
import pytest

import boost.math as bm

mp.mp.dps = 30
RTOL = 1e-12


def d(f, x, n=1):
    return mp.diff(f, x, n)


def fac2(n):
    return mp.fac2(n)


def tangent(n):
    # Boost's tangent_t2n(n) is the tangent number T_{2n-1}.
    if n == 0:
        return 0
    return (-1) ** (n - 1) * 2 ** (2 * n) * (2 ** (2 * n) - 1) * mp.bernoulli(2 * n) / (2 * n)


def heuman(k, phi):
    m = k * k
    mc = 1 - m
    K, E = mp.ellipk(m), mp.ellipe(m)
    return 2 / mp.pi * (E * mp.ellipf(phi, mc) + K * mp.ellipe(phi, mc) - K * mp.ellipf(phi, mc))


def owens_t(h, a):
    return mp.quad(lambda x: mp.exp(-h * h * (1 + x * x) / 2) / (1 + x * x), [0, a]) / (2 * mp.pi)


def sph(f, n, x):
    return mp.sqrt(mp.pi / (2 * x)) * f(n + mp.mpf(1) / 2, x)


PRIMES = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]

# (python name, args, mpmath reference)
CASES = [
    # gamma.cpp
    ("tgamma", (2.5,), mp.gamma),
    ("tgamma", (2.5, 1.5), lambda a, z: mp.gammainc(a, z)),
    ("tgamma1pm1", (0.25,), lambda z: mp.gamma(1 + z) - 1),
    ("lgamma", (7.5,), mp.loggamma),
    ("digamma", (3.5,), mp.digamma),
    ("trigamma", (3.5,), lambda z: mp.psi(1, z)),
    ("polygamma", (3, 2.5), mp.psi),
    ("tgamma_ratio", (5.5, 2.25), lambda a, b: mp.gamma(a) / mp.gamma(b)),
    ("tgamma_delta_ratio", (5.5, 0.25), lambda a, dd: mp.gamma(a) / mp.gamma(a + dd)),
    ("tgamma_lower", (2.5, 1.5), lambda a, z: mp.gammainc(a, 0, z)),
    ("gamma_p", (2.5, 1.5), lambda a, z: mp.gammainc(a, 0, z, regularized=True)),
    ("gamma_q", (2.5, 1.5), lambda a, z: mp.gammainc(a, z, regularized=True)),
    ("gamma_p_inv", (2.5, 0.3), None),
    ("gamma_q_inv", (2.5, 0.3), None),
    ("gamma_p_inva", (1.5, 0.3), None),
    ("gamma_q_inva", (1.5, 0.3), None),
    ("gamma_p_derivative", (2.5, 1.5), lambda a, z: z ** (a - 1) * mp.exp(-z) / mp.gamma(a)),
    ("lgamma_p", (2.5, 1.5), lambda a, z: mp.log(mp.gammainc(a, 0, z, regularized=True))),
    ("lgamma_q", (2.5, 1.5), lambda a, z: mp.log(mp.gammainc(a, z, regularized=True))),
    ("ligamma", (2.5, 1.5), lambda a, z: mp.log(mp.gammainc(a, z))),
    ("ligamma_lower", (2.5, 1.5), lambda a, z: mp.log(mp.gammainc(a, 0, z))),
    ("factorial", (10,), mp.factorial),
    ("unchecked_factorial", (10,), mp.factorial),
    ("double_factorial", (9,), fac2),
    ("rising_factorial", (2.5, 4), mp.rf),
    ("falling_factorial", (6.5, 3), mp.ff),
    ("binomial_coefficient", (10, 4), mp.binomial),
    ("bernoulli_b2n", (3,), lambda n: mp.bernoulli(2 * n)),
    ("tangent_t2n", (4,), tangent),
    # beta.cpp
    ("beta", (2.5, 1.5), mp.beta),
    ("beta", (2.5, 1.5, 0.3), lambda a, b, x: mp.betainc(a, b, 0, x)),
    ("betac", (2.5, 1.5, 0.3), lambda a, b, x: mp.betainc(a, b, x, 1)),
    ("lbeta", (25.5, 31.5), lambda a, b: mp.log(mp.beta(a, b))),
    ("ibeta", (2.5, 1.5, 0.3), lambda a, b, x: mp.betainc(a, b, 0, x, regularized=True)),
    ("ibetac", (2.5, 1.5, 0.3), lambda a, b, x: mp.betainc(a, b, x, 1, regularized=True)),
    ("libeta", (2.5, 1.5, 0.3), lambda a, b, x: mp.log(mp.betainc(a, b, 0, x, regularized=True))),
    ("libetac", (2.5, 1.5, 0.3), lambda a, b, x: mp.log(mp.betainc(a, b, x, 1, regularized=True))),
    ("ibeta_inv", (2.5, 1.5, 0.3), None),
    ("ibetac_inv", (2.5, 1.5, 0.3), None),
    ("ibeta_inva", (1.5, 0.3, 0.4), None),
    ("ibeta_invb", (2.5, 0.3, 0.4), None),
    ("ibetac_inva", (1.5, 0.3, 0.4), None),
    ("ibetac_invb", (2.5, 0.3, 0.4), None),
    ("ibeta_derivative", (2.5, 1.5, 0.3), lambda a, b, x: x ** (a - 1) * (1 - x) ** (b - 1) / mp.beta(a, b)),
    # erf_expint_zeta.cpp
    ("erf", (0.3,), mp.erf),
    ("erfc", (2.5,), mp.erfc),
    ("erf_inv", (0.3,), mp.erfinv),
    ("erfc_inv", (0.3,), lambda q: mp.erfinv(1 - mp.mpf(q))),
    ("expint", (1.5,), mp.ei),
    ("expint", (3, 1.5), mp.expint),
    ("zeta", (2.5,), mp.zeta),
    ("owens_t", (0.5, 0.75), owens_t),
    ("lambert_w0", (1.5,), lambda z: mp.lambertw(z, 0)),
    ("lambert_wm1", (-0.2,), lambda z: mp.lambertw(z, -1)),
    ("lambert_w0_prime", (1.5,), lambda z: (lambda w: w / (z * (1 + w)))(mp.lambertw(z, 0))),
    ("lambert_wm1_prime", (-0.2,), lambda z: (lambda w: w / (z * (1 + w)))(mp.lambertw(z, -1))),
    # elementary.cpp
    ("log1p", (1e-5,), lambda x: mp.log1p(x)),
    ("log1pmx", (0.25,), lambda x: mp.log1p(x) - x),
    ("expm1", (1e-5,), mp.expm1),
    ("cbrt", (5.5,), mp.cbrt),
    ("sqrt1pm1", (1e-5,), lambda x: mp.sqrt(1 + mp.mpf(x)) - 1),
    ("powm1", (1.5, 1e-3), lambda x, y: mp.power(x, y) - 1),
    ("pow1p", (1e-3, 2.5), lambda x, y: mp.power(1 + mp.mpf(x), y)),
    ("sin_pi", (0.3,), mp.sinpi),
    ("cos_pi", (0.3,), mp.cospi),
    ("sinc_pi", (0.3,), lambda x: mp.sin(x) / x),
    ("sinhc_pi", (0.3,), lambda x: mp.sinh(x) / x),
    ("asinh", (0.3,), mp.asinh),
    ("acosh", (1.3,), mp.acosh),
    ("atanh", (0.3,), mp.atanh),
    ("hypot", (3.0, 4.5), mp.hypot),
    ("rsqrt", (2.5,), lambda x: 1 / mp.sqrt(x)),
    ("logit", (0.3,), lambda p: mp.log(p / (1 - mp.mpf(p)))),
    ("logistic_sigmoid", (0.3,), lambda x: 1 / (1 + mp.exp(-x))),
    ("logaddexp", (0.3, 2.5), lambda x, y: mp.log(mp.exp(x) + mp.exp(y))),
    ("float_next", (1.0,), lambda x: np.nextafter(1.0, 2.0)),
    ("float_prior", (1.0,), lambda x: np.nextafter(1.0, 0.0)),
    ("float_distance", (1.0, 1.0 + 2**-50), lambda a, b: 4),
    ("float_advance", (1.0, 4), lambda x, n: 1 + 4 * 2**-52),
    ("nextafter", (1.0, 2.0), lambda x, y: np.nextafter(1.0, 2.0)),
    ("ulp", (1.0,), lambda x: 2**-52),
    ("relative_difference", (1.0, 1.5), lambda a, b: 0.5),
    ("epsilon_difference", (1.0, 1.0 + 2**-51), lambda a, b: 2),
    ("fibonacci", (30,), mp.fib),
    ("prime", (9,), lambda n: PRIMES[n]),
    # bessel.cpp
    ("cyl_bessel_j", (2.5, 1.5), mp.besselj),
    ("cyl_neumann", (2.5, 1.5), mp.bessely),
    ("cyl_bessel_i", (2.5, 1.5), mp.besseli),
    ("cyl_bessel_k", (2.5, 1.5), mp.besselk),
    ("cyl_bessel_j_prime", (2.5, 1.5), lambda v, x: mp.besselj(v, x, derivative=1)),
    ("cyl_neumann_prime", (2.5, 1.5), lambda v, x: d(lambda t: mp.bessely(v, t), x)),
    ("cyl_bessel_i_prime", (2.5, 1.5), lambda v, x: mp.besseli(v, x, derivative=1)),
    ("cyl_bessel_k_prime", (2.5, 1.5), lambda v, x: d(lambda t: mp.besselk(v, t), x)),
    ("sph_bessel", (3, 1.5), lambda n, x: sph(mp.besselj, n, x)),
    ("sph_neumann", (3, 1.5), lambda n, x: sph(mp.bessely, n, x)),
    ("sph_bessel_prime", (3, 1.5), lambda n, x: d(lambda t: sph(mp.besselj, n, t), x)),
    ("sph_neumann_prime", (3, 1.5), lambda n, x: d(lambda t: sph(mp.bessely, n, t), x)),
    ("cyl_bessel_j_zero", (2.5, 3), mp.besseljzero),
    ("cyl_neumann_zero", (2.5, 3), mp.besselyzero),
    ("cyl_hankel_1", (2.5, 1.5), mp.hankel1),
    ("cyl_hankel_2", (2.5, 1.5), mp.hankel2),
    ("sph_hankel_1", (2.5, 1.5), lambda v, x: sph(mp.hankel1, v, x)),
    ("sph_hankel_2", (2.5, 1.5), lambda v, x: sph(mp.hankel2, v, x)),
    ("airy_ai", (1.5,), mp.airyai),
    ("airy_bi", (1.5,), mp.airybi),
    ("airy_ai_prime", (1.5,), lambda x: mp.airyai(x, derivative=1)),
    ("airy_bi_prime", (1.5,), lambda x: mp.airybi(x, derivative=1)),
    ("airy_ai_zero", (3,), mp.airyaizero),
    ("airy_bi_zero", (3,), mp.airybizero),
    # elliptic.cpp (Boost takes the modulus k, mpmath the parameter m = k^2)
    ("ellint_1", (0.5,), lambda k: mp.ellipk(k * k)),
    ("ellint_1", (0.5, 0.7), lambda k, phi: mp.ellipf(phi, k * k)),
    ("ellint_2", (0.5,), lambda k: mp.ellipe(k * k)),
    ("ellint_2", (0.5, 0.7), lambda k, phi: mp.ellipe(phi, k * k)),
    ("ellint_3", (0.5, 0.25), lambda k, n: mp.ellippi(n, k * k)),
    ("ellint_3", (0.5, 0.25, 0.7), lambda k, n, phi: mp.ellippi(n, phi, k * k)),
    ("ellint_d", (0.5,), lambda k: (mp.ellipk(k * k) - mp.ellipe(k * k)) / (k * k)),
    ("ellint_d", (0.5, 0.7), lambda k, phi: (mp.ellipf(phi, k * k) - mp.ellipe(phi, k * k)) / (k * k)),
    ("ellint_rc", (1.5, 2.5), mp.elliprc),
    ("ellint_rd", (1.5, 2.5, 3.5), mp.elliprd),
    ("ellint_rf", (1.5, 2.5, 3.5), mp.elliprf),
    ("ellint_rg", (1.5, 2.5, 3.5), mp.elliprg),
    ("ellint_rj", (1.5, 2.5, 3.5, 4.5), mp.elliprj),
    ("jacobi_zeta", (0.5, 0.7), lambda k, phi: mp.ellipe(phi, k * k) - mp.ellipe(k * k) / mp.ellipk(k * k) * mp.ellipf(phi, k * k)),
    ("heuman_lambda", (0.5, 0.7), heuman),
] + [
    ("jacobi_" + f, (0.5, 0.7), (lambda f: lambda k, u: mp.ellipfun(f, u, m=k * k))(f))
    for f in ("sn", "cn", "dn", "cd", "dc", "ns", "sd", "ds", "nc", "nd", "sc", "cs")
] + [
    # jacobi_theta.cpp
    ("jacobi_theta%d" % n, (0.3, 0.2), (lambda n: lambda z, q: mp.jtheta(n, z, q))(n)) for n in (1, 2, 3, 4)
] + [
    ("jacobi_theta%dtau" % n, (0.3, 1.2), (lambda n: lambda z, t: mp.jtheta(n, z, mp.exp(-mp.pi * t)))(n)) for n in (1, 2, 3, 4)
] + [
    ("jacobi_theta3m1", (0.3, 0.2), lambda z, q: mp.jtheta(3, z, q) - 1),
    ("jacobi_theta4m1", (0.3, 0.2), lambda z, q: mp.jtheta(4, z, q) - 1),
    ("jacobi_theta3m1tau", (0.3, 1.2), lambda z, t: mp.jtheta(3, z, mp.exp(-mp.pi * t)) - 1),
    ("jacobi_theta4m1tau", (0.3, 1.2), lambda z, t: mp.jtheta(4, z, mp.exp(-mp.pi * t)) - 1),
    # polynomials.cpp
    ("legendre_p", (5, 0.3), mp.legendre),
    ("legendre_p", (5, 2, 0.3), mp.legenp),
    ("legendre_p_prime", (5, 0.3), lambda l, x: d(lambda t: mp.legendre(l, t), x)),
    ("legendre_q", (5, 0.3), lambda l, x: mp.legenq(l, 0, x)),
    ("laguerre", (5, 0.3), lambda n, x: mp.laguerre(n, 0, x)),
    ("laguerre", (5, 2, 0.3), mp.laguerre),
    ("hermite", (5, 0.3), mp.hermite),
    ("chebyshev_t", (5, 0.3), mp.chebyt),
    ("chebyshev_u", (5, 0.3), mp.chebyu),
    ("chebyshev_t_prime", (5, 0.3), lambda n, x: d(lambda t: mp.chebyt(n, t), x)),
    ("gegenbauer", (5, 1.5, 0.3), mp.gegenbauer),
    ("gegenbauer_prime", (5, 1.5, 0.3), lambda n, a, x: d(lambda t: mp.gegenbauer(n, a, t), x)),
    ("gegenbauer_derivative", (5, 1.5, 0.3, 2), lambda n, a, x, k: d(lambda t: mp.gegenbauer(n, a, t), x, k)),
    ("jacobi", (5, 1.5, 0.5, 0.3), mp.jacobi),
    ("jacobi_prime", (5, 1.5, 0.5, 0.3), lambda n, a, b, x: d(lambda t: mp.jacobi(n, a, b, t), x)),
    ("jacobi_derivative", (5, 1.5, 0.5, 0.3, 2), lambda n, a, b, x, k: d(lambda t: mp.jacobi(n, a, b, t), x, k)),
    ("spherical_harmonic", (3, 1, 0.3, 0.4), mp.spherharm),
    ("spherical_harmonic_r", (3, 1, 0.3, 0.4), lambda n, m, t, p: mp.re(mp.spherharm(n, m, t, p))),
    ("spherical_harmonic_i", (3, 1, 0.3, 0.4), lambda n, m, t, p: mp.im(mp.spherharm(n, m, t, p))),
    # hypergeometric.cpp
    ("hypergeometric_1F0", (2.5, 0.3), lambda a, z: mp.hyper([a], [], z)),
    ("hypergeometric_0F1", (2.5, 0.3), mp.hyp0f1),
    ("hypergeometric_2F0", (-3.0, 2.5, 0.3), mp.hyp2f0),
    ("hypergeometric_1F1", (2.5, 1.5, 0.3), mp.hyp1f1),
    ("hypergeometric_1F1_regularized", (2.5, 1.5, 0.3), lambda a, b, z: mp.hyp1f1(a, b, z) / mp.gamma(b)),
    ("log_hypergeometric_1F1", (2.5, 1.5, 0.3), lambda a, b, z: mp.log(mp.hyp1f1(a, b, z))),
    ("hypergeometric_pFq", ([1.5, 2.5], [3.5], 0.3), lambda a, b, z: mp.hyper(a, b, z)),
]

# Inverses are checked by round trip: (inverse, forward, which argument of the
# forward function the inverse returns, which argument it consumes as target).
INVERSES = {
    "gamma_p_inv": lambda a, p, x: bm.gamma_p(a, x) - p,
    "gamma_q_inv": lambda a, q, x: bm.gamma_q(a, x) - q,
    "gamma_p_inva": lambda z, p, a: bm.gamma_p(a, z) - p,
    "gamma_q_inva": lambda z, q, a: bm.gamma_q(a, z) - q,
    "ibeta_inv": lambda a, b, p, x: bm.ibeta(a, b, x) - p,
    "ibetac_inv": lambda a, b, q, x: bm.ibetac(a, b, x) - q,
    "ibeta_inva": lambda b, x, p, a: bm.ibeta(a, b, x) - p,
    "ibeta_invb": lambda a, x, p, b: bm.ibeta(a, b, x) - p,
    "ibetac_inva": lambda b, x, q, a: bm.ibetac(a, b, x) - q,
    "ibetac_invb": lambda a, x, q, b: bm.ibetac(a, b, x) - q,
}


def close(got, want, rtol=RTOL):
    want = complex(want) if isinstance(want, mp.mpc) or isinstance(want, complex) else float(want)
    return abs(got - want) <= rtol * max(abs(want), 1e-300)


def case_id(c):
    return "%s/%d" % (c[0], len(c[1]))


@pytest.mark.parametrize("case", CASES, ids=case_id)
def test_scalar_value(case):
    name, args, ref = case
    got = getattr(bm, name)(*args)
    if ref is None:
        assert abs(INVERSES[name](*args, got)) < 1e-12
    else:
        want = ref(*[mp.mpf(a) if isinstance(a, float) else a for a in args])
        assert close(got, want), (got, want)


@pytest.mark.parametrize("case", [c for c in CASES if c[0] != "hypergeometric_pFq"], ids=case_id)
def test_array_matches_scalar(case):
    """The ndarray overload broadcasts and agrees with the scalar overload."""
    name, args, _ = case
    f = getattr(bm, name)
    # Turn the last argument into a (2, 3) array of identical values.
    last = args[-1]
    arr = np.full((2, 3), last)
    got = f(*args[:-1], arr)
    assert got.shape == (2, 3)
    want = f(*args)
    assert np.all(got == want)


# Not elementwise, so tested on their own below.
SIEVE = {"prime_sieve", "prime_range", "prime_count", "polynomial_roots", "polynomial_roots_aberth", "polynomial_roots_companion_matrix"}
QUADRATURE = {"tanh_sinh", "exp_sinh", "sinh_sinh", "trapezoidal", "gauss", "gauss_kronrod", "gauss_hermite",
              "ooura_fourier_sin", "ooura_fourier_cos", "van_den_bos_unit_square"}
OPTIMIZATION = {"differential_evolution", "jso", "random_search", "cma_es"}


def test_every_function_has_a_case():
    names = {n for n in bm.__all__ if not inspect.isclass(getattr(bm, n))}
    covered = {c[0] for c in CASES} | SIEVE | QUADRATURE | OPTIMIZATION
    assert names - covered == set()


def test_optimization():
    def f(x):
        return (x[0] - 1) ** 2 + (x[1] + 2) ** 2
    lo, hi = [-5.0, -5.0], [5.0, 5.0]
    runs = [bm.differential_evolution(f, lo, hi, population=50, max_generations=200, threads=4), bm.jso(f, lo, hi),
            bm.random_search(f, lo, hi, max_function_calls=20000)]
    if hasattr(bm, "cma_es"):  # Only where Eigen was found.
        runs.append(bm.cma_es(f, lo, hi, max_generations=200))
    for x in runs:
        assert f(x) < 1e-2
    def bad(x):
        raise RuntimeError("boom")
    with pytest.raises(RuntimeError):
        bm.differential_evolution(bad, lo, hi, threads=4)
    with pytest.raises(ValueError):
        bm.jso(f, [0.0], [1.0, 2.0])


def test_quadrature():
    cases = [
        (bm.tanh_sinh(lambda x: 1 / math.sqrt(x), 0, 1), 2.0),
        (bm.exp_sinh(lambda x: math.exp(-x)), 1.0),
        (bm.sinh_sinh(lambda x: math.exp(-x * x)), math.sqrt(math.pi)),
        (bm.trapezoidal(lambda x: 1 / (5 - 4 * math.cos(x)), 0, 2 * math.pi), 2 * math.pi / 3),
        (bm.gauss(lambda x: x * x, 0, 1, points=7), 1 / 3),
        (bm.gauss_kronrod(math.sin, 0, math.pi, points=61), 2.0),
        (bm.gauss_hermite(lambda x: x * x), math.sqrt(math.pi) / 2),
        (bm.ooura_fourier_sin(lambda t: 1 / t, 1.0), math.pi / 2),
        (bm.ooura_fourier_cos(lambda t: 1 / (1 + t * t), 1.0), math.pi / (2 * math.e)),
        (bm.van_den_bos_unit_square(lambda x, y: x * y), 0.25),
    ]
    for result, exact in cases:
        assert math.isclose(result[0], exact, rel_tol=1e-12)
    with pytest.raises(ValueError):
        bm.gauss(math.sin, 0, 1, points=8)
    with pytest.raises(ZeroDivisionError):
        bm.tanh_sinh(lambda x: 1 / 0, 0, 1)


def test_sieve():
    want = [p for p in range(2, 100) if all(p % d for d in range(2, p))]
    got = bm.prime_sieve(100)
    assert got.dtype == np.uint64 and list(got) == want
    assert list(bm.prime_range(50, 100)) == [p for p in want if p >= 50]
    assert bm.prime_count(100) == len(want) and bm.prime_count(50, 100) == len([p for p in want if p >= 50])
    assert bm.prime_count(10**7) == 664579
    assert list(bm.prime_range(2**63 - 100, 2**63)) == [9223372036854775783]


def test_integer_results():
    assert type(bm.prime(3)) is int and bm.prime(3) == 7
    assert bm.prime(np.arange(4)).dtype == np.uint64
    assert bm.fibonacci(93) == 12200160415121876738 and bm.fibonacci(np.arange(3)).dtype == np.uint64
    with pytest.raises(ValueError):
        bm.prime(np.array([1, 20000]))


def test_python_scalars_give_python_floats():
    assert type(bm.libeta(1, 2, 0.1)) is float
    assert type(bm.tgamma(np.float32(3))) is np.float32
    assert bm.libeta(1, 2, 0) == -np.inf


# --- conventions ------------------------------------------------------------


def test_float32_is_evaluated_in_float():
    x = np.linspace(0.5, 3, 7, dtype=np.float32)
    y = bm.tgamma(x)
    assert y.dtype == np.float32
    assert isinstance(bm.tgamma(np.float32(0.5)), np.float32)
    # Python scalars do not upcast float32 arrays (NEP 50).
    assert bm.gamma_p(x, 1.5).dtype == np.float32
    assert bm.gamma_p(x, np.float64(1.5)).dtype == np.float64
    # float_distance is counted in float ulps.
    assert bm.float_distance(np.float32(1), np.nextafter(np.float32(1), np.float32(2))) == 1


def test_int_and_list_inputs():
    assert bm.tgamma(3) == 2.0
    assert np.allclose(bm.tgamma([1, 2, 3]), [1, 1, 2])
    assert bm.legendre_p(np.arange(3), 0.5).tolist() == [1.0, 0.5, -0.125]


def test_broadcasting_and_strides():
    a = np.array([[1.0], [2.0]])
    b = np.array([1.0, 2.0, 3.0])
    got = bm.ibeta(a, b, 0.5)
    assert got.shape == (2, 3)
    assert np.allclose(got, [[0.5, 0.75, 0.875], [0.25, 0.5, 0.6875]])
    base = np.arange(12.0).reshape(3, 4) + 1
    view = base[::2, ::-2]
    assert np.array_equal(bm.lgamma(view), bm.lgamma(np.ascontiguousarray(view)))
    with pytest.raises(ValueError):
        bm.beta(np.ones(2), np.ones(3))


def test_complex_results():
    got = bm.cyl_hankel_1(0.5, np.array([1.0, 2.0]))
    assert got.dtype == np.complex128
    assert close(got[1], mp.hankel1(0.5, 2))


def test_scalar_errors_raise():
    with pytest.raises(ValueError):
        bm.tgamma(-1.0)
    with pytest.raises(OverflowError):
        bm.tgamma(1e10)
    with pytest.raises(ValueError):
        bm.tgamma(np.array(-1.0))  # 0-d arrays behave like scalars
    with pytest.raises(TypeError):
        bm.legendre_p(1.5, 0.5)
    with pytest.raises(ValueError):
        bm.sph_bessel(-1, 0.5)


def test_array_errors_warn_and_give_nan():
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        got = bm.tgamma(np.array([-1.0, 1.5, -2.0]))
    assert np.isnan(got[0]) and np.isnan(got[2]) and close(got[1], mp.gamma(1.5))
    assert len(w) == 1 and issubclass(w[0].category, RuntimeWarning)
    assert "2 of 3" in str(w[0].message)


# --- distributions ----------------------------------------------------------

CONTINUOUS = [
    ("arcsine_distribution", ()),
    ("beta_distribution", (2.5, 1.5)),
    ("cauchy_distribution", ()),
    ("chi_squared_distribution", (3.5,)),
    ("exponential_distribution", (2.0,)),
    ("extreme_value_distribution", ()),
    ("generalized_extreme_value_distribution", (0.5, 2.0, 0.1)),
    ("fisher_f_distribution", (3.0, 7.0)),
    ("gamma_distribution", (2.5, 1.5)),
    ("inverse_chi_squared_distribution", (4.0,)),
    ("inverse_gamma_distribution", (3.0, 2.0)),
    ("inverse_gaussian_distribution", (1.5, 2.0)),
    ("laplace_distribution", ()),
    ("logistic_distribution", ()),
    ("lognormal_distribution", ()),
    ("normal_distribution", (1.0, 2.0)),
    ("pareto_distribution", (1.0, 3.0)),
    ("rayleigh_distribution", (1.5,)),
    ("skew_normal_distribution", (0.0, 1.0, 2.0)),
    ("students_t_distribution", (4.5,)),
    ("triangular_distribution", ()),
    ("uniform_distribution", (1.0, 3.0)),
    ("weibull_distribution", (2.0, 1.5)),
    ("holtsmark_distribution", ()),
    ("landau_distribution", ()),
    ("mapairy_distribution", ()),
    ("saspoint5_distribution", ()),
    ("von_mises_distribution", (0.5, 2.0)),
    ("kolmogorov_smirnov_distribution", (10.0,)),
    ("hyperexponential_distribution", ([0.3, 0.7], [1.0, 3.0])),
    ("non_central_beta_distribution", (2.5, 1.5, 1.0)),
    ("non_central_chi_squared_distribution", (3.0, 1.5)),
    ("non_central_f_distribution", (3.0, 7.0, 1.5)),
    ("non_central_t_distribution", (5.0, 1.0)),
]

DISCRETE = [
    ("bernoulli_distribution", (0.3,)),
    ("binomial_distribution", (20.0, 0.3)),
    ("geometric_distribution", (0.3,)),
    ("hypergeometric_distribution", (10, 20, 50)),
    ("negative_binomial_distribution", (5.0, 0.3)),
    ("poisson_distribution", (4.5,)),
]


# Classes that are not distributions.
TRANSFORMS = {"fft", "real_fft", "dct_ii"}


def test_every_distribution_has_a_case():
    names = {n for n in bm.__all__ if inspect.isclass(getattr(bm, n))}
    assert names - TRANSFORMS == {c[0] for c in CONTINUOUS + DISCRETE}
    assert TRANSFORMS <= names


@pytest.mark.parametrize("case", CONTINUOUS, ids=lambda c: c[0])
def test_continuous_distribution(case):
    name, params = case
    dist = getattr(bm, name)(*params)
    p = np.array([0.1, 0.5, 0.9])
    x = dist.quantile(p)
    assert np.allclose(dist.cdf(x), p, rtol=1e-9)
    assert np.allclose(dist.quantile_complement(1 - p), x, rtol=1e-9)
    assert np.allclose(dist.cdf(x) + dist.cdf_complement(x), 1, rtol=1e-12)
    pdf = dist.pdf(x)
    assert np.all(pdf > 0)
    assert np.allclose(dist.logpdf(x), np.log(pdf), rtol=1e-9)
    assert np.allclose(dist.logcdf(x), np.log(p), rtol=1e-9)
    assert np.allclose(dist.logcdf_complement(x), np.log1p(-p), rtol=1e-9)
    assert np.allclose(dist.hazard(x), pdf / (1 - p), rtol=1e-9)
    assert np.allclose(dist.chf(x), -np.log1p(-p), rtol=1e-9)
    lo, hi = dist.support()
    assert lo <= x[0] and x[-1] <= hi
    assert math.isfinite(dist.median())


@pytest.mark.parametrize("case", DISCRETE, ids=lambda c: c[0])
def test_discrete_distribution(case):
    name, params = case
    dist = getattr(bm, name)(*params)
    k = np.arange(0.0, 2.0 if name == "bernoulli_distribution" else 4.0)
    assert np.allclose(np.cumsum(dist.pdf(k)), dist.cdf(k), rtol=1e-12)
    assert np.allclose(dist.cdf(k) + dist.cdf_complement(k), 1, rtol=1e-12)
    q = dist.quantile(0.5)
    assert dist.cdf(q) >= 0.5 - 1e-12
    assert math.isfinite(dist.mean()) and math.isfinite(dist.variance())


@pytest.mark.parametrize("case", CONTINUOUS + DISCRETE, ids=lambda c: c[0])
def test_distribution_repr(case):
    # The repr names the constructor's arguments, each a property, and rebuilds the distribution.
    name, params = case
    dist = getattr(bm, name)(*params)
    r = repr(dist)
    assert r.startswith(name + "(")
    again = eval(r, vars(bm))
    assert repr(again) == r
    assert np.allclose(again.cdf(dist.quantile(0.3)), 0.3, rtol=1e-12) or name in {c[0] for c in DISCRETE}
    for i, value in enumerate(params):
        arg = r[len(name) + 1:-1].split(", ")[i].split("=")[0] if name != "hyperexponential_distribution" else None
        if arg is not None and not isinstance(value, (list, tuple)):
            got = getattr(dist, arg)
            assert (got() if callable(got) else got) == value


def test_normal_against_mpmath():
    dist = bm.normal_distribution(1.0, 2.0)
    for x in (-3.0, 0.5, 4.0):
        assert close(dist.cdf(x), mp.ncdf(x, 1, 2))
        assert close(dist.pdf(x), mp.npdf(x, 1, 2))
    assert dist.mean() == 1.0 and dist.standard_deviation() == 2.0


def test_undefined_moment_raises():
    with pytest.raises(ValueError):
        bm.cauchy_distribution().mean()
