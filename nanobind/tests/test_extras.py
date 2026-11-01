#  Copyright 2026 Nick Thompson
#  Distributed under the Boost Software License, Version 1.0.
#  https://www.boost.org/LICENSE_1_0.txt
"""Smoke tests: the special functions of other shapes, constants, differentiation, expansions, naive_monte_carlo."""

import cmath
import math

import mpmath as mp
import numpy as np
import pytest

import boost.math as bm
from boost.math import constants, differentiation, expansions, quadrature


def test_special_functions():
    assert bm.logsumexp(np.array([[1.0, 2.0], [3.0, 800.0]])) == pytest.approx(800) and bm.logsumexp([]) == -math.inf
    np.testing.assert_allclose(bm.legendre_p_zeros(5), np.polynomial.legendre.leggauss(5)[0][2:], atol=1e-15)
    np.testing.assert_allclose(bm.hermite_zeros(6), np.polynomial.hermite.hermgauss(6)[0][3:], rtol=1e-14)
    np.testing.assert_allclose(bm.laguerre_zeros(6), np.polynomial.laguerre.laggauss(6)[0], rtol=1e-14)
    np.testing.assert_array_equal(bm.laguerre_zeros(5, 0.5), quadrature.gauss_laguerre_rule(5, 0.5)[0])
    np.testing.assert_allclose(bm.jacobi_zeros(5, -0.5, -0.5), np.cos(np.arange(9, 0, -2) * np.pi / 10), atol=1e-15)
    with pytest.raises(ValueError):
        bm.jacobi_zeros(3, 0.5, -1.0)
    assert bm.jacobi_double_prime(4, 0.5, 0.25, 0.3) == pytest.approx(float(mp.diff(lambda t: mp.jacobi(4, 0.5, 0.25, t), 0.3, 2)))
    np.testing.assert_allclose(bm.cardinal_b_spline(3, [0.0, 1.0, 2.0]), [2 / 3, 1 / 6, 0], atol=1e-15)
    assert bm.cardinal_b_spline_prime(3, 1.0) == pytest.approx(-0.5) and bm.cardinal_b_spline_double_prime(3, 1.0) == pytest.approx(1)
    assert bm.forward_cardinal_b_spline(1, 1.0) == 1.0
    with pytest.raises(ValueError):
        bm.cardinal_b_spline(21, 0.0)
    haar = (1 - cmath.exp(-1.5j)) / 1.5j / math.sqrt(2 * math.pi)
    assert abs(bm.fourier_transform_daubechies_scaling(1, 1.5) - haar) < 1e-12
    assert isinstance(bm.fourier_transform_daubechies_wavelet(3, 2.0), complex)
    ref = np.polynomial.chebyshev.chebval(0.3, [1.0, 1.0, 0.5])
    assert bm.chebyshev_clenshaw_recurrence([2.0, 1.0, 0.5], 0.3) == pytest.approx(ref)
    assert bm.chebyshev_clenshaw_recurrence([2.0, 1.0, 0.5], 0.0, 2.0, 1.3) == pytest.approx(ref)


def test_chebyshev_transform_and_legendre_stieltjes():
    t = bm.chebyshev_transform(math.exp, 0.0, 1.0)
    assert t(0.5) == pytest.approx(math.exp(0.5)) and t.prime(0.5) == pytest.approx(math.exp(0.5))
    assert t.integrate() == pytest.approx(math.e - 1) and t.coefficients.ndim == 1 and "chebyshev_transform" in repr(t)
    with pytest.raises(ValueError):
        bm.chebyshev_transform(math.exp, 0.0, 1.0, max_refinements=0)
    E = bm.legendre_stieltjes(3)
    x = np.array([0.1, 0.7])
    np.testing.assert_allclose(E(x), np.polynomial.legendre.legval(x, [0, -9 / 14, 0, 1]), rtol=1e-13)
    assert E.zeros()[0] == 0 and len(E.zeros()) == 2 and E.norm_sq() > 0 and repr(E) == "legendre_stieltjes(3)"


@pytest.mark.parametrize("p", [2, 6])
def test_daubechies(p):
    phi = bm.daubechies_scaling(p, grid_refinements=10)
    x = np.linspace(0.1, 0.9, 5)
    # Partition of unity: the integer translates sum to 1.
    np.testing.assert_allclose(sum(phi(x + k) for k in range(2 * p)), 1, atol=1e-6)
    psi = bm.daubechies_wavelet(p, grid_refinements=10)
    assert psi.support == (1.0 - p, float(p)) and repr(psi) == f"daubechies_wavelet({p}, grid_refinements=10)"
    if p == 2:
        with pytest.raises(ValueError):
            phi.prime(1.0)
    else:
        assert np.isfinite(phi.double_prime(2.5))


def test_constants():
    assert all(isinstance(getattr(constants, n), float) for n in dir(constants) if not n.startswith("_"))
    assert constants.pi == math.pi and constants.e == math.e and constants.catalan == pytest.approx(float(mp.catalan), rel=1e-15)


def test_differentiation():
    d, err = differentiation.finite_difference_derivative(math.sin, 1.0)
    assert d == pytest.approx(math.cos(1.0), rel=1e-11) and err < 1e-10
    with pytest.raises(ValueError):
        differentiation.finite_difference_derivative(math.sin, 1.0, order=3)
    assert differentiation.complex_step_derivative(cmath.exp, 1.0) == pytest.approx(math.e, rel=1e-15)
    x = np.linspace(0, 1, 100)
    L = differentiation.discrete_lanczos_derivative(x[1] - x[0], n=5)
    np.testing.assert_allclose(L(x**2), 2 * x, atol=1e-10)
    assert L(x**2, 50) == pytest.approx(2 * x[50])
    with pytest.raises(ValueError):
        L(np.ones(10))


def test_expansions():
    s = expansions.simple_continued_fraction(math.pi)
    assert s.partial_denominators()[:5] == [3, 7, 15, 1, 292] and repr(s) == str(s) and str(s).startswith("[3; 7, 15, 1, 292,")
    assert expansions.centered_continued_fraction(math.pi).partial_denominators()[:4] == [3, 7, 16, -294]
    assert expansions.luroth_expansion(math.pi).digits()[:4] == [3, 7, 1, 1]
    assert expansions.engel_expansion(math.e).digits()[:5] == [1, 1, 2, 3, 4]
    with pytest.raises(ValueError):
        expansions.engel_expansion(0.0)


def test_naive_monte_carlo():
    r = quadrature.naive_monte_carlo(lambda v: v[0] * v[1], [(0.0, 1.0), (0.0, 2.0)], 1e-2, seed=7)
    assert type(r) is quadrature.naive_monte_carlo_result and r.value == pytest.approx(1, abs=0.05) and r.error < 1e-2 < r.calls
    with pytest.raises(KeyError):
        quadrature.naive_monte_carlo(lambda v: {}[0], [(0.0, 1.0)], 1e-2)
    with pytest.raises(ValueError):
        quadrature.naive_monte_carlo(lambda v: 1.0, [(1.0, 0.0)], 1e-2)
