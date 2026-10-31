#  Copyright 2026 Nick Thompson
#  Distributed under the Boost Software License, Version 1.0.
#  https://www.boost.org/LICENSE_1_0.txt
"""Smoke tests for boost.math.interpolators: each class once, against SciPy or data it must reproduce."""

import math
import warnings

import numpy as np
import pytest
from scipy import interpolate

import boost.math as bm
import boost.math.interpolators as bi

x = np.linspace(0, 3, 12)
y = np.sin(x)
h = float(x[1] - x[0])
inner = np.linspace(x[2], x[-3], 50)  # away from the end intervals, where the endpoint slopes differ from SciPy's


def test_against_scipy():
    assert np.allclose(bi.pchip(x, y)(inner), interpolate.PchipInterpolator(x, y)(inner), rtol=1e-12)
    assert np.allclose(bi.makima(x, y)(inner), interpolate.Akima1DInterpolator(x, y, method="makima")(inner), rtol=1e-12)
    ch, ref = bi.cubic_hermite(x, y, np.cos(x)), interpolate.CubicHermiteSpline(x, y, np.cos(x))
    assert np.allclose(ch(inner), ref(inner), rtol=1e-12) and np.allclose(ch.prime(inner), ref(inner, 1), rtol=1e-10)
    cch = bi.cardinal_cubic_hermite(y, np.cos(x), 0.0, h)
    assert np.allclose(cch(inner), ref(inner), rtol=1e-12) and cch.domain() == pytest.approx((0, 3))
    assert math.isnan(ch(math.nan)) and math.isnan(ch.prime(math.nan))


def test_reproduce_polynomials():
    p = np.polynomial.Polynomial([1, -2, 0.5, 0.25, -0.1, 0.03, 0.01, -0.002])
    t = np.linspace(0, 2.95, 37)
    # Each reproduces the polynomials up to its degree.
    for f, deg in [(bi.barycentric_rational(x, p.cutdeg(3)(x), approximation_order=3), 3),
                   (bi.quintic_hermite(x, *(p.cutdeg(5).deriv(k)(x) for k in range(3))), 5),
                   (bi.cardinal_quintic_hermite(*(p.cutdeg(5).deriv(k)(x) for k in range(3)), 0.0, h), 5),
                   (bi.septic_hermite(x, *(p.deriv(k)(x) for k in range(4))), 7),
                   (bi.cardinal_septic_hermite(*(p.deriv(k)(x) for k in range(4)), 0.0, h), 7)]:
        q = p.cutdeg(deg)
        assert np.allclose(f(t), q(t), rtol=1e-10, atol=1e-12)
        assert np.allclose(f.prime(t), q.deriv()(t), rtol=1e-8, atol=1e-10)
        if hasattr(f, "double_prime"):
            assert np.allclose(f.double_prime(t), q.deriv(2)(t), rtol=1e-6, atol=1e-8)


def test_cardinal():
    t = np.linspace(0.5, 2.5, 21)
    for f in [bi.cardinal_cubic_b_spline(y, 0.0, h, left_endpoint_derivative=1.0, right_endpoint_derivative=np.cos(3)),
              bi.cardinal_quadratic_b_spline(y, 0.0, h), bi.cardinal_quintic_b_spline(y, 0.0, h)]:
        assert np.allclose(f(t), np.sin(t), atol=2e-3) and np.allclose(f.prime(t), np.cos(t), atol=2e-2)
    # Band-limited and periodic data are reproduced:
    n, T = 16, 2.0
    s = np.arange(n) * (T / n)
    tp = bi.cardinal_trigonometric(np.cos(2 * np.pi * s / T), 0.0, T / n)
    assert np.allclose(tp(t), np.cos(2 * np.pi * t / T)) and tp.period() == pytest.approx(T)
    assert tp.integrate() == pytest.approx(0, abs=1e-12) and tp.squared_l2() == pytest.approx(T / 2)
    ws = bi.whittaker_shannon(y, 0.0, h)
    assert np.allclose(ws(x), y) and repr(ws) == "whittaker_shannon(n=12, t0=0.0, h=%r)" % h


def test_push_back_and_errors():
    m = bi.makima(x, y)
    m.push_back(3.5, np.sin(3.5))
    assert repr(m) == "makima(n=13)" and np.isfinite(m(3.4))
    for bad in [lambda: m.push_back(0.0, 1.0), lambda: m(10.0), lambda: bi.pchip(x[:3], y[:3])]:
        with pytest.raises(ValueError):
            bad()
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        r = bi.pchip(x, y)(np.array([1.0, 10.0]))
    assert np.isfinite(r[0]) and np.isnan(r[1]) and len(w) == 1
    with bm.errstate(domain="ignore"):
        assert np.isnan(bi.pchip(x, y)(10.0))
    big = np.linspace(0, 3, 50000)
    assert np.array_equal(bi.pchip(x, y)(big, threads=4), bi.pchip(x, y)(big))


def test_curves_and_bilinear():
    P = np.array([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0], [0.5, 2.0]])
    cr = bi.catmull_rom(P)
    assert repr(cr) == "catmull_rom(n=5, dim=2, closed=False, alpha=0.5)"
    assert np.allclose(cr(0.0), P[0]) and np.allclose(cr(cr.max_parameter()), P[-1]) and cr.prime(0.5).shape == (2,)
    assert cr(np.zeros((3, 4))).shape == (3, 4, 2)
    # A quadratic Bezier curve:
    Q = P[:3]
    bz, t = bi.bezier_polynomial(Q), 0.3
    assert np.allclose(bz(t), (1 - t) ** 2 * Q[0] + 2 * t * (1 - t) * Q[1] + t ** 2 * Q[2])
    assert np.allclose(bz.prime(t), 2 * (1 - t) * (Q[1] - Q[0]) + 2 * t * (Q[2] - Q[1]))
    with pytest.raises(ValueError):
        bz(1.5)
    s = np.linspace(0, 1, 8)
    vb = bi.vector_barycentric_rational(s, np.c_[s, s ** 2, 1 - s ** 3])
    assert np.allclose(vb([0.25, 0.6]), [[0.25, 0.0625, 1 - 0.25 ** 3], [0.6, 0.36, 1 - 0.6 ** 3]])
    with pytest.raises(ValueError):
        bi.vector_barycentric_rational(s[::-1], np.c_[s, s])
    xs, ys = np.arange(4) * 0.5, np.arange(3) * 2.0
    f = lambda u, v: 1 + 2 * u - v + 0.5 * u * v
    b = bi.bilinear_uniform(f(xs[None, :], ys[:, None]), dx=0.5, dy=2.0)
    u, v = np.array([0.1, 1.2, 1.5]), np.array([[0.5], [3.9]])
    assert np.allclose(b(u, v), f(u, v))
    with pytest.raises(ValueError):
        b(2.0, 0.0)
