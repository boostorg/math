#  Copyright 2026 Nick Thompson
#  Distributed under the Boost Software License, Version 1.0.
#  https://www.boost.org/LICENSE_1_0.txt
"""Smoke tests for boost.math.filters against NumPy."""

import numpy as np
import pytest

import boost.math as bm
from boost.math.filters import trend_filter

rng = np.random.default_rng(12345)
t = np.linspace(0, 1, 200)
Y = np.sin(6 * t) + 0.1 * rng.normal(size=t.size)


def test_lambda_zero_and_large_lambda():
    for k in range(3):
        z = trend_filter(Y, k, 0.0)
        assert type(z) is np.ndarray and z.dtype == np.float64 and z is not Y
        np.testing.assert_allclose(z, Y, rtol=0, atol=1e-12)
        # Past lambda_max the fit is the least-squares polynomial of degree k:
        np.testing.assert_allclose(trend_filter(Y, k, 1e8), np.polyval(np.polyfit(t, Y, k), t), rtol=0, atol=1e-6)


def test_cross_validation():
    step = np.where(t < 0.5, 0.0, 1.0)
    y = step + 0.1 * rng.normal(size=t.size)
    z = trend_filter(y)
    assert np.mean((z - step) ** 2) < 0.25 * np.mean((y - step) ** 2)
    assert np.array_equal(trend_filter(list(y)), z) and np.array_equal(trend_filter(Y, k=1, lambda_=1.0), trend_filter(Y, 1, 1.0))


def test_errors():
    for bad in [lambda: trend_filter(np.zeros((4, 4))), lambda: trend_filter(Y, k=1), lambda: trend_filter(Y, lambda_=1.0),
                lambda: trend_filter(Y, 3, 1.0), lambda: trend_filter(Y, 1, -1.0), lambda: trend_filter([1.0, 2.0])]:
        with pytest.raises(ValueError):
            bad()
    with pytest.raises(ValueError, match="trend_filter"):
        trend_filter(Y, -1, 1.0)
    with bm.errstate(domain="ignore"):
        assert np.isnan(trend_filter(Y, 1, -1.0)).all()
    with bm.errstate(domain="warn"), pytest.warns(RuntimeWarning, match="domain"):
        assert np.isnan(trend_filter([1.0, 2.0, 3.0])).all()
