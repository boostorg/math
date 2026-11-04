#  Copyright 2026 Nick Thompson
#  Distributed under the Boost Software License, Version 1.0.
#  https://www.boost.org/LICENSE_1_0.txt
"""Smoke tests for boost.math.filters against NumPy."""

import numpy as np
import pytest

import boost.math as bm
from boost.math.filters import trend_filter
from boost.math import filters

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


RULES = ["wiener", "universal_soft", "sure", "eb_posterior_median"]


def jump(n):
    t = np.arange(n) / n
    f = np.sin(4 * np.pi * t) + 2 * t**2 + np.where(t > 0.6, 0.5, 0.0)
    return t, f, f + 0.1 * rng.normal(size=n)


def mse(z, f):
    return np.mean((z - f) ** 2)


def test_wavelets():
    w = filters.symlet(8)
    assert isinstance(w, filters.OrthogonalWavelet) and w.vanishing_moments == 8 and repr(w) == "symlet(8)"
    assert repr(filters.daubechies(4)) == "daubechies(4)" and filters.daubechies(10).vanishing_moments == 10


@pytest.mark.parametrize("n", [2048, 2001])
def test_wavelet_denoise_rules(n):
    _, f, y = jump(n)
    noise = mse(y, f)
    for w in [filters.daubechies(4), filters.symlet(8)]:
        for rule in RULES:
            z = filters.wavelet_denoise(y, w, rule)
            assert type(z) is np.ndarray and z.dtype == np.float64 and z.shape == (n,)
            assert mse(z, f) < noise / 4
            assert mse(filters.wavelet_denoise(y, w, rule=rule, sigma=0.1), f) < noise / 4
    np.testing.assert_array_equal(filters.wavelet_denoise(y, w), filters.wavelet_denoise(y, w, "wiener"))
    np.testing.assert_array_equal(filters.wavelet_denoise(list(y), w), filters.wavelet_denoise(y, w))
    z = filters.wavelet_denoise(y)  # wavelet and rule by cross-validation
    assert z.shape == (n,) and mse(z, f) < noise / 10


@pytest.mark.parametrize("n", [2048, 2001])
def test_denoise(n):
    t, f, y = jump(n)
    blocks = (t > 0.1) - 2.0 * (t > 0.3) + 1.5 * (t > 0.55) + 0.7 * (t > 0.8)
    for g, z in [(f, y), (blocks, blocks + 0.1 * rng.normal(size=n))]:
        d = filters.denoise(z)
        assert d.dtype == np.float64 and d.shape == (n,)
        assert np.array_equal(d, filters.wavelet_denoise(z)) or np.array_equal(d, trend_filter(z))
        assert mse(d, g) < mse(z, g) / 10


def test_wavelet_errors():
    w = filters.daubechies(4)
    _, _, y = jump(2048)
    for bad in [lambda: filters.daubechies(1), lambda: filters.symlet(11), lambda: filters.wavelet_denoise(y, w, "median"),
                lambda: filters.wavelet_denoise(y, rule="sure"), lambda: filters.wavelet_denoise(y, sigma=0.1),
                lambda: filters.wavelet_denoise(y, rule="eb_posterior_median"),
                lambda: filters.wavelet_denoise(np.zeros((4, 4)), w), lambda: filters.wavelet_denoise(np.ones(14), w),
                lambda: filters.wavelet_denoise(y, w, "sure", -1.0), lambda: filters.wavelet_denoise(np.ones(13)),
                lambda: filters.denoise(np.ones(13))]:
        with pytest.raises(ValueError):
            bad()
    with pytest.raises(ValueError, match="rule and sigma need a wavelet"):
        filters.wavelet_denoise(y, rule="eb_posterior_median")
    with pytest.raises(ValueError, match="wavelet_denoise"):
        filters.wavelet_denoise(np.ones(14), w)
    with bm.errstate(domain="ignore"):
        assert np.isnan(filters.wavelet_denoise(np.ones(14), w)).all()
    with pytest.raises(TypeError):
        filters.OrthogonalWavelet()


def test_wavelet_transform():
    x = rng.normal(size=1001)
    w = filters.symlet(4)
    c = filters.wavelet_transform(x, w)
    assert type(c) is np.ndarray and c.dtype == np.float64 and c.shape == x.shape
    np.testing.assert_allclose(np.sum(c**2), np.sum(x**2), rtol=1e-13)
    np.testing.assert_allclose(filters.inverse_wavelet_transform(c, w), x, atol=1e-14)
    np.testing.assert_array_equal(filters.wavelet_transform(x, w, levels=0), x)
    c2 = filters.wavelet_transform(x, w, levels=2)  # [a_2 | d_2 | d_1]: d_1 = c[501:], d_2 = c[251:501]
    np.testing.assert_array_equal(c2[251:], c[251:])
    np.testing.assert_allclose(filters.inverse_wavelet_transform(c2, w, levels=2), x, atol=1e-14)
    t = np.arange(1001) / 1001
    d = filters.wavelet_transform(1 + 2 * t - t**2, w)[501:]  # degree < p: no finest details
    assert np.max(np.abs(d)) < 1e-12
    assert ((1001 - 1) // (4 * 4 - 2)).bit_length() == 7  # the documented count for sym4 at n = 1001
    filters.wavelet_transform(x, w, levels=7)
    with pytest.raises(ValueError):
        filters.wavelet_transform(x, w, levels=8)
    with pytest.raises(ValueError):
        filters.wavelet_transform(np.ones(14), w)


def test_signals():
    from boost.math.filters import signals
    for g in [signals.blocks, signals.bumps, signals.heavisine, signals.doppler]:
        f = g(1000)
        assert type(f) is np.ndarray and f.dtype == np.float64 and f.shape == (1000,)
        assert abs(np.std(f) - 7) < 1e-12
        y = f + rng.normal(size=f.size)
        assert mse(filters.denoise(y), f) < mse(y, f) / 3  # bumps at n = 1000: 0.26
    assert signals.doppler().shape == (2048,) and signals.heavisine(4)[0] == 0
    with pytest.raises(ValueError):
        signals.blocks(1)
