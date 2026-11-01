#  Copyright 2026 Nick Thompson
#  Distributed under the Boost Software License, Version 1.0.
#  https://www.boost.org/LICENSE_1_0.txt
"""The optimizers' OptimizeResult, the gradient optimizers, brent, and the quadrature result types."""

import math

import numpy as np
import pytest

import boost.math as bm

o = bm.optimize


def quadratic(x):
    return (x[0] - 1) ** 2 + (x[1] + 2) ** 2


def quadratic_gradient(x):
    return [2 * (x[0] - 1), 2 * (x[1] + 2)]


def test_derivative_free_results():
    lo, hi = [-5.0, -5.0], [5.0, 5.0]
    for r in [o.random_search(quadratic, lo, hi, max_function_calls=500, threads=2), o.turbo(quadratic, lo, hi, max_function_calls=40),
              o.dycors(quadratic, lo, hi, max_function_calls=40)]:
        assert type(r) is o.OptimizeResult and isinstance(r.x, np.ndarray)
        assert r.fun == quadratic(r.x) and r.nfev > 0 and r.nit is None and r.converged is None
    with pytest.raises(TypeError):
        o.turbo(f=quadratic, lower_bounds=lo, upper_bounds=hi)


@pytest.mark.parametrize("name", ["differential_evolution", "jso", "random_search", "multilevel_coordinate_search", "cma_es", "turbo", "dycors"])
def test_callback_errors_propagate(name):
    # Including from the worker threads, where threads= is taken.
    kw = {"threads": 2} if "threads=" in getattr(o, name).__doc__.splitlines()[0] else {}
    with pytest.raises(ZeroDivisionError):
        getattr(o, name)(lambda x: 1 / 0, [-1.0], [1.0], **kw)
    with pytest.raises(TypeError):
        getattr(o, name)(lambda x: "abc", [-1.0], [1.0], **kw)


def test_dycors():
    # The random sequence of a seed depends on the C++ standard library; over 100 seeds the error is below 2e-3.
    r = o.dycors(quadratic, [-5.0, -5.0], [5.0, 5.0], max_function_calls=100, seed=3)
    assert np.allclose(r.x, [1, -2], atol=0.1) and r.nfev == 100
    assert np.array_equal(o.dycors(quadratic, [-5.0, -5.0], [5.0, 5.0], max_function_calls=50, seed=3).x,
                          o.dycors(quadratic, [-5.0, -5.0], [5.0, 5.0], max_function_calls=50, seed=3).x)
    assert o.dycors(quadratic, [-5.0, -5.0], [5.0, 5.0], target_value=1.0).fun <= 1.0
    for bad in [{"initial_points": 3}, {"max_function_calls": 5}]:
        with pytest.raises(ValueError):
            o.dycors(quadratic, [-5.0, -5.0], [5.0, 5.0], **bad)


@pytest.mark.parametrize("name", ["gradient_descent", "nesterov", "lbfgs"])
def test_gradient_optimizers(name):
    fn = getattr(o, name)
    for r in [fn(quadratic, [0.0, 0.0], quadratic_gradient, tolerance=1e-8),
              fn(lambda x: (quadratic(x), np.array(quadratic_gradient(x))), np.zeros(2), tolerance=1e-8)]:
        assert np.allclose(r.x, [1, -2], atol=1e-8) and r.fun == quadratic(r.x) and r.converged and r.nfev > r.nit > 0
    assert not fn(quadratic, [0.0, 0.0], quadratic_gradient, max_iterations=1).converged
    with pytest.raises(ValueError):
        fn(quadratic, [0.0, 0.0], lambda x: [1.0])
    with pytest.raises(ZeroDivisionError):
        fn(lambda x: 1 / 0, [0.0, 0.0], quadratic_gradient)


def test_brent():
    r = o.brent(math.cos, 0, 6)
    assert abs(r.x - math.pi) < 1e-7 and r.fun == math.cos(r.x) and r.converged and r.nfev == r.nit + 1
    assert not o.brent(math.cos, 0, 6, max_iterations=2).converged
    with pytest.raises(ValueError):
        o.brent(math.cos, 6, 0)


def test_quadrature_result_types():
    q = bm.quadrature
    r = q.tanh_sinh(math.exp, 0, 1)
    assert type(r) is q.QuadratureResult and r.value == r[0] and r.l1_norm == r.value and r.error < 1e-10
    assert type(q.gauss_hermite(math.cos)) is q.GaussResult and q.gauss(math.exp, 0, 1).l1_norm > 0
    assert type(q.ooura_fourier_cos(lambda t: 1 / (1 + t * t), 1.0)) is q.OouraResult
    assert type(q.leja_legendre().integrate(math.cos)) is q.QuadratureResult
    assert q.QuadratureResult.__module__ == "boost.math.quadrature" and o.OptimizeResult.__module__ == "boost.math.optimize"
