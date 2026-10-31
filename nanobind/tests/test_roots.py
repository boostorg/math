#  Copyright 2026 Nick Thompson
#  Distributed under the Boost Software License, Version 1.0.
#  https://www.boost.org/LICENSE_1_0.txt
"""Smoke tests for boost.math.roots; the C++ is tested in test/."""

import math

import numpy as np
import pytest

import boost.math as bm
from boost.math import roots

SQRT2 = math.sqrt(2)


@pytest.mark.parametrize("solve", [roots.bisect, roots.toms748])
def test_bracketing(solve):
    r = solve(lambda x: x * x - 2, 0.0, 2.0)
    assert type(r) is roots.RootResult and r.lower <= SQRT2 <= r.upper and math.isclose(r.root, SQRT2, rel_tol=1e-15)
    assert solve(lambda x: x * x - 2, 0.0, 2.0, digits=10).iterations < r.iterations
    assert math.isclose(roots.bracket_and_solve_root(lambda x: math.log(x) - 3, 1.0, 2.0, True).root, math.exp(3))


@pytest.mark.parametrize("solve, f", [
    (roots.newton, lambda x: (math.cos(x) - x, -math.sin(x) - 1)),
    (roots.halley, lambda x: (math.cos(x) - x, -math.sin(x) - 1, -math.cos(x))),
    (roots.schroder, lambda x: (math.cos(x) - x, -math.sin(x) - 1, -math.cos(x))),
])
def test_derivative_based(solve, f):
    r = solve(f, 0.5, 0.0, 1.0)
    assert type(r) is roots.NewtonResult and math.isclose(r.root, 0.7390851332151607, rel_tol=1e-15)
    assert abs(roots.complex_newton(lambda z: (z * z + 1, 2 * z), 1 + 1j) - 1j) < 1e-15


def test_errors():
    with pytest.raises(ValueError, match="bracket"):
        roots.toms748(lambda x: x * x + 1, 0.0, 1.0)
    with bm.errstate(domain="ignore"):
        assert math.isnan(roots.toms748(lambda x: x * x + 1, 0.0, 1.0).root)
    with pytest.raises(RuntimeError, match="No change of sign"):
        roots.bisect(lambda x: x * x + 1, 0.0, 1.0)
    # An error in a boost.math function the callback calls propagates:
    with pytest.raises(ValueError, match="tgamma"):
        roots.toms748(lambda x: bm.tgamma(-1.0) if x > 0.5 else x - 0.3, 0.0, 1.0)
    with pytest.raises(ZeroDivisionError):
        roots.newton(lambda x: (1 / 0, 1.0), 0.5, 0.0, 1.0)
    with pytest.raises(TypeError):
        roots.toms748(math.sin, a=3.0, b=4.0)


def test_polynomials():
    assert np.allclose(roots.quadratic_roots(1, 0, -2), [-SQRT2, SQRT2]) and roots.quadratic_roots(1, 0, 2).size == 0
    assert np.allclose(roots.cubic_roots(1, -6, 11, -6), [1, 2, 3])
    assert np.allclose(roots.quartic_roots(1, 0, -5, 0, 4), [-2, -1, 1, 2])
    # x^2 - 2, coefficients from the constant term up.
    for solve in (roots.polynomial_roots, roots.polynomial_roots_aberth, roots.polynomial_roots_companion_matrix):
        assert np.allclose(sorted(solve([-2.0, 0.0, 1.0]).real), [-SQRT2, SQRT2])
