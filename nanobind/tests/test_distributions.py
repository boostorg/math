#  Copyright 2026 Nick Thompson
#  Distributed under the Boost Software License, Version 1.0.
#  https://www.boost.org/LICENSE_1_0.txt
"""rvs, interval, the parameter finders, find_location and find_scale, and the empirical cdf."""

import numpy as np
import pytest

import boost.math as bm

d = bm.distributions


class Fixed:
    """A stand-in for numpy.random.Generator whose random() returns the given uniforms."""

    def __init__(self, u):
        self.u = np.asarray(u, dtype=float)

    def random(self, size):
        return self.u


def test_rvs_and_interval():
    n = d.normal_distribution(1.0, 2.0)
    assert type(n.rvs()) is float and n.rvs((2, 3)).shape == (2, 3)
    assert np.array_equal(n.rvs(5, rng=np.random.default_rng(7)), n.rvs(5, rng=np.random.default_rng(7)))
    # The ends of the uniforms' range are inside the support, and the upper half is by complement.
    lo, mid, hi = n.rvs(rng=Fixed([0.0, 0.5, 1 - 2.0 ** -53]))
    assert np.isfinite(lo) and np.isfinite(hi) and lo - 1 == pytest.approx(1 - hi) and mid > 1
    lower, upper = n.interval(0.95)
    assert n.cdf(lower) == pytest.approx(0.025) and n.cdf_complement(upper) == pytest.approx(0.025)
    for bad in (-0.1, 1.5, np.nan):
        with pytest.raises(ValueError, match="confidence"):
            n.interval(bad)


@pytest.mark.parametrize("dist", [d.binomial_distribution(7.0, 0.5), d.poisson_distribution(4.5),
                                  d.hypergeometric_distribution(10, 20, 50), d.negative_binomial_distribution(5.0, 0.3)],
                         ids=repr)
def test_discrete_rvs_is_the_least_k_with_cdf_at_least_u(dist):
    assert type(dist.rvs()) is int
    u = np.arange(33) * (2 ** 53 // 33) * 2.0 ** -53  # 0 to 32/33, on the 2^-53 grid of Generator.random
    k = dist.rvs(rng=Fixed(u))
    assert k.dtype == np.int64
    p = u + 2.0 ** -54
    assert np.all(dist.cdf(k.astype(float)) >= p) and np.all((k == 0) | (dist.cdf(np.maximum(k - 1.0, 0)) < p))


def test_finders():
    t = d.students_t_distribution
    assert t(t.find_degrees_of_freedom(2.0, 0.9)).cdf(2.0) == pytest.approx(0.9)
    assert t.find_degrees_of_freedom(1.0, 0.05, 0.1, 2.0) > 1
    assert d.chi_squared_distribution.find_degrees_of_freedom(1.0, 0.05, 0.1, 2.0) > 1
    b = d.binomial_distribution
    assert b.find_lower_bound_on_p(20, 5, 0.05) < 0.25 < b.find_upper_bound_on_p(20, 5, 0.05)
    assert b.find_lower_bound_on_p(20, 5, 0.05, interval="jeffreys_prior") > b.find_lower_bound_on_p(20, 5, 0.05)
    with pytest.raises(ValueError):
        b.find_lower_bound_on_p(20, 5, 0.05, interval="wald")
    assert b.find_maximum_number_of_trials(5, 0.5, 0.05) < b.find_minimum_number_of_trials(5, 0.5, 0.05)
    for c in (d.negative_binomial_distribution, d.geometric_distribution):
        assert c.find_maximum_number_of_trials(10, 0.5, 0.05) < c.find_minimum_number_of_trials(10, 0.5, 0.05)
    assert d.geometric_distribution.find_lower_bound_on_p(20, 0.05) < d.geometric_distribution.find_upper_bound_on_p(20, 0.05)
    nc = d.non_central_chi_squared_distribution
    assert nc(nc.find_degrees_of_freedom(1.5, 4.0, 0.5), 1.5).cdf(4.0) == pytest.approx(0.5)
    assert nc(3.0, nc.find_non_centrality(3.0, 4.0, 0.5)).cdf(4.0) == pytest.approx(0.5)
    nt = d.non_central_t_distribution
    assert nt(5.0, nt.find_non_centrality(5.0, 2.0, 0.6)).cdf(2.0) == pytest.approx(0.6)
    nf = d.non_central_f_distribution
    assert nf(3.0, 7.0, nf.find_non_centrality(2.0, 3.0, 7.0, 0.5)).cdf(2.0) == pytest.approx(0.5)
    with pytest.raises(ValueError):
        t.find_degrees_of_freedom(2.0, 1.5)
    assert d.normal_distribution(d.find_location(d.normal_distribution, 2.0, 0.95, 1.5), 1.5).cdf(2.0) == pytest.approx(0.95)
    assert d.cauchy_distribution(0.5, d.find_scale(d.cauchy_distribution, 2.0, 0.95, 0.5)).cdf(2.0) == pytest.approx(0.95)
    with pytest.raises(TypeError):
        d.find_location(d.gamma_distribution, 2.0, 0.95, 1.0)


def test_empirical_cdf():
    F = d.empirical_cumulative_distribution_function(np.array([3.0, 1.0, 2.0, 2.0]))
    assert F(2.0) == 0.75 and np.array_equal(F([0, 1, 2.5, 3, 10]), [0, 0.25, 0.75, 1, 1])
    assert d.empirical_cumulative_distribution_function([1.0, 2.0], sorted=True)(1.5) == 0.5
    for bad in ([], [[1.0]], [1.0, np.nan]):
        with pytest.raises(ValueError):
            d.empirical_cumulative_distribution_function(bad)

