#  Copyright 2026 Nick Thompson
#  Distributed under the Boost Software License, Version 1.0.
#  https://www.boost.org/LICENSE_1_0.txt
"""Smoke tests for boost.math.statistics against NumPy, SciPy and small examples."""

import numpy as np
import pytest
import scipy.stats as ss

import boost.math.statistics as st

rng = np.random.default_rng(12345)
X = rng.normal(size=200)
Y = 0.5 * X + rng.normal(size=200)


def close(a, b, rtol=1e-10):
    np.testing.assert_allclose(a, b, rtol=rtol, atol=1e-13)


def test_univariate():
    x = X[::3]  # strided: converted, not mutated
    before = x.copy()
    close(st.mean(x), x.mean())
    close(st.variance(x), x.var())
    close(st.sample_variance(x), x.var(ddof=1))
    close(st.mean_and_sample_variance(x), (x.mean(), x.var(ddof=1)))
    m = st.first_four_moments(x)
    close(m, [x.mean()] + [ss.moment(x, k) for k in (2, 3, 4)])
    assert m.m2 == m[1]
    close(st.skewness(x), ss.skew(x))
    close(st.kurtosis(x), ss.kurtosis(x, fisher=False))
    close(st.excess_kurtosis(x), ss.kurtosis(x))
    close(st.median(x), np.median(x))
    close(st.median_absolute_deviation(x), ss.median_abs_deviation(x))
    np.testing.assert_array_equal(x, before)
    assert st.median_absolute_deviation([1.0, 2, 3, 4, 10], center=0.0) == 3
    assert st.interquartile_range([1.0, 2, 3, 4, 5, 6, 7, 8]) == 4
    assert st.gini_coefficient([0.0, 0, 0, 1]) == 0.75 and st.sample_gini_coefficient([0.0, 0, 0, 1]) == 1
    assert st.mode([3.0, 1.0, 2.0, 2.0, 3.0]) == [2.0, 3.0]


def test_bivariate_and_regression():
    close(st.covariance(X, Y), np.cov(X, Y, ddof=0)[0, 1])
    close(st.means_and_covariance(X, Y), (X.mean(), Y.mean(), np.cov(X, Y, ddof=0)[0, 1]))
    close(st.correlation_coefficient(X, Y), np.corrcoef(X, Y)[0, 1])
    r = ss.linregress(X, Y)
    close(st.simple_ordinary_least_squares(X, Y), (r.intercept, r.slope))
    close(st.simple_ordinary_least_squares_with_R_squared(X, Y), (r.intercept, r.slope, r.rvalue**2))
    with pytest.raises(ValueError):
        st.covariance(X, Y[:-1])


def test_signal():
    s, noisy = [1.0, -1.0], [1.5, -0.5]
    assert st.oracle_snr(s, noisy) == 4 and st.mean_invariant_oracle_snr(s, noisy) == 4
    assert st.oracle_snr_db(s, noisy) == pytest.approx(10 * np.log10(4)) == st.mean_invariant_oracle_snr_db(s, noisy)
    assert st.hoyer_sparsity([1.0, 0, 0, 0]) == 1 and st.hoyer_sparsity([1.0, 1, 1, 1]) == 0
    assert st.absolute_gini_coefficient([0.0, 0, 0, -1]) == 0.75 and st.sample_absolute_gini_coefficient([0.0, 0, 0, -1]) == 1
    x = [1.2, -0.9, 1.1, -1.3, 0.8, -1.0]
    assert st.m2m4_snr_estimator_db(x) == pytest.approx(10 * np.log10(st.m2m4_snr_estimator(x)))


def test_hypothesis_tests():
    r = st.one_sample_t_test(X, 0.1)
    close(r, ss.ttest_1samp(X, 0.1))
    assert r.pvalue == r[1]
    close(st.two_sample_t_test(X, Y[:150]), ss.ttest_ind(X, Y[:150]))
    close(st.two_sample_t_test(X, 3 * Y), ss.ttest_ind(X, 3 * Y, equal_var=False))  # Welch
    close(st.paired_samples_t_test(X, Y), ss.ttest_rel(X, Y))
    close(st.anderson_darling_normality_statistic(X), ss.anderson(X, method='interpolate').statistic)
    # A single sample needs the standard deviation:
    assert st.anderson_darling_normality_statistic([1.0], sd=1.0) > 0
    with pytest.raises(ValueError):
        st.anderson_darling_normality_statistic([1.0])
    # Lag-1 autocorrelation -3/4: Q = 4 * 6 * (9 / 16) / 3.
    close(st.ljung_box([1.0, -1, 1, -1], lags=1), (4.5, ss.chi2.sf(4.5, 1)))
    t, p = st.runs_above_and_below_median(X)
    assert -4 < t < 4 and 0 < p <= 1
    close(st.runs_above_and_below_threshold(X, float(np.median(X))), (t, p))
    x = np.array([1.5, 1.1, 1.4, 1.2, 1.7, 1.0, 1.3, 1.6, 0.9, 1.8])
    y = np.array([1.2, 0.8, 1.1, 0.9, 1.3, 1.0, 0.7, 1.4])
    z = (x.mean() - y.mean()) / np.sqrt(x.var(ddof=1) / x.size + y.var(ddof=1) / y.size)
    assert st.two_sample_z_test(x, y).statistic == pytest.approx(z, rel=1e-12) and len(st.one_sample_z_test(x, 1.0)) == 2


def test_chatterjee_and_errors():
    x = np.sort(X)
    close(st.chatterjee_correlation(x, x**2), ss.chatterjeexi(x, x**2).statistic)
    assert 0 < st.chatterjee_correlation(x, x**2, M=3) <= 1
    for bad in [lambda: st.chatterjee_correlation(X, Y), lambda: st.mean([]), lambda: st.sample_variance([1.0]),
                lambda: st.median(np.zeros((2, 2))), lambda: st.ljung_box(X, lags=1000)]:
        with pytest.raises(ValueError):
            bad()
