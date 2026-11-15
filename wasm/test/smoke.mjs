// Copyright 2026 Nick Thompson
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// Smoke test of the bindings; the functions themselves are tested in test/.

import { test } from 'node:test';
import assert from 'node:assert/strict';
import { pathToFileURL } from 'node:url';

const createModule = (await import(pathToFileURL(process.env.BOOST_MATH_WASM ?? 'build-wasm/boost_math.mjs'))).default;
const bm = await createModule();

function close(actual, expected, tol = 1e-13) {
   assert.ok(Math.abs(actual - expected) <= tol * Math.abs(expected), `${actual} != ${expected}`);
}

// References from mpmath at 40 digits.
test('special functions', () => {
   close(bm.tgamma(4.5), 11.631728396567452);
   close(bm.lgamma(100.5), 361.43554046777762);
   close(bm.digamma(2.5), 0.70315664064524319);
   close(bm.gamma_p(3, 2), 0.32332358381693654);
   close(bm.ibeta(2, 3, 0.4), 0.5248);
   close(bm.erf(0.5), 0.52049987781304654);
   close(bm.erfc(5), 1.5374597944280349e-12);
   close(bm.expint(1), 1.8951178163559368);
   close(bm.zeta(3), 1.2020569031595943);
   close(bm.cyl_bessel_j(2.5, 4), 0.44088497455734117);
   close(bm.cyl_bessel_k(1, 2), 0.13986588181652243);
   close(bm.sph_bessel(3, 2), 0.060722097662874828);
   close(bm.airy_ai(1), 0.13529241631288142);
   close(bm.ellint_1(0.5), 1.685750354812596);
   close(bm.ellint_2(0.5, 1), 0.96487645426862749);
   close(bm.lambert_w0(1), 0.56714329040978387);
   close(bm.owens_t(0.5, 0.5), 0.064488602847503757);
   // The recurrence's error terms come from two_product.
   close(bm.laguerre(100, 3.5), 0.38531832348347825);
   close(bm.hypergeometric_1F1(1, 2, 1), Math.E - 1);
   close(bm.cyl_bessel_j_zero(0, 1), 2.404825557695773);
});

test('errors return NaN or infinity', () => {
   assert.ok(Number.isNaN(bm.tgamma(-2)));
   assert.ok(Number.isNaN(bm.gamma_p(-1, 2)));
   assert.ok(Number.isNaN(bm.laguerre(-1, 2)));
   assert.ok(Number.isNaN(bm.legendre_p(2.7, 0.5)));
   assert.equal(bm.tgamma(200), Infinity);
   assert.ok(Number.isNaN(new bm.normal_distribution(0, -1).cdf(0)));
   assert.ok(Number.isNaN(new bm.cauchy_distribution(0, 1).mean()));
   assert.ok(Number.isNaN(new bm.geometric_distribution(2.2).mean()));
});

test('distributions', () => {
   close(new bm.normal_distribution(0, 1).cdf(1.5), 0.93319279873114193);
   close(new bm.students_t_distribution(5).ccdf(2), 0.050969739414929178);
   close(new bm.chi_squared_distribution(4).quantile(0.95), 9.4877290367811568);
   close(new bm.poisson_distribution(3).pdf(2), 0.22404180765538774);
   close(new bm.binomial_distribution(10, 0.3).cdf(3), 0.6496107184);
   const g = new bm.gamma_distribution(2.5, 1.5);
   close(g.quantile(g.cdf(4)), 4, 1e-12);
   close(g.mean(), 3.75);
   close(g.variance(), 5.625);
   close(new bm.normal_distribution(0, 1).entropy(), 1.4189385332046727);  // log(2 pi e) / 2
   assert.ok(Number.isNaN(new bm.non_central_t_distribution(5, 1).entropy()));
   assert.deepEqual(new bm.beta_distribution(2, 3).support(), [0, 1]);
   assert.ok(Number.isNaN(new bm.landau_distribution(0, 1).mean()));  // undefined
});

test('fft', () => {
   for (const n of [8, 7, 1000]) {  // radix 2, Bluestein, mixed radix
      const x = Float64Array.from({ length: 2 * n }, (_, i) => Math.sin(i * i));
      const f = new bm.fft(n), y = f.inverse_transform(f.transform(x));
      x.forEach((v, i) => assert.ok(Math.abs(y[i] - v) < 1e-13));
      f.delete();
   }
   // cos(2 pi 3 j / 16) has X[3] = 8 and nothing else.
   const r = new bm.real_fft(16), X = r.transform(Float64Array.from({ length: 16 }, (_, j) => Math.cos(2 * Math.PI * 3 * j / 16)));
   assert.equal(X.length, 18);
   X.forEach((v, i) => assert.ok(Math.abs(v - (i === 6 ? 8 : 0)) < 1e-13, `X[${i >> 1}]`));
   assert.throws(() => r.transform(new Float64Array(15)), RangeError);
   assert.throws(() => new bm.fft(0), RangeError);
   r.delete();
});

test('arrays', () => {
   const xs = Float64Array.from({ length: 101 }, (_, i) => -2 + i / 25);
   const ys = bm.erf_array(xs);
   assert.ok(ys instanceof Float64Array);
   xs.forEach((x, i) => assert.equal(ys[i], bm.erf(x)));
   assert.deepEqual([...bm.gamma_p_array(3, [1, 2])], [bm.gamma_p(3, 1), bm.gamma_p(3, 2)]);
   assert.deepEqual([...bm.ibeta_array(2, 3, [0.4])], [bm.ibeta(2, 3, 0.4)]);
   assert.deepEqual([...bm.laguerre_array(100, [3.5])], [bm.laguerre(100, 3.5)]);
   assert.ok(bm.laguerre_array(-1, [1, 2]).every(Number.isNaN));
   assert.equal(bm.tgamma_array([]).length, 0);
});

// Each of the calls throws a RangeError.
const rejects = (...calls) => calls.forEach((call) => assert.throws(call, RangeError, String(call)));
// Every number in expected, an array or object of them at any depth, is that of actual.
function same(actual, expected, tol = 1e-9) {
   if (typeof expected !== 'number') Object.keys(expected).forEach((k) => same(actual[k], expected[k], tol));
   else assert.ok(Object.is(actual, expected) || Math.abs(actual - expected) <= tol * (1 + Math.abs(expected)), `${actual} != ${expected}`);
}
const check = (rows, tol) => rows.forEach(([name, args, expected]) => same(bm[name](...args), expected, tol));  // rows [name, arguments, expected]

// name(arguments)=value for every other function, from mpmath and SciPy, or the library's own for heuman_lambda and the spheroidal
// angular and second radial functions. Each has an _array twin over its last argument, but for a few; those without a value are above.
const calls = `
   lgamma(100.5) digamma(2.5) erfc(5) expint(1) zeta(3) cyl_bessel_j(2.5,4) cyl_bessel_k(1,2) sph_bessel(3,2) airy_ai(1) ellint_1(0.5) ellint_2(0.5,1) lambert_w0(1)
   owens_t(0.5,0.5) hypergeometric_1F1(1,2,1) tgamma(2.5,1)=1.128803 trigamma(2.5)=0.4903578 polygamma(2,1.5)=-0.8287966 tgamma1pm1(0.1)=-0.04864923
   tgamma_lower(3,2)=0.6466472 tgamma_ratio(5.5,3)=26.17139 tgamma_delta_ratio(5,0.5)=0.458516 gamma_q(3,2)=0.6766764 gamma_p_inv(3,0.3)=1.913776
   gamma_q_inv(3,0.3)=3.615568 gamma_p_derivative(3,2)=0.2706706 beta(2,3)=0.08333333 beta(2,3,0.4)=0.04373333 betac(2,3,0.4)=0.0396 ibetac(2,3,0.4)=0.4752
   ibeta_inv(2,3,0.5)=0.3857276 ibetac_inv(2,3,0.3)=0.5084048 ibeta_derivative(2,3,0.4)=1.728 erf_inv(0.5)=0.4769363 erfc_inv(0.3)=0.7328691 expint(2,1.5)=0.07310079
   cyl_neumann(1.5,2)=-0.3956233 cyl_bessel_i(1.5,2)=1.099473 sph_neumann(2,1.5)=-1.345713 airy_bi(1)=1.207424 airy_ai_prime(1)=-0.1591474 airy_bi_prime(1)=0.9324359
   ellint_1(0.5,1)=1.037356 ellint_2(0.5)=1.467462 ellint_3(0.5,0.2)=1.892295 ellint_3(0.5,0.2,1)=1.101305 lambert_wm1(-0.2)=-2.542641 legendre_p(3,0.4)=-0.44
   hermite(3,0.4)=-4.288 chebyshev_t(3,0.4)=-0.944 log1p(0.3)=0.2623643 expm1(0.3)=0.3498588 cbrt(5)=1.709976 sqrt1pm1(0.3)=0.1401754 powm1(2,0.3)=0.2311444
   sinc_pi(0.7)=0.920311 sin_pi(0.3)=0.809017 cos_pi(0.3)=0.5877853 hypergeometric_0F1(2,0.7)=1.393301 hypergeometric_2F0(1,-2,0.3)=0.58 jacobi_sn(0.5,0.7)=0.6342933
   jacobi_cn(0.5,0.7)=0.7730925 jacobi_dn(0.5,0.7)=0.9483765 jacobi_zeta(0.5,0.7)=0.06561554 heuman_lambda(0.5,0.7)=0.6035963 ellint_d(0.5)=0.8731526
   ellint_d(0.5,0.7)=0.1071653 jacobi_theta1(0.7,0.3)=0.8381788 jacobi_theta2(0.7,0.3)=1.06383 jacobi_theta3(0.7,0.3)=1.086697 jacobi_theta4(0.7,0.3)=0.882775
   prolate_spheroidal_eigenvalue(2,1,1.5)=4.694785 oblate_spheroidal_eigenvalue(2,1,1.5)=7.265886 prolate_spheroidal_angular(2,1,1.5,0.4)=-1.146412
   prolate_spheroidal_angular_prime(2,1,1.5,0.4)=-2.174909 oblate_spheroidal_angular(2,1,1.5,0.4)=-1.051592 oblate_spheroidal_angular_prime(2,1,1.5,0.4)=-2.265388
   prolate_spheroidal_radial_1(2,1,1.5,1.7)=0.2379539 prolate_spheroidal_radial_1_prime(2,1,1.5,1.7)=0.2115861 prolate_spheroidal_radial_2(2,1,1.5,1.7)=-0.5825374
   prolate_spheroidal_radial_2_prime(2,1,1.5,1.7)=0.9643756 oblate_spheroidal_radial_1(2,1,1.5,0.6)=0.09031673 oblate_spheroidal_radial_1_prime(2,1,1.5,0.6)=0.1723524
   oblate_spheroidal_radial_2(2,1,1.5,0.6)=-1.496373 oblate_spheroidal_radial_2_prime(2,1,1.5,0.6)=2.571978 struve_h(1.5,2)=0.4506611 struve_h_prime(1.5,2)=0.4609795
   struve_l(1.5,2)=0.7028463 struve_l_prime(1.5,2)=1.031267 fresnel_c(1.5)=0.4452612 fresnel_s(1.5)=0.697505 ellint_rf(1,2,3)=0.7269459 ellint_rd(1,2,3)=0.2904603
   ellint_rc(1,2)=0.7853982 ellint_rg(1,2,3)=1.401847 cyl_bessel_j_prime(1.5,2)=0.1445458 cyl_neumann_prime(1.5,2)=0.5315032 cyl_bessel_i_prime(1.5,2)=1.221632
   cyl_bessel_k_prime(1.5,2)=-0.2548678 cyl_neumann_zero(0.5,2)=4.712389 airy_ai_zero(2)=-4.087949 airy_bi_zero(2)=-3.271093 gegenbauer(3,0.7,0.4)=-0.677824
   spherical_harmonic_r(3,2,0.7,0.3)=0.2677395 spherical_harmonic_i(3,2,0.7,0.3)=0.1831704`;

test('every other special function, and the _array twins', () => {
   const scalar_only = ['cyl_neumann_zero', 'airy_ai_zero', 'airy_bi_zero', 'gegenbauer', 'spherical_harmonic_r', 'spherical_harmonic_i'];
   for (const call of calls.trim().split(/\s+/)) {
      const [, name, list, expected] = call.match(/(\w+)\((.*)\)=?(.*)/), args = JSON.parse(`[${list}]`), x = args.pop(), y = bm[name](...args, x);
      if (expected) close(y, +expected, 1e-6);
      if (!scalar_only.includes(name)) assert.deepEqual([...bm[name + '_array'](...args, [x])], [y], name);
   }
   assert.deepEqual([...bm.fresnel_cis(1.5), ...bm.fresnel_cis_array([-2, 1.5])], [bm.fresnel_c(1.5), bm.fresnel_s(1.5), ...bm.fresnel_cis(-2), ...bm.fresnel_cis(1.5)]);
});

test('every distribution', () => {
   // [name, parameters, x, cdf(x)], from SciPy, or the library's own for landau, mapairy, kolmogorov_smirnov and non_central_beta.
   for (const [name, parameters, x, cdf] of [
      ['normal', [0, 1]], ['students_t', [5], 2, 0.9490303], ['chi_squared', [4], 9.487729, 0.95], ['poisson', [3], 2, 0.4231901], ['binomial', [10, 0.3]],
      ['beta', [2, 3], 0.4, 0.5248], ['gamma', [2.5, 1.5], 4, 0.6234324], ['exponential', [1.5], 0.4, 0.4511884], ['lognormal', [0.5, 0.8], 2, 0.5953906],
      ['weibull', [1.5, 2], 1.2, 0.3717129], ['cauchy', [1, 2], 2.5, 0.7048328], ['fisher_f', [4, 7], 1.3, 0.6434688], ['uniform', [1, 4], 2, 0.3333333],
      ['laplace', [1, 2], 2.5, 0.7638167], ['logistic', [1, 2], 2.5, 0.6791787], ['non_central_chi_squared', [3, 2], 4, 0.4838814], ['non_central_t', [5, 1], 2, 0.7780747],
      ['negative_binomial', [4, 0.3], 6, 0.3503893], ['geometric', [0.3], 2, 0.657], ['arcsine', [1, 4], 2, 0.3918266], ['bernoulli', [0.3], 0, 0.7],
      ['extreme_value', [1, 2], 2.5, 0.6235249], ['generalized_extreme_value', [1, 2, 0.3], 2.5, 0.6014519], ['holtsmark', [1, 2], 2.5, 0.7016819],
      ['hypergeometric', [6, 5, 20], 2, 0.8686791], ['inverse_chi_squared', [4, 0.5], 0.7, 0.5820097], ['inverse_gamma', [2.5, 1.5], 0.7, 0.5090552],
      ['inverse_gaussian', [1.5, 2], 1.2, 0.5430299], ['kolmogorov_smirnov', [10], 0.3, 0.6708952], ['landau', [1, 2], 2.5, 0.4417172], ['mapairy', [1, 2], 2.5, 0.7868677],
      ['non_central_beta', [2, 3, 1.5], 0.4, 0.3876983], ['non_central_f', [4, 7, 1.5], 1.3, 0.5026782], ['pareto', [1.5, 2.5], 2.5, 0.7211452],
      ['rayleigh', [1.5], 1.2, 0.273851], ['saspoint5', [1, 2], 2.5, 0.7037734], ['skew_normal', [1, 2, 3], 2.5, 0.5475253], ['triangular', [1, 2, 4], 2.5, 0.625],
      ['von_mises', [1, 2], 1.5, 0.7381922]]) {
      const d = new bm[name + '_distribution'](...parameters);
      if (cdf) close(d.cdf(x), cdf, 1e-6);  // those without are above
      for (const f of ['pdf', 'cdf', 'ccdf', 'quantile']) assert.deepEqual([...d[f + '_array'](new Float32Array([0.5, 0.25]))], [d[f](0.5), d[f](0.25)], name);
      d.delete();
   }
   same([new bm.gamma_distribution(2.5, 1.5).pdf(3), new bm.gamma_distribution(2.5, 1.5).ccdf(3)], [0.19196788093577977, 0.54941595135278032]);
   assert.ok(new bm.normal_distribution(0, -1).pdf_array([0, 1]).every(Number.isNaN));
});

test('transforms along an axis', () => {
   // Along an axis, of two rows of 8 complex or 16 real numbers. The power spectrum of 3 cos(2 pi 3 j / 16) has P[3] = 9/2.
   const z = Float64Array.from({ length: 32 }, (_, i) => Math.sin(i * i)), row = z.subarray(16), f = new bm.fft(8), r = new bm.real_fft(16), Z = bm.fft_axis(z, [2, 8], -1);
   assert.deepEqual([...Z.subarray(16)], [...f.transform(row)]);
   assert.deepEqual([...bm.inverse_fft_axis(Z, [2, 8], 1).subarray(16)], [...f.inverse_transform(Z.subarray(16))]);
   assert.equal(bm.fft_axis(z, [2, 8], 0)[16], z[0] - z[16]);
   assert.deepEqual([...bm.real_fft_axis(z, [2, 16], 1).subarray(18)], [...r.transform(row)]);
   assert.deepEqual([...bm.power_spectrum_axis(z, [2, 16], 1).subarray(9)], [...bm.power_spectrum(row)]);
   same(bm.power_spectrum(Float64Array.from({ length: 16 }, (_, j) => 3 * Math.cos(2 * Math.PI * 3 * j / 16))), [0, 0, 0, 4.5, 0, 0, 0, 0, 0]);
   rejects(() => bm.fft_axis(z, [2, 8], 2), () => bm.fft_axis(z, [2, 8], 0.5), () => bm.fft_axis(z, [2, 5], 0), () => bm.power_spectrum_axis(z, [6.4, 5], 0));
   f.delete();
   r.delete();
});

// Native Boost references, from clang++ -O2 on arm64.
test('statistics', () => {
   const x = [2.1, 3.7, 1.4, 5.9, 4.2, 3.3, 6.8, 2.6, 4.9, 3.1], y = new Float64Array([1.0, 2.9, 0.7, 5.2, 3.1, 2.2, 6.9, 1.8, 4.4, 2.0]);
   const t = (statistic, p_value) => ({ statistic, p_value }), runs = t(2.012461, 0.04417134);
   check([['mean', [x], 3.8], ['variance', [x], 2.562], ['sample_variance', [x], 2.846667], ['mean_and_sample_variance', [x], { mean: 3.8, sample_variance: 2.846667 }],
      ['skewness', [x], 0.407774], ['kurtosis', [x], 2.221927], ['excess_kurtosis', [x], -0.7780733], ['median', [x], 3.5], ['median_absolute_deviation', [x], 1.15],
      ['interquartile_range', [x], 2.3], ['gini_coefficient', [x], 0.2378947], ['sample_gini_coefficient', [x], 0.2643275], ['hoyer_sparsity', [x], 0.1146856],
      ['covariance', [x, y], 2.937], ['correlation_coefficient', [x, y], 0.9893738], ['anderson_darling_normality_statistic', [x], 0.1654371],
      ['simple_ordinary_least_squares', [x, y], { intercept: -1.336206, slope: 1.14637, r_squared: 0.9788606 }], ['runs_above_and_below_median', [x], runs],
      ['runs_above_and_below_threshold', [x, 3.5], runs], ['one_sample_t_test', [x, 3], t(1.499414, 0.1679996)], ['one_sample_z_test', [x, 3], t(1.499414, 0.1337662)],
      ['two_sample_t_test', [x, y], t(0.9551737, 0.352143)], ['two_sample_z_test', [x, y], t(0.9551737, 0.3394899)], ['ljung_box', [x, 3, 1], t(3.944437, 0.1391478)],
      ['paired_samples_t_test', [x, y], t(6.550742, 0.0001050808)], ['ljung_box', [x, -1, 0], t(3.944437, 0.2675268)]], 1e-6);
   assert.ok(x[0] === 2.1 && Number.isNaN(bm.median([1, NaN, 2])) && Number.isNaN(bm.runs_above_and_below_threshold([1, NaN, 2], 1.5).p_value));
   rejects(() => bm.mean([]), () => bm.covariance([1, 2], [1]), () => bm.simple_ordinary_least_squares([1, 1], [1, 2]), () => bm.ljung_box(x, 10, 0),
           () => bm.ljung_box(x, 3, 0.5));
});

test('filters', () => {
   const step = Float64Array.from({ length: 64 }, (_, i) => +(i >= 32)), y = step.map((v, i) => v + 0.1 * Math.sin(37 * i)), ramp = step.map((_, i) => 2 * i + 1);
   same(bm.trend_filter([0, 0, 0, 0, 1, 1, 1, 1, 1], 0, 0.1), [0.025, 0.025, 0.025, 0.025, 0.98, 0.98, 0.98, 0.98, 0.98]);
   // Without the order and lambda, or the wavelet and rule, they are those of the _select function; the rule alone is wiener.
   const t = bm.trend_filter_select(y), w = bm.wavelet_denoise_select(y), distance = (x) => x.reduce((s, v, i) => s + (v - step[i]) ** 2, 0);
   assert.deepEqual([bm.trend_filter(y), bm.wavelet_denoise(y), bm.wavelet_denoise(y, 'symlet', 4)],
                    [bm.trend_filter(y, t.k, t.lambda), bm.wavelet_denoise(y, w.family, w.p, w.rule), bm.wavelet_denoise(y, 'symlet', 4, 'wiener')]);
   same(bm.wavelet_denoise(y, 'symlet', 4, 'universal_soft', 1e-9), y, 1e-6);  // no noise to remove
   assert.ok(distance(bm.denoise(y)) < distance(y));
   // A ramp has no detail: the 64 / 2^levels scaling coefficients come first.
   const c = bm.wavelet_transform(ramp, 'symlet', 4), c1 = bm.wavelet_transform(y, 'daubechies', 2, 1), nonzero = (a) => a.findLastIndex((v) => Math.abs(v) > 1e-9);
   assert.deepEqual([bm.wavelet_levels(64, 'symlet', 4), nonzero(c), nonzero(bm.wavelet_transform(ramp, 'symlet', 4, 1))], [3, 7, 31]);
   same([bm.inverse_wavelet_transform(c, 'symlet', 4), bm.inverse_wavelet_transform(c1, 'daubechies', 2, 1)], [ramp, y]);
   rejects(() => bm.trend_filter(y, 1.5, 1), () => bm.trend_filter([1, 2, 3]), () => bm.trend_filter_select([1, 2, NaN, 4]), () => bm.wavelet_transform(y, 'haar', 4),
           () => bm.inverse_wavelet_transform(c, 'symlet', 4, 1.5), () => bm.wavelet_denoise(y, 'symlet', 4, 'hard'), () => bm.wavelet_denoise(y, 'symlet', 4, 'sure', -1));
});

test('interpolators and wavelets', () => {
   // makima and pchip from SciPy. A Catmull-Rom curve goes through its points at their parameters, the chord lengths to the power
   // alpha (1/2 unless given) added up. phi_2(1) = (1 + sqrt 3) / 2 = -psi_2(1), phi_2(2) = (1 - sqrt 3) / 2, and the nodes of every grid agree.
   const x = [0, 1, 2.5, 3, 4.5, 6, 7], y = new Float64Array([1, 3, 2, 5, 4, 4.5, 0]), cubes = x.map((v) => v ** 3), p = [[0, 0], [1, 1], [2, 0], [1, -1]];
   const a = (1 + Math.sqrt(3)) / 2, coarse = { grid_refinements: 5 }, closed = { alpha: 1, closed: true };
   check([['makima', [x, y, [2.8, 8]], [3.990877192982455, NaN]], ['pchip', [x, y, [2.8]], [3.944]], ['barycentric_rational', [x, cubes, [2.8]], [2.8 ** 3]],
      ['barycentric_rational', [x, cubes, [2.8], 1], [21.621059891941602]], ['cardinal_cubic_b_spline', [[0, 1, 4, 9, 16, 25, 36], 1, 0.5, [2.25]], [6.25]],
      ['catmull_rom_parameters', [p], [0, 2 ** 0.25, 2 * 2 ** 0.25, 3 * 2 ** 0.25]], ['catmull_rom_parameters', [p, closed], [0, 1, 2, 3, 4].map((k) => k * Math.SQRT2)],
      ['catmull_rom', [p, bm.catmull_rom_parameters(p)], p.flat()], ['catmull_rom', [p, [4 * Math.SQRT2], closed], [0, 0]],
      ['daubechies_scaling', [2, 1], a], ['daubechies_scaling', [2, 2, coarse], 1 - a], ['daubechies_wavelet', [2, 1], -a], ['daubechies_wavelet', [2, 1, coarse], -a],
      ['daubechies_scaling_array', [2, [1, 2, NaN]], [a, 1 - a, NaN]], ['daubechies_scaling_array', [2, [2], coarse], [1 - a]],
      ['daubechies_wavelet_array', [2, [1, -3.5]], [-a, 0]], ['daubechies_wavelet_array', [2, [1], coarse], [-a]]]);
   // The wavelet has p vanishing moments, so the transform of a line, whose samples reach 3000, is zero; that of a step is not.
   const ramp = Float64Array.from({ length: 1001 }, (_, i) => 3 * i - 7), step = Float64Array.from({ length: 2001 }, (_, i) => +(i >= 1000));
   same(bm.daubechies_wavelet_transform(ramp, 0.01, 4, [0.5, 1], [4, 5, 5.5]), [0, 0, 0, 0, 0, 0], 0.03);
   close(bm.daubechies_wavelet_transform(step, 0.001, 2, [0.25], [1], { max_refinements: 14, grid_refinements: 10 })[0], 0.1136, 0.02);
   rejects(() => bm.makima([0, 1, 1, 2], [1, 2, 3, 4], [1]), () => bm.pchip(x, [1, 2], [1]), () => bm.barycentric_rational(x, y, [1], 7), () => bm.catmull_rom(p, [-1]),
           () => bm.cardinal_cubic_b_spline([1, 2, 3, 4, 5], 0, 0, [1]), () => bm.catmull_rom(p.slice(1), [0]), () => bm.catmull_rom(p, [0], { alpha: 2 }),
           () => bm.daubechies_scaling(11, 0), () => bm.daubechies_scaling(10, 1, { grid_refinements: 21 }), () => bm.daubechies_wavelet_transform(ramp, 0, 4, [1], [1]));
});

const gaussian = (x) => Math.exp(-x * x), rootPi = Math.sqrt(Math.PI), one = () => 1, id = (x) => x, { exp, cos, E, PI } = Math;
const rules = [  // [rule, arguments, options, the integral with them]
   ['tanh_sinh', [exp, 0, 1], { tolerance: 1e-9, max_refinements: 8 }, E - 1], ['exp_sinh', [(x) => exp(-x) / Math.sqrt(x), 0, Infinity], {}, rootPi],
   ['sinh_sinh', [gaussian], {}, rootPi], ['gauss_kronrod', [gaussian, -Infinity, Infinity], { points: 61, max_depth: 5, tolerance: 1e-9 }, rootPi],
   ['gauss', [exp, 0, 1], { points: 15 }, E - 1], ['leja_legendre', [exp, 0, 1], {}, E - 1], ['leja_hermite', [one], {}, rootPi],
   ['leja_jacobi', [one, 0, 1], { alpha: 0, beta: -0.5 }, 2], ['leja_laguerre', [one], { alpha: -0.5, tolerance: 1e-10 }, rootPi],  // 1/sqrt(x); Gamma(1/2)
   ['gauss_hermite', [cos], { points: 16, averaged: true }, rootPi * exp(-0.25)], ['gauss_laguerre', [id], { points: 8, alpha: 1.5 }, 3.3233509704478426],  // Gamma(7/2)
   ['ooura_fourier_sin', [(x) => 1 / x, 1], {}, PI / 2], ['ooura_fourier_cos', [(x) => 1 / (1 + x * x), 1], { tolerance: 1e-12 }, PI / (2 * E)],
   ['gauss_jacobi', [one], { points: 4, alpha: 0.5, beta: 0.5 }, PI / 2], ['sparse_grid', [([x, y]) => x * y * y, [0, 0], [1, 2]], { max_calls: 999 }, 4 / 3]];

test('quadrature', () => {
   for (const [rule, args, options, expected] of rules) {  // each also runs without options
      close(bm[rule](...args, options).result, expected, 1e-9);
      assert.ok(Number.isFinite(bm[rule](...args).result), rule);
   }
   const r = bm.tanh_sinh((x) => Math.log(x) * Math.sqrt(x), 0, 1), h = bm.gauss_hermite(cos, { averaged: true });
   assert.ok(r.error < 1e-8 && r.L1 > 0.44 && r.levels >= 2 && h.error < 1e-8 && h.L1 > 1 && bm.ooura_fourier_sin((x) => 1 / x, 1).relative_error < 1e-8);
   rejects(() => bm.gauss_kronrod(exp, 0, 1, { points: 16 }), () => bm.leja_laguerre(one, { alpha: -2 }), () => bm.sparse_grid(one, [0, 0], [1]),
           () => bm.ooura_fourier_sin((x) => 1 / x, 1, { tolerance: 1e-20 }), () => bm.tanh_sinh(exp, 0, 1, { max_refinements: Infinity }), () => bm.tanh_sinh(7, 0, 1),
           () => bm.tanh_sinh(exp, 0, 1, { tolerance: '1e-4' }), () => bm.tanh_sinh(exp, 0, 1, { tolerence: 1e-4 }), () => bm.tanh_sinh(exp, 0, 1, 3));
   // An integrand that integrates with the same kind of rule and other exponents replaces the kept rule, which the outer one goes on using.
   const inner = (a) => bm.leja_jacobi(one, 0, 1, { beta: a }).result;
   close(bm.leja_jacobi(inner, 0, 1, { beta: 0.5, tolerance: 1e-10 }).result, bm.tanh_sinh((x) => Math.sqrt(x) / (x + 1), 0, 1).result, 1e-9);
   close(bm.sparse_grid(([x]) => bm.sparse_grid(([y]) => x * y, [0], [1]).result, [0], [1]).result, 0.25, 1e-12);
});

test('roots', () => {
   // x^2 + y^2 = 4 and x = y, from (1, 2): F, J and d2(v) = F''(x)[v, v]; or J flat, for Newton.
   const circle = (x) => [[x[0] ** 2 + x[1] ** 2 - 4, x[0] - x[1]], [[2 * x[0], 2 * x[1]], [1, -1]], (v) => [2 * v[0] ** 2 + 2 * v[1] ** 2, 0]];
   const cubic = (x) => [x * x * x - 2 * x - 5, 3 * x * x - 2, 6 * x], f = (x) => cubic(x)[0], root = 2.0945514815423265, bracket = { lower: root, upper: root };
   const flat = (x) => [circle(x)[0], circle(x)[1].flat()], both = [Math.SQRT2, Math.SQRT2], few = { max_iterations: 2 }, two = { iterations: 2 };
   check([['newton_raphson_iterate', [cubic, 3, 1, 4], { root }], ['halley_iterate', [cubic, 3, 1, 4], { root }], ['schroder_iterate', [cubic, 3, 1, 4], { root }],
      ['newton_raphson_iterate', [cubic, 3, 1, 4, few], two], ['halley_iterate', [cubic, 3, 1, 4, few], two], ['schroder_iterate', [cubic, 3, 1, 4, few], two],
      ['toms748_solve', [f, 1, 4], bracket], ['bisect', [f, 1, 4], bracket], ['toms748_solve', [f, 1, 4, few], two], ['bisect', [f, 1, 4, few], two],
      ['multivariate_newton', [circle, [1, 2]], { root: both }], ['multivariate_newton', [flat, [1, 2], 53], { root: both }],
      ['multivariate_newton', [circle, [1, 2], 53, 1], { root: [1.5, 1.5], iterations: 1 }], ['multivariate_halley', [circle, [1, 2]], { root: both }],
      ['multivariate_halley', [circle, [1, 2], 53], { root: both }], ['multivariate_halley', [circle, [1, 2], 53, 1], { root: [10 / 7, 1.4], iterations: 1 }]]);
   // 2 - 3 z + z^2 = (z - 1)(z - 2) as (re, im) pairs; z^3 = 1 from a 3 by 3 grid as (root, iterations) pairs: the centre is singular.
   for (const roots of [bm.polynomial_roots, bm.polynomial_roots_companion_matrix]) same(roots([2, -3, 1]).sort(), [0, 0, 1, 2]);
   const [newton, halley] = [false, true].map((h) => bm.polynomial_basins([-1, 0, 0, 1], h, 3, 2));
   same([newton.roots, [0, 4, 8, 12, 16].map((i) => halley.basin[i])], [[-0.5, -Math.sqrt(0.75), 1, 0, -0.5, Math.sqrt(0.75)], [2, 1, 255, 0, 1]]);
   assert.ok(halley.basin[1] < newton.basin[1] && bm.polynomial_basins([-1, 0, 0, 1], false, 2, 2, few).basin[1] === 2);
   rejects(() => bm.halley_iterate((x) => [x, 1], 3, 1, 4), () => bm.bisect(Math.cos, 0, 2, { digits: 100 }), () => bm.polynomial_roots([0, 0]),
           () => bm.polynomial_basins([1], false, 3, 2), () => bm.multivariate_halley(flat, [1, 2]), () => bm.multivariate_newton((x) => [[x[0]], [[1, 0]]], [1, 2]),
           () => bm.multivariate_newton(circle, [1, 2], 0), () => bm.multivariate_newton(circle, [1, 2], 53, Infinity));
});

const lower = [-2, -1], upper = [2, 3], rosenbrock = (x) => (1 - x[0]) * (1 - x[0]) + 100 * (x[1] - x[0] * x[0]) * (x[1] - x[0] * x[0]);
// The optimizers, and a count among the options of each.
const optimizers = { differential_evolution: 'max_generations', jso: 'max_function_evaluations', random_search: 'max_function_calls', cma_es: 'max_generations',
                     multilevel_coordinate_search: 'max_function_calls', dycors: 'max_function_calls', turbo: 'max_function_calls' };

test('test functions and optimizers', () => {
   check([['test_function', ['rosenbrock', [-1.5, 2]], 12.5], ['test_function', ['fonseca_fleming', [0, 0, 0]], [1 - Math.exp(-1), 1 - Math.exp(-1)]],
      ['test_function_array', ['himmelblau', [[3, 2], [0, 0]]], [0, 170]], ['test_function_array', ['zdt1', [new Array(30).fill(0)]], [[0, 1]]]]);
   for (const [name, key] of Object.entries(optimizers)) {
      const options = { seed: 1, [key]: 100 }, r = bm[name]('rosenbrock', lower, upper, options), best = Math.min(...r.queries.filter((_, i) => i % 3 === 2));
      assert.ok(r.x instanceof Float64Array && r.value === best && rosenbrock(r.x) === best && best < 25 && bm[name]('rosenbrock', lower, upper).value < 1, name);
      // The same arithmetic in a callback gives the same run, and so does the same seed. Counts are integers.
      assert.deepEqual(bm[name](rosenbrock, lower, upper, options), r, name);
      rejects(() => bm[name]('rosenbrock', lower, upper, { nope: 1 }), () => bm[name]('rosenbrock', upper, lower), () => bm[name](7, lower, upper),
              () => bm[name]('rosenbrock', lower, upper, { [key]: Infinity }), () => bm[name]('rosenbrock', lower, upper, { [key]: 1.5 }));
   }
   const { queries } = bm.jso(rosenbrock, lower, upper, { target_value: 30, initial_guess: [1, 1.5] });
   assert.deepEqual([queries.length, ...queries.subarray(0, 3)], [3 * 39, 1, 1.5, 25]);  // the guess first, and only its 39 initial points at the target
   rejects(() => bm.jso(rosenbrock, lower, upper, { max_function_evaluations: 38 }), () => bm.random_search(rosenbrock, lower, upper, { seed: -1 }),
           () => bm.cma_es(rosenbrock, lower, upper, { seed: 2 ** 60 }), () => bm.dycors(rosenbrock, lower, upper, { initial_guess: [9, 9] }),
           () => bm.test_function('nope', [1, 2]), () => bm.test_function('rosenbrock', [1, 2, 3]), () => bm.test_function_array('zdt1', [[1]]));
});

test('gradient optimizers', () => {
   const names = ['gradient_descent', 'nesterov', 'lbfgs', 'trust_region_newton'], start = [-1.5, 2];
   for (const [i, options] of [{ learning_rate: 0.001 }, { learning_rate: 0.001, momentum: 0.9 }, { history: 10 }, { initial_radius: 1 }].entries()) {
      const name = names[i], o = bm[name]('rosenbrock', start, options), newton = i === 3;
      o.step();
      const at = newton ? o.arguments() : start, [x, y] = at;  // the value and gradients are at the point stepped from, or to for the trust region
      same([o.objective_value(), ...o.gradients()], [rosenbrock(at), 2 * (x - 1) - 400 * x * (y - x * x), 200 * (y - x * x)]);
      const last = o.steps(20).subarray(95);
      assert.deepEqual([...last], [...o.arguments(), rosenbrock(o.arguments()), o.function_calls(), o.hessian_vector_products()]);
      assert.ok(last[2] < 12.5 && last[3] >= 20 && (last[4] > 0) === newton && (o.radius() > 0) === newton, name);
      // They share the tapes of the autodiff: another is refused while this one lives and, like a bad number of steps, leaves it as it was.
      for (const other of names) assert.throws(() => bm[other]('rosenbrock', start), /delete\(\)/);
      rejects(() => o.steps(Infinity), () => o.steps(-1), () => o.steps(1.5), () => o.steps(2 ** 31));
      assert.deepEqual([...o.arguments(), o.function_calls()], [last[0], last[1], last[3]]);
      o.delete();
      rejects(() => bm[name]((v) => v[0], [0, 0]), () => bm[name]('rosenbrock', [0, 0, 0]), () => bm[name]('rosenbrock', [0, 0], { nope: 1 }));
      bm[name]('rosenbrock', start).delete();
   }
});

const schaffer = (x) => [x[0] ** 2, (x[0] - 2) ** 2], fonseca = (x) => [...bm.test_function('fonseca_fleming', x)], plain = (rows) => rows.map((row) => [...row]);
// A front {x, f}: rows of points, and of the objectives at each.
const check_front = ({ x, f }, objectives) => assert.deepEqual([x.length > 0 && x[0] instanceof Float64Array, plain(f)], [true, x.map(objectives)]);
// Asks and tells twice, the objectives as rows and end to end: a wrong number of either leaves the optimizer as it was.
function ask_and_tell(o, points, objectives, first) {
   rejects(() => o.tell([[1, 2]]));  // nothing was asked for
   const x = o.ask();
   rejects(() => o.tell(x.slice(1).map(objectives)), () => o.tell(x.map(() => [])), () => o.tell(new Float64Array(15)), () => o.tell(7));
   assert.deepEqual([plain(x.slice(0, 2)), o.ask()], [first, x]);
   o.tell(x.map(objectives));
   o.tell(Float64Array.from(o.ask().flatMap(objectives)));
   check_front(o.front(), objectives);
   check_front({ x: o[points](), f: o.objectives() }, objectives);
   return o;
}

test('nsga2 and qlognehvi', () => {
   const one = [[-10], [10]], three = [[-4, -4, -4], [4, 4, 4]], few = { max_evaluations: 24, batch_size: 2, seed: 1 };
   for (const [name, objectives, [lo, hi], options] of [['nsga2', schaffer, one, { max_generations: 20, seed: 3 }], ['qlognehvi', fonseca, three, few]]) {
      const front = bm[name](objectives, lo, new Float64Array(hi), options), optimizer = bm[name + '_optimizer'];
      check_front(front, objectives);
      // An error caught inside the objectives leaves the run as it was, which is that of the seed.
      assert.deepEqual(bm[name]((x) => { assert.throws(() => bm.mean([]), RangeError); return objectives(x); }, lo, hi, options), front);
      rejects(() => bm[name](objectives, hi, lo), () => bm[name](objectives, [], [], { initial_population: [1] }), () => bm[name](objectives, lo, hi, 7),
              () => bm[name](7, lo, hi), () => bm[name](() => [], lo, hi), () => bm[name]((x) => (x[0] > 0 ? [1, 2] : [1]), lo, hi), () => new optimizer(hi, lo),
              () => new optimizer(lo, hi, { [Object.keys(options)[0]]: 10 }), () => bm[name](objectives, lo, hi, { initial_population: [lo.map((v) => 3 * v)] }));
   }
   assert.ok(bm.nsga2(schaffer, ...one).x.every((p) => p[0] > -0.1 && p[0] < 2.1));  // the Pareto set is 0 <= x <= 2
   assert.deepEqual(bm.qlognehvi('fonseca_fleming', ...three, few), bm.qlognehvi(fonseca, ...three, few));
   const first = [[1, 1, 1], [2, 2, 2]], q = new bm.qlognehvi_optimizer(...three, { initial_samples: 8, seed: 1, initial_population: first.flat() });
   const n = ask_and_tell(new bm.nsga2_optimizer(...one, { population_size: 8, seed: 3, initial_population: [[1], [2]] }), 'population', schaffer, [[1], [2]]);
   ask_and_tell(q, 'points', fonseca, first);
   // Results are copies, which outlive the memory they came from: it grows here.
   const ranks = n.ranks(), y = bm.erf_array([0.5]), size = bm.HEAPU8.length, defaults = [new bm.nsga2_optimizer(...one), new bm.qlognehvi_optimizer(...three)];
   assert.ok(bm.erf_array(new Float64Array(size / 8)).length === size / 8 && bm.HEAPU8.length > size && y[0] === bm.erf(0.5));
   assert.deepEqual([ranks instanceof Uint32Array, ranks.length, ranks[0], defaults[0].ask().length, q.ask(3).length, q.reference_point().length], [true, 8, 0, 100, 3, 2]);
   check([['hypervolume', [[[0, 1], [1, 0], [0.5, 0.5]], [1, 1]], 0.25], ['hypervolume', [[0, 1, 1, 0, 0.5, 0.5, 0.25, 0.25], [1, 1]], 0.5625]]);
   rejects(() => bm.hypervolume([[0, 0]], []), () => bm.hypervolume([[0, 0, 0]], [1, 1]), () => bm.nsga2(schaffer, ...one, { crossover_probability: 2 }),
           () => bm.qlognehvi(fonseca, ...three, { reference_point: [1] }), () => q.ask(0), () => q.ask(1.5), () => q.ask(Infinity));
   [n, q, ...defaults].forEach((o) => o.delete());
});

// The arguments after the callback of each function that takes one, without options.
const callers = Object.fromEntries([...rules.map(([rule, [, ...rest]]) => [rule, rest]), ['toms748_solve', [0, 1]], ['bisect', [0, 1]], ['multivariate_newton', [[1, 2]]],
   ['newton_raphson_iterate', [3, 1, 4]], ...[...Object.keys(optimizers), 'nsga2', 'qlognehvi'].map((name) => [name, [lower, upper]])]);

test('what a callback throws passes through', () => {
   const foreign = new WebAssembly.Exception(new WebAssembly.Tag({ parameters: [] }), []), error = new TypeError('mine');
   for (const [name, args] of Object.entries(callers)) {
      const call = (g) => bm[name](g, ...args);
      for (const thrown of [error, 'a string', { not: 'an error' }, foreign, undefined]) assert.throws(() => call(() => { throw thrown; }), (e) => e === thrown, name);
      assert.throws(() => call(() => bm.mean([])), RangeError, name);  // the error of a call into the module, and below one caught; no promises
      assert.throws(() => call(() => { try { bm.mean([]); } catch (e) { throw new EvalError(e.message); } }), EvalError, name);
      assert.throws(() => call(async () => 1), /synchronous/, name);
   }
});

// The stack pointer and the size of the memory, before and after 10^5 calls that fail: an exception that left the stack where
// it was thrown, or its allocations in place, would show. A callback that catches the errors of its own calls leaves the stack alone.
test('errors leave the stack and the heap as they were', { timeout: 600000 }, () => {
   const big = new Float64Array(10000), lo = new Array(10).fill(0), hi = new Array(10).fill(1), failing = new TypeError('mine'), fft = new bm.fft(8);
   let count = 0, depth;
   // A bad argument, a bad option, a constructor, a method, the library, a callback, an objective, objectives, and nested calls.
   const cases = [() => bm.mean([]), () => bm.nsga2(() => [0, 0], lo, hi, { population_size: 3, initial_population: big }), () => new bm.fft(0),
      () => fft.transform(big), () => bm.pchip(big, big, big), () => bm.tanh_sinh(() => { throw failing; }, 0, 1),
      () => bm.differential_evolution((x) => { if (++count % 7 === 0) throw failing; return x[0]; }, lo, hi, { NP: 8, initial_guess: lo }),
      () => bm.nsga2((x) => { if (++count % 5 === 0) throw 'a string'; return [x[0], x[1]]; }, lo, hi, { population_size: 8, initial_population: lo }),
      () => bm.gauss((x) => bm.gauss(() => bm.catmull_rom([[0, 0], [1, 1], [2, 0], [1, -1]], [NaN])[0], 0, x).result, 0, 1),
      () => bm.gauss((x) => { try { bm.mean([]); } catch (e) { new bm.fft(0); } return x; }, 0, 1)];
   const run = (n) => { for (const call of cases) for (let i = 0; i < n; ++i) assert.throws(call); };
   run(100);
   const sp = bm.stackSave(), memory = bm.HEAPU8.length;
   run(100000);
   assert.deepEqual([bm.stackSave(), bm.HEAPU8.length], [sp, memory]);
   bm.gauss((x) => { for (let i = 0; i < 100000; ++i) assert.throws(cases[0]); depth = bm.stackSave(); return x; }, 0, 1, { points: 7 });
   assert.ok(sp - depth < 4096 && bm.stackSave() === sp, `${sp} ${depth}`);
   fft.delete();
});
