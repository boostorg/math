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
