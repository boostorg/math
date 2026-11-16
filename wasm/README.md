# Boost.Math for WebAssembly

Special functions, distributions, FFTs, statistics, filters, interpolators, quadrature, root finding and optimizers
in double precision, built with [Emscripten](https://emscripten.org):

```sh
emcmake cmake -S wasm -B build-wasm && cmake --build build-wasm
BOOST_MATH_WASM=build-wasm/boost_math.mjs node --test wasm/test/smoke.mjs
```

```js
import createModule from './boost_math.mjs';
const bm = await createModule();
bm.tgamma(4.5);                      // 11.631728396567452
bm.tgamma_array([0.5, 1.5]);         // Float64Array [1.7724538509055159, 0.8862269254527579]
const t = new bm.students_t_distribution(5);
t.quantile(0.975);                   // 2.5705818356363146
t.delete();                          // embind objects live in the wasm heap
```

Arrays in (`Float64Array`, `Float32Array` or plain) give new `Float64Array`s out, which need no freeing. Most
functions have an `_array` twin over their last argument, and every distribution has `pdf_array`, `cdf_array`,
`ccdf_array` and `quantile_array`. Complex data are (re, im) pairs, and transforms of N-dimensional arrays, stored
row major, run along one axis:

```js
bm.power_spectrum(x);                         // P[0..n/2], with sum P = mean(x^2)
bm.real_fft_axis(image, [rows, cols], -1);    // each row; shape [rows, floor(cols/2) + 1]; fft_axis, inverse_fft_axis
bm.one_sample_t_test(x, 0);                   // {statistic, p_value}; also z, two sample, paired, runs, Ljung-Box
bm.simple_ordinary_least_squares(x, y);       // {intercept, slope, r_squared}
bm.trend_filter(y, 1, 10);                    // piecewise linear fit; trend_filter(y) picks k and lambda,
bm.trend_filter_select(y);                    // and this says which: {k, lambda, score}
bm.wavelet_denoise(y);                        // wavelet and rule by cross-validation (wavelet_denoise_select says
bm.wavelet_denoise(y, "symlet", 8, "sure");   // which); wavelet_transform(y, "daubechies", 4) and its inverse
bm.makima(x, y, t);                           // the interpolant of (x, y) at each t; pchip, barycentric_rational,
bm.catmull_rom(points, s);                    // cardinal_cubic_b_spline; a curve through the points, at each s
bm.daubechies_wavelet_array(4, t);            // daubechies_scaling; daubechies_wavelet_transform(y, dt, p, s, t)
```

The Daubechies functions interpolate a grid of `2^grid_refinements` points per unit, an option of each
(`bm.daubechies_scaling(2, x, { grid_refinements: 16 })`): 14 by default for p ≤ 5, 12 for p = 6 and 7, and 10 above.
The error against the exact function falls with the refinement, faster for the smoother functions of larger p:

| p  | 10 refinements | 14     | 16     |
|----|----------------|--------|--------|
| 2  | 1e-2           | 2e-3   | 9e-4   |
| 3  | 6e-4           | 3e-5   | 4e-6   |
| 4  | 2e-5           | 2e-7   | 2e-8   |
| 6  | 7e-8           | 1e-10  | 5e-12  |
| 10 | 1e-11          | 2e-15  | 2e-15  |

A grid holds 2 to 4 numbers at each of its (2p - 1) 2^refinements points, and one of more than 2^23 numbers is
refused: the C++ default, up to 2^21 points per unit, takes hundreds of megabytes.

A function to integrate, solve or minimize is a callback, and the last argument is an optional object of options.
Callbacks are synchronous: one that returns a promise is refused with a `TypeError`.

```js
bm.tanh_sinh(x => Math.log(x) * Math.sqrt(x), 0, 1);            // {result, error, L1, levels}; exp_sinh, sinh_sinh(f)
bm.gauss_kronrod(f, a, b, { points: 61, tolerance: 1e-10 });    // {result, error, L1}; gauss(f, a, b, { points: 30 })
bm.leja_jacobi(f, a, b, { alpha: 0, beta: -0.5 });              // nested Leja rules: leja_legendre, _hermite, _laguerre
bm.gauss_hermite(f, { points: 16, averaged: true });            // and gauss_laguerre, gauss_jacobi; averaged adds an error
bm.ooura_fourier_sin(x => 1 / x, omega);                        // {result, relative_error}; ooura_fourier_cos
bm.sparse_grid(([x, y]) => Math.exp(x * y), [-1, -1], [1, 1]);  // any dimension; { max_calls, tolerance }
bm.halley_iterate(x => [f(x), df(x), d2f(x)], guess, min, max); // {root, iterations}; newton_raphson_iterate, schroder_iterate
bm.toms748_solve(f, a, b);                                      // {lower, upper, iterations}; bisect
bm.polynomial_roots([-1, -1, 0, 0, 0, 1]);                      // x^5 - x - 1; polynomial_roots_companion_matrix
bm.multivariate_halley(x => [F(x), J(x), v => d2(x, v)], x0);   // {root, iterations}; multivariate_newton
bm.jso(x => f(x), lower, upper, { seed: 1 });                   // {x, value, queries}; differential_evolution, cma_es,
                                                                // random_search, multilevel_coordinate_search, turbo, dycors
bm.nsga2(x => [f1(x), f2(x)], lower, upper, { seed: 1 });       // {x, f}: the Pareto front; nsga2_optimizer has ask and tell
bm.qlognehvi(x => [f1(x), f2(x)], lower, upper, { seed: 1 });   // the same, for expensive objectives; qlognehvi_optimizer
bm.hypervolume(front.f, reference);                             // the volume the front dominates, below the reference point
bm.pareto_polish("zdt1", front, lower, upper);                  // {x, f, origin, cells, status, ...}: the front of a test problem,
                                                                // polished and refined with derivatives; { refinement_tolerance }
```

The options of an optimizer are the members of its parameters in C++ (`max_function_calls`, `NP`, ...) and `seed`,
`initial_guess` and `target_value`; `queries` holds every evaluation, as rows (x, f(x)) end to end. The module has
no threads: a long run blocks, so call it from a Worker. An objective can also be the name of a test problem
(`wasm/test_problems.hpp`), which runs without a callback: `rosenbrock`, `rastrigin`, `ackley`, `himmelblau` and
`styblinski_tang` of two variables, and `zdt1`, ..., `zdt4`, `zdt6`, `dtlz2`, `fonseca_fleming`, `kursawe`, `two_sites`
and `three_sites` of several objectives (`bm.test_function(name, x)` evaluates one). The gradient-based optimizers and
`pareto_polish` differentiate the objective by reverse-mode autodiff, so they take only those:
`bm.lbfgs("rosenbrock", [-1.5, 2])`, and `gradient_descent`,
`nesterov` and `trust_region_newton`, give an object with `step()`, `steps(n)`, `arguments()`, `objective_value()`
and `gradients()`. They share the tapes of the autodiff, so there is one at a time: making another before `delete()`
is called on the first throws.

Errors follow one rule:

- A function of numbers (the special functions, their `_array` twins, and the members of a distribution) returns NaN
  or infinity for a domain error, as `Math` does, and never throws.
- Every other call throws a `RangeError`, with a message, for an argument that is wrong whatever the data: a type,
  size, shape, count, option, name (of a wavelet family, a rule, a test problem), a wavelet order, the exponent of a
  weight, too few samples for a filter, and whatever the library throws in C++ (an optimizer's
  parameters, an interpolator's abscissas, a parameter outside a curve). Counts are integers: `Infinity`, `NaN`, 1.5
  and -1 are not. NaN in what such a call returns is a result, not a rejected argument: a solver that did not
  converge, an interpolant evaluated outside its data.
- What a callback throws reaches the caller unchanged, whatever it is, and the module can be used afterwards.

C++ exceptions are WebAssembly exceptions, which sets the oldest browsers the module loads in: Chrome 95, Firefox 100,
Safari 15.2, and Node 17.

`-DBOOST_MATH_WASM_RELAXED_SIMD=ON` gives the double-double arithmetic a hardware fma, but relaxed SIMD needs
Chrome or Edge 114+ or Firefox 146+, and Safari can't load it.
