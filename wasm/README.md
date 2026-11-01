# Boost.Math for WebAssembly

Special functions, distributions and FFTs in double precision, built with [Emscripten](https://emscripten.org):

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

Errors return NaN or infinity rather than throwing. `-DBOOST_MATH_WASM_RELAXED_SIMD=ON` gives the double-double
arithmetic a hardware fma, but relaxed SIMD needs Chrome or Edge 114+ or Firefox 146+, and Safari can't load it.
