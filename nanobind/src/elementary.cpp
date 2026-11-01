// Copyright 2026 Nick Thompson
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/math/special_functions/log1p.hpp>
#include <boost/math/special_functions/expm1.hpp>
#include <boost/math/special_functions/cbrt.hpp>
#include <boost/math/special_functions/sqrt1pm1.hpp>
#include <boost/math/special_functions/powm1.hpp>
#include <boost/math/special_functions/pow1p.hpp>
#include <boost/math/special_functions/sin_pi.hpp>
#include <boost/math/special_functions/cos_pi.hpp>
#include <boost/math/special_functions/sinc.hpp>
#include <boost/math/special_functions/sinhc.hpp>
#include <boost/math/special_functions/asinh.hpp>
#include <boost/math/special_functions/acosh.hpp>
#include <boost/math/special_functions/atanh.hpp>
#include <boost/math/special_functions/hypot.hpp>
#include <boost/math/special_functions/rsqrt.hpp>
#include <boost/math/special_functions/logit.hpp>
#include <boost/math/special_functions/logistic_sigmoid.hpp>
#include <boost/math/special_functions/logaddexp.hpp>
#include <boost/math/special_functions/next.hpp>
#include <boost/math/special_functions/ulp.hpp>
#include <boost/math/special_functions/relative_difference.hpp>
#include <boost/math/special_functions/fibonacci.hpp>
#include <boost/math/special_functions/prime.hpp>
#include <boost/math/special_functions/prime_sieve.hpp>
#include <cstdint>
#include <vector>

#include "bindings.hpp"
#include "macros.hpp"

void bmpy::bind_elementary(nb::module_& m) {
   BMPY_R1(log1p, "x", "log(1 + x), accurate for small x." BMPY_URL("powers/log1p"));
   BMPY_R1(log1pmx, "x", "log(1 + x) − x, accurate for small x." BMPY_URL("powers/log1p"));
   BMPY_R1(expm1, "x", "eˣ − 1, accurate for small x." BMPY_URL("powers/expm1"));
   BMPY_R1(cbrt, "x", "The real cube root ∛x, negative for x < 0." BMPY_URL("powers/cbrt"));
   BMPY_R1(sqrt1pm1, "x", "√(1 + x) − 1, accurate for small x." BMPY_URL("powers/sqrt1pm1"));
   BMPY_R2(powm1, "x", "y", "xʸ − 1, accurate when y log(x) is small." BMPY_URL("powers/powm1"));
   BMPY_R2(pow1p, "x", "y", "(1 + x)ʸ, accurate for small x." BMPY_URL("powers/pow1p"));
   BMPY_R1(sin_pi, "x", "sin(πx), exact at the integers." BMPY_URL("powers/sin_pi"));
   BMPY_R1(cos_pi, "x", "cos(πx), exact at the half-integers." BMPY_URL("powers/cos_pi"));
   BMPY_R1(sinc_pi, "x", "sin(x)/x, 1 at 0.\n\nUnnormalised: numpy.sinc(x) is sinc_pi(πx)." BMPY_URL("sinc/sinc_pi"));
   BMPY_R1(sinhc_pi, "x", "sinh(x)/x, 1 at 0." BMPY_URL("sinc/sinhc_pi"));
   BMPY_R1(asinh, "x", "The inverse hyperbolic sine, log(x + √(x² + 1))." BMPY_URL("inv_hyper/asinh"));
   BMPY_R1(acosh, "x", "The inverse hyperbolic cosine, log(x + √(x² − 1)), for x >= 1." BMPY_URL("inv_hyper/acosh"));
   BMPY_R1(atanh, "x", "The inverse hyperbolic tangent, ½ log((1 + x)/(1 − x)), for |x| < 1; ±1 overflow."
      BMPY_URL("inv_hyper/atanh"));
   BMPY_R2(hypot, "x", "y", "√(x² + y²), without intermediate overflow or underflow." BMPY_URL("powers/hypot"));
   BMPY_R1(rsqrt, "x", "1/√x." BMPY_URL("powers/rsqrt"));
   BMPY_R1(logit, "p", "The logit log(p/(1 − p)), for p in [0, 1]; the inverse of logistic_sigmoid." BMPY_URL("logistic/logit"));
   BMPY_R1(logistic_sigmoid, "x", "The logistic sigmoid 1/(1 + e^(−x)); as scipy.special.expit."
      BMPY_URL("logistic/logistic_sigmoid"));
   BMPY_R2(logaddexp, "x", "y", "log(eˣ + eʸ), without overflow; as numpy.logaddexp." BMPY_URL("powers/logaddexp"));

   // Floating-point navigation; float32 inputs work in float.
   BMPY_R1(float_next, "x", "The next representable value above x." BMPY_URL("next_float/float_next"));
   BMPY_R1(float_prior, "x", "The next representable value below x." BMPY_URL("next_float/float_prior"));
   BMPY_R2(float_distance, "a", "b", "The number of representable values from a to b, signed: positive when b > a."
      BMPY_URL("next_float/float_distance"));
   BMPY_R1I1(float_advance, int, "x", "n", "x moved n representable values (towards −∞ for n < 0)."
      BMPY_URL("next_float/float_advance"));
   BMPY_R2(nextafter, "x", "y", "The next representable value after x in the direction of y." BMPY_URL("next_float/nextafter"));
   BMPY_R1(ulp, "x", "The unit in the last place of x: the distance from |x| to the next value away from 0."
      BMPY_URL("next_float/ulp"));
   BMPY_R2(relative_difference, "a", "b", "The relative difference |a − b|/min(|a|, |b|), symmetric in a and b."
      BMPY_URL("float_comparison"));
   BMPY_R2(epsilon_difference, "a", "b", "relative_difference(a, b) in units of machine epsilon." BMPY_URL("float_comparison"));

   bmpy::def_vec<bmpy::Int<unsigned long long>>(
       m, "fibonacci", []<class T>(unsigned long long n) { return bm::fibonacci<std::uint64_t>(n); }, {"n"},
       "The Fibonacci number Fₙ, F₀ = 0, F₁ = 1, as an integer; n > 93 overflows uint64."
      BMPY_URL("number_series/fibonacci_numbers"));
   bmpy::def_vec<bmpy::Int<unsigned>>(
       m, "prime", []<class T>(unsigned n) { return static_cast<std::uint64_t>(bm::prime(n)); }, {"n"},
       "The nth prime, counting from prime(0) = 2; n < 10000." BMPY_URL("number_series/primes"));

   // The segmented sieve; results are uint64 arrays. prime_count needs an execution policy, and
   // libc++ may have no <execution>; without nvcc the cuda policy runs on the CPU.
   auto to_array = [](std::vector<std::uint64_t>&& primes) {
      auto* data = new std::vector<std::uint64_t>(std::move(primes));
      nb::capsule owner(data, [](void* p) noexcept { delete static_cast<std::vector<std::uint64_t>*>(p); });
      return nb::ndarray<nb::numpy, std::uint64_t, nb::ndim<1>>(data->data(), {data->size()}, owner);
   };
   m.def("prime_sieve", [to_array](std::uint64_t upper) {
      std::vector<std::uint64_t> primes;
      {
         nb::gil_scoped_release release;
         bm::prime_sieve(upper, primes);
      }
      return to_array(std::move(primes));
   }, nb::arg("upper"), "The primes in [2, upper), as a uint64 array, by a segmented sieve."
      BMPY_URL("number_series/prime_sieve"));
   m.def("prime_range", [to_array](std::uint64_t lower, std::uint64_t upper) {
      std::vector<std::uint64_t> primes;
      {
         nb::gil_scoped_release release;
         bm::prime_range(lower, upper, primes);
      }
      return to_array(std::move(primes));
   }, nb::arg("lower"), nb::arg("upper"), "The primes in [lower, upper), as a uint64 array."
      BMPY_URL("number_series/prime_sieve"));
   m.def("prime_count", [](std::uint64_t upper) {
      nb::gil_scoped_release release;
      return bm::prime_count(bm::execution::cuda, upper);
   }, nb::arg("upper"), "The number of primes in [2, upper), π(upper − 1)." BMPY_URL("number_series/prime_sieve"));
   m.def("prime_count", [](std::uint64_t lower, std::uint64_t upper) {
      nb::gil_scoped_release release;
      return bm::prime_count(bm::execution::cuda, lower, upper);
   }, nb::arg("lower"), nb::arg("upper"), "The number of primes in [lower, upper).");
}
