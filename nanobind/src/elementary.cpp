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
   BMPY_R1(log1p, "x");
   BMPY_R1(log1pmx, "x");
   BMPY_R1(expm1, "x");
   BMPY_R1(cbrt, "x");
   BMPY_R1(sqrt1pm1, "x");
   BMPY_R2(powm1, "x", "y");
   BMPY_R2(pow1p, "x", "y");
   BMPY_R1(sin_pi, "x");
   BMPY_R1(cos_pi, "x");
   BMPY_R1(sinc_pi, "x");
   BMPY_R1(sinhc_pi, "x");
   BMPY_R1(asinh, "x");
   BMPY_R1(acosh, "x");
   BMPY_R1(atanh, "x");
   BMPY_R2(hypot, "x", "y");
   BMPY_R1(rsqrt, "x");
   BMPY_R1(logit, "p");
   BMPY_R1(logistic_sigmoid, "x");
   BMPY_R2(logaddexp, "x", "y");

   // Floating-point navigation; float32 inputs work in float.
   BMPY_R1(float_next, "x");
   BMPY_R1(float_prior, "x");
   BMPY_R2(float_distance, "a", "b");
   BMPY_R1I1(float_advance, int, "x", "n");
   BMPY_R2(nextafter, "x", "y");
   BMPY_R1(ulp, "x");
   BMPY_R2(relative_difference, "a", "b");
   BMPY_R2(epsilon_difference, "a", "b");

   bmpy::def_vec<bmpy::Int<unsigned long long>>(
       m, "fibonacci", []<class T>(unsigned long long n) { return bm::fibonacci<std::uint64_t>(n); }, {"n"},
       "boost::math::fibonacci<std::uint64_t>");
   bmpy::def_vec<bmpy::Int<unsigned>>(
       m, "prime", []<class T>(unsigned n) { return static_cast<std::uint64_t>(bm::prime(n)); }, {"n"}, "boost::math::prime");

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
   }, nb::arg("upper"), "boost::math::prime_sieve: the primes in [2, upper)");
   m.def("prime_range", [to_array](std::uint64_t lower, std::uint64_t upper) {
      std::vector<std::uint64_t> primes;
      {
         nb::gil_scoped_release release;
         bm::prime_range(lower, upper, primes);
      }
      return to_array(std::move(primes));
   }, nb::arg("lower"), nb::arg("upper"), "boost::math::prime_range: the primes in [lower, upper)");
   m.def("prime_count", [](std::uint64_t upper) {
      nb::gil_scoped_release release;
      return bm::prime_count(bm::execution::cuda, upper);
   }, nb::arg("upper"), "boost::math::prime_count: the number of primes in [2, upper)");
   m.def("prime_count", [](std::uint64_t lower, std::uint64_t upper) {
      nb::gil_scoped_release release;
      return bm::prime_count(bm::execution::cuda, lower, upper);
   }, nb::arg("lower"), nb::arg("upper"), "boost::math::prime_count: the number of primes in [lower, upper)");
}
