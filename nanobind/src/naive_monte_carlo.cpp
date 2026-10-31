// Copyright 2026 Nick Thompson
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// quadrature.naive_monte_carlo. The C++ calls the integrand from its own threads, so each call takes the GIL;
// with the GIL serialising them, more than one thread gains nothing, and one is used.

#include <boost/math/quadrature/naive_monte_carlo.hpp>

#include <nanobind/stl/pair.h>
#include <nanobind/stl/vector.h>

#include <cmath>
#include <cstdint>
#include <exception>
#include <future>
#include <limits>
#include <mutex>
#include <utility>
#include <vector>

#include "bindings.hpp"

namespace {

namespace nb = nanobind;

// The constructor is noexcept and calls the integrand, so an exception from Python is kept here, NaN returned in
// its place, and the exception rethrown afterwards.
struct py_integrand {
   nb::callable f;
   std::exception_ptr* error;
   std::mutex* lock;
   double operator()(std::vector<double> const& x) const {
      nb::gil_scoped_acquire gil;
      try {
         return nb::cast<double>(f(bmpy::vector_to_array(x)));
      } catch (...) {
         std::lock_guard<std::mutex> g(*lock);
         if (!*error)
            *error = std::current_exception();
         return std::numeric_limits<double>::quiet_NaN();
      }
   }
};

}  // namespace

void bmpy::bind_naive_monte_carlo(nb::module_& m) {
   nb::object result_type = bmpy::named_tuple(m, "naive_monte_carlo_result", "value error calls");
   m.def("naive_monte_carlo", [result_type](nb::callable f, std::vector<std::pair<double, double>> const& bounds, double error_goal,
                                 bool singular, std::uint64_t seed) {
      for (auto const& [a, b] : bounds)
         if (!(a < b))
            throw nb::value_error("each bound (a, b) needs a < b");
      if (bounds.empty())
         throw nb::value_error("bounds is empty");
      if (!(error_goal > 0))
         throw nb::value_error("error_goal must be positive");
      std::exception_ptr error;
      std::mutex lock;
      py_integrand g{f, &error, &lock};
      double value = std::numeric_limits<double>::quiet_NaN(), estimate;
      std::uint64_t calls;
      {
         boost::math::quadrature::naive_monte_carlo<double, py_integrand> mc(g, bounds, error_goal, singular, 1, seed);
         if (!error) {
            nb::gil_scoped_release release;
            try {
               value = mc.integrate().get();
            } catch (...) {
               // A NaN from a failed call is reported as a domain error; prefer the Python exception behind it.
               if (!error)
                  error = std::current_exception();
            }
         }
         estimate = mc.current_error_estimate();
         calls = mc.calls();
      }
      if (error)
         std::rethrow_exception(error);
      return result_type(value, estimate, calls);
   }, nb::arg("f"), nb::arg("bounds"), nb::arg("error_goal"), nb::kw_only(), nb::arg("singular") = true, nb::arg("seed") = 0,
      nb::sig("def naive_monte_carlo(f, bounds, error_goal, /, *, singular=True, seed=0)"),
      "(value, error, calls): the Monte Carlo integral of f(x), x a 1D array, over the box of (a, b) bounds, either "
      "maybe infinite, until the error estimate is below error_goal; seed=0 seeds from std::random_device."
      BMPY_URL("naive_monte_carlo"));
}
