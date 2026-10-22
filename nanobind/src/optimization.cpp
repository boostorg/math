// Copyright 2026 Nick Thompson
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// The derivative-free optimizers (differential_evolution, jso, random_search, multilevel_coordinate_search, and cma_es
// where Eigen is found) on std::vector<double>,
// minimizing a Python callable f(x) -> float over the box [lower_bounds, upper_bounds]. Each returns the
// best argument found. The C++ calls f from its own worker threads, so the GIL is released for the run and
// taken for each call; an exception raised by f cancels the run and is re-raised when it returns.
// The gradient-based optimizers are left out: their objectives are written over reverse-mode autodiff
// variables, which a Python callable can't provide.

#include "bindings.hpp"

#include <boost/math/optimization/differential_evolution.hpp>
#include <boost/math/optimization/jso.hpp>
#include <boost/math/optimization/multilevel_coordinate_search.hpp>
#include <boost/math/optimization/random_search.hpp>
#if __has_include(<Eigen/Dense>)
#include <boost/math/optimization/cma_es.hpp>
#define BOOST_MATH_NB_HAS_CMA_ES
#endif

#include <nanobind/stl/optional.h>
#include <nanobind/stl/vector.h>

#include <atomic>
#include <cstdint>
#include <limits>
#include <mutex>
#include <optional>
#include <random>
#include <stdexcept>
#include <vector>

namespace bmpy {

namespace {

namespace bo = boost::math::optimization;
using vector = std::vector<double>;

// The cost function the optimizers call, from any of their threads.
class python_cost
{
public:
   explicit python_cost(nb::callable f) : f_(std::move(f)) {}

   double operator()(vector const& x) const
   {
      if (cancel_.load())
      {
         return std::numeric_limits<double>::infinity();
      }
      nb::gil_scoped_acquire acquire;
      try
      {
         return nb::cast<double>(f_(x));
      }
      catch (nb::python_error& e)
      {
         // The error indicator is per thread, so keep the exception itself to raise on the caller's thread.
         std::lock_guard<std::mutex> lock(mutex_);
         if (!error_)
         {
            error_.emplace(e);
         }
         cancel_.store(true);
         return std::numeric_limits<double>::infinity();
      }
   }

   std::atomic<bool>* cancellation() const { return &cancel_; }

   // Raises the first exception f raised, if any. Call with the GIL held.
   void rethrow() const
   {
      if (error_)
      {
         throw nb::python_error(*error_);
      }
   }

private:
   nb::callable f_;
   mutable std::atomic<bool> cancel_{false};
   mutable std::mutex mutex_;
   mutable std::optional<nb::python_error> error_;
};

void check_bounds(vector const& lower, vector const& upper, char const* name)
{
   if (lower.size() != upper.size() || lower.empty())
   {
      throw nb::value_error((std::string(name) + ": lower_bounds and upper_bounds must be non-empty and the same length").c_str());
   }
}

// Runs the optimizer with the GIL released, then re-raises the first exception from the cost function.
template <class Run>
vector run_released(python_cost const& cost, Run run)
{
   vector x;
   {
      nb::gil_scoped_release release;
      x = run();
   }
   cost.rethrow();
   return x;
}

}  // namespace

void bind_optimization(nb::module_& m)
{
   using nb::arg;
   const double no_target = std::numeric_limits<double>::quiet_NaN();

   m.def("differential_evolution", [](nb::callable f, vector lower_bounds, vector upper_bounds, std::optional<vector> initial_guess,
                                      std::size_t population, std::size_t max_generations, double mutation_factor,
                                      double crossover_probability, unsigned threads, double target_value, std::uint64_t seed) {
      check_bounds(lower_bounds, upper_bounds, "differential_evolution");
      bo::differential_evolution_parameters<vector> params;
      params.lower_bounds = std::move(lower_bounds);
      params.upper_bounds = std::move(upper_bounds);
      params.initial_guess = initial_guess ? &*initial_guess : nullptr;
      params.NP = population;
      params.max_generations = max_generations;
      params.mutation_factor = mutation_factor;
      params.crossover_probability = crossover_probability;
      params.threads = threads;
      python_cost cost(f);
      std::mt19937_64 gen(seed);
      return run_released(cost, [&] { return bo::differential_evolution(std::cref(cost), params, gen, target_value, cost.cancellation()); });
   }, arg("f"), arg("lower_bounds"), arg("upper_bounds"), arg("initial_guess") = nb::none(), arg("population") = 500,
      arg("max_generations") = 1000, arg("mutation_factor") = 0.65, arg("crossover_probability") = 0.5, arg("threads") = 1u,
      arg("target_value") = no_target, arg("seed") = 0,
      "optimization::differential_evolution: the argument in the box that minimizes f");

   m.def("jso", [](nb::callable f, vector lower_bounds, vector upper_bounds, std::optional<vector> initial_guess,
                   std::size_t initial_population_size, std::size_t max_function_evaluations, std::size_t threads,
                   double target_value, std::uint64_t seed) {
      check_bounds(lower_bounds, upper_bounds, "jso");
      bo::jso_parameters<vector> params;
      params.lower_bounds = std::move(lower_bounds);
      params.upper_bounds = std::move(upper_bounds);
      params.initial_guess = initial_guess ? &*initial_guess : nullptr;
      params.initial_population_size = initial_population_size;
      params.max_function_evaluations = max_function_evaluations;
      params.threads = threads;
      python_cost cost(f);
      std::mt19937_64 gen(seed);
      return run_released(cost, [&] { return bo::jso(std::cref(cost), params, gen, target_value, cost.cancellation()); });
   }, arg("f"), arg("lower_bounds"), arg("upper_bounds"), arg("initial_guess") = nb::none(), arg("initial_population_size") = 0,
      arg("max_function_evaluations") = 0, arg("threads") = 1, arg("target_value") = no_target, arg("seed") = 0,
      "optimization::jso: the argument in the box that minimizes f (0 means the defaults from the paper)");

   m.def("random_search", [](nb::callable f, vector lower_bounds, vector upper_bounds, std::optional<vector> initial_guess,
                             std::size_t max_function_calls, unsigned threads, double target_value, std::uint64_t seed) {
      check_bounds(lower_bounds, upper_bounds, "random_search");
      bo::random_search_parameters<vector> params;
      params.lower_bounds = std::move(lower_bounds);
      params.upper_bounds = std::move(upper_bounds);
      params.initial_guess = initial_guess ? &*initial_guess : nullptr;
      params.max_function_calls = max_function_calls;
      params.threads = threads;
      python_cost cost(f);
      std::mt19937_64 gen(seed);
      return run_released(cost, [&] { return bo::random_search(std::cref(cost), params, gen, target_value, cost.cancellation()); });
   }, arg("f"), arg("lower_bounds"), arg("upper_bounds"), arg("initial_guess") = nb::none(), arg("max_function_calls") = 10000,
      arg("threads") = 1u, arg("target_value") = no_target, arg("seed") = 0,
      "optimization::random_search: the best argument found in the box");

   m.def("multilevel_coordinate_search", [](nb::callable f, vector lower_bounds, vector upper_bounds,
                                            std::optional<vector> initial_guess, std::size_t max_function_calls,
                                            std::size_t max_level, bool local_search, std::size_t max_local_search_steps,
                                            std::size_t line_search_points, double target_value) {
      check_bounds(lower_bounds, upper_bounds, "multilevel_coordinate_search");
      bo::multilevel_coordinate_search_parameters<vector> params;
      params.lower_bounds = std::move(lower_bounds);
      params.upper_bounds = std::move(upper_bounds);
      params.initial_guess = initial_guess ? &*initial_guess : nullptr;
      params.max_function_calls = max_function_calls;
      params.max_level = max_level;
      params.local_search = local_search;
      params.max_local_search_steps = max_local_search_steps;
      params.line_search_points = line_search_points;
      python_cost cost(f);
      return run_released(cost, [&] { return bo::multilevel_coordinate_search(std::cref(cost), params, target_value, cost.cancellation()); });
   }, arg("f"), arg("lower_bounds"), arg("upper_bounds"), arg("initial_guess") = nb::none(), arg("max_function_calls") = 0,
      arg("max_level") = 0, arg("local_search") = true, arg("max_local_search_steps") = 50, arg("line_search_points") = 15,
      arg("target_value") = no_target,
      "optimization::multilevel_coordinate_search: the argument in the box that minimizes f (0 means the defaults from the paper)");

#ifdef BOOST_MATH_NB_HAS_CMA_ES
   m.def("cma_es", [](nb::callable f, vector lower_bounds, vector upper_bounds, std::optional<vector> initial_guess,
                      std::size_t max_generations, std::size_t population_size, double learning_rate, double target_value,
                      std::uint64_t seed) {
      check_bounds(lower_bounds, upper_bounds, "cma_es");
      bo::cma_es_parameters<vector> params;
      params.lower_bounds = std::move(lower_bounds);
      params.upper_bounds = std::move(upper_bounds);
      params.initial_guess = initial_guess ? &*initial_guess : nullptr;
      params.max_generations = max_generations;
      params.population_size = population_size;
      params.learning_rate = learning_rate;
      python_cost cost(f);
      std::mt19937_64 gen(seed);
      return run_released(cost, [&] { return bo::cma_es(std::cref(cost), params, gen, target_value, cost.cancellation()); });
   }, arg("f"), arg("lower_bounds"), arg("upper_bounds"), arg("initial_guess") = nb::none(), arg("max_generations") = 1000,
      arg("population_size") = 0, arg("learning_rate") = 1.0, arg("target_value") = no_target, arg("seed") = 0,
      "optimization::cma_es: the argument in the box that minimizes f");
#endif
}

}  // namespace bmpy
