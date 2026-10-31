// Copyright 2026 Nick Thompson
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// The optimizers, on std::vector<double>, minimizing a Python callable f(x) -> float of a 1D array, and Brent's
// one-dimensional minimizer. Each returns an OptimizeResult. The derivative-free ones search the box
// [lower_bounds, upper_bounds] and call f from their own worker threads, so the GIL is released for the run and
// taken for each call; an exception raised by f cancels the run and is re-raised when it returns. The gradient
// ones (gradient_descent, nesterov, lbfgs) take the gradient from Python, through their policy parameters, in
// place of the reverse-mode autodiff of the C++ defaults; they run with the GIL held.

#include "bindings.hpp"

#include <boost/math/optimization/cma_es.hpp>
#include <boost/math/optimization/differential_evolution.hpp>
#include <boost/math/optimization/jso.hpp>
#include <boost/math/optimization/minimizer.hpp>
#include <boost/math/optimization/multilevel_coordinate_search.hpp>
#include <boost/math/optimization/random_search.hpp>
#include <boost/math/optimization/turbo.hpp>
#include <boost/math/tools/minima.hpp>

#include <nanobind/stl/optional.h>
#include <nanobind/stl/vector.h>

#include <atomic>
#include <cstdint>
#include <exception>
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

// boost.math.optimize.OptimizeResult, made in bind_optimization; the module keeps it alive.
nb::handle optimize_result;

// The cost function the derivative-free optimizers call, from any of their threads. It counts the calls and keeps
// the least value, which is f at the argument returned: the optimizers return the best argument they evaluated.
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
         double v = nb::cast<double>(f_(vector_to_array(x)));
         // The GIL serializes these.
         ++calls_;
         best_ = v < best_ ? v : best_;
         return v;
      }
      catch (...)
      {
         // The error indicator is per thread, so keep the exception itself to raise on the caller's thread.
         std::lock_guard<std::mutex> lock(mutex_);
         if (!error_)
         {
            error_ = std::current_exception();
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
         std::rethrow_exception(error_);
      }
   }

   std::size_t calls() const { return calls_; }
   double best() const { return best_; }

private:
   nb::callable f_;
   mutable std::atomic<bool> cancel_{false};
   mutable std::mutex mutex_;
   mutable std::exception_ptr error_;
   mutable std::size_t calls_ = 0;
   mutable double best_ = std::numeric_limits<double>::infinity();
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
nb::object run_released(python_cost const& cost, Run run)
{
   vector x;
   {
      nb::gil_scoped_release release;
      x = run();
   }
   cost.rethrow();
   return optimize_result(vector_to_array(x), cost.best(), cost.calls());
}

// f with its gradient: grad(x), or, when grad is None, f(x) returns (value, gradient).
struct differentiable
{
   nb::callable f;
   nb::object grad;
   std::size_t nfev = 0;

   // f(x), and the gradient there in *g.
   double call(vector const& x, vector* g)
   {
      ++nfev;
      nb::object r = f(vector_to_array(x)), gr;
      if (grad.is_none())
      {
         gr = r[1];
         r = r[0];
      }
      else if (g)
      {
         gr = grad(vector_to_array(x));
      }
      double v = nb::cast<double>(r);
      if (!g)
      {
         return v;
      }
      nb::object a = numpy().attr("asarray")(gr, nb::arg("dtype") = numpy().attr("float64"));
      auto view = nb::cast<nb::ndarray<const double, nb::ndim<1>, nb::device::cpu>>(a);
      if (view.shape(0) != x.size())
      {
         throw nb::value_error("the gradient must be a 1D array of the length of x");
      }
      for (std::size_t i = 0; i < x.size(); ++i)
      {
         (*g)[i] = view(i);
      }
      return v;
   }
};

// The policies of the gradient optimizers, for a differentiable*.
struct no_initialization
{
   void operator()(vector&) const {}
};

struct value_policy
{
   double operator()(differentiable* d, vector const& x) const { return d->call(x, nullptr); }
};

struct gradient_policy
{
   template <class ValuePolicy>
   void operator()(differentiable* d, vector const& x, ValuePolicy&, double& v, vector& g) const
   {
      g.resize(x.size());
      v = d->call(x, &g);
   }
};

// minimize(opt) until the gradient's 2-norm is below tolerance, or max_iterations steps.
template <class Optimizer>
nb::object run_gradient(Optimizer& opt, differentiable& d, vector const& x, double tolerance, std::size_t max_iterations)
{
   auto r = bo::minimize(opt, bo::unconstrained_policy<vector>{}, bo::gradient_norm_convergence_policy<double>(tolerance),
                         bo::max_iter_termination_policy(max_iterations), false);
   // objective_value() is f before the last step.
   double fun = d.call(x, nullptr);
   return optimize_result(vector_to_array(x), fun, d.nfev, r.num_iter, r.converged);
}

#define BMPY_GRADIENT_DOC \
   " grad(x) is the gradient, or with grad None, f(x) returns (value, gradient). Stops when |gradient| < tolerance."

}  // namespace

void bind_optimization(nb::module_& m)
{
   using nb::arg;
   const double no_target = std::numeric_limits<double>::quiet_NaN();

   optimize_result = named_tuple(m, "OptimizeResult", "x fun nfev nit converged",
      "The best x found, fun = f(x), the calls of f, and the iterations and whether the stopping test was met (None "
      "for the derivative-free optimizers).", nb::make_tuple(nb::none(), nb::none()));

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
   }, arg(), arg(), arg(), nb::kw_only(), arg("initial_guess") = nb::none(), arg("population") = 500,
      arg("max_generations") = 1000, arg("mutation_factor") = 0.65, arg("crossover_probability") = 0.5, arg("threads") = 1u,
      arg("target_value") = no_target, arg("seed") = 0,
      nb::sig("def differential_evolution(f, lower_bounds, upper_bounds, /, *, initial_guess=None, population=500, max_generations=1000, "
              "mutation_factor=0.65, crossover_probability=0.5, threads=1, target_value=nan, seed=0)"),
      "Minimizes f(x), x a 1D array, over the box [lower_bounds, upper_bounds] by differential evolution."
      BMPY_URL("differential_evolution"));

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
   }, arg(), arg(), arg(), nb::kw_only(), arg("initial_guess") = nb::none(), arg("initial_population_size") = 0,
      arg("max_function_evaluations") = 0, arg("threads") = 1, arg("target_value") = no_target, arg("seed") = 0,
      nb::sig("def jso(f, lower_bounds, upper_bounds, /, *, initial_guess=None, initial_population_size=0, "
              "max_function_evaluations=0, threads=1, target_value=nan, seed=0)"),
      "Minimizes f over the box by jSO, a self-adapting differential evolution; 0 means the paper's defaults." BMPY_URL("jso"));

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
   }, arg(), arg(), arg(), nb::kw_only(), arg("initial_guess") = nb::none(), arg("max_function_calls") = 10000,
      arg("threads") = 1u, arg("target_value") = no_target, arg("seed") = 0,
      nb::sig("def random_search(f, lower_bounds, upper_bounds, /, *, initial_guess=None, max_function_calls=10000, threads=1, "
              "target_value=nan, seed=0)"),
      "Minimizes f over the box by evaluating it at uniformly random points." BMPY_URL("random_search"));

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
   }, arg(), arg(), arg(), nb::kw_only(), arg("initial_guess") = nb::none(), arg("max_function_calls") = 0,
      arg("max_level") = 0, arg("local_search") = true, arg("max_local_search_steps") = 50, arg("line_search_points") = 15,
      arg("target_value") = no_target,
      nb::sig("def multilevel_coordinate_search(f, lower_bounds, upper_bounds, /, *, initial_guess=None, max_function_calls=0, "
              "max_level=0, local_search=True, max_local_search_steps=50, line_search_points=15, target_value=nan)"),
      "Minimizes f over the box by multilevel coordinate search; 0 means the paper's defaults."
      BMPY_URL("multilevel_coordinate_search"));

   m.def("turbo", [](nb::callable f, vector lower_bounds, vector upper_bounds, std::optional<vector> initial_guess,
                     std::size_t max_function_calls, std::size_t initial_points, std::size_t batch_size, std::size_t candidates,
                     std::size_t max_training_points, std::size_t hyperparameter_steps, double target_value, std::uint64_t seed) {
      check_bounds(lower_bounds, upper_bounds, "turbo");
      bo::turbo_parameters<vector> params;
      params.lower_bounds = std::move(lower_bounds);
      params.upper_bounds = std::move(upper_bounds);
      params.initial_guess = initial_guess ? &*initial_guess : nullptr;
      params.max_function_calls = max_function_calls;
      params.initial_points = initial_points;
      params.batch_size = batch_size;
      params.candidates = candidates;
      params.max_training_points = max_training_points;
      params.hyperparameter_steps = hyperparameter_steps;
      python_cost cost(f);
      std::mt19937_64 gen(seed);
      return run_released(cost, [&] { return bo::turbo(std::cref(cost), params, gen, target_value, cost.cancellation()); });
   }, arg(), arg(), arg(), nb::kw_only(), arg("initial_guess") = nb::none(), arg("max_function_calls") = 0,
      arg("initial_points") = 0, arg("batch_size") = 1, arg("candidates") = 0, arg("max_training_points") = 200,
      arg("hyperparameter_steps") = 50, arg("target_value") = no_target, arg("seed") = 0,
      nb::sig("def turbo(f, lower_bounds, upper_bounds, /, *, initial_guess=None, max_function_calls=0, initial_points=0, "
              "batch_size=1, candidates=0, max_training_points=200, hyperparameter_steps=50, target_value=nan, seed=0)"),
      "Minimizes an expensive f over the box by TuRBO, trust-region Bayesian optimization; 0 means the defaults."
      BMPY_URL("turbo"));

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
   }, arg(), arg(), arg(), nb::kw_only(), arg("initial_guess") = nb::none(), arg("max_generations") = 1000,
      arg("population_size") = 0, arg("learning_rate") = 1.0, arg("target_value") = no_target, arg("seed") = 0,
      nb::sig("def cma_es(f, lower_bounds, upper_bounds, /, *, initial_guess=None, max_generations=1000, population_size=0, "
              "learning_rate=1.0, target_value=nan, seed=0)"),
      "Minimizes f over the box by the CMA evolution strategy; population_size 0 means 4 + floor(3 log(dimension))."
      BMPY_URL("cma_es"));

   m.def("gradient_descent", [](nb::callable f, vector x, nb::object grad, double learning_rate, double tolerance, std::size_t max_iterations) {
      differentiable d{f, grad};
      auto opt = bo::make_gradient_descent(&d, x, learning_rate, no_initialization{}, value_policy{}, gradient_policy{});
      return run_gradient(opt, d, x, tolerance, max_iterations);
   }, arg(), arg(), arg("grad") = nb::none(), nb::kw_only(), arg("learning_rate") = 0.01, arg("tolerance") = 1e-3,
      arg("max_iterations") = 100000,
      nb::sig("def gradient_descent(f, x0, /, grad=None, *, learning_rate=0.01, tolerance=0.001, max_iterations=100000)"),
      "Minimizes f from x0 by gradient descent." BMPY_GRADIENT_DOC
      BMPY_URL("gd_opt/gradient_descent"));

   m.def("nesterov", [](nb::callable f, vector x, nb::object grad, double learning_rate, double momentum, double tolerance,
                        std::size_t max_iterations) {
      differentiable d{f, grad};
      auto opt = bo::make_nag(&d, x, learning_rate, momentum, no_initialization{}, value_policy{}, gradient_policy{});
      return run_gradient(opt, d, x, tolerance, max_iterations);
   }, arg(), arg(), arg("grad") = nb::none(), nb::kw_only(), arg("learning_rate") = 0.01, arg("momentum") = 0.95,
      arg("tolerance") = 1e-3, arg("max_iterations") = 100000,
      nb::sig("def nesterov(f, x0, /, grad=None, *, learning_rate=0.01, momentum=0.95, tolerance=0.001, max_iterations=100000)"),
      "Minimizes f from x0 by Nesterov accelerated gradient descent." BMPY_GRADIENT_DOC
      BMPY_URL("gd_opt/nesterov"));

   m.def("lbfgs", [](nb::callable f, vector x, nb::object grad, std::size_t m, double tolerance, std::size_t max_iterations) {
      differentiable d{f, grad};
      auto opt = bo::make_lbfgs(&d, x, m, no_initialization{}, value_policy{}, gradient_policy{},
                                bo::strong_wolfe_line_search_policy<double>{});
      return run_gradient(opt, d, x, tolerance, max_iterations);
   }, arg(), arg(), arg("grad") = nb::none(), nb::kw_only(), arg("m") = 10, arg("tolerance") = 1e-3, arg("max_iterations") = 100000,
      nb::sig("def lbfgs(f, x0, /, grad=None, *, m=10, tolerance=0.001, max_iterations=100000)"),
      "Minimizes f from x0 by L-BFGS, with memory m and a strong Wolfe line search." BMPY_GRADIENT_DOC
      BMPY_URL("gd_opt/lbfgs"));

   m.def("brent", [](nb::callable f, double a, double b, int bits, std::uintmax_t max_iterations) {
      if (!(a < b))
      {
         throw nb::value_error("brent: a must be less than b");
      }
      std::size_t nfev = 0;
      auto g = [&](double x) {
         ++nfev;
         return nb::cast<double>(f(x));
      };
      std::uintmax_t it = max_iterations;
      auto [x, fx] = boost::math::tools::brent_find_minima(g, a, b, bits, it);
      return optimize_result(x, fx, nfev, it, it < max_iterations);
   }, arg(), arg(), arg(), nb::kw_only(), arg("bits") = 26, arg("max_iterations") = 500,
      nb::sig("def brent(f, a, b, /, *, bits=26, max_iterations=500)"),
      "Minimizes f of one float over [a, b] by Brent's method, to about `bits` bits (at most 26); x is a float."
      BMPY_URL("brent_minima"));
}

}  // namespace bmpy
