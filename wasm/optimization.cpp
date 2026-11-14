// Copyright 2026 Nick Thompson
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// Embind bindings for the optimizers. The module has no threads, so each runs on the caller's, and a long run blocks
// it: call it from a Worker.

#include "bindings.hpp"

#include <boost/math/optimization/cma_es.hpp>
#include <boost/math/optimization/detail/hypervolume.hpp>
#include <boost/math/optimization/differential_evolution.hpp>
#include <boost/math/optimization/dycors.hpp>
#include <boost/math/optimization/gradient_optimizers.hpp>
#include <boost/math/optimization/jso.hpp>
#include <boost/math/optimization/multilevel_coordinate_search.hpp>
#include <boost/math/optimization/nsga2.hpp>
#include <boost/math/optimization/qlognehvi.hpp>
#include <boost/math/optimization/random_search.hpp>
#include <boost/math/optimization/turbo.hpp>

#include "test_problems.hpp"

#include <emscripten/em_js.h>

#include <cstdint>
#include <functional>
#include <memory>
#include <random>
#include <utility>

using emscripten::EM_VAL;
namespace bo = bm::optimization;
namespace tp = test_problems;
using vec = std::vector<double>;

// The test problems of test_problems.hpp, by name: an objective can be one of these instead of a callback.
namespace named
{
template <class Names>
int find(const Names& names, const val& f)
{
   return f.isString() ? tp::find(names, f.as<std::string>()) : -1;
}

void check(int id, std::size_t d)
{
   if (d != 2)
      range_error(std::string(tp::single[id]) + " is a function of two variables.");
}

vec values(int id, const vec& x)
{
   if (x.size() < tp::variables(id))
      range_error(std::string(tp::multiple[id]) + " does not take " + std::to_string(x.size()) + " variables.");
   return tp::values(id, x);
}

// A number, or the array of the objectives.
val at(const val& name, const vec& x)
{
   const int id = find(tp::single, name);
   if (id >= 0)
      return check(id, x.size()), val(tp::value(id, x));
   if (find(tp::multiple, name) < 0)
      range_error("Unknown test function.");
   return to_js(values(find(tp::multiple, name), x));
}
} // namespace named

// f(x) for a callback f. Callbacks are synchronous: a promise, or anything else with a then, is refused. What f throws
// unwinds the optimizer and reaches the caller unchanged.
EM_JS_DEPS(boost_math_objective, "$Emval");
EM_JS(double, call_objective, (EM_VAL f, EM_VAL x), {
   const y = Emval.toValue(f)(Emval.toValue(x));
   if (typeof y === 'object' && y !== null && typeof y.then === 'function')
      throw new TypeError('The callback returned a promise: callbacks must be synchronous.');
   return y;
});

// An objective is a callback from a Float64Array to a number, or the name of a test function.
struct objective
{
   val f;
   int id;
   objective(const val& f_, std::size_t d) : f(f_), id(named::find(tp::single, f_))
   {
      if (id >= 0)
         named::check(id, d);
      else if (f.typeOf().as<std::string>() != "function")
         range_error("The objective must be a function or the name of a test function.");
   }
   double operator()(const vec& x) const { return id >= 0 ? tp::value(id, x) : call_objective(f.as_handle(), to_js(x).as_handle()); }
};

// Several objectives: a callback returning an array, the same length for every x, or a test function.
struct objectives
{
   val f;
   int id;
   std::shared_ptr<std::size_t> count = std::make_shared<std::size_t>(0);
   explicit objectives(const val& f_) : f(f_), id(named::find(tp::multiple, f_))
   {
      if (id < 0 && f.typeOf().as<std::string>() != "function")
         range_error("The objectives must be a function or the name of a test function.");
   }
   vec operator()(const vec& x) const
   {
      vec y = id < 0 ? numbers(sync(f(to_js(x)))) : named::values(id, x);
      if (y.empty() || (*count && y.size() != *count))
         range_error("f(x) must return the same number of objectives, at least one, for every x.");
      *count = y.size();
      return y;
   }
};

// f, the bounds and the options in; {x, value, queries} out, where queries holds every evaluation in order, as rows
// (x_1, ..., x_d, f(x)) end to end. Every optimizer takes the options seed (0), initial_guess and target_value, at or
// below which it stops; set(p, o) reads its own, named as the members of its parameters, and run calls it.
template <class Params, class Set, class Run>
val minimize(const val& f, const val& lower, const val& upper, const val& opts, Set set, Run run)
{
   options o(opts);
   Params p;
   p.lower_bounds = numbers(lower);
   p.upper_bounds = numbers(upper);
   vec guess;
   if (o.has("initial_guess"))
   {
      guess = numbers(o.get("initial_guess"));
      p.initial_guess = &guess;
   }
   std::mt19937_64 gen(o.seed());
   const double target = o.real("target_value", NAN);
   set(p, o);
   o.done();
   const objective cost(f, p.lower_bounds.size());
   std::vector<std::pair<vec, double>> queries;
   const vec x = run(cost, p, gen, target, &queries);
   vec flat;
   double best = NAN;
   for (const auto& q : queries)
   {
      flat.insert(flat.end(), q.first.begin(), q.first.end());
      flat.push_back(q.second);
      best = q.second < best || std::isnan(best) ? q.second : best;
   }
   val r = val::object();
   r.set("x", to_js(x));
   r.set("value", best);
   r.set("queries", to_js(flat));
   return r;
}

#define OPTIMIZER(name, set, call) \
   entry(#name, +[](const val& f, const val& l, const val& u) { return minimize<bo::name##_parameters<vec>>(f, l, u, val::undefined(), set, call); }, \
         +[](const val& f, const val& l, const val& u, const val& o) { return minimize<bo::name##_parameters<vec>>(f, l, u, o, set, call); })
#define RUN(name) [](const objective& f, auto& p, auto& gen, double target, auto* queries) { \
   return bo::name(f, p, gen, target, nullptr, nullptr, queries); }

// NSGA-II: {x, f}, the points of a front and their objectives, as arrays of rows.
namespace nsga2_js
{
val front(const std::vector<std::pair<vec, vec>>& members)
{
   std::vector<vec> x, f;
   for (const auto& p : members)
   {
      x.push_back(p.first);
      f.push_back(p.second);
   }
   val o = val::object();
   o.set("x", rows(x));
   o.set("f", rows(f));
   return o;
}

// Rows come in as arrays, or end to end in one array of numbers, which is then cut into rows of c > 0, or into r > 0
// rows.
std::vector<vec> matrix(const val& a, std::size_t r, std::size_t c)
{
   std::vector<vec> m;
   const std::size_t length = a["length"].isNumber() ? a["length"].as<std::size_t>() : 0;
   if (length > 0 && a[0].isNumber())
   {
      const vec v = numbers(a);
      c = c ? c : v.size() / r;
      for (std::size_t i = 0; c && i < v.size(); i += c)
         m.emplace_back(v.begin() + i, v.begin() + std::min(i + c, v.size()));
   }
   else
      for (std::size_t i = 0; i < length; ++i)
         m.push_back(numbers(a[i]));
   return m;
}

struct settings
{
   bo::nsga2_parameters<vec> params;
   std::vector<vec> initial;  // params.initial_population is set where the parameters are used
   std::uint64_t seed;
};

// max_generations is an option of the function only: the caller of ask and tell decides when to stop.
settings make(const val& lower, const val& upper, const val& opts, bool run)
{
   options o(opts);
   settings s;
   auto& p = s.params;
   p.threads = 1;
   p.lower_bounds = numbers(lower);
   p.upper_bounds = numbers(upper);
   // The rows of a flat initial population are cut by the number of variables.
   if (p.lower_bounds.empty() || p.lower_bounds.size() != p.upper_bounds.size())
      range_error("There must be as many lower as upper bounds, and at least one.");
   p.population_size = o.count("population_size", p.population_size);
   check_size(double(p.population_size) * double(p.lower_bounds.size()), "The population");
   if (run)
      p.max_generations = o.count("max_generations", p.max_generations);
   s.seed = o.seed();
   p.crossover_probability = o.real("crossover_probability", p.crossover_probability);
   p.mutation_probability = o.real("mutation_probability", p.mutation_probability);
   p.crossover_distribution_index = o.real("crossover_distribution_index", p.crossover_distribution_index);
   p.mutation_distribution_index = o.real("mutation_distribution_index", p.mutation_distribution_index);
   if (o.has("initial_population"))
      s.initial = matrix(o.get("initial_population"), 0, p.lower_bounds.size());
   o.done();
   return s;
}

val minimize(const val& f, const val& lower, const val& upper, const val& opts)
{
   settings s = make(lower, upper, opts, true);
   s.params.initial_population = &s.initial;
   std::mt19937_64 gen(s.seed);
   return front(bo::nsga2(objectives(f), s.params, gen));
}

struct optimizer
{
   bo::nsga2_optimizer<vec, vec> o;
   std::size_t n;
   optimizer(settings s) : o((s.params.initial_population = &s.initial, s.params), std::mt19937_64(s.seed)), n(s.params.population_size) {}
};
} // namespace nsga2_js

// qLogNEHVI, as NSGA-II: the options are the members of qlognehvi_parameters, and seed.
namespace qlognehvi_js
{
struct settings
{
   bo::qlognehvi_parameters<vec> params;
   std::vector<vec> initial;
   std::uint64_t seed;
};

// max_evaluations is an option of the function only: the caller of ask and tell decides when to stop.
settings make(const val& lower, const val& upper, const val& opts, bool run)
{
   options o(opts);
   settings s;
   auto& p = s.params;
   p.threads = 1;
   p.lower_bounds = numbers(lower);
   p.upper_bounds = numbers(upper);
   if (p.lower_bounds.empty() || p.lower_bounds.size() != p.upper_bounds.size())
      range_error("There must be as many lower as upper bounds, and at least one.");
   p.initial_samples = o.count("initial_samples", p.initial_samples);
   p.batch_size = o.count("batch_size", p.batch_size);
   if (run)
      p.max_evaluations = o.count("max_evaluations", p.max_evaluations);
   p.mc_samples = o.count("mc_samples", p.mc_samples);
   p.restarts = o.count("restarts", p.restarts);
   p.raw_samples = o.count("raw_samples", p.raw_samples);
   check_size(double(std::max({p.initial_samples, p.batch_size, p.raw_samples})) * double(p.lower_bounds.size()), "The points");
   if (o.has("reference_point"))
      p.reference_point = numbers(o.get("reference_point"));
   if (o.has("initial_population"))
      s.initial = nsga2_js::matrix(o.get("initial_population"), 0, p.lower_bounds.size());
   s.seed = o.seed();
   o.done();
   return s;
}

val minimize(const val& f, const val& lower, const val& upper, const val& opts)
{
   settings s = make(lower, upper, opts, true);
   s.params.initial_population = &s.initial;
   std::mt19937_64 gen(s.seed);
   return nsga2_js::front(bo::qlognehvi(objectives(f), s.params, gen));
}

struct optimizer
{
   bo::qlognehvi_optimizer<vec, vec> o;
   std::size_t d, asked = 0;  // the variables, and the points last asked for
   optimizer(settings s) : o((s.params.initial_population = &s.initial, s.params), std::mt19937_64(s.seed)), d(s.params.lower_bounds.size()) {}

   val ask(std::size_t n)
   {
      check_size(double(n) * double(d), "The points");
      const auto& points = n ? o.ask(n) : o.ask();
      asked = points.size();
      return rows(points);
   }
   void tell(const val& values)
   {
      if (!asked)
         range_error("tell takes the objectives of the points last asked for.");
      o.tell(nsga2_js::matrix(values, asked, 0));
      asked = 0;
   }
};

// The volume which the points dominate, below the reference point.
double hypervolume(const val& points, const val& reference)
{
   const vec r = numbers(reference);
   if (r.empty())
      range_error("The reference point must have a number for each objective.");
   const std::vector<vec> m = nsga2_js::matrix(points, 0, r.size());
   for (const vec& p : m)
      if (p.size() != r.size())
         range_error("Every point must have as many objectives as the reference point.");
   return bo::detail::hypervolume(m, r);
}
} // namespace qlognehvi_js

// The gradient-based optimizers, one step at a time: arguments(), objective_value() and gradients() are those of the
// last step, and function_calls() and hessian_vector_products() count from the start. steps(n) takes n steps and
// returns a row after each, (x_1, ..., x_d, f(x), function calls, Hessian-vector products), end to end. The tapes of
// the autodiff belong to the thread, not to an optimizer, and making an optimizer clears them: so there is one
// optimizer at a time, and making another before delete() is called on the first throws.
namespace gradient
{
namespace rd = bm::differentiation::reverse_mode;

struct counted
{
   int id;
   std::shared_ptr<std::size_t> calls;
   template <class X>
   typename X::value_type operator()(const X& x) const
   {
      ++*calls;
      return tp::value(id, x);
   }
};

bool& alive()
{
   static bool one = false;
   return one;
}

struct counted_products : bo::reverse_mode_hessian_vector_product_policy<double>
{
   std::shared_ptr<std::size_t> products;
   template <class... A>
   void operator()(A&&... a)
   {
      ++*products;
      bo::reverse_mode_hessian_vector_product_policy<double>::operator()(std::forward<A>(a)...);
   }
};

struct optimizer
{
   std::function<void()> step;
   std::function<vec()> arguments, gradients;
   std::function<double()> objective_value, radius;
   std::shared_ptr<std::size_t> calls = std::make_shared<std::size_t>(0), products = std::make_shared<std::size_t>(0);
   int id;
   std::shared_ptr<void> life;  // the last copy to go lets another optimizer be made

   val steps(double n)
   {
      const std::size_t taken = ::count(n, "The number of steps");
      check_size(double(taken) * double(arguments().size() + 3), "The path");
      vec path;
      for (std::size_t k = 0; k < taken; ++k)
      {
         step();
         const vec x = arguments();
         path.insert(path.end(), x.begin(), x.end());
         path.insert(path.end(), {tp::value(id, x), double(*calls), double(*products)});
      }
      return to_js(path);
   }
};

// Only trust_region_newton has a radius.
template <class O>
auto radius_of(const O& o, int) -> decltype(o.radius()) { return o.radius(); }
template <class O>
double radius_of(const O&, long) { return NAN; }

// make(f, x, o) gives the optimizer of f at x, a vector of rvar<double, Order>, reading its options from o.
template <std::size_t Order, class Make>
optimizer create(const val& f, const val& x0, const val& opts, Make make)
{
   const vec start = numbers(x0);
   const int id = named::find(tp::single, f);
   if (id < 0)
      range_error("The objective must be the name of a test function.");
   named::check(id, start.size());
   if (alive())
      range_error("Another gradient optimizer is alive, and they share the tapes of the autodiff: call delete() on it first.");
   bo::for_each_tape<double>([](auto& tape) { tape.clear(); }, std::integral_constant<std::size_t, 2>{});
   optimizer r;
   r.id = id;
   options o(opts);
   auto x = std::make_shared<std::vector<rd::rvar<double, Order>>>(start.begin(), start.end());
   auto opt = std::make_shared<decltype(make(counted{id, r.calls}, *x, o, r.products))>(make(counted{id, r.calls}, *x, o, r.products));
   o.done();
   r.step = [opt, x] { opt->step(); };
   r.arguments = [x] {
      vec v;
      for (const auto& xi : *x)
         v.push_back(xi.item());
      return v;
   };
   r.gradients = [opt] { return opt->gradients(); };
   r.objective_value = [opt] { return opt->objective_value(); };
   r.radius = [opt] { return radius_of(*opt, 0); };
   r.life = std::shared_ptr<void>(nullptr, [](void*) { alive() = false; });
   alive() = true;
   return r;
}

template <std::size_t Order, class Make>
void bind(const char* name, Make make)
{
   static const Make m = make;
   entry(name, +[](const val& f, const val& x0) { return create<Order>(f, x0, val::undefined(), m); },
         +[](const val& f, const val& x0, const val& o) { return create<Order>(f, x0, o, m); });
}
} // namespace gradient

EMSCRIPTEN_BINDINGS(boost_math_optimization)
{
   entry("test_function", +[](const val& name, const val& x) { return named::at(name, numbers(x)); });
   entry("test_function_array", +[](const val& name, const val& xs) {
      const std::size_t n = xs["length"].as<std::size_t>();
      if (named::find(tp::single, name) < 0)
      {
         val a = val::array();
         for (std::size_t i = 0; i < n; ++i)
            a.call<void>("push", named::at(name, numbers(xs[i])));
         return a;
      }
      vec y(n);
      for (std::size_t i = 0; i < n; ++i)
         y[i] = named::at(name, numbers(xs[i])).as<double>();
      return to_js(y);
   });

   OPTIMIZER(differential_evolution, [](auto& p, options& o) {
      p.threads = 1;
      p.NP = o.count("NP", p.NP);
      p.max_generations = o.count("max_generations", p.max_generations);
      p.mutation_factor = o.real("mutation_factor", p.mutation_factor);
      p.crossover_probability = o.real("crossover_probability", p.crossover_probability);
   }, RUN(differential_evolution));
   OPTIMIZER(jso, [](auto& p, options& o) {
      p.threads = 1;
      p.initial_population_size = o.count("initial_population_size", 0);
      p.max_function_evaluations = o.count("max_function_evaluations", 0);
   }, RUN(jso));
   OPTIMIZER(random_search, [](auto& p, options& o) {
      p.threads = 1;
      p.max_function_calls = o.count("max_function_calls", 10000);
   }, RUN(random_search));
   OPTIMIZER(cma_es, [](auto& p, options& o) {
      p.max_generations = o.count("max_generations", p.max_generations);
      p.population_size = o.count("population_size", 0);
      p.learning_rate = o.real("learning_rate", p.learning_rate);
   }, RUN(cma_es));
   OPTIMIZER(multilevel_coordinate_search, [](auto& p, options& o) {
      p.max_function_calls = o.count("max_function_calls", 0);
      p.max_level = o.count("max_level", 0);
      p.local_search = o.flag("local_search", true);
      p.max_local_search_steps = o.count("max_local_search_steps", p.max_local_search_steps);
      p.line_search_points = o.count("line_search_points", p.line_search_points);
   }, [](const objective& f, auto& p, auto&, double target, auto* queries) {
      return bo::multilevel_coordinate_search(f, p, target, nullptr, nullptr, queries);
   });
   OPTIMIZER(turbo, [](auto& p, options& o) {
      p.max_function_calls = o.count("max_function_calls", 0);
      p.initial_points = o.count("initial_points", 0);
      p.batch_size = o.count("batch_size", p.batch_size);
      p.candidates = o.count("candidates", 0);
      p.max_training_points = o.count("max_training_points", p.max_training_points);
      p.hyperparameter_steps = o.count("hyperparameter_steps", p.hyperparameter_steps);
   }, RUN(turbo));
   OPTIMIZER(dycors, [](auto& p, options& o) {
      p.max_function_calls = o.count("max_function_calls", 0);
      p.initial_points = o.count("initial_points", 0);
      p.candidates = o.count("candidates", 0);
   }, RUN(dycors));

   entry("nsga2", +[](const val& f, const val& l, const val& u) { return nsga2_js::minimize(f, l, u, val::undefined()); }, &nsga2_js::minimize);
   class_<nsga2_js::optimizer>("nsga2_optimizer")
      .constructor(+[](const val& l, const val& u) { return nsga2_js::optimizer(nsga2_js::make(l, u, val::undefined(), false)); })
      .constructor(+[](const val& l, const val& u, const val& o) { return nsga2_js::optimizer(nsga2_js::make(l, u, o, false)); })
      .function("ask", +[](nsga2_js::optimizer& self) { return rows(self.o.ask()); })
      .function("tell", +[](nsga2_js::optimizer& self, const val& values) { self.o.tell(nsga2_js::matrix(values, self.n, 0)); })
      .function("front", +[](const nsga2_js::optimizer& self) { return nsga2_js::front(self.o.front()); })
      .function("population", +[](const nsga2_js::optimizer& self) { return rows(self.o.population()); })
      .function("objectives", +[](const nsga2_js::optimizer& self) { return rows(self.o.objectives()); })
      .function("ranks", +[](const nsga2_js::optimizer& self) {
         return val(typed_memory_view(self.o.ranks().size(), self.o.ranks().data())).call<val>("slice"); });
   for (const char* name : {"new nsga2_optimizer", "nsga2_optimizer.ask", "nsga2_optimizer.tell"})
      throws(name);

   entry("qlognehvi", +[](const val& f, const val& l, const val& u) { return qlognehvi_js::minimize(f, l, u, val::undefined()); }, &qlognehvi_js::minimize);
   class_<qlognehvi_js::optimizer>("qlognehvi_optimizer")
      .constructor(+[](const val& l, const val& u) { return qlognehvi_js::optimizer(qlognehvi_js::make(l, u, val::undefined(), false)); })
      .constructor(+[](const val& l, const val& u, const val& o) { return qlognehvi_js::optimizer(qlognehvi_js::make(l, u, o, false)); })
      .function("ask", +[](qlognehvi_js::optimizer& self) { return self.ask(0); })
      .function("ask", +[](qlognehvi_js::optimizer& self, double n) { return self.ask(count(n, "The number of points", 1)); })
      .function("tell", &qlognehvi_js::optimizer::tell)
      .function("front", +[](const qlognehvi_js::optimizer& self) { return nsga2_js::front(self.o.front()); })
      .function("points", +[](const qlognehvi_js::optimizer& self) { return rows(self.o.points()); })
      .function("objectives", +[](const qlognehvi_js::optimizer& self) { return rows(self.o.objectives()); })
      .function("reference_point", +[](const qlognehvi_js::optimizer& self) { return to_js(self.o.reference_point()); });
   for (const char* name : {"new qlognehvi_optimizer", "qlognehvi_optimizer.ask", "qlognehvi_optimizer.tell"})
      throws(name);
   entry("hypervolume", &qlognehvi_js::hypervolume);

   using gradient::optimizer;
   class_<optimizer>("gradient_optimizer")
      .function("step", +[](optimizer& self) { self.step(); })
      .function("steps", &optimizer::steps)
      .function("arguments", +[](const optimizer& self) { return to_js(self.arguments()); })
      .function("gradients", +[](const optimizer& self) { return to_js(self.gradients()); })
      .function("objective_value", +[](const optimizer& self) { return self.objective_value(); })
      .function("radius", +[](const optimizer& self) { return self.radius(); })
      .function("function_calls", +[](const optimizer& self) { return double(*self.calls); })
      .function("hessian_vector_products", +[](const optimizer& self) { return double(*self.products); });
   for (const char* name : {"gradient_optimizer.step", "gradient_optimizer.steps"})
      throws(name);
   gradient::bind<1>("gradient_descent", [](auto f, auto& x, options& o, auto&) {
      return bo::make_gradient_descent(std::move(f), x, o.real("learning_rate", 0.01)); });
   gradient::bind<1>("nesterov", [](auto f, auto& x, options& o, auto&) {
      const double lr = o.real("learning_rate", 0.01);
      return bo::make_nag(std::move(f), x, lr, o.real("momentum", 0.95)); });
   gradient::bind<1>("lbfgs", [](auto f, auto& x, options& o, auto&) { return bo::make_lbfgs(std::move(f), x, o.count("history", 10)); });
   gradient::bind<2>("trust_region_newton", [](auto f, auto& x, options& o, auto& products) {
      return bo::make_trust_region_newton(std::move(f), x, o.real("initial_radius", 1), bo::tape_initializer_rvar<double, 2>{},
         bo::reverse_mode_function_eval_policy<double, 2>{}, bo::reverse_mode_gradient_evaluation_policy<double, 2>{},
         gradient::counted_products{{}, products}); });
}
