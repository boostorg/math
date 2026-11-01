// Copyright 2026 Nick Thompson
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// Discrete distributions use the default discrete_quantile policy.

#include <boost/math/distributions/bernoulli.hpp>
#include <boost/math/distributions/binomial.hpp>
#include <boost/math/distributions/geometric.hpp>
#include <boost/math/distributions/hypergeometric.hpp>
#include <boost/math/distributions/negative_binomial.hpp>
#include <boost/math/distributions/poisson.hpp>

#include <cstdint>
#include <string>

#include "bindings.hpp"
#include "distributions.hpp"

namespace bm = boost::math;
namespace nb = nanobind;
using nb::arg;
using bmpy::add_dist_methods;
using pol = bmpy::dist_policy;

void bmpy::bind_distributions_discrete(nb::module_& m) {
   {
      using D = bm::bernoulli_distribution<double, pol>;
      nb::class_<D> c(m, "bernoulli_distribution",
                     "Bernoulli distribution: 1 with probability p, else 0." BMPY_DIST_URL("bernoulli_dist"));
      c.def(bmpy::init<double>(), arg("p") = 0.5);
      add_dist_methods<D, true>(c);
      bmpy::add_dist_params(c, {BMPY_PARAM("p", success_fraction)});
   }
   {
      using D = bm::binomial_distribution<double, pol>;
      nb::class_<D> c(m, "binomial_distribution",
                     "Binomial distribution of the successes in n trials of success fraction p: pdf C(n, k) p^k (1 - p)^(n - k)." BMPY_DIST_URL("binomial_dist"));
      c.def(bmpy::init<double, double>(), arg("n") = 1.0, arg("p") = 0.5);
      add_dist_methods<D, true>(c);
      bmpy::add_dist_params(c, {BMPY_PARAM("n", trials), BMPY_PARAM("p", success_fraction)});
      // The bounds on p, by the Clopper-Pearson or the Jeffreys prior interval.
      auto bound = [](bool upper) {
         return [upper](double trials, double successes, double probability, std::string const& interval) {
            if (interval != "clopper_pearson" && interval != "jeffreys_prior")
               throw nb::value_error("interval must be 'clopper_pearson' or 'jeffreys_prior'");
            auto t = interval == "clopper_pearson" ? D::clopper_pearson_exact_interval : D::jeffreys_prior_interval;
            return bmpy::call_scalar("binomial_distribution.find_bound_on_p", [&] {
               return upper ? D::find_upper_bound_on_p(trials, successes, probability, t) : D::find_lower_bound_on_p(trials, successes, probability, t);
            });
         };
      };
      c.def_static("find_lower_bound_on_p", bound(false), arg(), arg(), arg(), nb::kw_only(), arg("interval") = "clopper_pearson",
         nb::sig("@staticmethod\ndef find_lower_bound_on_p(trials, successes, probability, /, *, interval='clopper_pearson')"),
         "The lower bound on the success fraction, given successes in trials, that the true p exceeds with confidence 1 - probability; "
         "interval is 'clopper_pearson' (exact) or 'jeffreys_prior'." BMPY_DIST_URL("binomial_dist"));
      c.def_static("find_upper_bound_on_p", bound(true), arg(), arg(), arg(), nb::kw_only(), arg("interval") = "clopper_pearson",
         nb::sig("@staticmethod\ndef find_upper_bound_on_p(trials, successes, probability, /, *, interval='clopper_pearson')"),
         "The upper bound on the success fraction, given successes in trials, that the true p is below with confidence 1 - probability; "
         "interval is 'clopper_pearson' (exact) or 'jeffreys_prior'." BMPY_DIST_URL("binomial_dist"));
      bmpy::def_static_scalar(c, "find_minimum_number_of_trials", "find_minimum_number_of_trials(k, p, alpha, /)",
         +[](double k, double p, double alpha) { return D::find_minimum_number_of_trials(k, p, alpha); },
         "The fewest trials of success fraction p for which more than k successes has probability at least 1 - alpha.");
      bmpy::def_static_scalar(c, "find_maximum_number_of_trials", "find_maximum_number_of_trials(k, p, alpha, /)",
         +[](double k, double p, double alpha) { return D::find_maximum_number_of_trials(k, p, alpha); },
         "The most trials of success fraction p for which k or fewer successes has probability at least 1 - alpha.");
   }
   {
      using D = bm::geometric_distribution<double, pol>;
      nb::class_<D> c(m, "geometric_distribution",
                     "Geometric distribution of the failures before the first success, of success fraction p: pdf p (1 - p)^k." BMPY_DIST_URL("geometric_dist"));
      c.def(bmpy::init<double>(), arg("p"));
      add_dist_methods<D, true>(c);
      bmpy::add_dist_params(c, {BMPY_PARAM("p", success_fraction)});
      bmpy::def_static_scalar(c, "find_lower_bound_on_p", "find_lower_bound_on_p(trials, alpha, /)",
         +[](double trials, double alpha) { return D::find_lower_bound_on_p(trials, alpha); },
         "The Clopper-Pearson lower bound on the success fraction, given trials, with confidence 1 - alpha.");
      bmpy::def_static_scalar(c, "find_upper_bound_on_p", "find_upper_bound_on_p(trials, alpha, /)",
         +[](double trials, double alpha) { return D::find_upper_bound_on_p(trials, alpha); },
         "The Clopper-Pearson upper bound on the success fraction, given trials, with confidence 1 - alpha.");
      bmpy::def_static_scalar(c, "find_minimum_number_of_trials", "find_minimum_number_of_trials(k, p, alpha, /)",
         +[](double k, double p, double alpha) { return D::find_minimum_number_of_trials(k, p, alpha); },
         "The fewest trials of success fraction p for which k or fewer failures has probability at most alpha.");
      bmpy::def_static_scalar(c, "find_maximum_number_of_trials", "find_maximum_number_of_trials(k, p, alpha, /)",
         +[](double k, double p, double alpha) { return D::find_maximum_number_of_trials(k, p, alpha); },
         "The most trials of success fraction p for which k or fewer failures has probability at least 1 - alpha.");
   }
   {
      using D = bm::hypergeometric_distribution<double, pol>;
      nb::class_<D> c(m, "hypergeometric_distribution",
                     "Hypergeometric distribution of the defective items in n drawn without replacement from N of which r are defective." BMPY_DIST_URL("hypergeometric_dist"));
      c.def(bmpy::init<std::uint64_t, std::uint64_t, std::uint64_t>(), arg("r"), arg("n"), arg("N"));
      add_dist_methods<D, true>(c);
      bmpy::add_dist_params(c, {BMPY_PARAM("r", defective), BMPY_PARAM("n", sample_count), BMPY_PARAM("N", total)});
   }
   {
      using D = bm::negative_binomial_distribution<double, pol>;
      nb::class_<D> c(m, "negative_binomial_distribution",
                     "Negative binomial distribution of the failures before the r-th success: pdf C(k + r - 1, k) p^r (1 - p)^k." BMPY_DIST_URL("negative_binomial_dist"));
      c.def(bmpy::init<double, double>(), arg("r"), arg("p"));
      add_dist_methods<D, true>(c);
      bmpy::add_dist_params(c, {BMPY_PARAM("r", successes), BMPY_PARAM("p", success_fraction)});
      bmpy::def_static_scalar(c, "find_lower_bound_on_p", "find_lower_bound_on_p(trials, successes, alpha, /)",
         +[](double trials, double successes, double alpha) { return D::find_lower_bound_on_p(trials, successes, alpha); },
         "The Clopper-Pearson lower bound on the success fraction, given successes in trials, with confidence 1 - alpha.");
      bmpy::def_static_scalar(c, "find_upper_bound_on_p", "find_upper_bound_on_p(trials, successes, alpha, /)",
         +[](double trials, double successes, double alpha) { return D::find_upper_bound_on_p(trials, successes, alpha); },
         "The Clopper-Pearson upper bound on the success fraction, given successes in trials, with confidence 1 - alpha.");
      bmpy::def_static_scalar(c, "find_minimum_number_of_trials", "find_minimum_number_of_trials(k, p, alpha, /)",
         +[](double k, double p, double alpha) { return D::find_minimum_number_of_trials(k, p, alpha); },
         "The fewest trials of success fraction p for which k or fewer failures has probability at most alpha.");
      bmpy::def_static_scalar(c, "find_maximum_number_of_trials", "find_maximum_number_of_trials(k, p, alpha, /)",
         +[](double k, double p, double alpha) { return D::find_maximum_number_of_trials(k, p, alpha); },
         "The most trials of success fraction p for which k or fewer failures has probability at least 1 - alpha.");
   }
   {
      using D = bm::poisson_distribution<double, pol>;
      nb::class_<D> c(m, "poisson_distribution",
                     "Poisson distribution of the given mean: pdf mean^k exp(-mean) / k!." BMPY_DIST_URL("poisson_dist"));
      c.def(bmpy::init<double>(), arg("mean") = 1.0);
      add_dist_methods<D, true>(c);
      bmpy::add_dist_params(c, {BMPY_PARAM("mean", mean)});
   }
}
