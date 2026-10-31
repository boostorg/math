// Copyright 2026 Nick Thompson
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include "bindings.hpp"

#include <string>

namespace {

namespace nb = nanobind;

// boost.math.<name>, importable as such: it goes in sys.modules beside this extension's package.
nb::module_ submodule(nb::module_& m, char const* name, char const* doc) {
   std::string package = nb::cast<std::string>(m.attr("__name__"));
   std::string full = package.substr(0, package.rfind('.')) + "." + name;
   nb::module_ sub = nb::steal<nb::module_>(PyModule_New(full.c_str()));
   sub.doc() = doc;
   m.attr(name) = sub;
   nb::module_::import_("sys").attr("modules")[full.c_str()] = sub;
   return sub;
}

}  // namespace

NB_MODULE(_ext, m) {
   m.doc() = "Development-only nanobind bindings for Boost.Math.";
   // A callback that returns the wrong type: TypeError, not RuntimeError("std::bad_cast").
   nb::register_exception_translator([](std::exception_ptr const& p, void*) {
      try {
         std::rethrow_exception(p);
      } catch (nb::cast_error const&) {
         PyErr_SetString(PyExc_TypeError, "a callback returned a value of the wrong type");
      }
   });
   bmpy::bind_errors(m);
   bmpy::bind_gamma(m);
   bmpy::bind_beta(m);
   bmpy::bind_erf_expint_zeta(m);
   bmpy::bind_elementary(m);
   bmpy::bind_bessel(m);
   bmpy::bind_elliptic(m);
   bmpy::bind_jacobi_theta(m);
   bmpy::bind_polynomials(m);
   bmpy::bind_hypergeometric(m);
   nb::module_ distributions = submodule(m, "distributions", "Probability distributions.");
   bmpy::bind_distributions_continuous(distributions);
   bmpy::bind_distributions_heavy(distributions);
   bmpy::bind_distributions_discrete(distributions);
   bmpy::bind_distributions_noncentral(distributions);
   nb::module_ fft = submodule(m, "fft", "Fourier and cosine transforms.");
   bmpy::bind_fft(fft);
   nb::module_ quadrature = submodule(m, "quadrature", "Numerical integration.");
   bmpy::bind_quadrature(quadrature);
   bmpy::bind_naive_monte_carlo(quadrature);
   nb::module_ optimize = submodule(m, "optimize", "Optimisation: derivative-free global, gradient-based, and one-dimensional (Brent).");
   bmpy::bind_optimization(optimize);
   nb::module_ roots = submodule(m, "roots",
      "Root finding: bracketing (bisect, toms748, bracket_and_solve_root) and Newton-type (newton, halley, schroder,\n"
      "complex_newton) solvers on a Python callable, and the roots of polynomials.\n\n"
      "https://www.boost.org/doc/libs/release/libs/math/doc/html/root_finding.html");
   bmpy::bind_polynomial_roots(roots);
   bmpy::bind_roots(roots);
   nb::module_ statistics = submodule(m, "statistics", "Descriptive statistics and hypothesis tests.");
   bmpy::bind_statistics(statistics);
   nb::module_ interpolators = submodule(m, "interpolators", "Interpolation.");
   bmpy::bind_interpolators(interpolators);
   bmpy::bind_special_extras(m);
   bmpy::bind_daubechies(m);
   nb::module_ constants = submodule(m, "constants", "Mathematical constants, as floats.\n\n"
      "https://www.boost.org/doc/libs/release/libs/math/doc/html/math_toolkit/constants.html");
   bmpy::bind_constants(constants);
   nb::module_ differentiation = submodule(m, "differentiation", "Numerical differentiation.");
   bmpy::bind_differentiation(differentiation);
   nb::module_ expansions = submodule(m, "expansions", "Continued fractions and other expansions of a float.");
   bmpy::bind_expansions(expansions);
#ifdef BMPY_HAVE_ULPS
   bmpy::bind_ulps(m);
#endif
}
