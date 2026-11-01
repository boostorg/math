// Copyright 2026 Nick Thompson
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// The expansions of a float in tools/: continued fractions, Luroth and Engel expansions, with int64 terms.
// repr and str are what operator<< prints.

#include <boost/math/tools/centered_continued_fraction.hpp>
#include <boost/math/tools/engel_expansion.hpp>
#include <boost/math/tools/luroth_expansion.hpp>
#include <boost/math/tools/simple_continued_fraction.hpp>

#include <nanobind/stl/vector.h>

#include <sstream>
#include <string>

#include "bindings.hpp"
#include "macros.hpp"

namespace {

namespace nb = nanobind;

namespace bt = boost::math::tools;

// operator<< takes some of these by non-const reference.
template <class E>
std::string to_string(E& e) {
   std::ostringstream os;
   os << e;
   return os.str();
}

template <class E>
nb::class_<E> bind_expansion(nb::module_& m, char const* name, char const* doc) {
   nb::class_<E> c(m, name, doc);
   c.def(nb::init<double>(), nb::arg("x"), nb::sig("def __init__(self, x, /) -> None"))
      .def("__repr__", &to_string<E>)
      .def("__str__", &to_string<E>);
   return c;
}

}  // namespace

void bmpy::bind_expansions(nb::module_& m) {
   using scf = bt::simple_continued_fraction<double>;
   using ccf = bt::centered_continued_fraction<double>;
   using luroth = bt::luroth_expansion<double>;
   using engel = bt::engel_expansion<double>;
   bind_expansion<scf>(m, "simple_continued_fraction",
      "The simple continued fraction of x, a0 + 1/(a1 + 1/(a2 + ...))." BMPY_URL("internals/simple_continued_fraction"))
      .def("partial_denominators", [](scf const& s) { return s.partial_denominators(); }, "[a0, a1, ...], a list.")
      .def("khinchin_geometric_mean", &scf::khinchin_geometric_mean,
           "(a1 a2 ... an)^(1/n).")
      .def("khinchin_harmonic_mean", &scf::khinchin_harmonic_mean, "n/(1/a1 + ... + 1/an).");
   bind_expansion<ccf>(m, "centered_continued_fraction",
      "The centered continued fraction of x, with nearest-integer a_i." BMPY_URL("internals/centered_continued_fraction"))
      .def("partial_denominators", [](ccf const& s) { return s.partial_denominators(); }, "[a0, a1, ...], a list.")
      .def("khinchin_geometric_mean", &ccf::khinchin_geometric_mean, "(|a1 a2 ... an|)^(1/n).");
   bind_expansion<luroth>(m, "luroth_expansion",
      "The Luroth expansion of x." BMPY_URL("internals/luroth_expansion"))
      .def("digits", [](luroth const& s) { return s.digits(); }, "[d0, d1, ...], a list.")
      .def("digit_geometric_mean", &luroth::digit_geometric_mean, "(d1 d2 ... dn)^(1/n).");
   bind_expansion<engel>(m, "engel_expansion",
      "The Engel expansion of x != 0, 1/a1 + 1/(a1 a2) + ...." BMPY_URL("internals/engel_expansion"))
      .def("digits", [](engel const& s) { return s.digits(); }, "[a1, a2, ...], a list.");
}
