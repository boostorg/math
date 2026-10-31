// Copyright 2026 Nick Thompson
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// The special functions of other shapes: of an array (logsumexp), with array results (legendre_p_zeros), with a
// compile-time order chosen at run time (cardinal B-splines, Fourier transforms of Daubechies functions), and the
// classes legendre_stieltjes and chebyshev_transform.

#include <boost/math/special_functions/cardinal_b_spline.hpp>
#include <boost/math/special_functions/chebyshev.hpp>
#include <boost/math/special_functions/chebyshev_transform.hpp>
#include <boost/math/special_functions/fourier_transform_daubechies.hpp>
#include <boost/math/special_functions/jacobi.hpp>
#include <boost/math/special_functions/legendre.hpp>
#include <boost/math/special_functions/legendre_stieltjes.hpp>
#include <boost/math/special_functions/logsumexp.hpp>

#include <nanobind/stl/pair.h>

#include <cmath>
#include <limits>
#include <string>
#include <vector>

#include "bindings.hpp"
#include "macros.hpp"

namespace {

namespace nb = nanobind;

struct py_legendre_stieltjes : bm::legendre_stieltjes<double> {
   std::size_t m;
   explicit py_legendre_stieltjes(std::size_t m_) : bm::legendre_stieltjes<double>(m_), m(m_) {}
};

struct py_chebyshev_transform : bm::chebyshev_transform<double> {
   double a, b;
   py_chebyshev_transform(bmpy::py_callable<double, double> const& f, double a_, double b_, double tol, std::size_t max_refinements)
       : bm::chebyshev_transform<double>(f, a_, b_, tol, max_refinements), a(a_), b(b_) {}
};

// B_n and its relatives for n in [Lo, 20].
template <int Lo, class Fn>
void def_b_spline(nb::module_& m, char const* name, Fn fn, char const* doc) {
   bmpy::def_vec<bmpy::Int<int>, bmpy::Real>(
       m, name, [fn]<class T>(int n, T x) { return bmpy::with_int<Lo, 20>("n", n, [&](auto N) { return fn(N, x); }); },
       {"n", "x"}, doc);
}

}  // namespace

void bmpy::bind_special_extras(nb::module_& m) {
   using bmpy::Int;
   using bmpy::Real;

   m.def("logsumexp", [](nb::handle x) {
      std::vector<double> v = array_to_vector(x);
      return v.empty() ? -std::numeric_limits<double>::infinity() : bm::logsumexp(v);
   }, nb::sig("def logsumexp(x, /) -> float"),
      "log(sum(exp(x))) over all the elements of x, without overflow; -inf for an empty array." BMPY_URL("powers/logaddexp"));

   m.def("legendre_p_zeros", [](int l) {
      std::vector<double> z;
      {
         nb::gil_scoped_release release;
         z = bm::legendre_p_zeros<double>(l);
      }
      return vector_to_array(std::move(z));
   }, nb::sig("def legendre_p_zeros(l, /)"),
      "The nonnegative zeros of the Legendre polynomial P_l, increasing; the others are their negatives." BMPY_URL("sf_poly/legendre"));

   bmpy::def_vec<Int<unsigned>, Real, Real, Real>(
       m, "jacobi_double_prime", []<class T>(unsigned n, T a, T b, T x) { return bm::jacobi_double_prime(n, a, b, x); },
       {"n", "alpha", "beta", "x"},
       "The second derivative in x of the Jacobi polynomial P_n^(alpha, beta)(x)." BMPY_URL("sf_poly/jacobi"));

   def_b_spline<0>(m, "cardinal_b_spline", []<class T>(auto N, T x) { return bm::cardinal_b_spline<N()>(x); },
       "The cardinal B-spline B_n(x), n in [0, 20], with support [-(n+1)/2, (n+1)/2]." BMPY_URL("sf_poly/cardinal_b_splines"));
   def_b_spline<0>(m, "cardinal_b_spline_prime", []<class T>(auto N, T x) { return bm::cardinal_b_spline_prime<N()>(x); },
       "The derivative B_n'(x) of the cardinal B-spline, n in [0, 20]." BMPY_URL("sf_poly/cardinal_b_splines"));
   def_b_spline<3>(m, "cardinal_b_spline_double_prime",
       []<class T>(auto N, T x) { return bm::cardinal_b_spline_double_prime<N()>(x); },
       "The second derivative B_n''(x) of the cardinal B-spline, n in [3, 20]." BMPY_URL("sf_poly/cardinal_b_splines"));
   def_b_spline<0>(m, "forward_cardinal_b_spline", []<class T>(auto N, T x) { return bm::forward_cardinal_b_spline<N()>(x); },
       "The cardinal B-spline shifted to support [0, n+1]: B_n(x - (n+1)/2), n in [0, 20]." BMPY_URL("sf_poly/cardinal_b_splines"));

   bmpy::def_vec<Int<int>, Real>(
       m, "fourier_transform_daubechies_scaling",
       []<class T>(int p, T omega) {
          return with_int<1, 10>("p", p, [&](auto P) { return bm::fourier_transform_daubechies_scaling<T, P()>(omega); });
       }, {"p", "omega"},
       "The Fourier transform (complex) of the Daubechies scaling function with p vanishing moments, p in [1, 10]."
       BMPY_URL("daubechies"));
   bmpy::def_vec<Int<int>, Real>(
       m, "fourier_transform_daubechies_wavelet",
       []<class T>(int p, T omega) {
          return with_int<1, 10>("p", p, [&](auto P) { return bm::fourier_transform_daubechies_wavelet<T, P()>(omega); });
       }, {"p", "omega"},
       "The Fourier transform (complex) of the Daubechies wavelet with p vanishing moments, p in [1, 10]." BMPY_URL("daubechies"));

   // The coefficients are an argument, so these are not def_vec functions of fixed arity.
   m.def("chebyshev_clenshaw_recurrence", [](nb::handle c, nb::handle x, int threads) {
      std::vector<double> cs = array_to_vector(c);
      auto f = [&cs]<class T>(double t) { return bm::chebyshev_clenshaw_recurrence(cs.data(), cs.size(), t); };
      return call_array<false, Real>("chebyshev_clenshaw_recurrence", f, {x}, {"x"}, threads);
   }, nb::arg(), nb::arg(), nb::kw_only(), nb::arg("threads") = 1,
      nb::sig("def chebyshev_clenshaw_recurrence(c, x, /, *, threads=1)"),
      "The Chebyshev series c[0]/2 + sum over k >= 1 of c[k] T_k(x), by Clenshaw's recurrence." BMPY_URL("sf_poly/chebyshev"));
   m.def("chebyshev_clenshaw_recurrence", [](nb::handle c, double a, double b, nb::handle x, int threads) {
      std::vector<double> cs = array_to_vector(c);
      auto f = [&]<class T>(double t) { return bm::chebyshev_clenshaw_recurrence(cs.data(), cs.size(), a, b, t); };
      return call_array<false, Real>("chebyshev_clenshaw_recurrence", f, {x}, {"x"}, threads);
   }, nb::arg(), nb::arg(), nb::arg(), nb::arg(), nb::kw_only(), nb::arg("threads") = 1,
      nb::sig("def chebyshev_clenshaw_recurrence(c, a, b, x, /, *, threads=1)"),
      "The same series on [a, b]: the T_k at (2x - a - b)/(b - a); x outside [a, b] is a domain error.");

   nb::class_<py_legendre_stieltjes> ls(m, "legendre_stieltjes",
      "The Legendre-Stieltjes polynomial E_m, m >= 1, whose zeros are the Kronrod nodes." BMPY_URL("sf_poly/legendre_stieltjes"));
   ls.def(nb::init<std::size_t>(), nb::arg("m"), nb::sig("def __init__(self, m, /) -> None"))
      .def_ro("m", &py_legendre_stieltjes::m)
      .def("norm_sq", [](py_legendre_stieltjes const& E) { return E.norm_sq(); }, "The integral of E_m(x)^2 over [-1, 1].")
      .def("zeros", [](py_legendre_stieltjes const& E) { return vector_to_array(E.zeros()); },
           "The nonnegative zeros of E_m, in increasing order (E_m is even or odd).")
      .def("__repr__", [](py_legendre_stieltjes const& E) { return "legendre_stieltjes(" + std::to_string(E.m) + ")"; });
   def_method_vec(ls, "__call__", "x", [](py_legendre_stieltjes const& E, double x) { return E(x); }, "E_m(x).");
   def_method_vec(ls, "prime", "x", [](py_legendre_stieltjes const& E, double x) { return E.prime(x); }, "E_m'(x).");

   nb::class_<py_chebyshev_transform> ct(m, "chebyshev_transform",
      "The Chebyshev series c[0]/2 + sum over k >= 1 of c[k] T_k(t), t = (2x - a - b)/(b - a), of f on [a, b]."
      BMPY_URL("sf_poly/chebyshev"));
   ct.def("__init__", [](py_chebyshev_transform* t, nb::callable f, double a, double b, double tol, std::size_t max_refinements) {
         // With none, the C++ leaves no coefficients, which integrate() then reads.
         if (max_refinements == 0)
            throw nb::value_error("chebyshev_transform: max_refinements must be positive");
         new (t) py_chebyshev_transform(py_callable<double, double>{f}, a, b, tol, max_refinements);
      }, nb::arg("f"), nb::arg("a"), nb::arg("b"), nb::kw_only(), nb::arg("tol") = 500 * std::numeric_limits<double>::epsilon(),
      nb::arg("max_refinements") = 16, nb::sig("def __init__(self, f, a, b, /, *, tol=1.1102230246251565e-13, max_refinements=16) -> None"))
      .def_ro("a", &py_chebyshev_transform::a)
      .def_ro("b", &py_chebyshev_transform::b)
      .def_prop_ro("coefficients", [](py_chebyshev_transform const& t) { return vector_to_array(t.coefficients()); },
           nb::rv_policy::move, "The c[k], a copy.")
      .def("integrate", [](py_chebyshev_transform const& t) { return t.integrate(); }, "The integral of the series over [a, b].")
      .def("__repr__", [](py_chebyshev_transform const& t) {
         return "<chebyshev_transform on [" + py_float_repr(t.a) + ", " + py_float_repr(t.b) + "], "
                + std::to_string(t.coefficients().size()) + " coefficients>";
      });
   def_method_vec(ct, "__call__", "x", [](py_chebyshev_transform const& t, double x) { return t(x); },
       "The series at x in [a, b], elementwise.");
   def_method_vec(ct, "prime", "x", [](py_chebyshev_transform const& t, double x) { return t.prime(x); },
       "The derivative of the series at x in [a, b], elementwise.");
}
