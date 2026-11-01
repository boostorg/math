// Copyright 2026 Nick Thompson
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// The interpolators, in double, each a class that owns a copy of the data it was built from. The scalar-valued
// ones evaluate elementwise (threads= as for the special functions); the curves (catmull_rom, bezier_polynomial,
// vector_barycentric_rational) take points with 1 to 4 coordinates, and map t of any shape to t.shape + (dim,).

#include "bindings.hpp"

#include <boost/math/interpolators/barycentric_rational.hpp>
#include <boost/math/interpolators/bezier_polynomial.hpp>
#include <boost/math/interpolators/bilinear_uniform.hpp>
#include <boost/math/interpolators/cardinal_cubic_b_spline.hpp>
#include <boost/math/interpolators/cardinal_quadratic_b_spline.hpp>
#include <boost/math/interpolators/cardinal_quintic_b_spline.hpp>
#include <boost/math/interpolators/cardinal_trigonometric.hpp>
#include <boost/math/interpolators/catmull_rom.hpp>
#include <boost/math/interpolators/cubic_hermite.hpp>
#include <boost/math/interpolators/makima.hpp>
#include <boost/math/interpolators/pchip.hpp>
#include <boost/math/interpolators/quintic_hermite.hpp>
#include <boost/math/interpolators/septic_hermite.hpp>
#include <boost/math/interpolators/vector_barycentric_rational.hpp>
#include <boost/math/interpolators/whittaker_shannon.hpp>

#include <nanobind/stl/optional.h>
#include <nanobind/stl/pair.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <variant>
#include <vector>

namespace bmpy {

namespace {

namespace bi = boost::math::interpolators;
using vec = std::vector<double>;
using opt = std::optional<double>;
constexpr double nan = std::numeric_limits<double>::quiet_NaN();

// An interpolator, with the number of points and the rest of its repr.
template <class I>
struct held {
   I f;
   std::size_t n;
   std::string params;
};

// A C-ordered copy of h as float64, which must have `ndim` dimensions; its shape goes to `shape`.
vec copy_array(nb::handle h, std::size_t ndim, char const* what, std::size_t* shape = nullptr) {
   nb::object a = numpy().attr("ascontiguousarray")(h, nb::arg("dtype") = numpy().attr("float64"));
   auto v = nb::cast<view<double>>(a);
   if (v.ndim() != ndim)
      throw nb::value_error((std::string(what) + " must be a " + std::to_string(ndim) + "-D array").c_str());
   for (std::size_t d = 0; shape && d < ndim; ++d) shape[d] = v.shape(d);
   return vec(v.data(), v.data() + v.size());
}
vec copy_vec(nb::handle h, char const* what) { return copy_array(h, 1, what); }

// Boost's constructors report bad data with std::logic_error and its subclasses: ValueError.
template <class F>
auto checked(F const& make) -> decltype(make()) {
   try {
      return make();
   } catch (std::logic_error const& e) {
      throw nb::value_error(e.what());
   }
}

template <class H>
void add_repr(nb::class_<H>& c) {
   std::string name = nb::cast<std::string>(c.attr("__name__"));
   c.def("__repr__", [name](H const& h) { return name + "(n=" + std::to_string(h.n) + h.params + ")"; });
}

// __call__, and prime, double_prime and domain where the C++ has them. The classes with push_back hold the GIL
// while they evaluate, so that push_back (which needs it) cannot change them meanwhile.
template <class H>
void add_calls(nb::class_<H>& c, bool release_gil = true) {
   // NaN passes the C++ range checks of some (the Hermite interpolators then read out of bounds): NaN out.
   def_method_vec(c, "__call__", "x", [](H const& h, double x) { return x == x ? h.f(x) : x; }, "The interpolant at x.",
                  release_gil);
   if constexpr (requires(H const& h) { h.f.prime(0.0); })
      def_method_vec(c, "prime", "x", [](H const& h, double x) { return x == x ? h.f.prime(x) : x; }, "Its derivative at x.",
                     release_gil);
   // septic_hermite's double_prime is a stub that returns NaN.
   if constexpr (requires(H const& h) { h.f.double_prime(0.0); } && !std::is_same_v<H, held<bi::septic_hermite<vec>>>)
      def_method_vec(c, "double_prime", "x", [](H const& h, double x) { return x == x ? h.f.double_prime(x) : x; },
                     "Its second derivative at x.", release_gil);
   if constexpr (requires(H const& h) { h.f.domain(); })
      c.def("domain", [](H const& h) { return h.f.domain(); }, "The interval (a, b) on which the interpolant is defined.");
   add_repr(c);
}

std::string cardinal(char const* t0_name, double t0, char const* h_name, double h) {
   return std::string(", ") + t0_name + "=" + py_float_repr(t0) + ", " + h_name + "=" + py_float_repr(h);
}

// The curves, with points in R^D.
template <std::size_t D>
using point = std::array<double, D>;
template <std::size_t D>
using catmull_rom_d = boost::math::catmull_rom<point<D>, std::vector<point<D>>>;
template <std::size_t D>
using bezier_d = bi::bezier_polynomial<std::vector<point<D>>>;
template <std::size_t D>
using vector_barycentric_d = bi::vector_barycentric_rational<vec, std::vector<point<D>>>;
template <template <std::size_t> class C>
using any_dim = std::variant<C<1>, C<2>, C<3>, C<4>>;

// make.template operator()<D>(points) for the points (rows of an (n, D) array) as std::array<double, D>.
template <template <std::size_t> class C, class Make>
held<any_dim<C>> make_curve(nb::handle points, char const* what, std::string const& params, Make const& make) {
   std::size_t shape[2];
   vec v = copy_array(points, 2, what, shape);
   auto build = [&]<std::size_t D>() -> any_dim<C> {
      std::vector<point<D>> p(shape[0]);
      for (std::size_t i = 0; i < p.size(); ++i) std::copy_n(v.data() + i * D, D, p[i].begin());
      return checked([&] { return make.template operator()<D>(std::move(p)); });
   };
   auto f = [&]() -> any_dim<C> {
      switch (shape[1]) {
      case 1: return build.template operator()<1>();
      case 2: return build.template operator()<2>();
      case 3: return build.template operator()<3>();
      case 4: return build.template operator()<4>();
      }
      throw nb::value_error((std::string(what) + " must have 1 to 4 coordinates").c_str());
   }();
   return {std::move(f), shape[0], ", dim=" + std::to_string(shape[1]) + params};
}

// g(curve, t) at each t, as an array of shape t.shape + (dim,); errors as for the elementwise functions.
template <class H, class G>
nb::object call_points(char const* name, H const& h, nb::handle t, G const& g) {
   nb::object a = numpy().attr("asarray")(t, nb::arg("dtype") = numpy().attr("float64"));
   nb::object flat = a.attr("ravel")();
   auto ts = nb::cast<view<double>>(flat);
   std::size_t n = ts.size(), dim = h.f.index() + 1;
   auto* out = new vec(n * dim);
   nb::capsule owner(out, [](void* p) noexcept { delete static_cast<vec*>(p); });
   error_actions actions = current_error_actions;
   error_record errors;
   {
      recording rec(errors);
      std::visit([&](auto const& f) {
         for (std::size_t i = 0; i < n; ++i) {
            try {
               if (std::isnan(ts.data()[i]))
                  throw std::domain_error(std::string(name) + ": t is NaN");
               auto p = g(f, ts.data()[i]);
               std::copy(p.begin(), p.end(), out->data() + i * dim);
            } catch (std::exception const& e) {
               record_exception(e);
               std::fill_n(out->data() + i * dim, dim, nan);
            }
            errors.next();
         }
      }, h.f);
   }
   if (errors.any())
      report(name, errors, n, nb::cast<std::size_t>(a.attr("ndim")) == 0, false, actions);
   nb::object shape = nb::steal(PySequence_Concat(a.attr("shape").ptr(), nb::make_tuple(dim).ptr()));
   return nb::cast(nb::ndarray<nb::numpy, double>(out->data(), {n * dim}, owner)).attr("reshape")(shape);
}

template <class H, class G>
void def_points(nb::class_<H>& c, char const* name, G g, char const* doc) {
   c.def(name, [name, g](H const& h, nb::handle t) { return call_points(name, h, t, g); }, nb::arg(),
         nb::sig(("def " + std::string(name) + "(self, t, /)").c_str()), doc);
}

auto curve_value = [](auto const& f, double t) { return f(t); };
auto curve_prime = [](auto const& f, double t) { return f.prime(t); };

// The C++ prints to std::cerr and returns NaN outside [0, 1].
template <class G>
auto bezier_checked(G g) {
   return [g](auto const& f, double t) {
      if (!(t >= 0 && t <= 1))
         throw std::domain_error("bezier_polynomial: t must be in [0, 1]");
      return g(f, t);
   };
}

using bary = held<bi::barycentric_rational<double>>;
using sinc = held<bi::whittaker_shannon<vec>>;
using cardinal_cubic_h = held<bi::cardinal_cubic_hermite<vec>>;
using cardinal_quintic_h = held<bi::cardinal_quintic_hermite<vec>>;
using septic_h = held<bi::septic_hermite<vec>>;
using cardinal_septic_h = held<bi::cardinal_septic_hermite<vec>>;
using cubic_b = held<bi::cardinal_cubic_b_spline<double>>;
using quadratic_b = held<bi::cardinal_quadratic_b_spline<double>>;
using quintic_b = held<bi::cardinal_quintic_b_spline<double>>;
using trig = held<bi::cardinal_trigonometric<vec>>;
using makima = held<bi::makima<vec>>;
using pchip = held<bi::pchip<vec>>;
using cubic_h = held<bi::cubic_hermite<vec>>;
using quintic_h = held<bi::quintic_hermite<vec>>;
using catmull = held<any_dim<catmull_rom_d>>;
using bezier = held<any_dim<bezier_d>>;
using vbary = held<any_dim<vector_barycentric_d>>;

// bilinear_uniform, with the rectangle it is defined on: the C++ prints to std::cerr outside it.
struct bilinear : held<bi::bilinear_uniform<vec>> {
   double x0, x1, y0, y1;
};

}  // namespace

void bind_interpolators(nb::module_& m) {
   nb::class_<bary> br(m, "barycentric_rational",
       "Floater-Hormann barycentric rational interpolant of (x[i], y[i]); approximation_order < len(x)." BMPY_URL("barycentric"));
   br.def("__init__", [](bary* self, nb::handle x, nb::handle y, std::size_t order) {
          vec xs = copy_vec(x, "x"), ys = copy_vec(y, "y");
          if (xs.size() != ys.size())
             throw nb::value_error("barycentric_rational: x and y must have the same length");
          new (self) bary{checked([&] { return bi::barycentric_rational<double>(xs.data(), ys.data(), xs.size(), order); }),
                                         xs.size(), ", approximation_order=" + std::to_string(order)};
       }, nb::arg(), nb::arg(), nb::arg("approximation_order") = 3,
       nb::sig("def __init__(self, x, y, /, approximation_order=3) -> None"));
   add_calls(br);

   nb::class_<cubic_b> cb(m, "cardinal_cubic_b_spline",
       "Cubic B-spline interpolant of y[i] = f(t0 + i*h); endpoint derivatives are estimated unless given."
       BMPY_URL("cardinal_cubic_b"));
   cb.def("__init__", [](cubic_b* self, nb::handle y, double t0, double h, opt d0, opt d1) {
          vec ys = copy_vec(y, "y");
          new (self) cubic_b{checked([&] { return bi::cardinal_cubic_b_spline<double>(ys.data(), ys.size(), t0, h, d0.value_or(nan), d1.value_or(nan)); }),
                             ys.size(), cardinal("t0", t0, "h", h)};
       }, nb::arg(), nb::arg(), nb::arg(), nb::arg("left_endpoint_derivative") = nb::none(), nb::arg("right_endpoint_derivative") = nb::none(),
       nb::sig("def __init__(self, y, t0, h, /, left_endpoint_derivative=None, right_endpoint_derivative=None) -> None"));
   add_calls(cb);

   nb::class_<quadratic_b> qb(m, "cardinal_quadratic_b_spline",
       "Quadratic B-spline interpolant of y[i] = f(t0 + i*h); endpoint derivatives are estimated unless given."
       BMPY_URL("cardinal_quadratic_b"));
   qb.def("__init__", [](quadratic_b* self, nb::handle y, double t0, double h, opt d0, opt d1) {
          vec ys = copy_vec(y, "y");
          new (self) quadratic_b{checked([&] { return bi::cardinal_quadratic_b_spline<double>(ys, t0, h, d0.value_or(nan), d1.value_or(nan)); }),
                                 ys.size(), cardinal("t0", t0, "h", h)};
       }, nb::arg(), nb::arg(), nb::arg(), nb::arg("left_endpoint_derivative") = nb::none(), nb::arg("right_endpoint_derivative") = nb::none(),
       nb::sig("def __init__(self, y, t0, h, /, left_endpoint_derivative=None, right_endpoint_derivative=None) -> None"));
   add_calls(qb);

   using pair = std::optional<std::pair<double, double>>;
   nb::class_<quintic_b> qnb(m, "cardinal_quintic_b_spline",
       "Quintic B-spline interpolant of y[i] = f(t0 + i*h); endpoint (first, second) derivatives are estimated unless given."
       BMPY_URL("cardinal_quintic_b"));
   qnb.def("__init__", [](quintic_b* self, nb::handle y, double t0, double h, pair d0, pair d1) {
          vec ys = copy_vec(y, "y");
          std::pair<double, double> none{nan, nan};
          new (self) quintic_b{checked([&] { return bi::cardinal_quintic_b_spline<double>(ys, t0, h, d0.value_or(none), d1.value_or(none)); }),
                               ys.size(), cardinal("t0", t0, "h", h)};
       }, nb::arg(), nb::arg(), nb::arg(), nb::arg("left_endpoint_derivatives") = nb::none(), nb::arg("right_endpoint_derivatives") = nb::none(),
       nb::sig("def __init__(self, y, t0, h, /, left_endpoint_derivatives=None, right_endpoint_derivatives=None) -> None"));
   add_calls(qnb);

   nb::class_<trig> ct(m, "cardinal_trigonometric",
       "Trigonometric interpolant of samples y[i] = f(t0 + i*h) of a function of period n*h." BMPY_URL("cardinal_trigonometric"));
   ct.def("__init__", [](trig* self, nb::handle y, double t0, double h) {
          vec ys = copy_vec(y, "y");
          new (self) trig{checked([&] { return bi::cardinal_trigonometric<vec>(ys, t0, h); }), ys.size(), cardinal("t0", t0, "h", h)};
       }, nb::arg(), nb::arg(), nb::arg(), nb::sig("def __init__(self, y, t0, h, /) -> None"))
     .def("period", [](trig const& t) { return t.f.period(); }, "The period n*h.")
     .def("integrate", [](trig const& t) { return t.f.integrate(); }, "The integral of the interpolant over one period.")
     .def("squared_l2", [](trig const& t) { return t.f.squared_l2(); }, "The integral of its square over one period.");
   add_calls(ct);

   nb::class_<sinc> ws(m, "whittaker_shannon",
       "Whittaker-Shannon (sinc) interpolant of samples y[i] = f(t0 + i*h)." BMPY_URL("whittaker_shannon"));
   ws.def("__init__", [](sinc* self, nb::handle y, double t0, double h) {
          vec ys = copy_vec(y, "y");
          std::size_t n = ys.size();
          new (self) sinc{checked([&] { return bi::whittaker_shannon<vec>(std::move(ys), t0, h); }), n, cardinal("t0", t0, "h", h)};
       }, nb::arg(), nb::arg(), nb::arg(), nb::sig("def __init__(self, y, t0, h, /) -> None"));
   add_calls(ws);

   nb::class_<makima> mk(m, "makima",
       "Modified Akima interpolant of (x[i], y[i]), x strictly increasing, at least 4 points." BMPY_URL("makima"));
   mk.def("__init__", [](makima* self, nb::handle x, nb::handle y, opt d0, opt d1) {
          vec xs = copy_vec(x, "x"), ys = copy_vec(y, "y");
          if (xs.size() != ys.size())  // Boost reads y before checking its length.
             throw nb::value_error("x and y must have the same length");
          std::size_t n = xs.size();
          new (self) makima{checked([&] { return bi::makima<vec>(std::move(xs), std::move(ys), d0.value_or(nan), d1.value_or(nan)); }), n, ""};
       }, nb::arg(), nb::arg(), nb::arg("left_endpoint_derivative") = nb::none(), nb::arg("right_endpoint_derivative") = nb::none(),
       nb::sig("def __init__(self, x, y, /, left_endpoint_derivative=None, right_endpoint_derivative=None) -> None"))
     .def("push_back", [](makima& h, double x, double y) { h.f.push_back(x, y); ++h.n; }, nb::arg(), nb::arg(),
          nb::sig("def push_back(self, x, y, /) -> None"), "Appends the point (x, y), with x beyond the last abscissa.");
   add_calls(mk, false);

   nb::class_<pchip> pc(m, "pchip",
       "Monotone piecewise cubic Hermite interpolant (PCHIP) of (x[i], y[i]), x strictly increasing, at least 4 points."
       BMPY_URL("pchip"));
   pc.def("__init__", [](pchip* self, nb::handle x, nb::handle y, opt d0, opt d1) {
          vec xs = copy_vec(x, "x"), ys = copy_vec(y, "y");
          if (xs.size() != ys.size())  // Boost reads y before checking its length.
             throw nb::value_error("x and y must have the same length");
          std::size_t n = xs.size();
          new (self) pchip{checked([&] { return bi::pchip<vec>(std::move(xs), std::move(ys), d0.value_or(nan), d1.value_or(nan)); }), n, ""};
       }, nb::arg(), nb::arg(), nb::arg("left_endpoint_derivative") = nb::none(), nb::arg("right_endpoint_derivative") = nb::none(),
       nb::sig("def __init__(self, x, y, /, left_endpoint_derivative=None, right_endpoint_derivative=None) -> None"))
     .def("push_back", [](pchip& h, double x, double y) { h.f.push_back(x, y); ++h.n; }, nb::arg(), nb::arg(),
          nb::sig("def push_back(self, x, y, /) -> None"), "Appends the point (x, y), with x beyond the last abscissa.");
   add_calls(pc, false);

   nb::class_<cubic_h> ch(m, "cubic_hermite",
       "Cubic Hermite interpolant of values y and derivatives dydx at strictly increasing x." BMPY_URL("cubic_hermite"));
   ch.def("__init__", [](cubic_h* self, nb::handle x, nb::handle y, nb::handle dydx) {
          vec xs = copy_vec(x, "x"), ys = copy_vec(y, "y"), ds = copy_vec(dydx, "dydx");
          std::size_t n = xs.size();
          new (self) cubic_h{checked([&] { return bi::cubic_hermite<vec>(std::move(xs), std::move(ys), std::move(ds)); }), n, ""};
       }, nb::arg(), nb::arg(), nb::arg(), nb::sig("def __init__(self, x, y, dydx, /) -> None"))
     .def("push_back", [](cubic_h& h, double x, double y, double dydx) { h.f.push_back(x, y, dydx); ++h.n; },
          nb::arg(), nb::arg(), nb::arg(), nb::sig("def push_back(self, x, y, dydx, /) -> None"),
          "Appends the point (x, y) with derivative dydx, with x beyond the last abscissa.");
   add_calls(ch, false);

   nb::class_<cardinal_cubic_h> cch(m, "cardinal_cubic_hermite",
       "Cubic Hermite interpolant of values y[i] and derivatives dydx[i] at the equispaced x0 + i*dx." BMPY_URL("cubic_hermite"));
   cch.def("__init__", [](cardinal_cubic_h* self, nb::handle y, nb::handle dydx, double x0, double dx) {
          vec ys = copy_vec(y, "y"), ds = copy_vec(dydx, "dydx");
          std::size_t n = ys.size();
          new (self) cardinal_cubic_h{checked([&] { return bi::cardinal_cubic_hermite<vec>(std::move(ys), std::move(ds), x0, dx); }),
                                          n, cardinal("x0", x0, "dx", dx)};
       }, nb::arg(), nb::arg(), nb::arg(), nb::arg(), nb::sig("def __init__(self, y, dydx, x0, dx, /) -> None"));
   add_calls(cch);

   nb::class_<quintic_h> qh(m, "quintic_hermite",
       "Quintic Hermite interpolant of y and its first two derivatives at strictly increasing x." BMPY_URL("quintic_hermite"));
   qh.def("__init__", [](quintic_h* self, nb::handle x, nb::handle y, nb::handle dydx, nb::handle d2ydx2) {
          vec xs = copy_vec(x, "x"), ys = copy_vec(y, "y"), ds = copy_vec(dydx, "dydx"), d2s = copy_vec(d2ydx2, "d2ydx2");
          std::size_t n = xs.size();
          new (self) quintic_h{checked([&] { return bi::quintic_hermite<vec>(std::move(xs), std::move(ys), std::move(ds), std::move(d2s)); }), n, ""};
       }, nb::arg(), nb::arg(), nb::arg(), nb::arg(), nb::sig("def __init__(self, x, y, dydx, d2ydx2, /) -> None"))
     .def("push_back", [](quintic_h& h, double x, double y, double dydx, double d2ydx2) { h.f.push_back(x, y, dydx, d2ydx2); ++h.n; },
          nb::arg(), nb::arg(), nb::arg(), nb::arg(), nb::sig("def push_back(self, x, y, dydx, d2ydx2, /) -> None"),
          "Appends the point (x, y) with derivatives dydx and d2ydx2, with x beyond the last abscissa.");
   add_calls(qh, false);

   nb::class_<cardinal_quintic_h> cqh(m, "cardinal_quintic_hermite",
       "Quintic Hermite interpolant of y[i], dydx[i] and d2ydx2[i] at the equispaced x0 + i*dx." BMPY_URL("quintic_hermite"));
   cqh.def("__init__", [](cardinal_quintic_h* self, nb::handle y, nb::handle dydx, nb::handle d2ydx2, double x0, double dx) {
          vec ys = copy_vec(y, "y"), ds = copy_vec(dydx, "dydx"), d2s = copy_vec(d2ydx2, "d2ydx2");
          std::size_t n = ys.size();
          new (self) cardinal_quintic_h{checked([&] { return bi::cardinal_quintic_hermite<vec>(std::move(ys), std::move(ds), std::move(d2s), x0, dx); }),
                                          n, cardinal("x0", x0, "dx", dx)};
       }, nb::arg(), nb::arg(), nb::arg(), nb::arg(), nb::arg(), nb::sig("def __init__(self, y, dydx, d2ydx2, x0, dx, /) -> None"));
   add_calls(cqh);

   nb::class_<septic_h> sh(m, "septic_hermite",
       "Septic Hermite interpolant of y and its first three derivatives at strictly increasing x."
       "\n\nhttps://www.boost.org/doc/libs/release/libs/math/doc/html/interpolation.html");
   sh.def("__init__", [](septic_h* self, nb::handle x, nb::handle y, nb::handle dydx, nb::handle d2ydx2, nb::handle d3ydx3) {
          vec xs = copy_vec(x, "x"), ys = copy_vec(y, "y"), ds = copy_vec(dydx, "dydx"), d2s = copy_vec(d2ydx2, "d2ydx2"), d3s = copy_vec(d3ydx3, "d3ydx3");
          std::size_t n = xs.size();
          new (self) septic_h{checked([&] { return bi::septic_hermite<vec>(std::move(xs), std::move(ys), std::move(ds), std::move(d2s), std::move(d3s)); }), n, ""};
       }, nb::arg(), nb::arg(), nb::arg(), nb::arg(), nb::arg(), nb::sig("def __init__(self, x, y, dydx, d2ydx2, d3ydx3, /) -> None"));
   add_calls(sh);

   nb::class_<cardinal_septic_h> csh(m, "cardinal_septic_hermite",
       "Septic Hermite interpolant of y[i] and its first three derivatives at the equispaced x0 + i*dx."
       "\n\nhttps://www.boost.org/doc/libs/release/libs/math/doc/html/interpolation.html");
   csh.def("__init__", [](cardinal_septic_h* self, nb::handle y, nb::handle dydx, nb::handle d2ydx2, nb::handle d3ydx3, double x0, double dx) {
          vec ys = copy_vec(y, "y"), ds = copy_vec(dydx, "dydx"), d2s = copy_vec(d2ydx2, "d2ydx2"), d3s = copy_vec(d3ydx3, "d3ydx3");
          std::size_t n = ys.size();
          new (self) cardinal_septic_h{checked([&] { return bi::cardinal_septic_hermite<vec>(std::move(ys), std::move(ds), std::move(d2s), std::move(d3s), x0, dx); }),
                                          n, cardinal("x0", x0, "dx", dx)};
       }, nb::arg(), nb::arg(), nb::arg(), nb::arg(), nb::arg(), nb::arg(),
       nb::sig("def __init__(self, y, dydx, d2ydx2, d3ydx3, x0, dx, /) -> None"));
   add_calls(csh);

   nb::class_<bilinear> bl(m, "bilinear_uniform",
       "Bilinear interpolant of a 2-D field, at least 2x2, field[i, j] = f(x0 + j*dx, y0 + i*dy): rows along y."
       BMPY_URL("bilinear_uniform"));
   bl.def("__init__", [](bilinear* self, nb::handle field, double dx, double dy, double x0, double y0) {
          std::size_t shape[2];
          vec v = copy_array(field, 2, "field", shape);
          std::string params = ", rows=" + std::to_string(shape[0]) + cardinal("dx", dx, "dy", dy) + cardinal("x0", x0, "y0", y0);
          auto f = checked([&] { return bi::bilinear_uniform<vec>(std::move(v), shape[0], shape[1], dx, dy, x0, y0); });
          new (self) bilinear{{std::move(f), shape[0] * shape[1], params}, x0, x0 + (shape[1] - 1) * dx, y0, y0 + (shape[0] - 1) * dy};
       }, nb::arg(), nb::arg("dx") = 1.0, nb::arg("dy") = 1.0, nb::arg("x0") = 0.0, nb::arg("y0") = 0.0,
       nb::sig("def __init__(self, field, /, dx=1.0, dy=1.0, x0=0.0, y0=0.0) -> None"))
     .def("__call__", [](bilinear const& b, nb::handle x, nb::handle y, int threads) {
          auto f = [&b]<class T>(double u, double v) {
             if (!(u >= b.x0 && u <= b.x1 && v >= b.y0 && v <= b.y1))
                throw std::domain_error("bilinear_uniform: (x, y) is outside the grid");
             return b.f(u, v);
          };
          return call_array<false, Real, Real>("__call__", f, {x, y}, {"x", "y"}, threads);
       }, nb::arg(), nb::arg(), nb::kw_only(), nb::arg("threads") = 1,
       nb::sig("def __call__(self, x, y, /, *, threads=1)"), "The interpolant at (x, y), broadcast elementwise.");
   add_repr(bl);

   nb::class_<catmull> cr(m, "catmull_rom",
       "Catmull-Rom spline through the rows of an (n, d) array, n >= 4, d = 1 to 4; t runs from 0 to max_parameter()."
       BMPY_URL("catmull_rom"));
   cr.def("__init__", [](catmull* self, nb::handle points, bool closed, double alpha) {
          std::string params = std::string(", closed=") + (closed ? "True" : "False") + ", alpha=" + py_float_repr(alpha);
          new (self) catmull(make_curve<catmull_rom_d>(points, "points", params, [&]<std::size_t D>(std::vector<point<D>>&& p) {
             return catmull_rom_d<D>(std::move(p), closed, alpha);
          }));
       }, nb::arg(), nb::arg("closed") = false, nb::arg("alpha") = 0.5,
       nb::sig("def __init__(self, points, /, closed=False, alpha=0.5) -> None"))
     .def("max_parameter", [](catmull const& c) { return std::visit([](auto const& f) { return f.max_parameter(); }, c.f); },
          "The parameter at the end of the curve; it starts at 0.");
   def_points(cr, "__call__", curve_value, "The point on the curve at parameter t; shape t.shape + (d,).");
   def_points(cr, "prime", curve_prime, "The curve's tangent dp/dt at t; shape t.shape + (d,).");
   add_repr(cr);

   nb::class_<bezier> bz(m, "bezier_polynomial",
       "Bezier curve, t in [0, 1], of the rows of an (n, d) array of control points, n >= 2, d = 1 to 4."
       BMPY_URL("bezier_polynomial"));
   bz.def("__init__", [](bezier* self, nb::handle points) {
          new (self) bezier(make_curve<bezier_d>(points, "control_points", "", []<std::size_t D>(std::vector<point<D>>&& p) {
             return bezier_d<D>(std::move(p));
          }));
       }, nb::arg(), nb::sig("def __init__(self, control_points, /) -> None"));
   def_points(bz, "__call__", bezier_checked(curve_value), "The point B(t) on the curve; shape t.shape + (d,).");
   def_points(bz, "prime", bezier_checked(curve_prime), "The tangent B'(t); shape t.shape + (d,).");
   add_repr(bz);

   nb::class_<vbary> vb(m, "vector_barycentric_rational",
       "barycentric_rational of the rows of a (len(t), d) array y, d = 1 to 4, at strictly increasing t."
       BMPY_URL("vector_barycentric"));
   vb.def("__init__", [](vbary* self, nb::handle t, nb::handle y, std::size_t order) {
          vec ts = copy_vec(t, "t");
          // The C++ only asserts these.
          if (std::adjacent_find(ts.begin(), ts.end(), std::greater_equal<>()) != ts.end())
             throw nb::value_error("vector_barycentric_rational: t must be strictly increasing");
          new (self) vbary(make_curve<vector_barycentric_d>(y, "y", ", approximation_order=" + std::to_string(order),
             [&]<std::size_t D>(std::vector<point<D>>&& p) {
                if (p.size() != ts.size() || order >= p.size())
                   throw std::domain_error("vector_barycentric_rational: need len(y) == len(t) > approximation_order");
                return vector_barycentric_d<D>(std::move(ts), std::move(p), order);
             }));
       }, nb::arg(), nb::arg(), nb::arg("approximation_order") = 3,
       nb::sig("def __init__(self, t, y, /, approximation_order=3) -> None"));
   def_points(vb, "__call__", curve_value, "The interpolant at t; shape t.shape + (d,).");
   def_points(vb, "prime", curve_prime, "Its derivative at t; shape t.shape + (d,).");
   add_repr(vb);
}

}  // namespace bmpy
