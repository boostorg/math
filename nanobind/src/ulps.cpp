// Copyright 2026 Nick Thompson
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// ulps plots and heatmaps for boost.math.ulps: the double implementation against the same
// function at cpp_bin_float_50, or against a Python reference.

#include <boost/math/tools/ulps_heatmap.hpp>
#include <boost/math/tools/ulps_plot.hpp>

#include <limits>
#include <memory>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <tuple>

#include <nanobind/stl/tuple.h>
#include <nanobind/stl/vector.h>

#include "bindings.hpp"
#include "ulps.hpp"

namespace nb = bmpy::nb;

namespace {

using bmpy::ulps::Precise;
using plot1 = boost::math::tools::ulps_plot<std::function<Precise(Precise const&)>, Precise, double>;
using heatmap2 = boost::math::tools::ulps_heatmap<std::function<Precise(Precise const&, Precise const&)>, Precise, double>;

bmpy::ulps::registry const& functions() {
   static bmpy::ulps::registry const r = [] {
      bmpy::ulps::registry r;
      bmpy::ulps::register_gamma(r);
      bmpy::ulps::register_special(r);
      bmpy::ulps::register_elliptic(r);
      return r;
   }();
   return r;
}

template <class Map>
auto const& lookup(Map const& m, std::string const& name) {
   auto it = m.find(name);
   if (it == m.end()) {
      throw nb::value_error(("no ulps reference for '" + name + "'").c_str());
   }
   return it->second;
}

// An RGB8 PNG, by boost.math.ulps._png with Python's zlib, so that the module links no zlib of its own.
std::vector<unsigned char> encode_png(std::vector<unsigned char> const& rgb, std::size_t width, std::size_t height) {
   nb::gil_scoped_acquire acquire;
   auto png = nb::cast<nb::bytes>(
       nb::module_::import_("boost.math.ulps").attr("_png")(nb::bytes(rgb.data(), rgb.size()), width, height));
   auto const* p = reinterpret_cast<unsigned char const*>(png.c_str());
   return {p, p + png.size()};
}

// A Python reference, called one at a time: not merely one holder of the GIL at a time, since the interpreter switches
// threads in the middle of a call, and mpmath keeps its working precision in a global, which its functions raise and
// restore. Arguments go as floats where they are doubles, and otherwise (the condition number's finite differences) as
// 55-digit strings, which boost.math.ulps turns into mpmath numbers. It returns a float or a decimal string; anything
// else, or an exception, is NaN. KeyboardInterrupt, SystemExit and a failed signal handler stop the run instead: every
// later call returns NaN at once, and rethrow() raises the error once the computation is over.
class python_reference {
public:
   explicit python_reference(nb::handle f) : f_(f) {}

   template <class... Args>
   Precise operator()(Args const&... args) {
      std::lock_guard<std::mutex> lock(mutex_);
      if (abort_) {
         return std::numeric_limits<Precise>::quiet_NaN();
      }
      nb::gil_scoped_acquire acquire;
      bool signal = false;
      try {
         signal = PyErr_CheckSignals() != 0;
         if (signal) {
            throw nb::python_error();
         }
         return from_python(f_(to_python(args)...));
      } catch (nb::python_error& e) {
         if (!signal && !e.matches(PyExc_KeyboardInterrupt) && !e.matches(PyExc_SystemExit)) {
            return std::numeric_limits<Precise>::quiet_NaN();
         }
         abort_ = true;
         error_ = std::make_unique<nb::python_error>(std::move(e));
         return std::numeric_limits<Precise>::quiet_NaN();
      } catch (std::exception const&) {
         return std::numeric_limits<Precise>::quiet_NaN();
      }
   }

   // With the GIL held.
   void rethrow() {
      if (error_) {
         nb::python_error e = std::move(*error_);
         error_.reset();
         throw e;
      }
   }

private:
   static nb::object to_python(Precise const& x) {
      double const d = static_cast<double>(x);
      if (Precise(d) == x) {
         return nb::float_(d);
      }
      return nb::str(x.str(55, std::ios_base::scientific).c_str());
   }

   static Precise from_python(nb::handle r) {
      if (nb::isinstance<nb::float_>(r) || nb::isinstance<nb::int_>(r)) {
         return Precise(nb::cast<double>(r));
      }
      if (nb::isinstance<nb::str>(r)) {
         return Precise(nb::cast<std::string>(r).c_str());
      }
      return std::numeric_limits<Precise>::quiet_NaN();
   }

   nb::handle f_;
   std::mutex mutex_;
   bool abort_ = false;
   std::unique_ptr<nb::python_error> error_;
};

std::tuple<std::string, std::string> plot(std::string const& name, double a, double b, std::size_t samples, bool log_x,
                                          std::string const& title, nb::handle reference, double clip, bool show_function,
                                          unsigned threads) {
   auto const& f = lookup(functions().unary, name);
   std::function<Precise(Precise const&)> precise = f.precise;
   python_reference py(reference);
   if (!reference.is_none()) {
      precise = [&py](Precise const& x) { return py(x); };
   }
   std::string svg;
   std::string summary;
   {
      nb::gil_scoped_release release;
      plot1 p(precise, a, b, samples, false, 20261013, log_x, threads);
      p.add_fn(f.coarse).title(title).show_function(show_function).png_encoder(encode_png);
      if (clip > 0) {
         p.clip(clip);
      }
      std::ostringstream os;
      os << p;
      svg = os.str();
      summary = p.summary();
   }
   py.rethrow();
   return {svg, summary};
}

std::tuple<std::string, std::string, std::string> heatmap(std::string const& name, double x_min, double x_max, double y_min,
                                                          double y_max, std::size_t width, std::size_t height, bool log_x,
                                                          bool log_y, std::string const& title, std::string const& x_label,
                                                          std::string const& y_label, bool log_scale, double clip,
                                                          nb::handle reference, bool show_function, unsigned threads) {
   auto const& f = lookup(functions().binary, name);
   std::function<Precise(Precise const&, Precise const&)> precise = f.precise;
   python_reference py(reference);
   if (!reference.is_none()) {
      precise = [&py](Precise const& x, Precise const& y) { return py(x, y); };
   }
   std::tuple<std::string, std::string, std::string> result;
   {
      nb::gil_scoped_release release;
      heatmap2 h(precise, x_min, x_max, y_min, y_max, width, height, log_x, log_y, threads);
      h.add_fn(f.coarse).title(title).x_label(x_label).y_label(y_label).log_scale(log_scale).show_function(show_function)
          .png_encoder(encode_png);
      if (clip > 0) {
         h.clip(clip);
      }
      result = {h.ulps_svg(), h.ulps_envelope_svg(), h.summary()};
   }
   py.rethrow();
   return result;
}

}  // namespace

void bmpy::bind_ulps(nb::module_& m) {
   m.def("_ulps_plot", &plot, nb::arg("name"), nb::arg("a"), nb::arg("b"), nb::arg("samples"), nb::arg("log_x"),
         nb::arg("title"), nb::arg("reference").none(), nb::arg("clip"), nb::arg("show_function"), nb::arg("threads"),
         "(svg, summary) of an ulps_plot; see boost.math.ulps.plot");
   m.def("_ulps_heatmap", &heatmap, nb::arg("name"), nb::arg("x_min"), nb::arg("x_max"), nb::arg("y_min"),
         nb::arg("y_max"), nb::arg("width"), nb::arg("height"), nb::arg("log_x"), nb::arg("log_y"), nb::arg("title"),
         nb::arg("x_label"), nb::arg("y_label"), nb::arg("log_scale"), nb::arg("clip"), nb::arg("reference").none(),
         nb::arg("show_function"), nb::arg("threads"),
         "(ulps svg, ulps/envelope svg, summary) of an ulps_heatmap; see boost.math.ulps.heatmap");
   m.def("_ulps_functions", [] {
      std::vector<std::string> unary, binary;
      for (auto const& kv : functions().unary) unary.push_back(kv.first);
      for (auto const& kv : functions().binary) binary.push_back(kv.first);
      return std::make_tuple(unary, binary);
   }, "The names boost.math.ulps.plot and heatmap accept.");
}
