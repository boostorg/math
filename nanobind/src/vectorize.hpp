// Copyright 2026 Nick Thompson
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// Elementwise (NumPy-broadcasting) wrappers for scalar Boost.Math functions.
//
// def_vec<Markers...>(m, "name", f, {"arg", ...}) registers two overloads:
//   1. a fast scalar overload taking Python floats/ints and returning a float,
//      which lets Boost.Math exceptions propagate (ValueError, OverflowError, ...);
//   2. a generic overload taking anything np.asarray understands. Arguments are
//      broadcast against each other; real arguments are evaluated in float when
//      np.result_type of the real arguments is float32 (or float16), otherwise in
//      double. 0-d inputs return a NumPy scalar and raise like the scalar path; for
//      arrays an element whose evaluation throws becomes NaN and one
//      RuntimeWarning summarises the failures.
//
// Each marker describes one argument: Real is the floating-point type T, Int<I>
// is an integer converted (with a range check) to I. f is a C++20 templated
// lambda `[]<class T>(...)` called as f.template operator()<T>(args...).

#ifndef BOOST_MATH_NANOBIND_VECTORIZE_HPP
#define BOOST_MATH_NANOBIND_VECTORIZE_HPP

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/complex.h>
#include <nanobind/stl/string.h>

#include <array>
#include <complex>
#include <cstdint>
#include <exception>
#include <limits>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

namespace bmpy {

namespace nb = nanobind;

struct Real {};
template <class I> struct Int {};

template <class M, class T> struct marker;
template <class T> struct marker<Real, T> {
   using storage = T;
   using arg = T;
   using scalar = double;
   static constexpr bool is_real = true;
};
template <class I, class T> struct marker<Int<I>, T> {
   using storage = std::int64_t;
   using arg = I;
   using scalar = std::int64_t;
   static constexpr bool is_real = false;
};

template <class I>
inline I checked_int(std::int64_t v) {
   if constexpr (std::is_signed_v<I>) {
      if (v < static_cast<std::int64_t>((std::numeric_limits<I>::min)()) || v > static_cast<std::int64_t>((std::numeric_limits<I>::max)()))
         throw std::domain_error("integer argument out of range");
   } else {
      if (v < 0 || static_cast<std::uint64_t>(v) > static_cast<std::uint64_t>((std::numeric_limits<I>::max)()))
         throw std::domain_error("integer argument out of range");
   }
   return static_cast<I>(v);
}

template <class M, class T>
inline typename marker<M, T>::arg convert_scalar(typename marker<M, T>::storage v) {
   if constexpr (marker<M, T>::is_real)
      return v;
   else
      return checked_int<typename marker<M, T>::arg>(v);
}

template <class R>
inline R nan_of() {
   if constexpr (std::is_floating_point_v<R>)
      return std::numeric_limits<R>::quiet_NaN();
   else {
      using V = typename R::value_type;
      return R(std::numeric_limits<V>::quiet_NaN(), std::numeric_limits<V>::quiet_NaN());
   }
}

inline nb::module_& numpy() {
   static nb::module_ np = nb::module_::import_("numpy");
   return np;
}

// Converts obj to an ndarray of the storage type of marker M.
template <class M, class T>
nb::object as_array(nb::handle obj, const char* argname) {
   nb::module_& np = numpy();
   if constexpr (marker<M, T>::is_real) {
      nb::object dt = std::is_same_v<T, float> ? np.attr("float32") : np.attr("float64");
      return np.attr("asarray")(obj, nb::arg("dtype") = dt);
   } else {
      nb::object a = np.attr("asarray")(obj);
      std::string kind = nb::cast<std::string>(a.attr("dtype").attr("kind"));
      if (kind != "i" && kind != "u" && kind != "b")
         throw nb::type_error((std::string("argument '") + argname + "' must be an integer or an integer array").c_str());
      return a.attr("astype")(np.attr("int64"));
   }
}

template <class S>
using view = nb::ndarray<const S, nb::device::cpu>;

template <class T, class... Ms, class F, std::size_t... Is>
nb::object call_array_impl(const char* name, const F& f, const std::array<nb::handle, sizeof...(Ms)>& objs,
                           const std::array<const char*, sizeof...(Ms)>& names, std::index_sequence<Is...>) {
   constexpr std::size_t N = sizeof...(Ms);
   std::array<nb::object, N> arrays{as_array<Ms, T>(objs[Is], names[Is])...};
   std::tuple<view<typename marker<Ms, T>::storage>...> views{nb::cast<view<typename marker<Ms, T>::storage>>(arrays[Is])...};

   // Broadcast shape (NumPy rules).
   std::size_t ndim = 0;
   std::array<std::size_t, N> nd{std::get<Is>(views).ndim()...};
   for (std::size_t k : nd) ndim = k > ndim ? k : ndim;
   std::vector<std::size_t> shape(ndim, 1);
   std::array<std::vector<std::int64_t>, N> strides;
   auto fold_shape = [&](std::size_t a, const auto& v) {
      std::size_t off = ndim - v.ndim();
      for (std::size_t d = 0; d < v.ndim(); ++d) {
         std::size_t s = v.shape(d);
         if (s != 1) {
            if (shape[off + d] != 1 && shape[off + d] != s)
               throw nb::value_error((std::string(name) + ": operands could not be broadcast together").c_str());
            shape[off + d] = s;
         }
      }
      (void)a;
   };
   (fold_shape(Is, std::get<Is>(views)), ...);
   auto fold_strides = [&](std::size_t a, const auto& v) {
      strides[a].assign(ndim, 0);
      std::size_t off = ndim - v.ndim();
      for (std::size_t d = 0; d < v.ndim(); ++d)
         strides[a][off + d] = v.shape(d) == 1 ? 0 : v.stride(d);
   };
   (fold_strides(Is, std::get<Is>(views)), ...);
   std::size_t size = 1;
   for (std::size_t s : shape) size *= s;

   using R = decltype(f.template operator()<T>(std::declval<typename marker<Ms, T>::arg>()...));
   std::tuple<const typename marker<Ms, T>::storage*...> base{std::get<Is>(views).data()...};

   R* out = new R[size ? size : 1];
   nb::capsule owner(out, [](void* p) noexcept { delete[] static_cast<R*>(p); });

   if (ndim == 0) {
      // Scalar semantics: let exceptions propagate.
      out[0] = f.template operator()<T>(convert_scalar<Ms, T>(*std::get<Is>(base))...);
   } else {
      std::size_t nerr = 0;
      std::string first;
      {
         nb::gil_scoped_release release;
         std::vector<std::size_t> idx(ndim, 0);
         std::array<std::int64_t, N> pos{};
         for (std::size_t i = 0; i < size; ++i) {
            try {
               out[i] = f.template operator()<T>(convert_scalar<Ms, T>(std::get<Is>(base)[pos[Is]])...);
            } catch (const std::exception& e) {
               if constexpr (!std::is_integral_v<R>)
                  out[i] = nan_of<R>();
               if (nerr++ == 0) first = e.what();
            }
            // Advance the odometer.
            for (std::size_t d = ndim; d-- > 0;) {
               if (++idx[d] < shape[d]) {
                  for (std::size_t a = 0; a < N; ++a) pos[a] += strides[a][d];
                  break;
               }
               for (std::size_t a = 0; a < N; ++a) pos[a] -= strides[a][d] * static_cast<std::int64_t>(shape[d] - 1);
               idx[d] = 0;
            }
         }
      }
      // Integer results have no NaN to mark a failed element, so raise instead.
      if constexpr (std::is_integral_v<R>) {
         if (nerr)
            throw std::domain_error(std::string(name) + ": " + first);
      }
      if (nerr) {
         std::string msg = std::string(name) + ": " + std::to_string(nerr) + " of " + std::to_string(size) +
                           " evaluations raised and were set to NaN; first error: " + first;
         if (PyErr_WarnEx(PyExc_RuntimeWarning, msg.c_str(), 1) != 0)
            throw nb::python_error();
      }
   }
   nb::object result = nb::cast(nb::ndarray<nb::numpy, R>(out, ndim, shape.data(), owner));
   if (ndim == 0)
      return result.attr("__getitem__")(nb::tuple());
   return result;
}

template <bool AllowFloat, class... Ms, class F>
nb::object call_array(const char* name, const F& f, const std::array<nb::handle, sizeof...(Ms)>& objs,
                      const std::array<const char*, sizeof...(Ms)>& names) {
   // Choose the evaluation type from the real arguments only.
   bool use_float = false;
   constexpr bool any_real = (marker<Ms, double>::is_real || ...);
   if constexpr (any_real && AllowFloat) {
      nb::module_& np = numpy();
      nb::list reals;
      std::size_t a = 0;
      // Python scalars stay weakly typed (NEP 50); anything else goes through
      // np.asarray, since result_type would read a list as a dtype spec.
      auto add = [&](nb::handle h) {
         if (PyFloat_Check(h.ptr()) || PyLong_Check(h.ptr()) || PyComplex_Check(h.ptr()))
            reals.append(h);
         else
            reals.append(np.attr("asarray")(h));
      };
      ((marker<Ms, double>::is_real ? add(objs[a++]) : (void)a++), ...);
      nb::object rt = np.attr("result_type")(*reals);
      std::string kind = nb::cast<std::string>(rt.attr("kind"));
      if (kind == "c")
         throw nb::type_error((std::string(name) + ": complex arguments are not supported").c_str());
      if (kind == "f" && nb::cast<int>(rt.attr("itemsize")) <= 4)
         use_float = true;
   }
   nb::object result = use_float
      ? call_array_impl<float, Ms...>(name, f, objs, names, std::index_sequence_for<Ms...>{})
      : call_array_impl<double, Ms...>(name, f, objs, names, std::index_sequence_for<Ms...>{});
   // Plain Python numbers in, plain Python number out, as from the math module.
   bool python_scalars = true;
   for (nb::handle h : objs)
      python_scalars = python_scalars && (PyFloat_Check(h.ptr()) || PyLong_Check(h.ptr()));
   if (python_scalars)
      return result.attr("item")();
   return result;
}

template <class... Ms, class F, std::size_t... Is>
void def_vec_impl(nb::module_& m, const char* name, F f, std::array<const char*, sizeof...(Ms)> names,
                  const char* doc, std::index_sequence<Is...>) {
   auto scalar = [f](typename marker<Ms, double>::scalar... xs) {
      return f.template operator()<double>(convert_scalar<Ms, double>(xs)...);
   };
   m.def(name, scalar, nb::arg(names[Is])..., doc);
   auto generic = [f, names, name](std::conditional_t<true, nb::object, decltype(Is)>... xs) {
      return call_array<true, Ms...>(name, f, std::array<nb::handle, sizeof...(Ms)>{xs...}, names);
   };
   m.def(name, generic, nb::arg(names[Is])...);
}

// Registers `name` in m; see the comment at the top of this file.
template <class... Ms, class F>
void def_vec(nb::module_& m, const char* name, F f, std::array<const char*, sizeof...(Ms)> names, const char* doc = "") {
   def_vec_impl<Ms...>(m, name, f, names, doc, std::index_sequence_for<Ms...>{});
}

}  // namespace bmpy

#endif
