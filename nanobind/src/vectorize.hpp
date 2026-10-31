// Copyright 2026 Nick Thompson
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// Elementwise (NumPy-broadcasting) wrappers for scalar Boost.Math functions.
//
// def_vec<Markers...>(m, "name", f, {"arg", ...}, doc) registers name(arg, ..., /, *, threads=1), taking
// anything np.asarray understands. Arguments are broadcast against each other; real arguments are
// evaluated in float when np.result_type of the real arguments is float32 (or float16), otherwise in
// double. Python numbers in give a Python number out, as from the math module; NumPy scalars and 0-d
// arrays give a NumPy scalar, and arrays an array. Python and NumPy scalars skip the array machinery.
// Errors follow boost.math.errstate (errors.hpp): by default scalars raise, and for arrays a failed
// element becomes NaN (an overflow, an infinity) and one RuntimeWarning summarises the failures.
// threads > 1 (0: all cores) splits large arrays across that many threads.
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

#include <algorithm>
#include <array>
#include <complex>
#include <cstdint>
#include <exception>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <thread>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include "errors.hpp"

// Ends a docstring with the link to the docs page math_toolkit/page.html.
#define BMPY_URL(page) "\n\nhttps://www.boost.org/doc/libs/release/libs/math/doc/html/math_toolkit/" page ".html"

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

// The result of a failed evaluation.
template <class R>
inline R nan_of() {
   if constexpr (std::is_floating_point_v<R>)
      return std::numeric_limits<R>::quiet_NaN();
   else if constexpr (requires { std::numeric_limits<typename R::value_type>::quiet_NaN(); }) {
      using V = typename R::value_type;
      return R(std::numeric_limits<V>::quiet_NaN(), std::numeric_limits<V>::quiet_NaN());
   } else
      return R{};
}

inline nb::module_& numpy() {
   static nb::module_ np = nb::module_::import_("numpy");
   return np;
}

// f() with the errors handled as for a scalar argument; see errors.hpp.
template <class F>
auto call_scalar(char const* name, F const& f) {
   using R = decltype(f());
   error_actions actions = current_error_actions;
   error_record errors;
   R r = [&] {
      recording rec(errors);
      try {
         return f();
      } catch (std::exception const& e) {
         record_exception(e);
         return nan_of<R>();
      }
   }();
   if (errors.any())
      report(name, errors, 1, true, std::is_integral_v<R>, actions);
   return r;
}

// The NumPy scalar type of R (np.float64, ...).
template <class R>
nb::handle numpy_type() {
   static nb::handle t = [] {
      std::string name;
      if constexpr (std::is_integral_v<R>)
         name = (std::is_signed_v<R> ? "int" : "uint") + std::to_string(8 * sizeof(R));
      else if constexpr (std::is_floating_point_v<R>)
         name = "float" + std::to_string(8 * sizeof(R));
      else
         name = "complex" + std::to_string(8 * sizeof(R));
      return nb::object(numpy().attr(name.c_str())).release();
   }();
   return t;
}

inline bool is_numpy_integer(nb::handle h) {
   static nb::handle t = nb::object(numpy().attr("integer")).release();
   return PyObject_IsInstance(h.ptr(), t.ptr()) == 1;
}

// The scalar fast path. What the arguments were: NumPy scalars give a NumPy scalar out, and a float32
// among the real arguments makes the evaluation type float.
struct scalar_kinds {
   bool numpy = false, float32 = false, float64 = false, integer = false;
};

// Reads h into v if it is a Python or NumPy scalar of the kind marker M takes.
template <class M>
bool read_scalar(nb::handle h, typename marker<M, double>::scalar& v, scalar_kinds& k) {
   PyObject* o = h.ptr();
   PyTypeObject* t = Py_TYPE(o);
   if constexpr (marker<M, double>::is_real) {
      if (t == &PyFloat_Type) {
         v = PyFloat_AsDouble(o);
         return true;
      }
      if (t == &PyLong_Type) {
         v = PyLong_AsDouble(o);
      } else {
         bool f64 = h.type().is(numpy_type<double>()), f32 = !f64 && h.type().is(numpy_type<float>());
         if (!f64 && !f32 && !is_numpy_integer(h))
            return false;
         k.numpy = true;
         (f64 ? k.float64 : f32 ? k.float32 : k.integer) = true;
         v = PyFloat_AsDouble(o);
      }
   } else {
      nb::object i;
      if (t != &PyLong_Type) {
         if (!is_numpy_integer(h))
            return false;
         k.numpy = true;
         i = nb::steal(PyNumber_Index(o));
         o = i.ptr();
      }
      v = o ? PyLong_AsLongLong(o) : -1;
   }
   if (v == -1 && PyErr_Occurred()) {
      PyErr_Clear();
      return false;
   }
   return true;
}

// Evaluates f at Python or NumPy scalars, storing the result in `result`; false for anything else.
template <bool AllowFloat, class... Ms, class F, std::size_t... Is>
bool call_scalars(nb::object& result, char const* name, F const& f, std::array<nb::handle, sizeof...(Ms)> const& objs,
                  std::index_sequence<Is...>) {
   scalar_kinds k;
   std::tuple<typename marker<Ms, double>::scalar...> v;
   // NumPy promotes float32 with some integer types to float64 and not with others; result_type decides.
   if (!(read_scalar<Ms>(objs[Is], std::get<Is>(v), k) && ...) || (k.float32 && k.integer))
      return false;
   auto finish = [&](auto r) { result = k.numpy ? numpy_type<decltype(r)>()(r) : nb::cast(r); };
   if (AllowFloat && k.float32 && !k.float64)
      finish(call_scalar(name, [&] { return f.template operator()<float>(convert_scalar<Ms, float>(std::get<Is>(v))...); }));
   else
      finish(call_scalar(name, [&] { return f.template operator()<double>(convert_scalar<Ms, double>(std::get<Is>(v))...); }));
   return true;
}

// The number of threads to split `size` evaluations across.
inline std::size_t thread_count(int threads, std::size_t size) {
   std::size_t n = threads > 0 ? threads : (std::max)(1u, std::thread::hardware_concurrency());
   // Fewer than this many evaluations per thread are not worth a thread.
   constexpr std::size_t min_chunk = 4096;
   return (std::max<std::size_t>)(1, (std::min)(n, size / min_chunk));
}

// "def name(a, b, /, *, threads=1)", or with self first for a method.
template <std::size_t N>
std::string vec_signature(char const* name, std::array<char const*, N> const& names, bool method = false) {
   std::string s = std::string("def ") + name + (method ? "(self, " : "(");
   for (char const* n : names)
      s += std::string(n) + ", ";
   return s + "/, *, threads=1)";
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
                           const std::array<const char*, sizeof...(Ms)>& names, int threads, bool release_gil,
                           std::index_sequence<Is...>) {
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

   auto eval = [&](std::size_t i, std::array<std::int64_t, N> const& pos) {
      try {
         out[i] = f.template operator()<T>(convert_scalar<Ms, T>(std::get<Is>(base)[pos[Is]])...);
      } catch (std::exception const& e) {
         record_exception(e);
         out[i] = nan_of<R>();
      }
      active_record->next();
   };
   // Evaluates elements [lo, hi).
   auto run = [&](std::size_t lo, std::size_t hi) {
      std::vector<std::size_t> idx(ndim, 0);
      std::array<std::int64_t, N> pos{};
      for (std::size_t d = ndim, r = lo; d-- > 0; r /= shape[d]) {
         idx[d] = r % shape[d];
         for (std::size_t a = 0; a < N; ++a) pos[a] += strides[a][d] * static_cast<std::int64_t>(idx[d]);
      }
      for (std::size_t i = lo; i < hi; ++i) {
         eval(i, pos);
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
   };

   error_actions actions = current_error_actions;
   error_record errors;
   if (size != 0) {
      std::optional<nb::gil_scoped_release> release;
      if (release_gil)
         release.emplace();
      std::size_t n = thread_count(threads, size);
      std::vector<error_record> records(n);
      std::vector<std::thread> pool;
      auto work = [&](std::size_t k) {
         recording rec(records[k]);
         run(size * k / n, size * (k + 1) / n);
      };
      for (std::size_t k = 1; k < n; ++k) pool.emplace_back(work, k);
      work(0);
      for (std::thread& t : pool) t.join();
      for (error_record const& r : records) errors.merge(r);
   }
   // 0-d arrays behave like scalars.
   if (errors.any())
      report(name, errors, size, ndim == 0, std::is_integral_v<R>, actions);
   nb::object result = nb::cast(nb::ndarray<nb::numpy, R>(out, ndim, shape.data(), owner));
   if (ndim == 0)
      return result.attr("__getitem__")(nb::tuple());
   return result;
}

template <bool AllowFloat, class... Ms, class F>
nb::object call_array(const char* name, const F& f, const std::array<nb::handle, sizeof...(Ms)>& objs,
                      const std::array<const char*, sizeof...(Ms)>& names, int threads = 1, bool release_gil = true) {
   if (threads < 0)
      throw nb::value_error((std::string(name) + ": threads must be non-negative").c_str());
   nb::object scalar;
   if (call_scalars<AllowFloat, Ms...>(scalar, name, f, objs, std::index_sequence_for<Ms...>{}))
      return scalar;
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
      ? call_array_impl<float, Ms...>(name, f, objs, names, threads, release_gil, std::index_sequence_for<Ms...>{})
      : call_array_impl<double, Ms...>(name, f, objs, names, threads, release_gil, std::index_sequence_for<Ms...>{});
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
   auto fn = [f, names, name](std::conditional_t<true, nb::handle, decltype(Is)>... xs, int threads) {
      return call_array<true, Ms...>(name, f, std::array<nb::handle, sizeof...(Ms)>{xs...}, names, threads);
   };
   m.def(name, fn, std::conditional_t<true, nb::arg, decltype(Is)>()..., nb::kw_only(), nb::arg("threads") = 1,
         nb::sig(vec_signature(name, names).c_str()), doc);
}

// Registers `name` in m; see the comment at the top of this file.
template <class... Ms, class F>
void def_vec(nb::module_& m, const char* name, F f, std::array<const char*, sizeof...(Ms)> names, const char* doc = "") {
   def_vec_impl<Ms...>(m, name, f, names, doc, std::index_sequence_for<Ms...>{});
}

// cls.name(argname, /, *, threads=1): g(self, x) elementwise in double. With release_gil false, the GIL is held
// throughout, for a class whose methods may change it.
template <class C, class G>
void def_method_vec(nb::class_<C>& cls, const char* name, const char* argname, G g, const char* doc = "",
                    bool release_gil = true) {
   cls.def(
       name,
       [g, name, argname, release_gil](const C& c, nb::handle x, int threads) {
          auto f = [&c, &g]<class T>(double v) { return g(c, v); };
          return call_array<false, Real>(name, f, {x}, {argname}, threads, release_gil);
       },
       nb::arg(), nb::kw_only(), nb::arg("threads") = 1,
       nb::sig(vec_signature(name, std::array<const char*, 1>{argname}, true).c_str()), doc);
}

// A NumPy array that owns v.
template <class V>
nb::ndarray<nb::numpy, typename V::value_type, nb::ndim<1>> vector_to_array(V v) {
   auto* data = new V(std::move(v));
   nb::capsule owner(data, [](void* p) noexcept { delete static_cast<V*>(p); });
   return {data->data(), {data->size()}, owner};
}

// x as a contiguous 1D double vector (np.asarray(x, float).ravel()).
inline std::vector<double> array_to_vector(nb::handle x) {
   nb::object a = numpy().attr("ascontiguousarray")(x, nb::arg("dtype") = numpy().attr("float64")).attr("reshape")(-1);
   auto v = nb::cast<nb::ndarray<const double, nb::ndim<1>, nb::c_contig, nb::device::cpu>>(a);
   return std::vector<double>(v.data(), v.data() + v.shape(0));
}

// A Python callable as a C++ function of X... returning R, called with the GIL held. Boost.Math errors in the
// Python code it runs raise, even inside a recording.
template <class R, class... X>
struct py_callable {
   nb::callable f;
   R operator()(X... x) const {
      recording suspended(nullptr);
      return nb::cast<R>(f(x...));
   }
};

// The named tuple type m.name, with the given fields (a space-separated string), docstring and defaults.
inline nb::object named_tuple(nb::module_& m, char const* name, char const* fields, char const* doc = nullptr,
                              nb::object defaults = nb::none()) {
   nb::object t = nb::module_::import_("collections").attr("namedtuple")(name, fields, nb::arg("defaults") = defaults,
                                                                          nb::arg("module") = m.attr("__name__"));
   if (doc)
      t.attr("__doc__") = doc;
   m.attr(name) = t;
   return t;
}

// fn(std::integral_constant<int, I>()) for the I in [Lo, Hi] equal to n; a domain_error for any other n.
template <int Lo, int Hi, class Fn>
auto with_int(char const* what, int n, Fn&& fn) {
   return [&]<int... Is>(std::integer_sequence<int, Is...>) {
      using R = decltype(fn(std::integral_constant<int, Lo>()));
      if (n < Lo || n > Hi)
         throw std::domain_error(std::string(what) + " must be in [" + std::to_string(Lo) + ", " + std::to_string(Hi)
                                 + "], got " + std::to_string(n));
      R r{};
      (void)((n == Lo + Is ? (r = fn(std::integral_constant<int, Lo + Is>()), true) : false) || ...);
      return r;
   }(std::make_integer_sequence<int, Hi - Lo + 1>());
}

inline std::string py_float_repr(double x) { return nb::cast<std::string>(nb::repr(nb::float_(x))); }

}  // namespace bmpy

#endif
