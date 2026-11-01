// Copyright 2026 Nick Thompson
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// tools::fft<std::complex<double>>, tools::real_fft<std::complex<double>> and the DCT-II behind
// chebyshev_transform (detail::dct_ii<double>), as classes. Arrays must be 1-D, C-contiguous NumPy arrays of
// exactly the right dtype and length (writable for the in-place methods); anything else raises ValueError,
// except that fft.transform and inverse_transform convert real and integer arrays to complex128.
// Transforms release the GIL. Each object owns its workspace, so a mutex keeps two threads from using one
// object at the same time.

#include "bindings.hpp"

#include <boost/math/special_functions/chebyshev_transform.hpp>
#include <boost/math/tools/fft.hpp>

#include <mutex>
#include <span>
#include <stdexcept>
#include <string>

namespace bmpy {

namespace {

using C = std::complex<double>;

template <class T>
using array1d = nb::ndarray<nb::numpy, T, nb::ndim<1>, nb::c_contig, nb::device::cpu>;

template <class T>
char const* dtype_name()
{
   return std::is_same_v<T, C> ? "complex128" : "float64";
}

// Checks that a is a 1-D, C-contiguous NumPy array of dtype T and length n (and writable if asked),
// raising ValueError otherwise, and returns a pointer to its data.
template <class T>
T* checked_data(nb::handle a, std::size_t n, bool writable, char const* what)
{
   nb::module_ np = nb::module_::import_("numpy");
   std::string prefix = std::string(what) + ": ";
   if (!nb::isinstance(a, np.attr("ndarray")))
   {
      throw nb::value_error((prefix + "expected a NumPy array").c_str());
   }
   if (!a.attr("dtype").equal(np.attr("dtype")(dtype_name<T>())))
   {
      throw nb::value_error((prefix + "expected dtype " + dtype_name<T>() + ", got " +
                             nb::cast<std::string>(a.attr("dtype").attr("name"))).c_str());
   }
   if (nb::cast<int>(a.attr("ndim")) != 1)
   {
      throw nb::value_error((prefix + "expected a 1-D array").c_str());
   }
   nb::handle flags = a.attr("flags");
   if (!nb::cast<bool>(flags.attr("c_contiguous")))
   {
      throw nb::value_error((prefix + "expected a contiguous array").c_str());
   }
   if (writable && !nb::cast<bool>(flags.attr("writeable")))
   {
      throw nb::value_error((prefix + "expected a writable array").c_str());
   }
   std::size_t len = nb::cast<std::size_t>(a.attr("shape")[0]);
   if (len != n)
   {
      throw nb::value_error((prefix + "expected length " + std::to_string(n) + ", got " + std::to_string(len)).c_str());
   }
   if (writable)
   {
      return nb::cast<array1d<T>>(a).data();
   }
   return const_cast<T*>(nb::cast<nb::ndarray<nb::numpy, T const, nb::ndim<1>, nb::c_contig, nb::device::cpu>>(a).data());
}

template <class T>
array1d<T> new_array(std::size_t n, T** data)
{
   T* buf = new T[n == 0 ? 1 : n];
   nb::capsule owner(buf, [](void* p) noexcept { delete[] static_cast<T*>(p); });
   *data = buf;
   std::size_t shape[1] = {n};
   return array1d<T>(buf, 1, shape, owner);
}

struct py_fft
{
   explicit py_fft(std::size_t n) : f(n) {}
   boost::math::tools::fft<C> f;
   std::mutex mutex;

   array1d<C> out_of_place(nb::handle x, bool inverse)
   {
      std::size_t n = f.size();
      // Real and integer arrays are transformed as complex128, as by np.fft.
      nb::object converted;
      nb::module_ np = nb::module_::import_("numpy");
      std::string kind = nb::isinstance(x, np.attr("ndarray")) ? nb::cast<std::string>(x.attr("dtype").attr("kind")) : "";
      if (kind.size() == 1 && std::string("biuf").find(kind) != std::string::npos)
      {
         converted = np.attr("asarray")(x, np.attr("complex128"));
         x = converted;
      }
      C const* in = checked_data<C>(x, n, false, inverse ? "inverse_transform" : "transform");
      C* out;
      array1d<C> result = new_array<C>(n, &out);
      {
         nb::gil_scoped_release release;
         std::lock_guard<std::mutex> lock(mutex);
         std::span<C const> s_in(in, n);
         std::span<C> s_out(out, n);
         if (inverse)
            f.inverse_transform(s_in, s_out);
         else
            f.transform(s_in, s_out);
      }
      return result;
   }

   void in_place(nb::handle x, bool inverse)
   {
      std::size_t n = f.size();
      C* v = checked_data<C>(x, n, true, inverse ? "inverse_transform_in_place" : "transform_in_place");
      nb::gil_scoped_release release;
      std::lock_guard<std::mutex> lock(mutex);
      std::span<C> s(v, n);
      if (inverse)
         f.inverse_transform(s);
      else
         f.transform(s);
   }
};

struct py_real_fft
{
   explicit py_real_fft(std::size_t n) : f(n), n(n) {}
   boost::math::tools::real_fft<C> f;
   std::size_t n;
   std::mutex mutex;

   array1d<C> transform(nb::handle x)
   {
      double const* in = checked_data<double>(x, n, false, "transform");
      C* out;
      array1d<C> result = new_array<C>(n / 2 + 1, &out);
      {
         nb::gil_scoped_release release;
         std::lock_guard<std::mutex> lock(mutex);
         std::span<double const> s_in(in, n);
         std::span<C> s_out(out, n / 2 + 1);
         f.transform(s_in, s_out);
      }
      return result;
   }
};

struct py_dct_ii
{
   explicit py_dct_ii(std::size_t n) : d(n), n(n) {}
   boost::math::detail::dct_ii<double> d;
   std::size_t n;
   std::mutex mutex;

   array1d<double> transform(nb::handle x)
   {
      double const* in = checked_data<double>(x, n, false, "transform");
      double* out;
      array1d<double> result = new_array<double>(n, &out);
      {
         nb::gil_scoped_release release;
         std::lock_guard<std::mutex> lock(mutex);
         d.execute(in, out);
      }
      return result;
   }

   void transform_in_place(nb::handle x)
   {
      double* v = checked_data<double>(x, n, true, "transform_in_place");
      nb::gil_scoped_release release;
      std::lock_guard<std::mutex> lock(mutex);
      d.execute(v);
   }
};

void check_length(std::size_t n)
{
   if (n == 0)
   {
      throw nb::value_error("The length of the transform must be positive.");
   }
}

}  // namespace

void bind_fft(nb::module_& m)
{
   nb::class_<py_fft>(m, "fft",
                      "Discrete Fourier transform of length n, X[k] = sum_j x[j] exp(-2 pi i j k/n), on complex128 arrays "
                      "(real input is converted for the out-of-place transforms); the inverse divides by n.\n\n"
                      "As numpy.fft.fft, but planned once for the length n and reusable.\n\n"
                      "https://www.boost.org/doc/libs/release/libs/math/doc/html/math_toolkit/fft.html")
      .def("__init__", [](py_fft* self, std::size_t n) { check_length(n); new (self) py_fft(n); }, nb::arg("n"))
      .def_prop_ro("size", [](py_fft& self) { return self.f.size(); })
      .def("__len__", [](py_fft& self) { return self.f.size(); })
      .def("__repr__", [](py_fft& self) { return "fft(n=" + std::to_string(self.f.size()) + ")"; })
      .def("transform", [](py_fft& self, nb::handle x) { return self.out_of_place(x, false); }, nb::arg("x"),
           "Returns the transform of x as a new array.")
      .def("inverse_transform", [](py_fft& self, nb::handle x) { return self.out_of_place(x, true); }, nb::arg("x"),
           "Returns the inverse transform of x as a new array.")
      .def("transform_in_place", [](py_fft& self, nb::handle x) { self.in_place(x, false); }, nb::arg("x"),
           "Replaces x (writable, contiguous) with its transform.")
      .def("inverse_transform_in_place", [](py_fft& self, nb::handle x) { self.in_place(x, true); }, nb::arg("x"),
           "Replaces x (writable, contiguous) with its inverse transform.");

   nb::class_<py_real_fft>(m, "real_fft",
                           "Transform of n float64 samples: returns X[0], ..., X[n/2] (complex128); the rest are conj(X[n-k]).\n\n"
                           "As numpy.fft.rfft.\n\n"
                           "https://www.boost.org/doc/libs/release/libs/math/doc/html/math_toolkit/fft.html")
      .def("__init__", [](py_real_fft* self, std::size_t n) { check_length(n); new (self) py_real_fft(n); }, nb::arg("n"))
      .def_prop_ro("size", [](py_real_fft& self) { return self.n; })
      .def("__len__", [](py_real_fft& self) { return self.n; })
      .def("__repr__", [](py_real_fft& self) { return "real_fft(n=" + std::to_string(self.n) + ")"; })
      .def("transform", &py_real_fft::transform, nb::arg("x"));

   nb::class_<py_dct_ii>(m, "dct_ii",
                         "DCT-II with FFTW's REDFT10 normalisation, y[k] = 2 sum_j x[j] cos(pi k (2j+1)/(2n)), on float64 arrays "
                         "(the transform behind chebyshev_transform).\n\n"
                         "As scipy.fft.dct(x, type=2), which uses the same scaling.\n\n"
                         "https://www.boost.org/doc/libs/release/libs/math/doc/html/math_toolkit/fft.html")
      .def("__init__", [](py_dct_ii* self, std::size_t n) { check_length(n); new (self) py_dct_ii(n); }, nb::arg("n"))
      .def_prop_ro("size", [](py_dct_ii& self) { return self.n; })
      .def("__len__", [](py_dct_ii& self) { return self.n; })
      .def("__repr__", [](py_dct_ii& self) { return "dct_ii(n=" + std::to_string(self.n) + ")"; })
      .def("transform", &py_dct_ii::transform, nb::arg("x"))
      .def("transform_in_place", &py_dct_ii::transform_in_place, nb::arg("x"));
}

}  // namespace bmpy
