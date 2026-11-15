// Copyright 2026 Nick Thompson
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// Shared by the files of the WebAssembly bindings: include it first. Doubles stay doubles: long double is software
// quad here. The error convention has two halves:
//
// 1. A function of numbers (the special functions, their _array twins, and the members of a distribution) returns NaN
//    or infinity for a domain error, as JavaScript's Math does, and never throws. These are bound with function().
// 2. Everything else is bound with entry(), and throws a RangeError for an argument that is wrong whatever the data: a
//    type, size, shape, count, option or name, and whatever the library throws in C++. NaN in its result is a number
//    the method computed (no convergence, a point outside the data), not a rejected argument. What a callback throws
//    passes through unchanged.

#ifndef BOOST_MATH_WASM_BINDINGS_HPP
#define BOOST_MATH_WASM_BINDINGS_HPP

#define BOOST_MATH_DOMAIN_ERROR_POLICY errno_on_error
#define BOOST_MATH_POLE_ERROR_POLICY ignore_error
#define BOOST_MATH_OVERFLOW_ERROR_POLICY ignore_error
#define BOOST_MATH_EVALUATION_ERROR_POLICY ignore_error
#define BOOST_MATH_ROUNDING_ERROR_POLICY ignore_error
#define BOOST_MATH_ASSERT_UNDEFINED_POLICY false
#define BOOST_MATH_PROMOTE_DOUBLE_POLICY false

#include <boost/math/tools/config.hpp>

#include <emscripten/bind.h>
#include <emscripten/em_asm.h>

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

namespace bm = boost::math;
using emscripten::class_;
using emscripten::function;
using emscripten::typed_memory_view;
using emscripten::val;
using cplx = std::complex<double>;

[[noreturn]] inline void range_error(const std::string& what)
{
   throw std::range_error(what);
}

// Marks an export as one that throws: "f", "C.method", or "new C" for the constructor of a class. post.js wraps
// each, so that a C++ exception reaches JavaScript as a RangeError and the stack pointer is put back.
inline void throws(const char* name)
{
   EM_ASM({ (Module.throwing ||= []).push(UTF8ToString($0)); }, name);
}

// function(), for each overload f of an export that throws. Their arguments are numbers, booleans and vals: embind
// frees a std::string argument only when the call returns.
template <class... F>
void entry(const char* name, F... f)
{
   (function(name, f), ...);
   throws(name);
}

inline std::vector<double> numbers(const val& x)
{
   return emscripten::convertJSArrayToNumberVector<double>(x);
}

inline std::string text(const val& x, const char* what)
{
   if (!x.isString())
      range_error(std::string(what) + " must be a string.");
   return x.as<std::string>();
}

// Arrays go out as new Float64Arrays; complex data are interleaved (re, im).
template <class T>
val to_js(const std::vector<T>& v)
{
   return val(typed_memory_view(v.size() * sizeof(T) / sizeof(double), reinterpret_cast<const double*>(v.data()))).call<val>("slice");
}

// An array of rows, each a Float64Array.
inline val rows(const std::vector<std::vector<double>>& m)
{
   val a = val::array();
   for (const auto& r : m)
      a.call<void>("push", to_js(r));
   return a;
}

// Callbacks are synchronous: the result of an async function is refused, not coerced.
inline val sync(val result)
{
   static const val promise = val::global("Promise");
   if (result.instanceof(promise))
      val::global("TypeError").new_(std::string("The callback returned a promise: callbacks must be synchronous.")).throw_();
   return result;
}

// v, if it is an integer in [least, most]: NaN and infinity are not. A count is at most 2^31, which a std::size_t holds.
inline double whole(double v, const char* name, double least, double most)
{
   if (!(v == std::trunc(v) && v >= least && v <= most))
      range_error(std::string(name) + " must be an integer in [" + std::to_string(std::int64_t(least)) + ", " + std::to_string(std::int64_t(most)) + "].");
   return v;
}
inline std::size_t count(double v, const char* name, double least = 0, double most = 0x1p31)
{
   return static_cast<std::size_t>(whole(v, name, least, most));
}

// The numbers a result may hold, 2^27 (a gigabyte): products of counts are checked against it before anything is
// allocated.
inline void check_size(double count, const char* what)
{
   if (!(count <= 0x1p27))
      range_error(std::string(what) + " would hold more than 2^27 numbers.");
}

// An options object: every option has a default, and one that was not asked for by the time done() runs is an error.
struct options
{
   val o;
   std::vector<std::string> known;

   explicit options(const val& v) : o(v.isUndefined() || v.isNull() ? val::object() : v)
   {
      if (o.typeOf().as<std::string>() != "object")
         range_error("The options must be an object.");
   }
   val get(const char* name)
   {
      known.push_back(name);
      val v = o[name];
      return v.isNull() ? val::undefined() : v;
   }
   bool has(const char* name) { return !get(name).isUndefined(); }
   double real(const char* name, double fallback)
   {
      val v = get(name);
      if (!v.isUndefined() && !v.isNumber())
         range_error(std::string(name) + " must be a number.");
      return v.isUndefined() ? fallback : v.as<double>();
   }
   bool flag(const char* name, bool fallback = false)
   {
      val v = get(name);
      return v.isUndefined() ? fallback : v.as<bool>();
   }
   // An integer in [least, most].
   double whole(const char* name, double fallback, double least, double most) { return ::whole(real(name, fallback), name, least, most); }
   std::size_t count(const char* name, std::size_t fallback, double least = 0, double most = 0x1p31)
   {
      return ::count(real(name, double(fallback)), name, least, most);
   }
   std::uint64_t seed() { return std::uint64_t(whole("seed", 0, 0, 0x1p53)); }
   void done() const
   {
      const val keys = val::global("Object").call<val>("keys", o);
      for (std::size_t i = 0; i < keys["length"].as<std::size_t>(); ++i)
         if (std::find(known.begin(), known.end(), keys[i].as<std::string>()) == known.end())
            range_error("Unknown option " + keys[i].as<std::string>() + ".");
   }
};

#endif
