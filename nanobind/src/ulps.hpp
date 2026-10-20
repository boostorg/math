// Copyright 2026 Nick Thompson
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// The functions boost.math.ulps can plot: each at double (the implementation) and at
// cpp_bin_float_50 (the reference), type-erased so that ulps_plot and ulps_heatmap are
// instantiated once per arity.

#ifndef BOOST_MATH_NANOBIND_ULPS_HPP
#define BOOST_MATH_NANOBIND_ULPS_HPP

#include <boost/multiprecision/cpp_bin_float.hpp>

#include <functional>
#include <map>
#include <string>

namespace bmpy::ulps {

using Precise = boost::multiprecision::cpp_bin_float_50;

struct fn1 {
   std::function<double(double)> coarse;
   std::function<Precise(Precise const&)> precise;
};

struct fn2 {
   std::function<double(double, double)> coarse;
   std::function<Precise(Precise const&, Precise const&)> precise;
};

struct registry {
   std::map<std::string, fn1> unary;
   std::map<std::string, fn2> binary;
};

template <class F>
fn1 make1(F f) {
   return {[f](double x) { return static_cast<double>(f(x)); },
           [f](Precise const& x) { return static_cast<Precise>(f(x)); }};
}

template <class F>
fn2 make2(F f) {
   return {[f](double x, double y) { return static_cast<double>(f(x, y)); },
           [f](Precise const& x, Precise const& y) { return static_cast<Precise>(f(x, y)); }};
}

// One entry per function; fn is called through bm::.
#define BMPY_ULPS1(fn) r.unary.emplace(#fn, bmpy::ulps::make1([](auto x) { return bm::fn(x); }))
#define BMPY_ULPS2(fn) r.binary.emplace(#fn, bmpy::ulps::make2([](auto x, auto y) { return bm::fn(x, y); }))

void register_gamma(registry& r);
void register_special(registry& r);
void register_elliptic(registry& r);

}  // namespace bmpy::ulps

#endif
