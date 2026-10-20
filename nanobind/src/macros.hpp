// Copyright 2026 Nick Thompson
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// Shorthands for the common all-real signatures. The function name is the
// Python name and the Boost.Math function called through bm::.

#ifndef BOOST_MATH_NANOBIND_MACROS_HPP
#define BOOST_MATH_NANOBIND_MACROS_HPP

#define BMPY_R1(fn, a) \
   bmpy::def_vec<bmpy::Real>(m, #fn, []<class T>(T x0) { return bm::fn(x0); }, {a}, "boost::math::" #fn)
#define BMPY_R2(fn, a, b) \
   bmpy::def_vec<bmpy::Real, bmpy::Real>(m, #fn, []<class T>(T x0, T x1) { return bm::fn(x0, x1); }, {a, b}, "boost::math::" #fn)
#define BMPY_R3(fn, a, b, c)                                                                                          \
   bmpy::def_vec<bmpy::Real, bmpy::Real, bmpy::Real>(m, #fn, []<class T>(T x0, T x1, T x2) { return bm::fn(x0, x1, x2); }, \
                                                     {a, b, c}, "boost::math::" #fn)
#define BMPY_R4(fn, a, b, c, d)                                                                               \
   bmpy::def_vec<bmpy::Real, bmpy::Real, bmpy::Real, bmpy::Real>(                                             \
       m, #fn, []<class T>(T x0, T x1, T x2, T x3) { return bm::fn(x0, x1, x2, x3); }, {a, b, c, d}, \
       "boost::math::" #fn)
// One integer argument of type I followed by one real argument.
#define BMPY_I1R1(fn, I, a, b) \
   bmpy::def_vec<bmpy::Int<I>, bmpy::Real>(m, #fn, []<class T>(I n, T x) { return bm::fn(n, x); }, {a, b}, "boost::math::" #fn)
// One real argument followed by one integer argument of type I.
#define BMPY_R1I1(fn, I, a, b) \
   bmpy::def_vec<bmpy::Real, bmpy::Int<I>>(m, #fn, []<class T>(T x, I n) { return bm::fn(x, n); }, {a, b}, "boost::math::" #fn)
// Functions of one integer whose result type is the template argument (factorial<T>(n), ...); double only.
#define BMPY_T_I1(fn, I, a) \
   bmpy::def_vec<bmpy::Int<I>>(m, #fn, []<class T>(I n) { return bm::fn<double>(n); }, {a}, "boost::math::" #fn "<double>")

namespace bm = boost::math;

#endif
