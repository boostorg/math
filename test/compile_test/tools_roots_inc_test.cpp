//  Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Basic sanity check that header <boost/math/tools/roots.hpp>
// #includes all the files that it needs to.
//
#include <boost/math/tools/roots.hpp>

typedef double (*F)(double);
typedef std::pair<double, double> (*F2)(double);
typedef std::tr1::tuple<double, double, double> (*F3)(double);
#define T double
typedef boost::math::tools::eps_tolerance<double> Tol;

template std::pair<T, T> boost::math::tools::bisect<F, T, Tol>(F f, T min, T max, Tol tol, boost::uintmax_t& max_iter);
template std::pair<T, T> boost::math::tools::bisect<F, T, Tol>(F f, T min, T max, Tol tol);
template T boost::math::tools::newton_raphson_iterate<F2, T>(F2 f, T guess, T min, T max, int digits, boost::uintmax_t& max_iter);
template T boost::math::tools::halley_iterate<F3, T>(F3 f, T guess, T min, T max, int digits, boost::uintmax_t& max_iter);
template T boost::math::tools::schroeder_iterate<F3, T>(F3 f, T guess, T min, T max, int digits, boost::uintmax_t& max_iter);


