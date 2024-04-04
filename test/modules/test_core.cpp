// Copyright John Maddock 2022.

// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <iomanip>
#include <limits>
#include <tuple>
#include <complex>
#include <array>
#include <cstdint>

import boost.math.core;

int main()
{
   std::cout << std::setprecision(std::numeric_limits<float>::max_digits10) << boost::math::sign(-1.245f) << std::endl;
   std::cout << std::setprecision(std::numeric_limits<double>::max_digits10) << boost::math::sign(-1.245) << std::endl;
   std::cout << std::setprecision(std::numeric_limits<long double>::max_digits10) << boost::math::sign(-1.245L) << std::endl;
   std::cout << std::setprecision(std::numeric_limits<float>::max_digits10) << boost::math::signbit(-1.245f) << std::endl;
   std::cout << std::setprecision(std::numeric_limits<double>::max_digits10) << boost::math::signbit(-1.245) << std::endl;
   std::cout << std::setprecision(std::numeric_limits<long double>::max_digits10) << boost::math::signbit(-1.245L) << std::endl;
   std::cout << std::setprecision(std::numeric_limits<float>::max_digits10) << boost::math::changesign(-1.245f) << std::endl;
   std::cout << std::setprecision(std::numeric_limits<double>::max_digits10) << boost::math::changesign(-1.245) << std::endl;
   std::cout << std::setprecision(std::numeric_limits<long double>::max_digits10) << boost::math::changesign(-1.245L) << std::endl;
   std::cout << std::setprecision(std::numeric_limits<float>::max_digits10) << boost::math::copysign(1.0f, -1.245f) << std::endl;
   std::cout << std::setprecision(std::numeric_limits<double>::max_digits10) << boost::math::copysign(1.0, -1.245) << std::endl;
   std::cout << std::setprecision(std::numeric_limits<long double>::max_digits10) << boost::math::copysign(1.0L, -1.245L) << std::endl;
   
   std::cout << boost::math::fpclassify(1.0f) << std::endl;
   std::cout << boost::math::isnormal(1.0f) << std::endl;
   std::cout << boost::math::isnan(1.0f) << std::endl;
   std::cout << boost::math::isinf(1.0f) << std::endl;
   std::cout << boost::math::isfinite(1.0f) << std::endl;
   std::cout << boost::math::fpclassify(1.0) << std::endl;
   std::cout << boost::math::isnormal(1.0) << std::endl;
   std::cout << boost::math::isnan(1.0) << std::endl;
   std::cout << boost::math::isinf(1.0) << std::endl;
   std::cout << boost::math::isfinite(1.0) << std::endl;
   std::cout << boost::math::fpclassify(1.0L) << std::endl;
   std::cout << boost::math::isnormal(1.0L) << std::endl;
   std::cout << boost::math::isnan(1.0L) << std::endl;
   std::cout << boost::math::isinf(1.0L) << std::endl;
   std::cout << boost::math::isfinite(1.0L) << std::endl;

   std::cout << std::setprecision(std::numeric_limits<double>::max_digits10) << boost::math::trunc(3.14) << std::endl;
   std::cout << boost::math::itrunc(3.14) << std::endl;
   std::cout << boost::math::ltrunc(3.14) << std::endl;
   std::cout << boost::math::lltrunc(3.14) << std::endl;
   std::cout << boost::math::iconvert(3.14, boost::math::policies::make_policy()) << std::endl;
   std::cout << boost::math::lconvert(3.14, boost::math::policies::make_policy()) << std::endl;
   std::cout << boost::math::llconvert(3.14, boost::math::policies::make_policy()) << std::endl;

   std::cout << std::setprecision(std::numeric_limits<double>::max_digits10) << boost::math::float_next(3.1456789876) << std::endl;
   std::cout << std::setprecision(std::numeric_limits<double>::max_digits10) << boost::math::float_prior(3.1456789876) << std::endl;
   std::cout << boost::math::float_distance(3.1456789876, boost::math::float_prior(3.1456789876)) << std::endl;
   std::cout << boost::math::float_distance(3.1456789876, boost::math::float_advance(3.1456789876, 2)) << std::endl;

   auto f = [](double x) { return 3 * x - 2; };
   std::uintmax_t max_iter = 300;
   std::cout << boost::math::tools::toms748_solve(f, 0.0, 4.0, f(0.0), f(4.0), boost::math::tools::eps_tolerance<double>(), max_iter, boost::math::policies::make_policy()).first << std::endl;
   max_iter = 300;
   std::cout << boost::math::tools::toms748_solve(f, 0.0, 4.0, boost::math::tools::eps_tolerance<double>(), max_iter, boost::math::policies::make_policy()).first << std::endl;
   max_iter = 300;
   std::cout << boost::math::tools::toms748_solve(f, 0.0, 4.0, f(0.0), f(4.0), boost::math::tools::eps_tolerance<double>(), max_iter).first << std::endl;
   max_iter = 300;
   std::cout << boost::math::tools::toms748_solve(f, 0.0, 4.0, boost::math::tools::eps_tolerance<double>(), max_iter).first << std::endl;
   max_iter = 300;
   std::cout << boost::math::tools::bracket_and_solve_root(f, 1.0, 2.0, true, boost::math::tools::eps_tolerance<double>(), max_iter, boost::math::policies::make_policy()).first << std::endl;
   max_iter = 300;
   std::cout << boost::math::tools::bracket_and_solve_root(f, 1.0, 2.0, true, boost::math::tools::eps_tolerance<double>(), max_iter).first << std::endl;

   max_iter = 300;
   std::cout << boost::math::tools::bisect(f, 0.0, 3.0, boost::math::tools::eps_tolerance<double>(), max_iter, boost::math::policies::make_policy()).first << std::endl;
   max_iter = 300;
   std::cout << boost::math::tools::bisect(f, 0.0, 3.0, boost::math::tools::eps_tolerance<double>(), max_iter).first << std::endl;
   max_iter = 300;
   std::cout << boost::math::tools::bisect(f, 0.0, 3.0, boost::math::tools::eps_tolerance<double>()).first << std::endl;
   max_iter = 300;
   auto f2 = [](double x) { return std::make_tuple(3 * x - 2, 3.0); };
   std::cout << boost::math::tools::newton_raphson_iterate(f2, 2.0, 0.0, 4.0, 50, max_iter) << std::endl;
   std::cout << boost::math::tools::newton_raphson_iterate(f2, 2.0, 0.0, 4.0, 50) << std::endl;
   auto f3 = [](double x) { return std::make_tuple(3 * x * x - 5 * x - 3, 6 * x - 5, 6.0); };
   max_iter = 200;
   std::cout << boost::math::tools::halley_iterate(f3, 2.0, 1.0, 5.0, 50, max_iter) << std::endl;
   max_iter = 200;
   std::cout << boost::math::tools::halley_iterate(f3, 2.0, 1.0, 5.0, 50) << std::endl;
   max_iter = 200;
   std::cout << boost::math::tools::schroder_iterate(f3, 2.0, 1.0, 5.0, 50, max_iter) << std::endl;
   max_iter = 200;
   std::cout << boost::math::tools::schroder_iterate(f3, 2.0, 1.0, 5.0, 50) << std::endl;
   max_iter = 200;
   std::cout << boost::math::tools::schroeder_iterate(f3, 2.0, 1.0, 5.0, 50, max_iter) << std::endl;
   max_iter = 200;
   std::cout << boost::math::tools::schroeder_iterate(f3, 2.0, 1.0, 5.0, 50) << std::endl;
   auto f4 = [](std::complex<double> z) { return std::make_pair(z * z + 1.0, 2.0 * z); };
   std::cout << boost::math::tools::complex_newton(f4, std::complex<double>(0)) << std::endl;
   std::cout << boost::math::tools::quadratic_roots(1.0, 0.0, -1.0).first << std::endl;

   double coef[] = { 1, 4, -3, 5, -2 };
   double coef2[] = { 5, -7, 1, -2, 0 };
   std::cout << boost::math::tools::evaluate_polynomial(coef, 2.0) << std::endl;
   std::cout << boost::math::tools::evaluate_even_polynomial(coef, 2.0) << std::endl;
   std::cout << boost::math::tools::evaluate_odd_polynomial(coef, 2.0) << std::endl;
   std::cout << boost::math::tools::evaluate_rational(coef, coef2, 2.0) << std::endl;
}
