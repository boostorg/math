// Copyright 2026 Nick Thompson
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// Test problems for the optimizers of the WebAssembly bindings, by name: they run without a callback, and the
// gradient-based optimizers, which differentiate their objective, take nothing else. Every function is generic in the
// scalar type of its argument, a container X of reals or of autodiff variables.

#ifndef BOOST_MATH_WASM_TEST_PROBLEMS_HPP
#define BOOST_MATH_WASM_TEST_PROBLEMS_HPP

#include <boost/math/constants/constants.hpp>

#include <cmath>
#include <cstddef>
#include <iterator>
#include <string>
#include <vector>

namespace test_problems
{
// One objective, of two variables. Each is 0 at its minimum, but styblinski_tang: -78.33233140754282.
const char* const single[] = {"rosenbrock", "rastrigin", "ackley", "himmelblau", "styblinski_tang"};
// Several objectives: two, and dtlz2 three. The zdt and dtlz2 are posed on [0, 1]^n, fonseca_fleming on [-4, 4]^n
// and kursawe on [-5, 5]^3.
const char* const multiple[] = {"zdt1", "zdt2", "zdt3", "zdt4", "zdt6", "dtlz2", "fonseca_fleming", "kursawe"};

// The index of a name, or -1.
template <class Names>
int find(const Names& names, const std::string& name)
{
   for (int i = 0; i < int(std::size(names)); ++i)
      if (name == names[i])
         return i;
   return -1;
}

template <class X>
typename X::value_type value(int id, const X& v)
{
   using std::cos;
   using std::exp;
   using std::sqrt;
   const double two_pi = 2 * boost::math::constants::pi<double>(), e = boost::math::constants::e<double>();
   const auto& x = v[0];
   const auto& y = v[1];
   switch (id)
   {
   case 0:
      return (1 - x) * (1 - x) + 100 * (y - x * x) * (y - x * x);
   case 1:
      return 20 + x * x - 10 * cos(two_pi * x) + y * y - 10 * cos(two_pi * y);
   case 2:
      return -20 * exp(-0.2 * sqrt(0.5 * (x * x + y * y))) - exp(0.5 * (cos(two_pi * x) + cos(two_pi * y))) + e + 20;
   case 3:
      return (x * x + y - 11) * (x * x + y - 11) + (x + y * y - 7) * (x + y * y - 7);
   default:
      return (x * x * x * x - 16 * x * x + 5 * x + y * y * y * y - 16 * y * y + 5 * y) / 2;
   }
}

// Zitzler, Deb and Thiele (2000), problem 1, 2, 3, 4 or 6: the Pareto front is at x[i] = 0 for i > 0.
template <class X>
std::vector<typename X::value_type> zdt(const X& x, int problem)
{
   using Real = typename X::value_type;
   using std::cos;
   using std::exp;
   using std::pow;
   using std::sin;
   using std::sqrt;
   const double pi = boost::math::constants::pi<double>();
   Real s = 0;
   for (std::size_t i = 1; i < x.size(); ++i)
      s += problem == 4 ? Real(x[i] * x[i] - 10 * cos(4 * pi * x[i])) : Real(x[i]);
   const double n = double(x.size() - 1);
   Real f1 = x[0];
   Real g = problem == 4 ? Real(1 + 10 * n + s) : Real(1 + 9 * s / n);
   if (problem == 6)
   {
      f1 = 1 - exp(-4 * x[0]) * pow(sin(6 * pi * x[0]), 6);
      g = 1 + 9 * pow(s / n, 0.25);
   }
   Real h = 1 - sqrt(f1 / g);
   if (problem == 2 || problem == 6)
      h = 1 - (f1 / g) * (f1 / g);
   else if (problem == 3)
      h -= (f1 / g) * sin(10 * pi * f1);
   return {f1, g * h};
}

// Deb, Thiele, Laumanns and Zitzler (2005): the Pareto front is the unit sphere in the positive octant.
template <class X>
std::vector<typename X::value_type> dtlz2(const X& x)
{
   using Real = typename X::value_type;
   using std::sin;
   const double half_pi = boost::math::constants::half_pi<double>();
   Real g = 1;
   for (std::size_t i = 2; i < x.size(); ++i)
      g += (x[i] - 0.5) * (x[i] - 0.5);
   // cos(t) as sin(pi/2 - t), so that no objective is negative at x = 1.
   const Real c0 = sin(half_pi * (1 - x[0])), c1 = sin(half_pi * (1 - x[1]));
   return {g * c0 * c1, g * c0 * sin(half_pi * x[1]), g * sin(half_pi * x[0])};
}

// Fonseca and Fleming (1995).
template <class X>
std::vector<typename X::value_type> fonseca_fleming(const X& x)
{
   using Real = typename X::value_type;
   using std::exp;
   const double s = 1 / std::sqrt(double(x.size()));
   Real d1 = 0, d2 = 0;
   for (const auto& xi : x)
   {
      d1 += (xi - s) * (xi - s);
      d2 += (xi + s) * (xi + s);
   }
   return {1 - exp(-d1), 1 - exp(-d2)};
}

// Kursawe (1991): the Pareto front is in three pieces, and is not convex.
template <class X>
std::vector<typename X::value_type> kursawe(const X& x)
{
   using Real = typename X::value_type;
   using std::abs;
   using std::exp;
   using std::pow;
   using std::sin;
   using std::sqrt;
   Real f1 = 0, f2 = 0;
   for (std::size_t i = 0; i < x.size(); ++i)
   {
      if (i + 1 < x.size())
         f1 -= 10 * exp(-sqrt(x[i] * x[i] + x[i + 1] * x[i + 1]) / 5);
      f2 += pow(abs(x[i]), 0.8) + 5 * sin(x[i] * x[i] * x[i]);
   }
   return {f1, f2};
}

// The objectives of multiple[id] at x, which needs at least variables(id) variables.
inline std::size_t variables(int id) { return id == 5 ? 3 : 2; }

template <class X>
std::vector<typename X::value_type> values(int id, const X& x)
{
   static const int zdts[] = {1, 2, 3, 4, 6};
   return id < 5 ? zdt(x, zdts[id]) : id == 5 ? dtlz2(x) : id == 6 ? fonseca_fleming(x) : kursawe(x);
}
} // namespace test_problems

#endif
