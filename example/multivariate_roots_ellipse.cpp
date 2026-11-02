// Copyright Nick Thompson 2026.
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0. (See accompanying file
// LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// The point of the ellipse x^2/a^2 + y^2/b^2 = 1 closest to q, by Newton's method and Cuyt's multivariate Halley
// method, with the Jacobian and the second directional derivative from forward-mode autodiff.

#include <boost/math/differentiation/autodiff.hpp>
#include <boost/math/tools/multivariate_roots.hpp>
#include <array>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <tuple>

//[multivariate_roots_ellipse_system
// At the closest point p, p is on the ellipse and p - q is along the normal (x/a^2, y/b^2).
struct ellipse_closest_point
{
   double a;
   double b;
   std::array<double, 2> q;

   template <class T>
   std::array<T, 2> operator()(std::array<T, 2> const& p) const
   {
      return {(p[0] - q[0]) * p[1] / (b * b) - (p[1] - q[1]) * p[0] / (a * a),
              p[0] * p[0] / (a * a) + p[1] * p[1] / (b * b) - 1};
   }
};
//]

//[multivariate_roots_autodiff
using boost::math::differentiation::make_fvar;

// F(x + t v) with its first Order derivatives in t at t = 0: F, J v and F''[v, v]. One forward pass in the single
// variable t, however many variables F has.
template <int Order, class System, std::size_t N>
auto along(System const& F, std::array<double, N> const& x, std::array<double, N> const& v)
{
   auto const t = make_fvar<double, Order>(0);
   std::array<boost::math::differentiation::autodiff_fvar<double, Order>, N> p;
   for (std::size_t k = 0; k < N; ++k)
      p[k] = x[k] + v[k] * t;
   return F(p);
}

// (F, J, d2) for multivariate_halley: J a column at a time, d2(v) = F''[v, v] when the solver asks for it.
template <class System, std::size_t N>
auto with_derivatives(System const& F, std::array<double, N> const& x)
{
   std::array<double, N> value{};
   std::array<std::array<double, N>, N> J{};
   for (std::size_t j = 0; j < N; ++j)
   {
      std::array<double, N> e{};
      e[j] = 1;
      auto const Fj = along<1>(F, x, e);
      for (std::size_t i = 0; i < N; ++i)
      {
         value[i] = Fj[i].derivative(0);
         J[i][j] = Fj[i].derivative(1);
      }
   }
   auto d2 = [F, x](std::array<double, N> const& v) {
      auto const Fv = along<2>(F, x, v);
      std::array<double, N> w;
      for (std::size_t i = 0; i < N; ++i)
         w[i] = Fv[i].derivative(2);
      return w;
   };
   return std::make_tuple(value, J, d2);
}
//]

int main()
{
   using boost::math::tools::multivariate_halley;
   using boost::math::tools::multivariate_newton;
   double const a = 3;
   double const b = 1;
   std::cout << std::setprecision(17);
   for (std::array<double, 2> q : {std::array<double, 2>{4, 2}, {1, 0.25}, {-0.5, 3}, {2.9, -0.1}})
   {
      //[multivariate_roots_ellipse_solve
      ellipse_closest_point F{a, b, q};
      auto f = [&F](std::array<double, 2> const& p) { return with_derivatives(F, p); };
      // Start from the point of the ellipse at the angle of q, stretched to the ellipse's shape:
      double const theta = std::atan2(a * q[1], b * q[0]);
      std::array<double, 2> const p0{a * std::cos(theta), b * std::sin(theta)};
      std::uintmax_t halley_iterations = 50;
      auto p = multivariate_halley(f, p0, std::numeric_limits<double>::digits, halley_iterations);
      //]
      std::uintmax_t newton_iterations = 50;
      auto pn = multivariate_newton(f, p0, std::numeric_limits<double>::digits, newton_iterations);
      std::cout << "q = (" << q[0] << ", " << q[1] << "): closest point (" << p[0] << ", " << p[1]
                << "), distance " << std::hypot(p[0] - q[0], p[1] - q[1]) << "\n"
                << "   Halley iterations: " << halley_iterations << ", Newton iterations: " << newton_iterations
                << ", |Halley - Newton| = " << std::hypot(p[0] - pn[0], p[1] - pn[1]) << "\n";
   }
}
