//  Copyright Nick Thompson 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_QUADRATURE_SPARSE_GRID_HPP
#define BOOST_MATH_QUADRATURE_SPARSE_GRID_HPP

#ifdef _MSC_VER
#pragma once
#endif

#include <boost/math/quadrature/leja.hpp>

#ifndef BOOST_MATH_BUILD_MODULE
#include <cstddef>
#include <functional>
#include <map>
#include <set>
#include <type_traits>
#include <utility>
#include <vector>
#endif

BOOST_MATH_NAMESPACE_BEGIN namespace quadrature {

namespace detail {
// The point f is called with: std::vector<Real> if f takes one, otherwise the argument of f's only call operator, such as
// std::array<Real, N> or Eigen::VectorXd.
template <class R, class A>
std::decay_t<A> leja_argument(std::function<R(A)>);
template <class F, class Real, class = void>
struct leja_point { using type = std::vector<Real>; };
template <class F, class Real>
struct leja_point<F, Real, std::enable_if_t<!std::is_invocable_v<F&, const std::vector<Real>&>,
                                            std::void_t<decltype(leja_argument(std::function{std::declval<F&>()}))>>>
{
   using type = decltype(leja_argument(std::function{std::declval<F&>()}));
};

template <class Point, class = void>
struct leja_resizable : std::false_type {};
template <class Point>
struct leja_resizable<Point, std::void_t<decltype(std::declval<Point&>().resize(0u))>> : std::true_type {};
}

// Dimension-adaptive sparse grid (Gerstner and Griebel) for the integral of f over the product of the rule's weight,
// as [-1, 1]^d for leja_legendre and exp(-|x|^2) over R^d for leja_hermite. Level k of a coordinate is the 2k - 1 point
// Leja rule U_k, and the multi-index k contributes Delta_k f = (U_{k_1} - U_{k_1 - 1}) x ... x (U_{k_d} - U_{k_d - 1}) f,
// with U_0 = 0. The levels grow by two points, not one: for even f, as x^2, the rules with 1 and 2 points can agree,
// and a zero contribution would stop the refinement in that coordinate. f takes a container of dimension Reals, resized if it can be.
template <class F, class Rule, class Point = typename detail::leja_point<F, typename Rule::value_type>::type>
auto sparse_grid(F f, const Rule& rule, unsigned dimension, std::size_t max_calls,
                 typename Rule::value_type tolerance = tools::root_epsilon<typename Rule::value_type>(),
                 typename Rule::value_type* error = nullptr, typename Rule::value_type* pL1 = nullptr)
   -> decltype(f(std::declval<const Point&>()))
{
   using Real = typename Rule::value_type;
   using K = decltype(f(std::declval<const Point&>()));
   using index = std::vector<unsigned>;
   static_assert(!std::is_integral<K>::value,
                 "The return type cannot be integral, it must be a real or complex floating point type, or a vector of them.");
   constexpr const char* function = "boost::math::quadrature::sparse_grid<%1%>";
   Point x{};
   if constexpr (detail::leja_resizable<Point>::value)
   {
      if (x.size() == 0)
         x.resize(dimension);
   }
   if (dimension == 0 || max_calls == 0 || rule.weights(1).empty() || static_cast<std::size_t>(x.size()) != dimension)
   {
      Real nan = policies::raise_domain_error(function, "The dimension (%1%) and max_calls must be > 0, and the dimension the size of f's argument.",
                                              Real(dimension), typename Rule::policy_type());
      if (error)
         *error = nan;
      if (pL1)
         *pL1 = nan;
      return detail::leja_nan<K, Real>();
   }

   // delta[k][j], j < 2k - 1, are the weights of U_k - U_{k-1}; the rule computes them as they are needed.
   std::vector<std::vector<Real>> delta(1);
   auto possible = [&](unsigned level) {
      return 2 * level - 1 <= rule.max_points();
   };
   auto grow = [&](unsigned level) {
      for (unsigned k = static_cast<unsigned>(delta.size()); k <= level; ++k)
      {
         delta.push_back(rule.weights(2 * k - 1));
         for (unsigned j = 0; k > 1 && j < 2 * k - 3; ++j)
            delta[k][j] -= rule.weights(2 * k - 3)[j];
      }
   };
   grow(1);

   // The points, in a trie by their indices j_i into rule.abscissa(): kids[n][j] is the child of node n, at depth i, for
   // j_i = j; at depth d - 1 it indexes points, which hold f and the weight of the point in the sparse grid.
   constexpr std::size_t none = static_cast<std::size_t>(-1);
   std::vector<std::vector<std::size_t>> kids(1);
   std::vector<std::pair<K, Real>> points;
   Real L1 = 0;
   // Adds c times the product of delta[k_l][j_l], l >= i, times f(x_j) to sum over the box j_l < 2 k_l - 1.
   auto walk = [&](auto& self, const index& k, std::size_t n, unsigned i, const Real& c, K& sum, bool& first) -> void {
      const unsigned m = 2 * k[i] - 1;
      if (kids[n].size() < m)
         kids[n].resize(m, none);
      for (unsigned j = 0; j < m; ++j)
      {
         const Real cj = c * delta[k[i]][j];
         x[i] = rule.abscissa()[j];
         std::size_t next = kids[n][j];
         if (next == none)
         {
            next = i + 1 < dimension ? kids.size() : points.size();
            kids[n][j] = next;
            if (i + 1 < dimension)
               kids.emplace_back();
            else
               points.emplace_back(f(x), Real(0));
         }
         if (i + 1 < dimension)
         {
            self(self, k, next, i + 1, cj, sum, first);
            continue;
         }
         detail::leja_accumulate(sum, cj, points[next].first, first);
         first = false;
         points[next].second += cj;
         L1 += cj * detail::leja_distance<Real>(points[next].first);
      }
   };
   // Adds Delta_k f to result, evaluating f at the new points, and returns its max norm. The sum of the Delta_k rounds far
   // less than the sum over the points, whose weights are large and of both signs.
   K result{};
   bool result_first = true;
   auto add = [&](const index& k) {
      K sum{};
      bool first = true;
      walk(walk, k, 0, 0, Real(1), sum, first);
      detail::leja_accumulate(result, Real(1), sum, result_first);
      result_first = false;
      return detail::leja_distance<Real>(sum);
   };
   // The new points of k are those with j_i in {2 k_i - 3, 2 k_i - 2} where k_i > 1, and j_i = 0 elsewhere.
   auto cost = [&](const index& k) {
      std::size_t n = 1;
      for (unsigned i = 0; i < dimension; ++i)
         n *= k[i] > 1 ? 2 : 1;
      return n;
   };

   // The active indices, with the norms of their contributions, and the old ones, whose forward neighbours are added.
   std::map<index, Real> active;
   std::set<index> old;
   // The active indices by decreasing norm, ties to the smallest, and the running sum of their norms; added, the sum of all
   // the norms ever added, bounds its rounding error.
   auto by_norm = [](const std::pair<Real, index>& a, const std::pair<Real, index>& b) {
      return a.first > b.first || (a.first == b.first && a.second < b.second);
   };
   std::set<std::pair<Real, index>, decltype(by_norm)> queue(by_norm);
   Real active_sum = 0, added = 0;
   std::size_t updates = 0;
   auto activate = [&](const index& m) {
      const Real v = add(m);
      active[m] = v;
      if (v == v)
         queue.emplace(v, m);
      active_sum += v;
      added += v;
      ++updates;
   };
   auto deactivate = [&](const index& k) {
      auto a = active.find(k);
      if (a->second == a->second)
         queue.erase({a->second, k});
      active_sum -= a->second;
      active.erase(a);
      ++updates;
   };
   // The contributions of old indices that the rule has too few points to refine stay in the error estimate.
   Real unrefined = 0;
   activate(index(dimension, 1));
   std::size_t calls = 1;
   Real err = 0;
   // When the estimate passes, every active index is refined once more, and the result is accepted only if the estimate
   // still passes, since a contribution can vanish by symmetry or by chance.
   bool checked = false;
   for (;;)
   {
      err = unrefined + active_sum;
      // Near the tolerance the running sum's rounding matters: sum again.
      if (err <= tolerance * L1 + (updates + 1) * tools::epsilon<Real>() * (added + unrefined))
      {
         active_sum = 0;
         for (const auto& a : active)
            active_sum += a.second;
         added = active_sum;
         updates = 0;
         err = unrefined + active_sum;
      }
      std::vector<index> batch;
      // A NaN from f: no refinement can help.
      if (active.empty() || (err <= tolerance * L1 && checked) || err != err || L1 != L1)
         break;
      if (err <= tolerance * L1)
      {
         for (const auto& a : active)
            batch.push_back(a.first);
         checked = true;
      }
      else
      {
         batch.push_back(queue.begin()->second);
         checked = false;
      }
      // The admissible forward neighbours, whose backward neighbours are all old once the batch is.
      std::set<index> next;
      std::size_t new_calls = 0;
      Real capped = 0;
      for (const index& k : batch)
         old.insert(k);
      for (const index& k : batch)
      {
         bool full = false;
         for (unsigned i = 0; i < dimension; ++i)
         {
            index m = k;
            if (!possible(++m[i]))
            {
               full = true;
               continue;
            }
            bool admissible = !active.count(m) && !old.count(m);
            for (unsigned l = 0; l < dimension && admissible; ++l)
            {
               if (m[l] > 1)
               {
                  --m[l];
                  admissible = old.count(m) > 0;
                  ++m[l];
               }
            }
            if (admissible && next.insert(m).second)
               new_calls += cost(m);
         }
         if (full)
            capped += active[k];
      }
      if (calls + new_calls > max_calls)
      {
         // Out of calls before the check is done: the estimate is not to be trusted.
         if (checked)
            err = std::numeric_limits<Real>::quiet_NaN();
         break;
      }
      unrefined += capped;
      for (const index& m : next)
         grow(*std::max_element(m.begin(), m.end()));
      for (const index& k : batch)
         deactivate(k);
      for (const index& m : next)
         activate(m);
      calls += new_calls;
   }

   L1 = 0;
   for (const auto& p : points)
      L1 += p.second * detail::leja_distance<Real>(p.first);
   if (error)
      *error = err;
   if (pL1)
      *pL1 = L1;
   return result;
}

} // namespace quadrature
BOOST_MATH_NAMESPACE_END

#endif // BOOST_MATH_QUADRATURE_SPARSE_GRID_HPP
