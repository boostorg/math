//  Copyright Jacob Hass 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_QUADRATURE_GAUSS_JACOBI_HPP
#define BOOST_MATH_QUADRATURE_GAUSS_JACOBI_HPP

#ifdef _MSC_VER
#pragma once
#endif

#include <boost/math/tools/config.hpp>

#ifdef BOOST_MATH_NO_CXX17_IF_CONSTEXPR
#error "The header <boost/math/quadrature/gauss_jacobi.hpp> requires C++17 or later."
#endif

#ifndef BOOST_MATH_BUILD_MODULE
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <type_traits>
#include <utility>
#include <vector>
#endif
#include <boost/math/constants/constants.hpp>
#include <boost/math/policies/error_handling.hpp>
#include <boost/math/policies/policy.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/math/tools/big_constant.hpp>
#include <boost/math/tools/precision.hpp>
#include <boost/math/quadrature/detail/quadrature_constant.hpp>
#include <boost/math/special_functions/detail/orthogonal_polynomial.hpp>
#include <boost/math/special_functions/jacobi.hpp>

BOOST_MATH_NAMESPACE_BEGIN namespace quadrature { namespace detail {

template <class Real, unsigned Category>
class gauss_jacobi_detail
{
   unsigned N;
   Real alpha;
   Real beta;
   std::vector<Real> abscissa_vals;
   std::vector<Real> weight_vals;
   bool anti_gauss_internal;
   std::vector<Real> anti_abscissa_vals;
   std::vector<Real> anti_weight_vals;
   using polynomial = BOOST_MATH_NAMESPACE::detail::orthogonal_polynomial<Real, BOOST_MATH_NAMESPACE::detail::jacobi_family<Real> >;
public:
   gauss_jacobi_detail(unsigned N_, Real alpha_, Real beta_) : N(N_), alpha(alpha_), beta(beta_)
   {
      BOOST_MATH_NAMESPACE::detail::jacobi_family<Real> family(N, alpha, beta);
      polynomial value(N, family);
      abscissa_vals = value.abscissa();
      weight_vals = value.weights();
      // Laurie 1996, Theorem 3: the anti-Gauss rule has a node outside (-1, 1) when alpha or beta is below -1/2
      // and N is large, and for some N when it is near -1/2.
      anti_gauss_internal = polynomial::anti_gauss_is_internal(N + 1, family, Real(1)) && polynomial::anti_gauss_is_internal(N + 1, family, Real(-1));
    }
   const std::vector<Real>& abscissa()
   {
      return abscissa_vals;
   }

   const std::vector<Real>& weights()
   {
      return weight_vals;
   }

   bool anti_gauss_is_internal() const
   {
      return anti_gauss_internal;
   }

   // Laurie's anti-Gauss rule with N + 1 points, computed the first time an error estimate is requested.
   const std::vector<Real>& anti_gauss_abscissa()
   {
      if (anti_abscissa_vals.empty())
      {
         polynomial value(N + 1, BOOST_MATH_NAMESPACE::detail::jacobi_family<Real>(N + 1, alpha, beta), true);
         anti_abscissa_vals = value.abscissa();
         anti_weight_vals = value.weights();
      }
      return anti_abscissa_vals;
   }

   const std::vector<Real>& anti_gauss_weights()
   {
      anti_gauss_abscissa();
      return anti_weight_vals;
   }
};

} // namespace detail

template <class Real, class Policy = BOOST_MATH_NAMESPACE::policies::policy<> >
class gauss_jacobi : public detail::gauss_jacobi_detail<Real, detail::quadrature_constant_category<Real>::value>
{

public:
   gauss_jacobi(unsigned N_, Real alpha_, Real beta_) : detail::gauss_jacobi_detail<Real, detail::quadrature_constant_category<Real>::value>(N_, alpha_, beta_) {}

   // Integrates f(x) (1 - x)^alpha * (1 + x)^beta over [-1, 1].
   template <class F>
   auto integrate(F f, Real* pL1 = nullptr)->decltype(std::declval<F>()(std::declval<Real>()))
   {
      // In many math texts, K represents the field of real or complex numbers.
      // Too bad we can't put blackboard bold into C++ source!
      using K = decltype(f(Real(0)));
      static_assert(!std::is_integral<K>::value,
                    "The return type cannot be integral, it must be either a real or complex floating point type.");
      using std::abs;
      const char* function = "boost::math::quadrature::gauss_jacobi<%1>::integrate";
      Real error;
      if (detail::overflowed(this->abscissa(), function, pL1, &error, Policy()))
      {
         return static_cast<K>(error);
      }
      K result = Real(0);
      Real L1 = abs(result);
      for (unsigned i = 0; i < this->abscissa().size(); ++i)
      {
         K fp = f(static_cast<Real>(this->abscissa()[i]));
         result += fp * static_cast<Real>(this->weights()[i]);
         L1 += abs(fp) * static_cast<Real>(this->weights()[i]);
      }
      if (pL1)
         *pL1 = L1;
      return result;
   }

   // Also sets *error to half the difference from Laurie's anti-Gauss rule, which takes N + 1 more evaluations of f.
   // If that rule has a node on or outside [-1, 1], there is no estimate, which is an evaluation error.
   template <class F>
   auto integrate(F f, Real* error, Real* pL1)->decltype(std::declval<F>()(std::declval<Real>()))
   {
      using K = decltype(f(Real(0)));
      using std::abs;
      if (error && !this->anti_gauss_is_internal())
      {
         *error = policies::raise_evaluation_error(
            "boost::math::quadrature::gauss_jacobi<%1%>::integrate",
            "The anti-Gauss rule for the error estimate has a node on or outside [-1, 1] for these alpha, beta and points.",
            std::numeric_limits<Real>::quiet_NaN(), Policy());
         error = nullptr;
      }
      Real L1;
      K result = integrate(f, &L1);
      if (pL1)
         *pL1 = L1;
      if (!error)
         return result;
      const char* function = "boost::math::quadrature::gauss_jacobi<%1%>::integrate";
      if (detail::overflowed(this->abscissa(), function, nullptr, error, Policy()) || detail::overflowed(this->anti_gauss_abscissa(), function, nullptr, error, Policy()))
         return result;
      K anti_result = Real(0);
      for (unsigned i = 0; i < this->anti_gauss_abscissa().size(); ++i)
      {
         anti_result += f(this->anti_gauss_abscissa()[i]) * this->anti_gauss_weights()[i];
      }
      // Rounding makes the difference meaningless below a few ulps of L1, and aliasing can make it exactly zero.
      *error = (std::max)(static_cast<Real>(abs(anti_result - result)) / 2, 2 * tools::epsilon<Real>() * L1);
      return result;
   }

   // Explicit zero and norm for vector- and matrix-valued integrands.
   template <class F, class Norm>
   auto integrate(F f, const decltype(std::declval<F>()(std::declval<Real>()))& zero, Norm norm, Real* pL1 = nullptr)
      ->decltype(static_cast<Real>(norm(f(Real(0)))), f(Real(0)))
   {
      using K = decltype(f(Real(0)));
      static_assert(!std::is_integral<K>::value, "The return type cannot be integral.");
      const char* function = "boost::math::quadrature::gauss_jacobi<%1>::integrate";
      Real error;
      if (detail::overflowed(this->abscissa(), function, pL1, &error, Policy()))
      {
         K result = zero * error;
         return result;
      }
      K result = zero;
      Real L1 = static_cast<Real>(norm(result));
      for (unsigned i = 0; i < this->abscissa().size(); ++i)
      {
         K fp = f(static_cast<Real>(this->abscissa()[i]));
         result += fp * static_cast<Real>(this->weights()[i]);
         L1 += static_cast<Real>(norm(fp)) * static_cast<Real>(this->weights()[i]);
      }
      if (pL1)
         *pL1 = L1;
      return result;
   }
};
} // namespace quadrature

BOOST_MATH_NAMESPACE_END

#endif // BOOST_MATH_QUADRATURE_GAUSS_JACOBI_HPP
