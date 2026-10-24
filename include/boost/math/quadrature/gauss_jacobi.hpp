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
      // gauss_jacobi raises the domain error.
      if (!in_domain())
      {
         abscissa_vals.assign(N, std::numeric_limits<Real>::quiet_NaN());
         weight_vals.assign(N, std::numeric_limits<Real>::quiet_NaN());
         anti_gauss_internal = false;
         return;
      }
      BOOST_MATH_NAMESPACE::detail::jacobi_family<Real> family(N, alpha, beta);
      polynomial value(N, family);
      abscissa_vals = value.abscissa();
      weight_vals = value.weights();
      // Laurie 1996, Theorem 3: the anti-Gauss rule has a node outside [-1, 1] when alpha or beta is below -1/2
      // and N is large, and for some N when it is near -1/2. At -1/2 it can have nodes on -1 and 1.
      anti_gauss_internal = polynomial::anti_gauss_is_internal(N + 1, family, Real(1)) && polynomial::anti_gauss_is_internal(N + 1, family, Real(-1));
      // Laurie's anti-Gauss rule with N + 1 points, for the error estimate.
      if (anti_gauss_internal)
      {
         polynomial anti(N + 1, BOOST_MATH_NAMESPACE::detail::jacobi_family<Real>(N + 1, alpha, beta, true), true);
         anti_abscissa_vals = anti.abscissa();
         anti_weight_vals = anti.weights();
         // Nodes on the endpoints come out a few eps either side.
         for (Real& x : anti_abscissa_vals)
            x = (std::min)((std::max)(x, Real(-1)), Real(1));
      }
    }
   const std::vector<Real>& abscissa() const
   {
      return abscissa_vals;
   }

   const std::vector<Real>& weights() const
   {
      return weight_vals;
   }

   bool anti_gauss_is_internal() const
   {
      return anti_gauss_internal;
   }

   const std::vector<Real>& anti_gauss_abscissa() const
   {
      return anti_abscissa_vals;
   }

   const std::vector<Real>& anti_gauss_weights() const
   {
      return anti_weight_vals;
   }

protected:
   bool in_domain() const
   {
      return alpha > -1 && beta > -1 && (BOOST_MATH_NAMESPACE::isfinite)(alpha) && (BOOST_MATH_NAMESPACE::isfinite)(beta);
   }
};

} // namespace detail

template <class Real, class Policy = BOOST_MATH_NAMESPACE::policies::policy<> >
class gauss_jacobi : public detail::gauss_jacobi_detail<Real, detail::quadrature_constant_category<Real>::value>
{
   // When the anti-Gauss rule has a node outside [-1, 1], the error estimate is NaN rather than evaluate f there;
   // the integral itself is unaffected.
   Real* estimable(Real* error) const
   {
      if (error && !this->anti_gauss_is_internal())
      {
         *error = std::numeric_limits<Real>::quiet_NaN();
         return nullptr;
      }
      return error;
   }

   template <class F, class K, class Norm>
   Real anti_gauss_error(F& f, K anti_result, const K& result, Real L1, Norm norm) const
   {
      return detail::anti_gauss_error(f, this->anti_gauss_abscissa(), this->anti_gauss_weights(), false, anti_result, result, L1, norm, "boost::math::quadrature::gauss_jacobi<%1%>::integrate", Policy());
   }

public:
   gauss_jacobi(unsigned N_, Real alpha_, Real beta_) : detail::gauss_jacobi_detail<Real, detail::quadrature_constant_category<Real>::value>(N_, alpha_, beta_)
   {
      if (!(alpha_ > -1) || !(BOOST_MATH_NAMESPACE::isfinite)(alpha_))
      {
         policies::raise_domain_error("boost::math::quadrature::gauss_jacobi<%1%>::gauss_jacobi",
            "The parameter alpha must be finite and > -1, but got %1%.", alpha_, Policy());
      }
      if (!(beta_ > -1) || !(BOOST_MATH_NAMESPACE::isfinite)(beta_))
      {
         policies::raise_domain_error("boost::math::quadrature::gauss_jacobi<%1%>::gauss_jacobi",
            "The parameter beta must be finite and > -1, but got %1%.", beta_, Policy());
      }
   }

   // Integrates f(x) (1 - x)^alpha * (1 + x)^beta over [-1, 1].
   template <class F>
   auto integrate(F f, Real* error = nullptr, Real* pL1 = nullptr) const->decltype(std::declval<F>()(std::declval<Real>()))
   {
      // In many math texts, K represents the field of real or complex numbers.
      // Too bad we can't put blackboard bold into C++ source!
      using K = decltype(f(Real(0)));
      static_assert(!std::is_integral<K>::value,
                    "The return type cannot be integral, it must be either a real or complex floating point type.");
      using std::abs;
      error = estimable(error);
      Real nan = std::numeric_limits<Real>::quiet_NaN();
      if (!this->in_domain() || detail::overflowed(this->abscissa(), "boost::math::quadrature::gauss_jacobi<%1%>::integrate", pL1, &nan, Policy()))
      {
         if (pL1)
            *pL1 = nan;
         if (error)
            *error = nan;
         return static_cast<K>(nan);
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
      if (error)
         *error = anti_gauss_error(f, K(Real(0)), result, L1, [](const K& z) { using std::abs; return abs(z); });
      return result;
   }

   // Explicit zero and norm for vector- and matrix-valued integrands.
   template <class F, class Norm>
   auto integrate(F f, const decltype(std::declval<F>()(std::declval<Real>()))& zero, Norm norm, Real* error = nullptr, Real* pL1 = nullptr) const
      ->decltype(static_cast<Real>(norm(f(Real(0)))), f(Real(0)))
   {
      using K = decltype(f(Real(0)));
      static_assert(!std::is_integral<K>::value, "The return type cannot be integral.");
      error = estimable(error);
      Real nan = std::numeric_limits<Real>::quiet_NaN();
      if (!this->in_domain() || detail::overflowed(this->abscissa(), "boost::math::quadrature::gauss_jacobi<%1%>::integrate", pL1, &nan, Policy()))
      {
         if (pL1)
            *pL1 = nan;
         if (error)
            *error = nan;
         K result = zero * nan;
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
      if (error)
         *error = anti_gauss_error(f, zero, result, L1, norm);
      return result;
   }
};
} // namespace quadrature

BOOST_MATH_NAMESPACE_END

#endif // BOOST_MATH_QUADRATURE_GAUSS_JACOBI_HPP
