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
public:
   gauss_jacobi_detail(unsigned N_, Real alpha_, Real beta_) : N(N_), alpha(alpha_), beta(beta_)
   {
      BOOST_MATH_NAMESPACE::detail::jacobi_family<Real> family(N, alpha, beta);
      BOOST_MATH_NAMESPACE::detail::orthogonal_polynomial<Real, BOOST_MATH_NAMESPACE::detail::jacobi_family<Real> > value(N, family);
      abscissa_vals = value.abscissa();
      weight_vals = value.weights();
    }
   const std::vector<Real>& abscissa()
   {
      return abscissa_vals;
   }

   const std::vector<Real>& weights()
   {
      return weight_vals;
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
