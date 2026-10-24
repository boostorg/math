//  Copyright John Maddock 2015.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)


#ifndef BOOST_MATH_QUADRATURE_DETAIL_QUADRATURE_CONSTANT_HPP
#define BOOST_MATH_QUADRATURE_DETAIL_QUADRATURE_CONSTANT_HPP

#ifndef BOOST_MATH_BUILD_MODULE
#include <algorithm>
#include <limits>
#include <type_traits>
#endif

BOOST_MATH_NAMESPACE_BEGIN namespace quadrature{ namespace detail {

template <class T>
struct quadrature_constant_category
{
   static const unsigned value =
      (std::numeric_limits<T>::is_specialized == 0) ? 999 :
      (std::numeric_limits<T>::radix == 2) ?
      (
#ifdef BOOST_HAS_FLOAT128
         (std::numeric_limits<T>::digits <= 113) && std::is_constructible<T, __float128>::value ? 0 :
#else
         (std::numeric_limits<T>::digits <= std::numeric_limits<long double>::digits) && std::is_constructible<T, long double>::value ? 0 :
#endif
         (std::numeric_limits<T>::digits10 <= 110) && std::is_constructible<T, const char*>::value ? 4 : 999
      ) : (std::numeric_limits<T>::digits10 <= 110) && std::is_constructible<T, const char*>::value ? 4 : 999;

   using storage_type =
      std::conditional_t<(std::numeric_limits<T>::is_specialized == 0), T,
         std::conditional_t<(std::numeric_limits<T>::radix == 2),
            std::conditional_t< ((std::numeric_limits<T>::digits <= std::numeric_limits<float>::digits) && std::is_constructible<T, float>::value),
               float,
               std::conditional_t<((std::numeric_limits<T>::digits <= std::numeric_limits<double>::digits) && std::is_constructible<T, double>::value),
                  double,
                  std::conditional_t<((std::numeric_limits<T>::digits <= std::numeric_limits<long double>::digits) && std::is_constructible<T, long double>::value),
                     long double,
#ifdef BOOST_HAS_FLOAT128
                     std::conditional_t<((std::numeric_limits<T>::digits <= 113) && std::is_constructible<T, __float128>::value),
                        __float128,
                        T
                     >
                  >
#else
                     T
                  >
#endif
               >
            >, T
         >
      >;
};

template<typename RandomAccessContainer, typename Real, typename Policy>
bool overflowed(const RandomAccessContainer& abscissa, const char* function, Real* pL1, Real* error, const Policy& pol)
{
   if (!(BOOST_MATH_NAMESPACE::isfinite)(static_cast<Real>(abscissa.back())))
   {
      *error = policies::raise_evaluation_error(
               function,
               "The recurrence relation overflows with this many points; use fewer points or a type with a wider exponent range.",
               std::numeric_limits<Real>::quiet_NaN(), pol);
      if (pL1)
         *pL1 = *error;
      return true;
   }
   return false;
}

// The averaged rule (G + A) / 2, exact for polynomials of degree 2N + 1. Half the difference is set as its error;
// it bounds the error when G and A bracket the integral. L1 goes in as that of G and comes out as that of the average.
template <class F, class RandomAccessContainer, class K, class Real, class Norm, class Policy>
K anti_gauss_average(F& f, const RandomAccessContainer& abscissa, const RandomAccessContainer& weights, bool symmetric,
                     K anti_result, const K& result, Real& L1, Norm norm, Real* error, const char* function, const Policy& pol)
{
   Real nan;
   if (!(BOOST_MATH_NAMESPACE::isfinite)(L1) || overflowed(abscissa, function, static_cast<Real*>(nullptr), &nan, pol))
   {
      L1 = std::numeric_limits<Real>::quiet_NaN();
      if (error)
         *error = L1;
      return anti_result * L1;
   }
   Real anti_L1 = 0;
   for (unsigned i = 0; i < abscissa.size(); ++i)
   {
      Real x = abscissa[i];
      K fx = f(x);
      Real n = static_cast<Real>(norm(fx));
      if (symmetric && x != 0)
      {
         K fmx = f(-x);
         n += static_cast<Real>(norm(fmx));
         fx += fmx;
      }
      anti_result += fx * static_cast<Real>(weights[i]);
      anti_L1 += n * static_cast<Real>(weights[i]);
   }
   L1 = (L1 + anti_L1) / 2;
   if (error)
      *error = (std::max)(static_cast<Real>(norm(anti_result - result)) / 2, 2 * tools::epsilon<Real>() * L1);
   anti_result += result;
   return anti_result * Real(0.5);
}

} // namespace detail
} // namespace quadrature
BOOST_MATH_NAMESPACE_END

#endif // BOOST_MATH_QUADRATURE_DETAIL_QUADRATURE_CONSTANT_HPP