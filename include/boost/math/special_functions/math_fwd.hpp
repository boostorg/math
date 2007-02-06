// math_fwd.hpp

// TODO revise completely for new distribution classes.

// Copyright Paul A. Bristow 2006.
// Copyright John Maddock 2006.

// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

// Omnibus list of forward declarations of math special functions.

// IT = Integer type.
// RT = Real type (built-in floating-point types, float, double, long double) & User Defined Types
// AT = Integer or Real type 

#ifndef BOOST_MATH_SPECIAL_MATH_FWD_HPP
#define BOOST_MATH_SPECIAL_MATH_FWD_HPP

#include <boost/math/tools/promotion.hpp> // for argument promotion.
#include <complex>

namespace boost
{
	namespace math
	{ // Math functions (in roughly alphabetic order).

   // Beta functions.
   template <class RT1, class RT2>
   typename tools::promote_args<RT1, RT2>::type 
         beta(RT1 a, RT2 b); // Beta function (2 arguments).

   template <class RT1, class RT2, class RT3>
   typename tools::promote_args<RT1, RT2, RT3>::type 
         beta(RT1 a, RT2 b, RT3 x); // Beta function (3 arguments).

   template <class RT1, class RT2, class RT3>
   typename tools::promote_args<RT1, RT2, RT3>::type 
         betac(RT1 a, RT2 b, RT3 x);

   template <class RT1, class RT2, class RT3>
   typename tools::promote_args<RT1, RT2, RT3>::type 
         ibeta(RT1 a, RT2 b, RT3 x); // Incomplete beta function.

   template <class RT1, class RT2, class RT3>
   typename tools::promote_args<RT1, RT2, RT3>::type 
         ibetac(RT1 a, RT2 b, RT3 x); // Incomplete beta complement function.

   template <class T1, class T2, class T3, class T4>
   typename tools::promote_args<T1, T2, T3, T4>::type  
         ibeta_inv(T1 a, T2 b, T3 p, T4* py);

   template <class RT1, class RT2, class RT3>
   typename tools::promote_args<RT1, RT2, RT3>::type 
         ibeta_inv(RT1 a, RT2 b, RT3 p); // Incomplete beta inverse function.

   template <class T1, class T2, class T3, class T4>
   typename tools::promote_args<T1, T2, T3, T4>::type 
         ibetac_inv(T1 a, T2 b, T3 q, T4* py);

   template <class RT1, class RT2, class RT3>
   typename tools::promote_args<RT1, RT2, RT3>::type 
         ibetac_inv(RT1 a, RT2 b, RT3 q); // Incomplete beta complement inverse function.

   template <class RT1, class RT2, class RT3>
   typename tools::promote_args<RT1, RT2, RT3>::type 
         ibeta_derivative(RT1 a, RT2 b, RT3 x);  // derivative of incomplete beta

   // erf & erfc error functions.
   template <class RT> // Error function.
   typename tools::promote_args<RT>::type erf(RT z);

   template <class RT>// Error function complement.
   typename tools::promote_args<RT>::type erfc(RT z);

   template <class RT>// Error function inverse.
   typename tools::promote_args<RT>::type erf_inv(RT z);

   template <class RT>// Error function complement inverse.
   typename tools::promote_args<RT>::type erfc_inv(RT z);

   // Polynomials:
   template <class T1, class T2, class T3>
   typename tools::promote_args<T1, T2, T3>::type 
         legendre_next(unsigned l, T1 x, T2 Pl, T3 Plm1);

   template <class T>
   typename tools::promote_args<T>::type 
         legendre_p(int l, T x);

   template <class T>
   typename tools::promote_args<T>::type 
         legendre_q(unsigned l, T x);

   template <class T1, class T2, class T3>
   typename tools::promote_args<T1, T2, T3>::type 
         legendre_next(unsigned l, unsigned m, T1 x, T2 Pl, T3 Plm1);

   template <class T>
   typename tools::promote_args<T>::type 
         legendre_p(int l, int m, T x);

   template <class T1, class T2, class T3>
   typename tools::promote_args<T1, T2, T3>::type  
         laguerre_next(unsigned n, T1 x, T2 Ln, T3 Lnm1);

   template <class T1, class T2, class T3>
   typename tools::promote_args<T1, T2, T3>::type  
      laguerre_next(unsigned n, unsigned l, T1 x, T2 Pl, T3 Plm1);

   template <class T>
   typename tools::promote_args<T>::type 
      hermite(unsigned n, T x);

   template <class T1, class T2, class T3>
   typename tools::promote_args<T1, T2, T3>::type 
      hermite_next(unsigned n, T1 x, T2 Hn, T3 Hnm1);

   template <class T1, class T2>
   std::complex<typename tools::promote_args<T1, T2>::type> 
         spherical_harmonic(unsigned n, int m, T1 theta, T2 phi);

   template <class T1, class T2>
   typename tools::promote_args<T1, T2>::type 
         spherical_harmonic_r(unsigned n, int m, T1 theta, T2 phi);

   template <class T1, class T2>
   typename tools::promote_args<T1, T2>::type 
         spherical_harmonic_i(unsigned n, int m, T1 theta, T2 phi);

   // Elliptic integrals:
   template <class T1, class T2, class T3>
   typename tools::promote_args<T1, T2, T3>::type 
         ellint_rf(T1 x, T2 y, T3 z);

   template <class T1, class T2, class T3>
   typename tools::promote_args<T1, T2, T3>::type 
         ellint_rd(T1 x, T2 y, T3 z);

   template <class T1, class T2>
   typename tools::promote_args<T1, T2>::type 
         ellint_rc(T1 x, T2 y);

   template <class T1, class T2, class T3, class T4>
   typename tools::promote_args<T1, T2, T3, T4>::type 
         ellint_rj(T1 x, T2 y, T3 z, T4 p);

   template <typename T>
   typename tools::promote_args<T>::type ellint_2(T k);

   template <class T1, class T2>
   typename tools::promote_args<T1, T2>::type ellint_2(T1 k, T2 phi);

   template <typename T>
   typename tools::promote_args<T>::type ellint_1(T k);

   template <class T1, class T2>
   typename tools::promote_args<T1, T2>::type ellint_1(T1 k, T2 phi);

   template <class T1, class T2, class T3>
   typename tools::promote_args<T1, T2, T3>::type ellint_3(T1 k, T2 v, T3 phi);

   template <class T1, class T2>
   typename tools::promote_args<T1, T2>::type ellint_3(T1 k, T2 v);

   // Factorial functions.
   // Note: not for integral types, at present.
   template <class RT>
   struct max_factorial;
   template <class RT>
   RT factorial(unsigned int);
   template <class RT>
   RT unchecked_factorial(unsigned int); 
   template <class RT>
   RT double_factorial(unsigned i);

   template <class RT>
   typename tools::promote_args<RT>::type falling_factorial(RT x, unsigned n);

   template <class RT>
   typename tools::promote_args<RT>::type rising_factorial(RT x, unsigned n);


   // Fpclassify - classify floating-point as NaN or infinity...
   template <class T>
   int fpclassify (T);

   // Gamma functions.
   template <class RT>
   typename tools::promote_args<RT>::type tgamma(RT z);

   template <class RT>
   typename tools::promote_args<RT>::type tgamma1pm1(RT z);

   template <class RT1, class RT2>
   typename tools::promote_args<RT1, RT2>::type tgamma(RT1 a, RT2 z);

   template <class RT>
   typename tools::promote_args<RT>::type lgamma(RT z, int* sign);

   template <class RT>
   typename tools::promote_args<RT>::type lgamma(RT x);

   template <class RT1, class RT2>
   typename tools::promote_args<RT1, RT2>::type tgamma_lower(RT1 a, RT2 z);

   template <class RT1, class RT2>
   typename tools::promote_args<RT1, RT2>::type gamma_q(RT1 a, RT2 z);

   template <class RT1, class RT2>
   typename tools::promote_args<RT1, RT2>::type gamma_p(RT1 a, RT2 z);

   template <class T1, class T2>
   typename tools::promote_args<T1, T2>::type tgamma_delta_ratio(T1 z, T2 delta);

   template <class T1, class T2>
   typename tools::promote_args<T1, T2>::type tgamma_ratio(T1 a, T2 b);

   template <class T1, class T2>
   typename tools::promote_args<T1, T2>::type gamma_p_derivative(T1 a, T2 x);

   // gamma inverse.
   template <class T1, class T2>
   typename tools::promote_args<T1, T2>::type gamma_p_inv(T1 a, T2 p);

   template <class T1, class T2>
   typename tools::promote_args<T1, T2>::type gamma_q_inv(T1 a, T2 q);

   // digamma:
   template <class T>
   typename tools::promote_args<T>::type digamma(T x);

   // Hypotenuse function sqrt(x ^ 2 + y ^ 2).
   template <class T1, class T2>
   typename tools::promote_args<T1, T2>::type 
         hypot(T1 x, T2 y);

   // cbrt - cube root.
   template <class RT>
   typename tools::promote_args<RT>::type cbrt(RT z);

   // log1p is log(x + 1)
   template <class T>
   typename tools::promote_args<T>::type log1p(T);

   // Exp (x) minus 1 functions.
   template <class T>
   typename tools::promote_args<T>::type expm1(T);

   // Power - 1
   template <class T1, class T2>
   inline typename tools::promote_args<T1, T2>::type 
         powm1(const T1 a, const T2 z);

   // sqrt(1+x) - 1
   template <class T>
   typename tools::promote_args<T>::type sqrt1pm1(const T& val);

   // sinus cardinals:
   template <class T>
   typename tools::promote_args<T>::type sinc_pi(T x);

   template <class T>
   typename tools::promote_args<T>::type sinhc_pi(T x);

   // inverse hyperbolics:
   template<typename T>
   typename tools::promote_args<T>::type asinh(const T x);

   template<typename T>
   typename tools::promote_args<T>::type acosh(const T x);

   template<typename T>
   typename tools::promote_args<T>::type atanh(const T x);

 	} // namespace math
} // namespace boost

#endif // BOOST_MATH_SPECIAL_MATH_FWD_HPP
