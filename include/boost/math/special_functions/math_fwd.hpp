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

//#include <boost/math/tools/promotion.hpp> // for argument promotion.

namespace boost
{
	namespace math
	{ // Math functions (in roughly alphabetic order).

   // Beta functions.
   template <class RT>
   RT beta(RT a, RT b); // Beta function (2 arguments).
   template <class RT>
   RT beta(RT a, RT b, RT x);// Beta function (3 arguments).
   template <class RT>
   RT betac(RT a, RT b, RT x);
   template <class RT> 
   RT ibeta(RT a, RT b, RT x); // Incomplete beta function.
   template <class RT>
   RT ibetac(RT a, RT b, RT x); // Incomplete beta complement function.
   template <class RT>
   RT ibeta_inv(RT a, RT b, RT p); // Incomplete beta inverse function.
   template <class RT>
   RT ibetac_inv(RT a, RT b, RT q); // Incomplete beta complement inverse function.


   // cbrt - cube root.
   template <class RT>
   RT cbrt(RT z);


   // erf & erfc error functions.
   template <class RT> // Error function.
   RT erf(RT z);
   template <class RT>// Error function complement.
   RT erfc(RT z);
   template <class RT>// Error function inverse.
   RT erf_inv(RT z);
   template <class RT>// Error function complement inverse.
   RT erfc_inv(RT z);

   // Exp (x minus 1) functions.
   template <class T>
   T expm1(T);

   // Factorial functions.
   // Note: not for integral types, at present.
   template <class RT>
   struct max_factorial;
   template <class RT>
   RT factorial(unsigned int);
   template <class RT>
   RT unchecked_factorial(unsigned int); 


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
   template <class T>
   T hypot(T, T);

   // log1p is log(x + 1)
   template <class T>
   T log1p(T);

   // Power
   template <class T>
   T powm1(const T, const T);

   // sqrt
   template <class T>
   T sqrt1pm1(const T&);

 	} // namespace math
} // namespace boost

#endif // BOOST_MATH_SPECIAL_MATH_FWD_HPP
