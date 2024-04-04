//  (C) Copyright John Maddock 2022.
//  (C) Copyright Matt Borland 2024.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

module;

#define BOOST_MATH_SF_AS_MODULE
#define BOOST_MATH_AS_MODULE
#define BOOST_MATH_STANDALONE

#if defined(__cpp_lib_modules) && __cpp_lib_modules >= 202207L

import std;
import std.compat;

#else

#include <cmath>
#include <climits>
#include <algorithm>
#include <type_traits>
#include <cstdint>
#include <limits>
#include <mutex>
#include <vector>
#include <string>

#endif

#include <boost/math/tools/promotion.hpp>


export module boost.math.special_functions;

import boost.math.core;
import boost.math.constants;

//
// Forward declarations in math_fwd.hpp are not available, so forward declare what we need to first:
//
export namespace boost::math {

   template <class RT>
   struct max_factorial;
   template <class RT>
   RT factorial(unsigned int);
   template <class RT, class Policy>
   RT factorial(unsigned int, const Policy& pol);
   template <class RT>
   RT unchecked_factorial(unsigned int BOOST_MATH_APPEND_EXPLICIT_TEMPLATE_TYPE(RT));
   template <class RT>
   RT double_factorial(unsigned i);
   template <class RT, class Policy>
   RT double_factorial(unsigned i, const Policy& pol);

   template <class RT>
   typename tools::promote_args<RT>::type falling_factorial(RT x, unsigned n);

   template <class RT, class Policy>
   typename tools::promote_args<RT>::type falling_factorial(RT x, unsigned n, const Policy& pol);

   template <class RT>
   typename tools::promote_args<RT>::type rising_factorial(RT x, int n);

   template <class RT, class Policy>
   typename tools::promote_args<RT>::type rising_factorial(RT x, int n, const Policy& pol);

   // Gamma functions.
   template <class RT>
   typename tools::promote_args<RT>::type tgamma(RT z);

   template <class RT>
   typename tools::promote_args<RT>::type tgamma1pm1(RT z);

   template <class RT, class Policy>
   typename tools::promote_args<RT>::type tgamma1pm1(RT z, const Policy& pol);

   template <class RT1, class RT2>
   typename tools::promote_args<RT1, RT2>::type tgamma(RT1 a, RT2 z);

   template <class RT1, class RT2, class Policy>
   typename tools::promote_args<RT1, RT2>::type tgamma(RT1 a, RT2 z, const Policy& pol);

   template <class RT>
   typename tools::promote_args<RT>::type lgamma(RT z, int* sign);

   template <class RT, class Policy>
   typename tools::promote_args<RT>::type lgamma(RT z, int* sign, const Policy& pol);

   template <class RT>
   typename tools::promote_args<RT>::type lgamma(RT x);

   template <class RT, class Policy>
   typename tools::promote_args<RT>::type lgamma(RT x, const Policy& pol);

   template <class RT1, class RT2>
   typename tools::promote_args<RT1, RT2>::type tgamma_lower(RT1 a, RT2 z);

   template <class RT1, class RT2, class Policy>
   typename tools::promote_args<RT1, RT2>::type tgamma_lower(RT1 a, RT2 z, const Policy&);

   template <class RT1, class RT2>
   typename tools::promote_args<RT1, RT2>::type gamma_q(RT1 a, RT2 z);

   template <class RT1, class RT2, class Policy>
   typename tools::promote_args<RT1, RT2>::type gamma_q(RT1 a, RT2 z, const Policy&);

   template <class RT1, class RT2>
   typename tools::promote_args<RT1, RT2>::type gamma_p(RT1 a, RT2 z);

   template <class T1, class T2>
   typename tools::promote_args<T1, T2>::type tgamma_delta_ratio(T1 z, T2 delta);

   template <class T1, class T2, class Policy>
   typename tools::promote_args<T1, T2>::type tgamma_delta_ratio(T1 z, T2 delta, const Policy&);

   template <class T1, class T2>
   typename tools::promote_args<T1, T2>::type tgamma_ratio(T1 a, T2 b);

   template <class T1, class T2, class Policy>
   typename tools::promote_args<T1, T2>::type tgamma_ratio(T1 a, T2 b, const Policy&);

   template <class T1, class T2>
   typename tools::promote_args<T1, T2>::type gamma_p_derivative(T1 a, T2 x);

   template <class T1, class T2, class Policy>
   typename tools::promote_args<T1, T2>::type gamma_p_derivative(T1 a, T2 x, const Policy&);

   // gamma inverse.
   template <class T1, class T2>
   typename tools::promote_args<T1, T2>::type gamma_p_inv(T1 a, T2 p);

   template <class T1, class T2, class Policy>
   typename tools::promote_args<T1, T2>::type gamma_p_inva(T1 a, T2 p, const Policy&);

   template <class T1, class T2>
   typename tools::promote_args<T1, T2>::type gamma_p_inva(T1 a, T2 p);

   template <class T1, class T2, class Policy>
   typename tools::promote_args<T1, T2>::type gamma_p_inv(T1 a, T2 p, const Policy&);

   template <class T1, class T2>
   typename tools::promote_args<T1, T2>::type gamma_q_inv(T1 a, T2 q);

   template <class T1, class T2, class Policy>
   typename tools::promote_args<T1, T2>::type gamma_q_inv(T1 a, T2 q, const Policy&);

   template <class T1, class T2>
   typename tools::promote_args<T1, T2>::type gamma_q_inva(T1 a, T2 q);

   template <class T1, class T2, class Policy>
   typename tools::promote_args<T1, T2>::type gamma_q_inva(T1 a, T2 q, const Policy&);

   // digamma:
   template <class T>
   typename tools::promote_args<T>::type digamma(T x);

   template <class T, class Policy>
   typename tools::promote_args<T>::type digamma(T x, const Policy&);

   // trigamma:
   template <class T>
   typename tools::promote_args<T>::type trigamma(T x);

   template <class T, class Policy>
   typename tools::promote_args<T>::type trigamma(T x, const Policy&);

   // polygamma:
   template <class T>
   typename tools::promote_args<T>::type polygamma(int n, T x);

   template <class T, class Policy>
   typename tools::promote_args<T>::type polygamma(int n, T x, const Policy&);

   // erf & erfc error functions.
   template <class RT> // Error function.
   typename tools::promote_args<RT>::type erf(RT z);
   template <class RT, class Policy> // Error function.
   typename tools::promote_args<RT>::type erf(RT z, const Policy&);

   template <class RT>// Error function complement.
   typename tools::promote_args<RT>::type erfc(RT z);
   template <class RT, class Policy>// Error function complement.
   typename tools::promote_args<RT>::type erfc(RT z, const Policy&);

   template <class RT>// Error function inverse.
   typename tools::promote_args<RT>::type erf_inv(RT z);
   template <class RT, class Policy>// Error function inverse.
   typename tools::promote_args<RT>::type erf_inv(RT z, const Policy& pol);

   template <class RT>// Error function complement inverse.
   typename tools::promote_args<RT>::type erfc_inv(RT z);
   template <class RT, class Policy>// Error function complement inverse.
   typename tools::promote_args<RT>::type erfc_inv(RT z, const Policy& pol);

   // Zeta:
   template <class T, class Policy>
   typename tools::promote_args<T>::type zeta(T s, const Policy&);
}

#include <boost/math/special_functions/log1p.hpp>
#include <boost/math/special_functions/expm1.hpp>
#include <boost/math/special_functions/powm1.hpp>
#include <boost/math/special_functions/sqrt1pm1.hpp>
#include <boost/math/special_functions/detail/unchecked_factorial.hpp>
#include <boost/math/special_functions/bernoulli.hpp>
#include <boost/math/special_functions/gamma.hpp>
