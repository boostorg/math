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

   // Binomial distribution.

   //template <class AT, class RT> // Binomial distribution (k, n, x)
   //// Probability of number of events between 0 and k-1 inclusive, if expected probability of success events is success_fraction.
   //typename tools::promote_arg3<AT, AT, RT>::type
   //// Return type is the wider of the two (perhaps promoted) floating-point types.
   //binomial(AT k, AT n, RT success_fraction);

   //template <class AT, class RT> // Binomial distribution complement (k, n, x)
   //// Probability of number of events between 0 and k-1 inclusive, if expected mean is x.
   //typename tools::promote_arg3<AT, AT, RT>::type
   //binomial_c(AT k, AT n, RT success_fraction);

   //template <class AT, class RT>
   //// success_fraction if number of events between 0 and k-1 inclusive, and probability p.
   //typename tools::promote_arg3<AT, AT, RT>::type
   //binomial_inv(AT k,  AT n, RT probability); // Binomial distribution inverse (k, n, p)

   // cbrt - cube root.
   template <class RT>
   RT cbrt(RT z);

   // chi_sqr functions.
   //template <class AT, class RT> // probability chisqr(degrees_of_freedom, chisqr)
   //typename tools::promote_arg2<RT, AT>::type
   //// return type is the wider of the two (?promoted) floating point types.
   //chisqr(AT degrees_of_freedom, RT chisqr);

   //template <class AT, class RT> // complement probability chisqr_c(degrees_of_freedom, chisqr)
   //typename tools::promote_arg2<RT, AT>::type
   //chisqr_c(AT degrees_of_freedom, RT chisqr);

   //template <class AT, class RT> // chisqr = chisqr_inv(degrees_of_freedom,  probability)
   //typename tools::promote_arg2<RT, AT>::type
   //chisqr_inv(AT degrees_of_freedom, RT probability);

   //template <class AT, class RT> // degrees_of_freedom = chisqr_inv_df(chisqr, probability)
   //typename tools::promote_arg2<RT, AT>::type // <RT, AT> but both RT.
   //chisqr_df_inv(RT chisqr, RT probability);

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

   // Fisher-Snedecor functions.

   //template <class AT, class RT>
   //typename tools::promote_arg3<AT, AT, RT>::type
   //// return type is the wider of the three (possibly promoted) floating point types.
   //fisher(AT degrees_of_freedom1, AT degrees_of_freedom2, RT fisher);

   //template <class AT, class RT>
   //typename tools::promote_arg3<AT, AT, RT>::type
   //fisher_c(AT degrees_of_freedom1, AT degrees_of_freedom2, RT fisher);

   //template <class AT, class RT>
   //typename tools::promote_arg3<AT, AT, RT>::type
   //fisher_inv(AT degrees_of_freedom1, AT degrees_of_freedom2, RT fisher);

   //template <class AT, class RT>
   //typename tools::promote_arg3<AT, AT, RT>::type
   //fisher_c_inv(AT degrees_of_freedom1, AT degrees_of_freedom2, RT fisher);

   // Fpclassify - classify floating-point as NaN or infinity...
   template <class T>
   int fpclassify (T);

   // Gamma functions.
   template <class RT>
   RT tgamma(RT z);
   template <class RT>
   RT tgamma(RT a, RT z);
   template <class RT>
   RT lgamma(RT z, int* sign);
   template <class RT>
   RT lgamma(RT x);
   template <class RT>
   RT tgamma_lower(RT a, RT z);
   template <class RT>
   RT gamma_Q(RT a, RT z);
   template <class RT>
   RT gamma_P(RT a, RT z);
   // gamma inverse.
   template <class RT>
   RT gamma_P_inv(RT a, RT p);
   template <class RT>
   RT gamma_Q_inv(RT a, RT q);

   // Hypotenuse function sqrt(x ^ 2 + y ^ 2).
   template <class T>
   T hypot(T, T);

   // log1p is log(x + 1)
   template <class T>
   T log1p(T);

   // Power
   template <class T>
   T powm1(const T, const T);

   //template <class AT, class RT>
   //typename tools::promote_arg2<RT, AT>::type
   //students_t(AT degrees_of_freedom1, RT probability); //  probability from t.

   //template <class AT, class RT>
   //typename tools::promote_arg2<RT, AT>::type
   //students_t_inv(AT degrees_of_freedom1, RT probability); // t from probability.

   // sqrt
   template <class T>
   T sqrtp1m1(const T&);

 	} // namespace math
} // namespace boost

#endif // BOOST_MATH_SPECIAL_MATH_FWD_HPP
