//  Copyright (c) 2007 John Maddock
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// This header just defines the function entry points, and adds dispatch
// to the right implementation method.  Most of the implementation details
// are in separate headers and copyright Xiaogang Zhang.
//
#ifndef BOOST_MATH_BESSEL_HPP
#define BOOST_MATH_BESSEL_HPP

#include <boost/math/special_functions/detail/bessel_jy.hpp>
#include <boost/math/special_functions/detail/bessel_jn.hpp>
#include <boost/math/special_functions/detail/bessel_yn.hpp>
#include <boost/math/special_functions/detail/bessel_ik.hpp>
#include <boost/math/special_functions/detail/bessel_i0.hpp>
#include <boost/math/special_functions/detail/bessel_i1.hpp>
#include <boost/math/special_functions/detail/bessel_kn.hpp>
#include <boost/math/special_functions/sin_pi.hpp>
#include <boost/math/special_functions/cos_pi.hpp>
#include <boost/math/special_functions/sinc.hpp>
#include <boost/math/tools/evaluation_type.hpp>
#include <boost/math/tools/rational.hpp>
#include <boost/math/tools/promotion.hpp>

namespace boost{ namespace math{

namespace detail{

typedef mpl::int_<0> bessel_no_int_tag;      // No integer optimisation possible.
typedef mpl::int_<1> bessel_maybe_int_tag;   // Maybe integer optimisation.
typedef mpl::int_<2> bessel_int_tag;         // Definite integer optimistaion.

template <class T1, class T2>
struct bessel_traits
{
   typedef typename tools::promote_args<
      T1, T2
   >::type result_type;

   typedef typename mpl::if_c<
      (std::numeric_limits<result_type>::is_specialized == 0)
      || (std::numeric_limits<result_type>::digits == 0)
      || (std::numeric_limits<result_type>::digits > 64),
      bessel_no_int_tag,
      typename mpl::if_<
         is_integral<T1>,
         bessel_int_tag,
         bessel_maybe_int_tag
      >::type
   >::type optimisation_tag;
};

template <class T>
struct bessel_j_small_z_series_term
{
   typedef T result_type;

   bessel_j_small_z_series_term(T v_, T x)
      : N(0), v(v_)
   {
      using namespace std;
      mult = x / 2;
      term = pow(mult, v) / tgamma(v+1);
      mult *= -mult;
   }
   T operator()()
   {
      T r = term;
      ++N;
      term *= mult / (N * (N + v));
      return r;
   }
private:
   unsigned N;
   T v;
   T mult;
   T term;
};

template <class T>
struct sph_bessel_j_small_z_series_term
{
   typedef T result_type;

   sph_bessel_j_small_z_series_term(unsigned v_, T x)
      : N(0), v(v_)
   {
      using namespace std;
      mult = x / 2;
      term = pow(mult, T(v)) / tgamma(v+1+T(0.5f));
      mult *= -mult;
   }
   T operator()()
   {
      T r = term;
      ++N;
      term *= mult / (N * T(N + v + 0.5f));
      return r;
   }
private:
   unsigned N;
   unsigned v;
   T mult;
   T term;
};

template <class T>
inline T bessel_j_small_z_series(T v, T x)
{
   bessel_j_small_z_series_term<T> s(v, x);
   boost::uintmax_t max_iter = BOOST_MATH_MAX_ITER;
   T result = boost::math::tools::sum_series(s, boost::math::tools::digits<T>(), max_iter);
   tools::check_series_iterations(BOOST_CURRENT_FUNCTION, max_iter);
   return result;
}

template <class T>
inline T sph_bessel_j_small_z_series(unsigned v, T x)
{
   using namespace std; // ADL of std names
   sph_bessel_j_small_z_series_term<T> s(v, x);
   boost::uintmax_t max_iter = BOOST_MATH_MAX_ITER;
   T result = boost::math::tools::sum_series(s, boost::math::tools::digits<T>(), max_iter);
   tools::check_series_iterations(BOOST_CURRENT_FUNCTION, max_iter);
   return result * sqrt(constants::pi<T>() / 4);
}

template <class T>
T cyl_bessel_j_imp(T v, T x, const bessel_no_int_tag& t)
{
   using namespace std;
   if(x < 0)
   {
      // better have integer v:
      if(floor(v) == v)
      {
         T r = cyl_bessel_j_imp(v, -x, t);
         if(tools::real_cast<int>(v) & 1)
            r = -r;
         return r;
      }
      else
         return tools::domain_error<T>(
            BOOST_CURRENT_FUNCTION,
            "Got x = %1%, but we need x >= 0", x);
   }
   if(x == 0)
      return (v == 0) ? 1 : (v > 0) ? 0 : 
         tools::domain_error<T>(
            BOOST_CURRENT_FUNCTION, 
            "Got v = %1%, but require v >= 0 or a negative integer: the result would be complex.", v);
   
   
   if((v >= 0) && ((x < 1) || (v > x * x / 4)))
   {
      return bessel_j_small_z_series(v, x);
   }
   
   T j, y;
   bessel_jy(v, x, &j, &y, need_j);
   return j;
}

template <class T>
inline T cyl_bessel_j_imp(T v, T x, const bessel_maybe_int_tag&)
{
   using namespace std;  // ADL of std names.
   typedef typename bessel_asymptotic_tag<T>::type tag_type;
   if((fabs(v) < 200) && (floor(v) == v))
   {
      if(fabs(x) > asymptotic_bessel_j_limit<T>(v, tag_type()))
         return asymptotic_bessel_j_large_x_2(v, x);
      else
         return bessel_jn(tools::real_cast<int>(v), x);
   }
   return cyl_bessel_j_imp(v, x, bessel_no_int_tag());
}

template <class T>
inline T cyl_bessel_j_imp(int v, T x, const bessel_int_tag&)
{
   using namespace std;
   typedef typename bessel_asymptotic_tag<T>::type tag_type;
   if(fabs(x) > asymptotic_bessel_j_limit<T>(abs(v), tag_type()))
   {
      T r = asymptotic_bessel_j_large_x_2(static_cast<T>(abs(v)), x);
      if((v < 0) && (v & 1))
         r = -r;
      return r;
   }
   else
      return bessel_jn(v, x);
}

template <class T>
T sph_bessel_j_imp(unsigned n, T x)
{
   using namespace std; // ADL of std names
   if(x < 0)
      return tools::domain_error<T>(
         BOOST_CURRENT_FUNCTION,
         "Got x = %1%, but function requires x > 0.", x);
   //
   // Special case, n == 0 resolves down to the sinus cardinal of x:
   //
   if(n == 0)
      return boost::math::sinc_pi(x);
   //
   // When x is small we may end up with 0/0, use series evaluation
   // instead, especially as it converges rapidly:
   //
   if(x < 1)
      return sph_bessel_j_small_z_series(n, x);
   //
   // Default case is just a naive evaluation of the definition:
   //
   return sqrt(constants::pi<T>() / (2 * x)) 
      * cyl_bessel_j_imp(T(n)+T(0.5f), x, bessel_no_int_tag());
}

template <class T>
T cyl_bessel_i_imp(T v, T x)
{
   //
   // This handles all the bessel I functions, note that we don't optimise
   // for integer v, other than the v = 0 or 1 special cases, as Millers
   // algorithm is at least as inefficient as the general case (the general
   // case has better error handling too).
   //
   using namespace std;
   if(x < 0)
   {
      // better have integer v:
      if(floor(v) == v)
      {
         T r = cyl_bessel_i_imp(v, -x);
         if(tools::real_cast<int>(v) & 1)
            r = -r;
         return r;
      }
      else
         return tools::domain_error<T>(
            BOOST_CURRENT_FUNCTION,
            "Got x = %1%, but we need x >= 0", x);
   }
   if(x == 0)
   {
      return (v == 0) ? 1 : 0;
   }
   if(v == 0.5f)
   {
      // common special case, note try and avoid overflow in exp(x):
      T e = exp(x / 2);
      return e * (e / sqrt(2 * x * constants::pi<T>()));
   }
   if(tools::digits<T>() <= 64)
   {
      if(v == 0)
      {
         return bessel_i0(x);
      }
      if(v == 1)
      {
         return bessel_i1(x);
      }
   }
   T I, K;
   bessel_ik(v, x, &I, &K, need_i);
   return I;
}

template <class T>
T cyl_bessel_k_imp(T v, T x, const bessel_no_int_tag& t)
{
   using namespace std;
   if(x < 0)
   {
      return tools::domain_error<T>(
         BOOST_CURRENT_FUNCTION,
         "Got x = %1%, but we need x > 0", x);
   }
   if(x == 0)
   {
      return (v == 0) ? tools::overflow_error<T>(BOOST_CURRENT_FUNCTION)
         : tools::domain_error<T>(
         BOOST_CURRENT_FUNCTION,
         "Got x = %1%, but we need x > 0", x);
   }
   T I, K;
   bessel_ik(v, x, &I, &K, need_k);
   return K;
}

template <class T>
inline T cyl_bessel_k_imp(T v, T x, const bessel_maybe_int_tag&)
{
   using namespace std;
   if((floor(v) == v))
   {
      return bessel_kn(tools::real_cast<int>(v), x);
   }
   return cyl_bessel_k_imp(v, x, bessel_no_int_tag());
}

template <class T>
inline T cyl_bessel_k_imp(int v, T x, const bessel_int_tag&)
{
   return bessel_kn(v, x);
}

template <class T>
inline T cyl_neumann_imp(T v, T x, const bessel_no_int_tag&)
{
   if(x <= 0)
   {
      return (v == 0) && (x == 0) ?
         tools::overflow_error<T>(BOOST_CURRENT_FUNCTION)
         : tools::domain_error<T>(
               BOOST_CURRENT_FUNCTION,
               "Got x = %1%, but result is complex for x <= 0", x);
   }
   T j, y;
   bessel_jy(v, x, &j, &y, need_y);
   // 
   // Post evaluation check for internal overflow during evaluation,
   // can occur when x is small and v is large, in which case the result
   // is -INF:
   //
   if(!(boost::math::isfinite)(y))
      return -tools::overflow_error<T>(BOOST_CURRENT_FUNCTION);
   return y;
}

template <class T>
inline T cyl_neumann_imp(T v, T x, const bessel_maybe_int_tag&)
{
   using namespace std;
   typedef typename bessel_asymptotic_tag<T>::type tag_type;
   if(floor(v) == v)
   {
      if((fabs(x) > asymptotic_bessel_y_limit<T>(tag_type())) && (fabs(x) > 5 * abs(v)))
      {
         T r = asymptotic_bessel_y_large_x_2(static_cast<T>(abs(v)), x);
         if((v < 0) && (tools::real_cast<int>(v) & 1))
            r = -r;
         return r;
      }
      else
         return bessel_yn(tools::real_cast<int>(v), x);
   }
   return cyl_neumann_imp<T>(v, x, bessel_no_int_tag());
}

template <class T>
inline T cyl_neumann_imp(int v, T x, const bessel_int_tag&)
{
   using namespace std;
   typedef typename bessel_asymptotic_tag<T>::type tag_type;
   if((fabs(x) > asymptotic_bessel_y_limit<T>(tag_type())) && (fabs(x) > 5 * abs(v)))
   {
      T r = asymptotic_bessel_y_large_x_2(static_cast<T>(abs(v)), x);
      if((v < 0) && (v & 1))
         r = -r;
      return r;
   }
   else
      return bessel_yn(tools::real_cast<int>(v), x);
}

template <class T>
T sph_neumann_imp(unsigned v, T x)
{
   using namespace std; // ADL of std names
   //
   // Nothing much to do here but check for errors, and
   // evaluate the function's definition directly:
   //
   if(x < 0)
      return tools::domain_error<T>(
         BOOST_CURRENT_FUNCTION,
         "Got x = %1%, but function requires x > 0.", x);

   if(x < 2 * tools::min_value<T>())
      return -tools::overflow_error<T>(BOOST_CURRENT_FUNCTION);

   T result = cyl_neumann_imp(T(v)+0.5f, x, bessel_no_int_tag());
   T tx = sqrt(constants::pi<T>() / (2 * x));

   if((tx > 1) && (tools::max_value<T>() / tx < result))
      return -tools::overflow_error<T>(BOOST_CURRENT_FUNCTION);

   return result * tx;
}

} // namespace detail

template <class T1, class T2>
inline typename detail::bessel_traits<T1, T2>::result_type cyl_bessel_j(T1 v, T2 x)
{
   BOOST_FPU_EXCEPTION_GUARD
   typedef typename detail::bessel_traits<T1, T2>::result_type result_type;
   typedef typename detail::bessel_traits<T1, T2>::optimisation_tag tag_type;
   typedef typename tools::evaluation<result_type>::type value_type;
   return tools::checked_narrowing_cast<result_type>(detail::cyl_bessel_j_imp<value_type>(v, static_cast<value_type>(x), tag_type()), BOOST_CURRENT_FUNCTION);
}

template <class T>
inline typename detail::bessel_traits<T, T>::result_type sph_bessel(unsigned v, T x)
{
   BOOST_FPU_EXCEPTION_GUARD
   typedef typename detail::bessel_traits<T, T>::result_type result_type;
   typedef typename tools::evaluation<result_type>::type value_type;
   return tools::checked_narrowing_cast<result_type>(detail::sph_bessel_j_imp<value_type>(v, static_cast<value_type>(x)), BOOST_CURRENT_FUNCTION);
}

template <class T1, class T2>
inline typename detail::bessel_traits<T1, T2>::result_type cyl_bessel_i(T1 v, T2 x)
{
   BOOST_FPU_EXCEPTION_GUARD
   typedef typename detail::bessel_traits<T1, T2>::result_type result_type;
   typedef typename tools::evaluation<result_type>::type value_type;
   return tools::checked_narrowing_cast<result_type>(detail::cyl_bessel_i_imp<value_type>(v, static_cast<value_type>(x)), BOOST_CURRENT_FUNCTION);
}

template <class T1, class T2>
inline typename detail::bessel_traits<T1, T2>::result_type cyl_bessel_k(T1 v, T2 x)
{
   BOOST_FPU_EXCEPTION_GUARD
   typedef typename detail::bessel_traits<T1, T2>::result_type result_type;
   typedef typename detail::bessel_traits<T1, T2>::optimisation_tag tag_type;
   typedef typename tools::evaluation<result_type>::type value_type;
   return tools::checked_narrowing_cast<result_type>(detail::cyl_bessel_k_imp<value_type>(v, static_cast<value_type>(x), tag_type()), BOOST_CURRENT_FUNCTION);
}

template <class T1, class T2>
inline typename detail::bessel_traits<T1, T2>::result_type cyl_neumann(T1 v, T2 x)
{
   BOOST_FPU_EXCEPTION_GUARD
   typedef typename detail::bessel_traits<T1, T2>::result_type result_type;
   typedef typename detail::bessel_traits<T1, T2>::optimisation_tag tag_type;
   typedef typename tools::evaluation<result_type>::type value_type;
   return tools::checked_narrowing_cast<result_type>(detail::cyl_neumann_imp<value_type>(v, static_cast<value_type>(x), tag_type()), BOOST_CURRENT_FUNCTION);
}

template <class T>
inline typename detail::bessel_traits<T, T>::result_type sph_neumann(unsigned v, T x)
{
   BOOST_FPU_EXCEPTION_GUARD
   typedef typename detail::bessel_traits<T, T>::result_type result_type;
   typedef typename tools::evaluation<result_type>::type value_type;
   return tools::checked_narrowing_cast<result_type>(detail::sph_neumann_imp<value_type>(v, static_cast<value_type>(x)), BOOST_CURRENT_FUNCTION);
}

} // namespace math
} // namespace boost

#endif // BOOST_MATH_BESSEL_HPP
