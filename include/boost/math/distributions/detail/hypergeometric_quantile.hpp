// Copyright 2008 John Maddock
//
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_DISTRIBUTIONS_DETAIL_HG_QUANTILE_HPP
#define BOOST_MATH_DISTRIBUTIONS_DETAIL_HG_QUANTILE_HPP

#include <boost/math/policies/error_handling.hpp>
#include <boost/math/distributions/detail/hypergeometric_pdf.hpp>

namespace boost{ namespace math{ namespace detail{

template <class T>
inline unsigned round_x_from_p(unsigned x, T p, T cum, T fudge_factor, const policies::discrete_quantile<policies::integer_round_down>&)
{
   if(p * fudge_factor <= cum)
      return x;
   return --x;
}

template <class T>
inline unsigned round_x_from_p(unsigned x, T p, T cum, T /*fudge_factor*/, const policies::discrete_quantile<policies::integer_round_up>&)
{
   if(cum < p)
      return ++x;
   return x;
}

template <class T>
inline unsigned round_x_from_p(unsigned x, T p, T cum, T fudge_factor, const policies::discrete_quantile<policies::integer_round_inwards>&)
{
   if((p > 0.5) && (p * fudge_factor <= cum))
      return --x;
   return x;
}

template <class T>
inline unsigned round_x_from_p(unsigned x, T p, T cum, T fudge_factor, const policies::discrete_quantile<policies::integer_round_outwards>&)
{
   if((p < 0.5) && (p * fudge_factor <= cum))
      return --x;
   return x;
}

template <class T>
inline unsigned round_x_from_q(unsigned x, T /*q*/, T /*cum*/, T /*fudge_factor*/, const policies::discrete_quantile<policies::integer_round_down>&)
{
   return x;
}

template <class T>
inline unsigned round_x_from_q(unsigned x, T q, T cum, T fudge_factor, const policies::discrete_quantile<policies::integer_round_up>&)
{
   if((q * fudge_factor <= cum) && (q < 0.5))
      return ++x;
   return x;
}

template <class T>
inline unsigned round_x_from_q(unsigned x, T q, T cum, T fudge_factor, const policies::discrete_quantile<policies::integer_round_inwards>&)
{
   if((q > 0.5) && (q * fudge_factor <= cum))
      return ++x;
   return x;
}

template <class T>
inline unsigned round_x_from_q(unsigned x, T q, T cum, T fudge_factor, const policies::discrete_quantile<policies::integer_round_outwards>&)
{
   if((q < 0.5) && (q * fudge_factor <= cum))
      return ++x;
   return x;
}

template <class T, class Policy>
unsigned hypergeometric_quantile_imp(T p, T q, unsigned r, unsigned n, unsigned N, const Policy& pol)
{
#ifdef BOOST_MSVC
#  pragma warning(push)
#  pragma warning(disable:4267)
#endif
   typedef typename Policy::discrete_quantile_type discrete_quantile_type;
   BOOST_MATH_STD_USING
   T result;
   T fudge_factor = 1 + tools::epsilon<T>() * 50;
   if(p <= 0.5)
   {
      unsigned base = static_cast<unsigned>((std::max)(0, (int)(n + r) - (int)(N)));
      unsigned x = base;
      result = hypergeometric_pdf<T>(x, r, n, N, pol);
      T diff = result;
      while(result * fudge_factor < p)
      {
         diff = (n - x) * (r - x) * diff / ((x + 1) * (N + x + 1 - n - r));
         ++x;
         result += diff;
      }
      return x == base ? x : round_x_from_p(x, p, fudge_factor, result, discrete_quantile_type());
   }
   else
   {
      unsigned lim = (std::min)(r, n);
      unsigned x = lim;
      result = 0;
      T diff = hypergeometric_pdf<T>(x, r, n, N, pol);
      while(result * fudge_factor <= q)
      {
         result += diff;
         diff = x * (N + x - n - r) * diff / ((1 + n - x) * (1 + r - x));
         --x;
      }
      return x == lim ? x : round_x_from_q(x, q, result, fudge_factor, discrete_quantile_type());
   }
#ifdef BOOST_MSVC
#  pragma warning(pop)
#endif
}

template <class T, class Policy>
inline unsigned hypergeometric_quantile(T p, T q, unsigned r, unsigned n, unsigned N, const Policy&)
{
   BOOST_FPU_EXCEPTION_GUARD
   typedef typename tools::promote_args<T>::type result_type;
   typedef typename policies::evaluation<result_type, Policy>::type value_type;
   typedef typename policies::normalise<
      Policy, 
      policies::promote_float<false>, 
      policies::promote_double<false>, 
      policies::assert_undefined<> >::type forwarding_policy;

   return detail::hypergeometric_quantile_imp<value_type>(p, q, r, n, N, forwarding_policy());
}

}}} // namespaces

#endif

