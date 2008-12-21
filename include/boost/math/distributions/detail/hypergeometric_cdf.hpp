// Copyright 2008 John Maddock
//
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_DISTRIBUTIONS_DETAIL_HG_CDF_HPP
#define BOOST_MATH_DISTRIBUTIONS_DETAIL_HG_CDF_HPP

#include <boost/math/policies/error_handling.hpp>
#include <boost/math/distributions/detail/hypergeometric_pdf.hpp>

namespace boost{ namespace math{ namespace detail{

template <class T, class Policy>
T hypergeometric_cdf_imp(unsigned x, unsigned r, unsigned n, unsigned N, bool invert, const Policy& pol)
{
#ifdef BOOST_MSVC
#  pragma warning(push)
#  pragma warning(disable:4267)
#endif
   BOOST_MATH_STD_USING
   T result = 0;
   T mean = T(r * n) / N;
   if(x < mean)
   {
      result = hypergeometric_pdf<T>(x, r, n, N, pol);
      T diff = result;
      unsigned lower_limit = static_cast<unsigned>((std::max)(0, (int)(n + r) - (int)(N)));
      while(result / diff > tools::epsilon<T>())
      {
         diff = x * (N + x - n - r) * diff / ((1 + n - x) * (1 + r - x));
         if(x == lower_limit)
            break;
         --x;
         result += diff;
      }
   }
   else
   {
      invert = !invert;
      unsigned upper_limit = (std::min)(r, n);
      if(x != upper_limit)
      {
         ++x;
         result = hypergeometric_pdf<T>(x, r, n, N, pol);
         T diff = result;
         while((x <= upper_limit) && (result / diff > tools::epsilon<T>()))
         {
            diff = (n - x) * (r - x) * diff / ((x + 1) * (N + x + 1 - n - r));
            ++x;
            result += diff;
         }
      }
   }
   if(invert)
      result = 1 - result;
   return result;
#ifdef BOOST_MSVC
#  pragma warning(pop)
#endif
}

template <class T, class Policy>
inline T hypergeometric_cdf(unsigned x, unsigned r, unsigned n, unsigned N, bool invert, const Policy&)
{
   BOOST_FPU_EXCEPTION_GUARD
   typedef typename tools::promote_args<T>::type result_type;
   typedef typename policies::evaluation<result_type, Policy>::type value_type;
   typedef typename policies::normalise<
      Policy, 
      policies::promote_float<false>, 
      policies::promote_double<false>, 
      policies::discrete_quantile<>,
      policies::assert_undefined<> >::type forwarding_policy;

   value_type result;
   result = detail::hypergeometric_cdf_imp<value_type>(x, r, n, N, invert, forwarding_policy());
   return policies::checked_narrowing_cast<result_type, forwarding_policy>(result, "boost::math::hypergeometric_cdf<%1%>(%1%,%1%,%1%,%1%)");
}

}}} // namespaces

#endif

