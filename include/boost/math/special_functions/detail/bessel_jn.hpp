//  Copyright (c) 2006 Xiaogang Zhang
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_BESSEL_JN_HPP
#define BOOST_MATH_BESSEL_JN_HPP

#ifdef _MSC_VER
#pragma once
#endif

#include <boost/math/special_functions/detail/bessel_j0.hpp>
#include <boost/math/special_functions/detail/bessel_j1.hpp>
#include <boost/math/special_functions/detail/bessel_jy.hpp>
#include <boost/math/special_functions/detail/bessel_jy_asym.hpp>

// Bessel function of the first kind of integer order
// J_n(z) is the minimal solution
// n < abs(z), forward recurrence stable and usable
// n >= abs(z), forward recurrence unstable, use Miller's algorithm

namespace boost { namespace math { namespace detail{

template <class T, class Policy>
struct bessel_j_small_z_series_term
{
   typedef T result_type;

   bessel_j_small_z_series_term(T v_, T x)
      : N(0), v(v_)
   {
      BOOST_MATH_STD_USING
      mult = x / 2;
      mult *= -mult;
      term = 1;
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

template <class T, class Policy>
inline T bessel_j_small_z_series(T v, T x, const Policy& pol)
{
   BOOST_MATH_STD_USING
   T prefix;
   if(v < max_factorial<T>::value)
   {
      prefix = pow(x / 2, v) / boost::math::tgamma(v+1, pol);
   }
   else
   {
      prefix = v * log(x / 2) - boost::math::lgamma(v+1, pol);
      prefix = exp(prefix);
   }
   if(0 == prefix)
      return prefix;

   bessel_j_small_z_series_term<T, Policy> s(v, x);
   boost::uintmax_t max_iter = policies::get_max_series_iterations<Policy>();
#if BOOST_WORKAROUND(__BORLANDC__, BOOST_TESTED_AT(0x582))
   T zero = 0;
   T result = boost::math::tools::sum_series(s, boost::math::policies::get_epsilon<T, Policy>(), max_iter, zero);
#else
   T result = boost::math::tools::sum_series(s, boost::math::policies::get_epsilon<T, Policy>(), max_iter);
#endif
   policies::check_series_iterations("boost::math::bessel_j_small_z_series<%1%>(%1%,%1%)", max_iter, pol);
   return prefix * result;
}

template <typename T, typename Policy>
T bessel_jn(int n, T x, const Policy& pol)
{
    T value(0), factor, current, prev, next;

    BOOST_MATH_STD_USING

    //
    // Reflection has to come first:
    //
    if (n < 0)
    {
        factor = (n & 0x1) ? -1 : 1;  // J_{-n}(z) = (-1)^n J_n(z)
        n = -n;
    }
    else
    {
        factor = 1;
    }
    //
    // Special cases:
    //
    if (n == 0)
    {
        return factor * bessel_j0(x);
    }
    if (n == 1)
    {
        return factor * bessel_j1(x);
    }

    if (x == 0)                             // n >= 2
    {
        return static_cast<T>(0);
    }

    typedef typename bessel_asymptotic_tag<T, Policy>::type tag_type;
    if(fabs(x) > asymptotic_bessel_j_limit<T>(n, tag_type()))
      return factor * asymptotic_bessel_j_large_x_2<T>(n, x);

    BOOST_ASSERT(n > 1);
    if (n < abs(x))                         // forward recurrence
    {
        prev = bessel_j0(x);
        current = bessel_j1(x);
        for (int k = 1; k < n; k++)
        {
            value = 2 * k * current / x - prev;
            prev = current;
            current = value;
        }
    }
    else if(x < 1)
    {
       return factor * bessel_j_small_z_series(T(n), x, pol);
    }
    else                                    // backward recurrence
    {
        T fn; int s;                        // fn = J_(n+1) / J_n
        // |x| <= n, fast convergence for continued fraction CF1
        boost::math::detail::CF1_jy(static_cast<T>(n), x, &fn, &s, pol);
        // tiny initial value to prevent overflow
        T init = sqrt(tools::min_value<T>());
        prev = fn * init;
        current = init;
        for (int k = n; k > 0; k--)
        {
            next = 2 * k * current / x - prev;
            prev = current;
            current = next;
        }
        T ratio = init / current;           // scaling ratio
        value = bessel_j0(x) * ratio;       // normalization
    }
    value *= factor;

    return value;
}

}}} // namespaces

#endif // BOOST_MATH_BESSEL_JN_HPP

