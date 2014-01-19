
///////////////////////////////////////////////////////////////////////////////
//  Copyright 2013 Nikhar Agrawal
//  Copyright 2013 Christopher Kormanyos
//  Copyright 2013 John Maddock
//  Copyright 2013 Paul Bristow
//  Distributed under the Boost
//  Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef _BOOST_POLYGAMMA_DETAIL_2013_07_30_HPP_
  #define _BOOST_POLYGAMMA_DETAIL_2013_07_30_HPP_

  #include <cmath>
  #include <limits>
  #include <boost/cstdint.hpp>
  #include <boost/math/policies/policy.hpp>
  #include <boost/math/special_functions/bernoulli.hpp>
  #include <boost/math/special_functions/trunc.hpp>
  #include <boost/math/special_functions/zeta.hpp>
  #include <boost/mpl/if.hpp>
  #include <boost/mpl/int.hpp>
  #include <boost/static_assert.hpp>
  #include <boost/type_traits/is_convertible.hpp>

  namespace boost { namespace math { namespace detail {

  template<class T>
  struct max_iteration
  {
    //TODO Derive a suitable formula based on the precision of T
    static const int value=2500;
  };

  template<class T>
  bool factorial_overflow(const int n)
  {
    // Use Stirling's approximation to check if n! would overflow when data type is T.
    static const long int max_precision = std::numeric_limits<T>::max_exponent10;

    T nn                   = T(n);
    T log_n                = log(nn);
    T n_log_n              = n * log_n;
    T n_log_n_minus_n      = n_log_n - n;
    T base_10              = n_log_n_minus_n/log(10);
    long int base_10_ceil  = boost::math::ltrunc(base_10) + 1;

    // Since nlogn - n < log(n!) by a small margin, we add 10 as safety measure.
    return (((base_10_ceil + 10) > max_precision )? 1 : 0);
  }

  template<class T>
  int possible_factorial_overflow_index()
  {
    // we use binary search to determine a good approximation for an index that might overflow

    int upper_limit = max_iteration<T>::value;
    int lower_limit = 8;

    if(factorial_overflow<T>(upper_limit) == 0)
    {
      return upper_limit;
    }

    while(upper_limit > (lower_limit + 4))
    {
      const int mid = (upper_limit + lower_limit) / 2;

      if(factorial_overflow<T>(mid) == 0)
      {
        lower_limit = mid;
      }
      else
      {
        upper_limit = mid;
      }
    }

    return lower_limit;
  }

  template<class T, class Policy>
  T digamma_atinfinityplus(const int, const T &x, const Policy&)
  {
    BOOST_MATH_STD_USING

    // calculate a high bernoulli number upfront to make use of cache
    unsigned int bernoulli_index = 100;
    boost::math::bernoulli_b2n<T>(bernoulli_index);

    T z(x);
    T log_z(log(z));
    T one_over_2z= T(1) / (2 * z);
    T sum(0);

    for(int two_k = 2; two_k < max_iteration<T>::value; two_k += 2)
    {
      if(two_k/2 > static_cast<boost::int32_t>(bernoulli_index))
      {
        try
        {
          int temp = static_cast<int>(bernoulli_index * 1.5F);
          boost::math::bernoulli_b2n<T>(temp);
          bernoulli_index = temp;
        }
        catch(...)
        {
          break;
        }
      }

      T term(1);
      T one_over_two_k       = T(1) / two_k;
      T z_pow_two_k          = pow(z, static_cast<boost::int32_t>(two_k));
      T one_over_z_pow_two_k = T(1) / z_pow_two_k;
      T bernoulli_term       = boost::math::bernoulli_b2n<T>(two_k / 2);

      term = (bernoulli_term * one_over_two_k) * one_over_z_pow_two_k;

      if(term == 0)
      {
        continue;
      }

      sum += term;

      T term_base_10_exp = ((term < 0) ? -term : term);
      T sum_base_10_exp  = ((sum  < 0) ? -sum  : sum);

      int exponent_value;

      static_cast<void>(frexp(term_base_10_exp, &exponent_value));
      term_base_10_exp = T(exponent_value) * 0.303F;

      static_cast<void>(frexp(sum_base_10_exp, &exponent_value));
      sum_base_10_exp = T(exponent_value) * 0.303F;

      long int order_check =  boost::math::ltrunc(term_base_10_exp) - boost::math::ltrunc(sum_base_10_exp);
      long int tol         =  std::numeric_limits<T>::digits10;


      if((two_k > 24) && (order_check < -tol))
      {
        break;
      }
    }

    return (log_z - one_over_2z) - sum;
  }

  template<class T, class Policy>
  T polygamma_atinfinityplus(const int n, const T& x, const Policy& pol) // for large values of x such as for x> 400
  {
     BOOST_MATH_STD_USING

     if(n == 0)
     {
       return digamma_atinfinityplus(n, x, pol);
    }

     //TODO try calculating for bernoulli_index= max_iteration, if error then set bernoulli_index=100
     unsigned int bernoulli_index = 100;
     boost::math::bernoulli_b2n<T>(bernoulli_index);

     const bool b_negate = ((n % 2) == 0);

     const T n_minus_one_fact            = boost::math::factorial<T>(n - 1);
     const T nn                          = T(n);
     const T n_fact                      = n_minus_one_fact * nn;
     const T one_over_z                  = T(1) / x;
     const T one_over_z2                 = one_over_z * one_over_z;
     const T one_over_z_pow_n            = T(1) / pow(x, n);
           T one_over_x_pow_two_k_plus_n = one_over_z_pow_n * one_over_z2;
           T two_k_plus_n_minus_one      = nn + T(1);
           T two_k_plus_n_minus_one_fact = n_fact * (n + 1); //(n+3)! ?
           T one_over_two_k_fact         = T(1) / 2;
           T sum                         = (  (boost::math::bernoulli_b2n<T>(1) * two_k_plus_n_minus_one_fact)
                                            * (one_over_two_k_fact * one_over_x_pow_two_k_plus_n));

     // Perform the Bernoulli series expansion.
     for(int two_k = 4; two_k < max_iteration<T>::value; two_k += 2)
     {
       if((two_k / 2) > static_cast<int>(bernoulli_index))
       {
         try
         {
           //TODO the multiplication factor should depend upon T, small precision, smaller multiplication factor
           int temp = static_cast<int>(bernoulli_index * 2.0F);
           boost::math::bernoulli_b2n<T>(temp);
           bernoulli_index = temp;
         }
         catch(...)
         {
           break;
         }
       }

       one_over_x_pow_two_k_plus_n *= one_over_z2;
       two_k_plus_n_minus_one_fact *= ++two_k_plus_n_minus_one;
       two_k_plus_n_minus_one_fact *= ++two_k_plus_n_minus_one;
       one_over_two_k_fact         /= static_cast<boost::int32_t>(two_k * static_cast<boost::int32_t>(two_k - static_cast<boost::int32_t>(1)));

       const T term = (  (boost::math::bernoulli_b2n<T>(two_k/2) * two_k_plus_n_minus_one_fact)
                       * (one_over_two_k_fact * one_over_x_pow_two_k_plus_n));

        if(term == 0)
        {
          continue;
        }

        sum += term;

        T term_base_10_exp = ((term < 0) ? -term : term);
        T sum_base_10_exp  = ((sum  < 0) ? -sum  : sum);

        int exponent_value;

        static_cast<void>(frexp(term_base_10_exp, &exponent_value));
        term_base_10_exp = T(exponent_value) * 0.303F;

        static_cast<void>(frexp(sum_base_10_exp, &exponent_value));
        sum_base_10_exp = T(exponent_value) * 0.303F;

        long int order_check =  boost::math::ltrunc(term_base_10_exp) - boost::math::ltrunc(sum_base_10_exp);
        long int tol         =  std::numeric_limits<T>::digits10;

        if((two_k > 24) && (order_check < -tol))
        {
          break;
        }
     }

     sum += ((((n_minus_one_fact * (nn + (x * static_cast<boost::int32_t>(2)))) * one_over_z_pow_n) * one_over_z) / 2);

     return ((!b_negate) ? sum : -sum);
  }

  template<class T, class Policy>
  T polygamma_attransitionplus(const int n, const T& x, const Policy&)
  {
    // this doesn't work for digamma either

    // Use Euler-Maclaurin summation.

    // Use N = (0.4 * digits) + (4 * n)
    BOOST_MATH_STD_USING
    const int d4d  = static_cast<boost::int32_t>(0.4F * std::numeric_limits<T>::digits10);
    const int N4dn = static_cast<boost::int32_t>(d4d + (4 * n));
    const int N    = static_cast<boost::int32_t>((std::min)(N4dn, (std::numeric_limits<int>::max)()));
    const int m    = n;

    const int minus_m_minus_one = -m - 1;

    T z(x);
    T sum0(0);
    T z_plus_k_pow_minus_m_minus_one(0);

    for(int k = 1; k <= N; ++k)
    {
      z_plus_k_pow_minus_m_minus_one = pow(z, minus_m_minus_one);
      sum0 += z_plus_k_pow_minus_m_minus_one;
      ++z;
    }

    const T one_over_z_plus_N_pow_minus_m           = pow(z, -m);
    const T one_over_z_plus_N_pow_minus_m_minus_one = one_over_z_plus_N_pow_minus_m / z;

    const T term0 = one_over_z_plus_N_pow_minus_m_minus_one / 2;
    const T term1 = one_over_z_plus_N_pow_minus_m           / m;

          T   sum1                                      = T(0);
          T   one_over_two_k_fact                       = T(1) / 2;
          int mk                                        = m + 1;
          T   am                                        = T(mk);
    const T   one_over_z_plus_N_squared                 = T(1) / (z * z);
          T   one_over_z_plus_N_pow_minus_m_minus_two_k = one_over_z_plus_N_pow_minus_m * one_over_z_plus_N_squared;

    for(int k = 1; k < max_iteration<T>::value; ++k)
    {
      const int two_k = 2 * k; // k << 1

      const T term = ((boost::math::bernoulli_b2n<T>(two_k / 2) * am) * one_over_two_k_fact) * one_over_z_plus_N_pow_minus_m_minus_two_k;

      T term_base_10_exp = ((term < 0) ? -term : term);
      T sum_base_10_exp  = ((sum1 < 0) ? -sum1 : sum1);

      int exponent_value;

      static_cast<void>(frexp(term_base_10_exp, &exponent_value));
      term_base_10_exp = T(exponent_value) * 0.303F;

      static_cast<void>(frexp(sum_base_10_exp, &exponent_value));
      sum_base_10_exp = T(exponent_value) * 0.303F;

      long int order_check =  boost::math::ltrunc(term_base_10_exp) - boost::math::ltrunc(sum_base_10_exp);
      long int tol         =  std::numeric_limits<T>::digits10;

      if((two_k > 24) && (order_check < -tol))
      {
        break;
      }

      sum1 += term;

      one_over_two_k_fact /= (two_k + 1);
      one_over_two_k_fact /= (two_k + 2);

      ++mk;
      am *= mk;
      ++mk;
      am *= mk;

      one_over_z_plus_N_pow_minus_m_minus_two_k *= one_over_z_plus_N_squared;
    }

    const T pg = (((sum0 + term0) + term1) + sum1) * factorial<T>(m);

    const bool b_negate = ((m % 2) == 0);

    return ((!b_negate) ? pg : -pg);
  }

  template<class T, class Policy>
  T polygamma_nearzero(const int n, const T& x, const Policy&)
  {
    BOOST_MATH_STD_USING
    // not defined for digamma

    // Use a series expansion for x near zero which uses poly_gamma(m, 1) which,
    // in turn, uses the Riemann zeta function for integer arguments.
    // http://functions.wolfram.com/GammaBetaErf/PolyGamma2/06/01/03/01/02/
    const bool b_negate = (( n % 2 ) == 0 ) ;

    const T n_fact               =  boost::math::factorial<T>(n);
    const T z_pow_n_plus_one     =  pow(x, static_cast<boost::int64_t>(n + 1));
    const T n_fact_over_pow_term =  n_fact / z_pow_n_plus_one;
    const T term0                =  !b_negate ? n_fact_over_pow_term : -n_fact_over_pow_term;

          T one_over_k_fact      =  T(1);
          T z_pow_k              =  T(1);
          T k_plus_n_fact        =  boost::math::factorial<T>(n);
          T k_plus_n_plus_one    =  T(n + 1);
    const T pg_kn                =  k_plus_n_fact * boost::math::zeta<T>(k_plus_n_plus_one);
          bool    b_neg_term     =  ((n % 2) == 0);
          T sum                  =  ((!b_neg_term) ? pg_kn : -pg_kn);

    for(int k = 1; k < max_iteration<T>::value; k++)
    {
      k_plus_n_fact   *= k_plus_n_plus_one++;
      one_over_k_fact /= k;
      z_pow_k         *= x;

      const T pg = k_plus_n_fact * boost::math::zeta<T>(k_plus_n_plus_one);

      const T term = (pg * z_pow_k) * one_over_k_fact;

      T term_base_10_exp = ((term < 0) ? -term : term);
      T sum_base_10_exp  = ((sum  < 0) ? -sum  : sum);

      int exponent_value;

      static_cast<void>(frexp(term_base_10_exp, &exponent_value));
      term_base_10_exp = T(exponent_value) * 0.303F;

      static_cast<void>(frexp(sum_base_10_exp, &exponent_value));
      sum_base_10_exp = T(exponent_value) * 0.303F;

      long int order_check =  boost::math::ltrunc(term_base_10_exp) - boost::math::ltrunc(sum_base_10_exp);
      long int tol         =  std::numeric_limits<T>::digits10;

      if((k > 12) && (order_check < -tol))
      {
        break;
      }

      b_neg_term = !b_neg_term;

      ((!b_neg_term) ? sum += term : sum -= term);
    }

    return term0 + sum;

  }

  template<class T, class Policy>
  inline T polygamma_imp(const int n, T x, const Policy &pol)
  {
    if(x < 0.5F)
    {
      return polygamma_nearzero(n, x, pol);
    }
    else if (x > 400.0F)
    {
      return polygamma_atinfinityplus(n, x, pol); //just a test return value
    }
    else
    {
      return polygamma_attransitionplus(n, x, pol);
    }
  }

} } } // namespace boost::math::detail

#endif // _BOOST_POLYGAMMA_DETAIL_2013_07_30_HPP_

