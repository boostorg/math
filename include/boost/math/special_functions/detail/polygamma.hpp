
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

  namespace boost { namespace math { namespace detail{

#if 0

  template<class T, class Policy>
  T digamma_atinfinityplus(const int, const T &x, const Policy&)
  {
    BOOST_MATH_STD_USING

    T z(x);
    T log_z(log(z));
    T one_over_2z= T(1) / (2 * z);
    T sum(0);

    for(int two_k = 2; two_k < max_iteration<T>::value; two_k += 2)
    {


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
      long int tol         =  std::numeric_limits<T>::digits_base10;


      if((two_k > 24) && (order_check < -tol))
      {
        break;
      }
    }

    return (log_z - one_over_2z) - sum;
  }
#endif

  template<class T, class Policy>
  T polygamma_atinfinityplus(const int n, const T& x, const Policy& pol) // for large values of x such as for x> 400
  {
     // See http://functions.wolfram.com/GammaBetaErf/PolyGamma2/06/02/0001/
     BOOST_MATH_STD_USING
     //
     // sum       == current value of accumulated sum.
     // term      == value of current term to be added to sum.
     // part_term == value of current term excluding the Bernoulli number part
     //
     T term, sum, part_term;
     T x_squared = x * x;
     //
     // Start by setting part_term to:
     //
     // (n-1)! / x^(n+1)
     //
     // which is common to both the first term of the series (with k = 1)
     // and to the leading part.  
     // We can then get to the leading term by:
     //
     // part_term * (n + 2 * x) / 2
     //
     // and to the first term in the series 
     // (excluding the Bernoulli number) by:
     //
     // part_term n * (n + 1) / (2x)
     //
     // If either the factorial would overflow,
     // or the power term underflows, this just gets set to 0 and then we
     // know that we have to use logs for the initial terms:
     //
     part_term = ((n > boost::math::max_factorial<T>::value) && (n * n > tools::log_max_value<T>())) 
        ? 0 : boost::math::factorial<T>(n - 1, pol) * pow(x, -n - 1);
     if(part_term == 0)
     {
        // Either n is very large, or the power term underflows,
        // set the initial values of part_term, term and sum via logs:
        part_term = boost::math::lgamma(n, pol) - (n + 1) * log(x);
        sum = exp(part_term + log(n + 2 * x) - boost::math::constants::ln_two<T>());
        part_term += log(n * (n + 1)) - boost::math::constants::ln_two<T>() - log(x);
        part_term = exp(part_term);
     }
     else
     {
        sum = part_term * (n + 2 * x) / 2;
        part_term *= n * (n + 1) / 2;
        part_term /= x;
     }
     //
     // If the leading term is 0, so is the result:
     //
     if(sum == 0)
        return sum;

     for(unsigned k = 1;;)
     {
        term = part_term * boost::math::bernoulli_b2n<T>(k, pol);
        sum += term;
        //
        // Normal termination condition:
        //
        if(fabs(term / sum) < tools::epsilon<T>())
           break;
        //
        // Increment our counter, and move part_term on to the next value:
        //
        ++k;
        part_term *= (n + 2 * k - 2) * (n - 1 + 2 * k);
        part_term /= (2 * k - 1) * 2 * k;
        part_term /= x_squared;
        //
        // Emergency get out termination condition:
        //
        if(k > policies::get_max_series_iterations<Policy>())
        {
           policies::raise_evaluation_error("polygamma<%1%>(int, %1%)", "Series did not converge, closest value was %1%", sum, pol);
           break;
        }
     }
     
     if((n - 1) & 1)
        sum = -sum;

     return sum;
  }

  template<class T, class Policy>
  T polygamma_attransitionplus(const int n, const T& x, const Policy& pol)
  {
    // See: http://functions.wolfram.com/GammaBetaErf/PolyGamma2/16/01/01/0017/

    // Use N = (0.4 * digits) + (4 * n) for target value for x:
    BOOST_MATH_STD_USING
    const int d4d  = static_cast<boost::int32_t>(0.4F * policies::digits_base10<T, Policy>());
    const int N4dn = static_cast<boost::int32_t>(d4d + (4 * n));
    const int N    = static_cast<boost::int32_t>((std::min)(N4dn, (std::numeric_limits<int>::max)()));
    const int m    = n;
    const int iter = N - itrunc(x);

    const int minus_m_minus_one = -m - 1;

    T z(x);
    T sum0(0);
    T z_plus_k_pow_minus_m_minus_one(0);

    // Forward recursion to larger x:
    for(int k = 1; k <= iter; ++k)
    {
      z_plus_k_pow_minus_m_minus_one = pow(z, minus_m_minus_one);
      sum0 += z_plus_k_pow_minus_m_minus_one;
      ++z;
    }
    sum0 *= boost::math::factorial<T>(n);
    if((n - 1) & 1)
       sum0 = -sum0;

    return sum0 + polygamma_atinfinityplus(n, z, pol);
  }

  template<class T, class Policy>
  T polygamma_nearzero(const int n, const T& x, const Policy& pol)
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

    for(unsigned k = 1;; k++)
    {
      k_plus_n_fact   *= k_plus_n_plus_one;
      k_plus_n_plus_one += 1;
      one_over_k_fact /= k;
      z_pow_k         *= x;

      const T pg = k_plus_n_fact * boost::math::zeta<T>(k_plus_n_plus_one);

      const T term = (pg * z_pow_k) * one_over_k_fact;

      if(fabs(term / sum) < tools::epsilon<T>())
      {
        break;
      }

      b_neg_term = !b_neg_term;

      ((!b_neg_term) ? sum += term : sum -= term);

      if(k > policies::get_max_series_iterations<Policy>())
      {
         policies::raise_evaluation_error("polygamma<%1%>(int, %1%)", "Series did not converge, closest value was %1%", sum, pol);
         break;
      }
    }

    return term0 + sum;

  }

  template<class T, class Policy>
  inline T polygamma_imp(const int n, T x, const Policy &pol)
  {
    if(n == 0)
       return boost::math::digamma(x);
    if(n < 0)
       return policies::raise_domain_error<T>("boost::math::polygamma<%1%>(int, %1%)", "Order must be >= 0, but got %1%", n, pol);
    if(x < 0.5F)
    {
      return polygamma_nearzero(n, x, pol);
    }
    else if(x > 0.4F * policies::digits_base10<T, Policy>() + 4 * n)
    {
      return polygamma_atinfinityplus(n, x, pol);
    }
    else
    {
      return polygamma_attransitionplus(n, x, pol);
    }
  }

} } } // namespace boost::math::detail

#endif // _BOOST_POLYGAMMA_DETAIL_2013_07_30_HPP_

