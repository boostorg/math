
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
  #include <boost/math/special_functions/digamma.hpp>
  #include <boost/math/special_functions/sin_pi.hpp>
  #include <boost/math/special_functions/cos_pi.hpp>
  #include <boost/math/special_functions/pow.hpp>
  #include <boost/mpl/if.hpp>
  #include <boost/mpl/int.hpp>
  #include <boost/static_assert.hpp>
  #include <boost/type_traits/is_convertible.hpp>

  namespace boost { namespace math { namespace detail{

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
        ? T(0) : boost::math::factorial<T>(n - 1, pol) * pow(x, -n - 1);
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
      z += 1;
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
    const T z_pow_n_plus_one     =  pow(x, n + 1);
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

  template <class T, class Policy>
  T poly_cot_pi(int n, T x, const Policy& pol)
  {
     // Return n'th derivative of cot(pi*x) at x, these are simply
     // tabulated for up to n = 9, beyond that it is possible to
     // calculate coefficients as follows:
     //
     // The general form of each derivative is:
     //
     // pi^n * SUM{k=0, n} C[k,n] * cos(k * pi * x) * csc^(n+1)(pi * x)
     //
     // With constant C[0,1] = -1 and all other C[k,n] = 0;
     // Then for each k < n+1:
     // C[|1 - k|, n+1]  += -(k + n + 2) * C[k, n] / 2;
     // C[k + 1, n+1]    += -(n + 2 - k) * C[k, n] / 2;
     //
     // It's worth noting however, that as well as requiring quite a bit
     // of storage space, this method has no better accuracy than recursion
     // on x to x > 0 when computing polygamma :-(
     //
     T s = boost::math::sin_pi(x);
     switch(n)
     {
     case 1:
        return -constants::pi<T, Policy>() / (s * s);
     case 2:
     {
        T c = boost::math::cos_pi(x);
        return 2 * constants::pi<T, Policy>() * constants::pi<T, Policy>() * c / boost::math::pow<3>(s);
     }
     case 3:
     {
        T c = boost::math::cos_pi(2 * x);
        return -2 * boost::math::pow<3>(constants::pi<T, Policy>()) * (c + 2) / boost::math::pow<4>(s);
     }
     case 4:
     {
        T c = boost::math::cos_pi(x);
        T c2 = boost::math::cos_pi(2 * x);
        return 4 * boost::math::pow<4>(constants::pi<T, Policy>()) * (c2 + 5) * c / boost::math::pow<5>(s);
     }
     case 5:
     {
        T c2 = boost::math::cos_pi(2 * x);
        T c4 = boost::math::cos_pi(4 * x);
        return -2 * boost::math::pow<5>(constants::pi<T, Policy>()) *(26 * c2 + c4 + 33) / boost::math::pow<6>(s);
     }
     case 6:
     {
        T c = boost::math::cos_pi(x);
        T c2 = boost::math::cos_pi(2 * x);
        T c4 = boost::math::cos_pi(4 * x);
        return 4 * boost::math::pow<6>(constants::pi<T, Policy>()) * (56 * c2 + c4 + 123) * c / boost::math::pow<7>(s);
     }
     case 7:
     {
        T c2 = boost::math::cos_pi(2 * x);
        T c4 = boost::math::cos_pi(4 * x);
        T c6 = boost::math::cos_pi(6 * x);
        return -2 * boost::math::pow<7>(constants::pi<T, Policy>()) * (1191 * c2 + 120 * c4 + c6 + 1208) / boost::math::pow<8>(s);
     }
     case 8:
     {
        T c = boost::math::cos_pi(x);
        T c3 = boost::math::cos_pi(3 * x);
        T c5 = boost::math::cos_pi(5 * x);
        T c7 = boost::math::cos_pi(7 * x);
        return 2 * boost::math::pow<8>(constants::pi<T, Policy>()) * (15619 * c + 4293 * c3 + 247 * c5 + c7) / boost::math::pow<9>(s);
     }
     case 9:
     {
        T c2 = boost::math::cos_pi(2 * x);
        T c4 = boost::math::cos_pi(4 * x);
        T c6 = boost::math::cos_pi(6 * x);
        T c8 = boost::math::cos_pi(8 * x);
        return -2 * boost::math::pow<9>(constants::pi<T, Policy>()) * (88234 * c2 + 14608 * c4 + 502 * c6 + c8 + 78095) / boost::math::pow<10>(s);
     }
     }
     return policies::raise_domain_error<T>("boost::math::polygamma<%1%>(int, %1%)", "Derivative of cotangent not implemented at n = %1%", n, pol);
  }

  template<class T, class Policy>
  inline T polygamma_imp(const int n, T x, const Policy &pol)
  {
    BOOST_MATH_STD_USING
    static const char* function = "boost::math::polygamma<%1%>(int, %1%)";
    if(n == 0)
       return boost::math::digamma(x);
    if(n < 0)
       return policies::raise_domain_error<T>(function, "Order must be >= 0, but got %1%", n, pol);
    if(x < 0)
    {
       if(floor(x) == x)
       {
          //
          // Result is infinity if x is odd, and a pole error if x is even.
          //
          if(lltrunc(x) & 1)
             return policies::raise_overflow_error<T>(function, 0, pol);
          else
             return policies::raise_pole_error<T>(function, "Evaluation at negative integer %1%", x, pol);
       }
       if(n < 10)
       {
          //
          // We have tabulated the derivatives of cot(x) up to the 9th derivative, which
          // allows us to use: http://functions.wolfram.com/06.15.16.0001.01
          T z = 1 - x;
          T result = polygamma_imp(n, z, pol) + constants::pi<T, Policy>() * poly_cot_pi(n, z, pol);
          return n & 1 ? -result : result;
       }
       //
       // Try http://functions.wolfram.com/06.15.16.0007.01
       //
       if(x <= -static_cast<int>(policies::get_max_series_iterations<Policy>()))
          return policies::raise_evaluation_error<T>(function, "Argument is outside the bounds for which we can reasonably calculate polygamma (got x = %1%)", x, pol);
       int m = boost::math::itrunc(ceil(-x));
       T z = x + m;
       T sum = 0;
       for(int k = 1; k <= m; ++k)
       {
          sum += pow(z - k, -n - 1);
       }
       sum *= boost::math::factorial<T>(n);
       if(n & 1)
          sum = -sum;
       return polygamma_imp(n, z, pol) - sum;
    }

    if(x < 0.125F)
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

