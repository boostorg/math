
///////////////////////////////////////////////////////////////////////////////
//  Copyright 2014 Anton Bikineev
//  Copyright 2014 Christopher Kormanyos
//  Copyright 2014 John Maddock
//  Copyright 2014 Paul Bristow
//  Distributed under the Boost
//  Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_HYPERGEOMETRIC_1F1_RECURRENCE_HPP_
#define BOOST_HYPERGEOMETRIC_1F1_RECURRENCE_HPP_

#include <boost/math/special_functions/modf.hpp>
#include <boost/math/special_functions/next.hpp>

#include <boost/math/tools/recurrence.hpp>

  namespace boost { namespace math { namespace detail {

  // forward declaration for initial values
  template <class T, class Policy>
  inline T hypergeometric_1F1_imp(const T& a, const T& b, const T& z, const Policy& pol);

  template <class T, class Policy>
  inline T hypergeometric_1F1_imp(const T& a, const T& b, const T& z, const Policy& pol, int& log_scaling);

  template <class T>
  struct hypergeometric_1F1_recurrence_a_coefficients
  {
    typedef boost::math::tuple<T, T, T> result_type;

    hypergeometric_1F1_recurrence_a_coefficients(const T& a, const T& b, const T& z):
    a(a), b(b), z(z)
    {
    }

    result_type operator()(boost::intmax_t i) const
    {
      const T ai = a + i;

      const T an = b - ai;
      const T bn = (2 * ai - b + z);
      const T cn = -ai;

      return boost::math::make_tuple(an, bn, cn);
    }

  private:
    const T a, b, z;
  };

  template <class T>
  struct hypergeometric_1F1_recurrence_b_coefficients
  {
    typedef boost::math::tuple<T, T, T> result_type;

    hypergeometric_1F1_recurrence_b_coefficients(const T& a, const T& b, const T& z):
    a(a), b(b), z(z)
    {
    }

    result_type operator()(boost::intmax_t i) const
    {
      const T bi = b + i;

      const T an = bi * (bi - 1);
      const T bn = bi * (1 - bi - z);
      const T cn = z * (bi - a);

      return boost::math::make_tuple(an, bn, cn);
    }

  private:
    const T a, b, z;
  };
  //
  // for use when we're recursing to a small b:
  //
  template <class T>
  struct hypergeometric_1F1_recurrence_small_b_coefficients
  {
     typedef boost::math::tuple<T, T, T> result_type;

     hypergeometric_1F1_recurrence_small_b_coefficients(const T& a, const T& b, const T& z, int N) :
        a(a), b(b), z(z), N(N)
     {
     }

     result_type operator()(boost::intmax_t i) const
     {
        const T bi = b + (i + N);
        const T bi_minus_1 = b + (i + N - 1);

        const T an = bi * bi_minus_1;
        const T bn = bi * (-bi_minus_1 - z);
        const T cn = z * (bi - a);

        return boost::math::make_tuple(an, bn, cn);
     }

  private:
     const T a, b, z;
     int N;
  };

  template <class T>
  struct hypergeometric_1F1_recurrence_a_and_b_coefficients
  {
    typedef boost::math::tuple<T, T, T> result_type;

    hypergeometric_1F1_recurrence_a_and_b_coefficients(const T& a, const T& b, const T& z):
    a(a), b(b), z(z)
    {
    }

    result_type operator()(boost::intmax_t i) const
    {
      const T ai = a + i;
      const T bi = b + i;

      const T an = bi * (bi - 1);
      const T bn = bi * ((1 - bi) + z);
      const T cn = -ai * z;

      return boost::math::make_tuple(an, bn, cn);
    }

  private:
    const T a, b, z;
  };

  template <class T, class Policy>
  inline T hypergeometric_1F1_backward_recurrence_for_negative_a(const T& a, const T& b, const T& z, const Policy& pol, const char* function, int& log_scaling)
  {
    BOOST_MATH_STD_USING // modf, frexp, fabs, pow

    boost::intmax_t integer_part = 0;
    T ak = modf(a, &integer_part);
    //
    // We need ak-1 positive to avoid infinite recursion below:
    //
    if (0 != ak)
    {
       ak += 2;
       integer_part -= 2;
    }

    if (-integer_part > static_cast<boost::intmax_t>(policies::get_max_series_iterations<Policy>()))
       return policies::raise_evaluation_error<T>(function, "1F1 arguments sit in a range with a so negative that we have no evaluation method, got a = %1%", std::numeric_limits<T>::quiet_NaN(), pol);

    T first, second;
    if(ak == 0)
    { 
       first = 1;
       ak -= 1;
       second = 1 - z / b;
    }
    else
    {
       int scaling1(0), scaling2(0);
       first = detail::hypergeometric_1F1_imp(ak, b, z, pol, scaling1);
       ak -= 1;
       second = detail::hypergeometric_1F1_imp(ak, b, z, pol, scaling2);
       if (scaling1 != scaling2)
       {
          second *= exp(T(scaling2 - scaling1));
       }
       log_scaling += scaling1;
    }
    ++integer_part;

    detail::hypergeometric_1F1_recurrence_a_coefficients<T> s(ak, b, z);

    return tools::apply_recurrence_relation_backward(s,
                                                     static_cast<unsigned int>(std::abs(integer_part)),
                                                     first,
                                                     second, &log_scaling);
  }


  template <class T, class Policy>
  T hypergeometric_1F1_backwards_recursion_on_b_for_negative_a(const T& a, const T& b, const T& z, const Policy& pol, const char* function, int& log_scaling)
  {
     BOOST_MATH_STD_USING // modf, frexp, fabs, pow
     //
     // We compute 
     //
     // M[a + a_shift, b + b_shift; z] 
     //
     // and recurse backwards on a and b down to
     //
     // M[a, b, z]
     //
     // With a + a_shift > 1 and b + b_shift > z
     // 
     // There are 3 distinct regions to ensure stability during the recursions:
     //
     // a > 0         :  stable for backwards on a
     // a < 0, b > 0  :  stable for backwards on a and b
     // a < 0, b < 0  :  stable for backwards on b (as long as |b| is small). 
     // 
     // We could simplify things by ignoring the middle region, but it's more efficient
     // to recurse on a and b together when we can.
     //
     int b_shift = itrunc(z - b) + 2;

     int a_shift = itrunc(-a);
     if (a + a_shift != 0)
     {
        a_shift += 2;
     }

     if (b_shift > static_cast<boost::intmax_t>(boost::math::policies::get_max_series_iterations<Policy>()))
        return boost::math::policies::raise_evaluation_error<T>(function, "1F1 arguments sit in a range with z - b so large that we have no evaluation method, got z - b = %1%", z - b, pol);

     if (a_shift > static_cast<boost::intmax_t>(boost::math::policies::get_max_series_iterations<Policy>()))
        return boost::math::policies::raise_evaluation_error<T>(function, "1F1 arguments sit in a range with a so negative that we have no evaluation method, got a = %1%", a, pol);

     int a_b_shift = b < 0 ? itrunc(b + b_shift) : b_shift;   // The max we can shift on a and b together
     int leading_a_shift = 3;               // Just enough to make a negative
     if (a_b_shift > a_shift - 3)
     {
        a_b_shift = a_shift - 3;
     }
     else
     {
        // Need to ensure that leading_a_shift is large enough that a will reach it's target
        // after the first 2 phases (-,0) and (-,-) are over:
        leading_a_shift = a_shift - a_b_shift;
     }
     int trailing_b_shift = b_shift - a_b_shift;
     //
     // We could probably remove these constraints but it would greatly complexify things
     // and we simply do not use this routine in any area where these asserts will fail:
     //
     BOOST_ASSERT(leading_a_shift > 1);
     BOOST_ASSERT(a_b_shift > 1);

     T first, second;
     first = boost::math::detail::hypergeometric_1F1_imp(T(a + a_shift), T(b + b_shift), z, pol, log_scaling);
     //
     // It would be good to compute "second" from first and the ratio - unfortunately we are right on the cusp
     // recursion on a switching from stable backwards to stable forwards behaviour and so this is not possible here.
     //
     second = boost::math::detail::hypergeometric_1F1_imp(T(a + a_shift - 1), T(b + b_shift), z, pol, log_scaling);

     //
     // Now we have [a + a_shift, b + b_shift, z] and [a + a_shift - 1, b + b_shift, z]
     // and want to recurse until [a + a_shift - leading_a_shift, b + b_shift, z] and [a + a_shift - leadng_a_shift - 1, b + b_shift, z]
     // which is leading_a_shift -1 steps.
     //
     hypergeometric_1F1_recurrence_a_coefficients<T> coef_a(a + a_shift - 1, b + b_shift, z);
     second = boost::math::tools::apply_recurrence_relation_backward(coef_a, leading_a_shift, first, second, &log_scaling, &first);
     //
     // Now we need to switch to an a+b shift so that we have:
     // [a + a_shift - leading_a_shift, b + b_shift, z] and [a + a_shift - leadng_a_shift - 1, b + b_shift - 1, z]
     // A&S 13.4.3 gives us what we need:
     //
     {
        // local a's and b's:
        T la = a + a_shift - leading_a_shift - 1;
        T lb = b + b_shift;
        second = ((1 + la - lb) * second - la * first) / (1 - lb);
     }
     //
     // Now apply a_b_shift - 1 recursions to get down to
     // [a + 1, b + trailing_b_shift + 1, z] and [a, b + trailing_b_shift, z]
     //
     hypergeometric_1F1_recurrence_a_and_b_coefficients<T> coef_ab(a + a_shift - leading_a_shift - 1, b + b_shift - 1, z);
     second = boost::math::tools::apply_recurrence_relation_backward(
        coef_ab,
        a_b_shift - 1, first, second, &log_scaling, &first);
     //
     // Now we need to switch to a b shift, a different application of A&S 13.4.3
     // will get us there, we leave "second" where it is, and move "first" sideways:
     //
     {
        T lb = b + trailing_b_shift + 1;
        first = (second * (lb - 1) - a * first) / -(1 + a - lb);
     }
     //
     // Finish off by applying trailing_b_shift recursions:
     //
     if (trailing_b_shift)
     {
        hypergeometric_1F1_recurrence_small_b_coefficients<T> coef_b(a, b, z, trailing_b_shift);
        second = boost::math::tools::apply_recurrence_relation_backward(coef_b, trailing_b_shift, first, second, &log_scaling);
     }
     return second;
  }



  } } } // namespaces

#endif // BOOST_HYPERGEOMETRIC_1F1_RECURRENCE_HPP_
