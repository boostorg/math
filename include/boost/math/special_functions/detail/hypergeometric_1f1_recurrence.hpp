
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

  template <class T>
  struct hypergeometric_1f1_recurrence_a_coefficients
  {
    typedef boost::math::tuple<T, T, T> result_type;

    hypergeometric_1f1_recurrence_a_coefficients(const T& a, const T& b, const T& z):
    a(a), b(b), z(z)
    {
    }

    result_type operator()(boost::intmax_t i) const
    {
      const T ai = a + i;

      const T an = -ai;
      const T bn = (b - (2 * ai)) - z;
      const T cn = b - ai;

      return boost::math::make_tuple(an, bn, cn);
    }

  private:
    const T a, b, z;
  };

  template <class T>
  struct hypergeometric_1f1_recurrence_b_coefficients
  {
    typedef boost::math::tuple<T, T, T> result_type;

    hypergeometric_1f1_recurrence_b_coefficients(const T& a, const T& b, const T& z):
    a(a), b(b), z(z)
    {
    }

    result_type operator()(boost::intmax_t i) const
    {
      const T bi = b + i;

      const T an = z * (bi - a);
      const T bn = bi * ((z + bi) - 1);
      const T cn = bi * (bi - 1);

      return boost::math::make_tuple(an, bn, cn);
    }

  private:
    const T a, b, z;
  };

  template <class T>
  struct hypergeometric_1f1_recurrence_a_and_b_coefficients
  {
    typedef boost::math::tuple<T, T, T> result_type;

    hypergeometric_1f1_recurrence_a_and_b_coefficients(const T& a, const T& b, const T& z):
    a(a), b(b), z(z)
    {
    }

    result_type operator()(boost::intmax_t i) const
    {
      const T ai = a + i;
      const T bi = b + i;

      const T an = ai * z;
      const T bn = bi * ((1 - bi) + z);
      const T cn = bi * (1 - bi);

      return boost::math::make_tuple(an, bn, cn);
    }

  private:
    const T a, b, z;
  };

  // forward declaration for initial values
  template <class T, class Policy>
  inline T hypergeometric_1f1_imp(const T& a, const T& b, const T& z, const Policy& pol);

  template <class T, class Policy>
  inline T hypergeometric_1f1_backward_recurrence_for_negative_a(const T& a, const T& b, const T& z, const Policy& pol)
  {
    BOOST_MATH_STD_USING // modf, frexp, fabs, pow

    boost::intmax_t integer_part = 0;
    const T bk = modf(b, &integer_part);
    T ak = modf(a, &integer_part);

    int exp_of_a = 0; frexp(a, &exp_of_a);
    int exp_of_b = 0; frexp(b, &exp_of_b);

    const bool are_fractional_parts_close_enough =
      fabs(boost::math::float_distance(ak, bk)) <= pow(2, (std::max)(exp_of_a, exp_of_b));

    if ((a < b) && (b < 0) && (are_fractional_parts_close_enough)) // TODO: has to be researched deeper
    {
      ak = b - 1;
      integer_part -= (boost::math::lltrunc(ceil(b)) - 1);
    }

    T first = detail::hypergeometric_1f1_imp(ak, b, z, pol);
    --ak;
    T second = detail::hypergeometric_1f1_imp(ak, b, z, pol);

    detail::hypergeometric_1f1_recurrence_a_coefficients<T> s(ak, b, z);

    return tools::solve_recurrence_relation_backward(s,
                                                     static_cast<unsigned int>(std::abs(integer_part)),
                                                     first,
                                                     second);
  }

  template <class T, class Policy>
  inline T hypergeometric_1f1_forward_recurrence_for_positive_a(const T& a, const T& b, const T& z, const Policy& pol)
  {
    BOOST_MATH_STD_USING // modf, fabs

    boost::intmax_t integer_part = 0;
    T ak = modf(a, &integer_part);

    T first = detail::hypergeometric_1f1_imp(ak, b, z, pol);
    ++ak;
    T second = detail::hypergeometric_1f1_imp(ak, b, z, pol);

    detail::hypergeometric_1f1_recurrence_a_coefficients<T> s(ak, b, z);

    return tools::solve_recurrence_relation_forward(s, integer_part, first, second);
  }

  template <class T, class Policy>
  inline T hypergeometric_1f1_backward_recurrence_for_negative_b(const T& a, const T& b, const T& z, const Policy& pol)
  {
    BOOST_MATH_STD_USING // modf, fabs

    boost::intmax_t integer_part = 0;
    T bk = modf(b, &integer_part);

    T first = detail::hypergeometric_1f1_imp(a, bk, z, pol);
    --bk;
    T second = detail::hypergeometric_1f1_imp(a, bk, z, pol);

    detail::hypergeometric_1f1_recurrence_b_coefficients<T> s(a, bk, z);

    return tools::solve_recurrence_relation_backward(s,
                                                     static_cast<unsigned int>(std::abs(integer_part)),
                                                     first,
                                                     second);
  }

  // this method works provided that integer part of a is the same as integer part of b
  // we don't make this check inside it.
  template <class T, class Policy>
  inline T hypergeometric_1f1_backward_recurrence_for_negative_a_and_b(const T& a, const T& b, const T& z, const Policy& pol)
  {
    BOOST_MATH_STD_USING // modf, fabs

    boost::intmax_t integer_part = 0;
    T ak = modf(a, &integer_part);
    T bk = modf(b, &integer_part);

    T first = detail::hypergeometric_1f1_imp(ak, bk, z, pol);
    --ak; --bk;
    T second = detail::hypergeometric_1f1_imp(ak, bk, z, pol);

    detail::hypergeometric_1f1_recurrence_a_and_b_coefficients<T> s(ak, bk, z);

    return tools::solve_recurrence_relation_backward(s, fabs(integer_part), first, second);
  }

  // ranges
  template <class T>
  inline bool hypergeometric_1f1_is_a_small_enough(const T& a)
  {
    return a < -10; // TODO: make dependent on precision
  }

  } } } // namespaces

#endif // BOOST_HYPERGEOMETRIC_1F1_RECURRENCE_HPP_
