
///////////////////////////////////////////////////////////////////////////////
//  Copyright 2014 Anton Bikineev
//  Copyright 2014 Christopher Kormanyos
//  Copyright 2014 John Maddock
//  Copyright 2014 Paul Bristow
//  Distributed under the Boost
//  Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
#ifndef BOOST_MATH_HYPERGEOMETRIC_ASYM_HPP
  #define BOOST_MATH_HYPERGEOMETRIC_ASYM_HPP

  #include <boost/math/special_functions/gamma.hpp>

  namespace boost { namespace math {

  // forward declaration of used function 2f0
  template <class T1, class T2, class T3, class Policy>
  inline typename tools::promote_args<T1, T2, T3>::type hypergeometric_2f0(T1 a1, T2 a2, T3 z, const Policy& /* pol */);

  namespace detail {

  // assumes a and b are not non-positive integers
  template <class T, class Policy>
  inline T hypergeometric_1f1_asym_positive_series(const T& a, const T& b, const T& z, const Policy& pol)
  {
    BOOST_MATH_STD_USING
    static const char* const function = "boost::math::hypergeometric_1f1_asym_positive_series<%1%>(%1%,%1%,%1%)";

    const bool is_a_integer = (a == floor(a));
    const bool is_b_integer = (b == floor(b));

    BOOST_ASSERT(!(is_a_integer && a <= 0) && !(is_b_integer && b <= 0));

    const T gamma_ratio = (b > 0 && a > 0) ?
      boost::math::tgamma_ratio(b, a, pol) :
      T(boost::math::tgamma(b, pol) / boost::math::tgamma(a, pol));
    const T prefix_a = (exp(z) * gamma_ratio) * pow(z, (a - b));

    return prefix_a * boost::math::hypergeometric_2f0(b - a, 1 - a, 1 / z, pol);
  }

  // assumes b and (b - a) are not non-positive integers
  template <class T, class Policy>
  inline T hypergeometric_1f1_asym_negative_series(const T& a, const T& b, const T& z, const Policy& pol)
  {
    BOOST_MATH_STD_USING
    static const char* const function = "boost::math::hypergeometric_1f1_asym_negative_series<%1%>(%1%,%1%,%1%)";

    const T b_minus_a = b - a;

    const bool is_a_integer = (a == floor(a));
    const bool is_b_minus_a_integer = (b_minus_a == floor(b_minus_a));

    BOOST_ASSERT(!(is_a_integer && a <= 0) && !(is_b_minus_a_integer && b_minus_a <= 0));

    const T gamma_ratio = (b > 0 && b_minus_a > 0) ?
      boost::math::tgamma_ratio(b, b_minus_a, pol) :
      boost::math::tgamma(b) / boost::math::tgamma((b_minus_a), pol);
    const T prefix_b = gamma_ratio / pow(-z, a);

    return prefix_b * boost::math::hypergeometric_2f0(a, 1 - b_minus_a, -1 / z, pol);
  }

  // experimental range
  template <class T>
  inline bool hypergeometric_1f1_asym_region(const T& a, const T& b, const T& z)
  {
    BOOST_MATH_STD_USING

    const T the_max_of_one_and_b_minus_a  ((std::max)(T(1), fabs(b - a)));
    const T the_max_of_one_and_one_minus_a((std::max)(T(1), fabs(1 - a)));

    const T the_product_of_these_maxima(the_max_of_one_and_b_minus_a * the_max_of_one_and_one_minus_a);

    const T abs_of_z(fabs(z));

    if ((abs_of_z > 100) && (the_product_of_these_maxima < (abs_of_z / 2))) // TODO: not a good way
      return true;

    return false;
  }

  } } } // namespaces

#endif // BOOST_MATH_HYPERGEOMETRIC_ASYM_HPP
