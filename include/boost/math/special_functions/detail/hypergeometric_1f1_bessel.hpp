
///////////////////////////////////////////////////////////////////////////////
//  Copyright 2014 Anton Bikineev
//  Copyright 2014 Christopher Kormanyos
//  Copyright 2014 John Maddock
//  Copyright 2014 Paul Bristow
//  Distributed under the Boost
//  Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
#ifndef BOOST_MATH_HYPERGEOMETRIC_1F1_BESSEL_HPP
  #define BOOST_MATH_HYPERGEOMETRIC_1F1_BESSEL_HPP

  #include <boost/math/tools/series.hpp>
  #include <boost/math/special_functions/bessel.hpp>
  #include <boost/math/special_functions/laguerre.hpp>

  namespace boost { namespace math { namespace detail {

  // declarations of helpers for 13_3_7 and 13_3_8:
  // bessel_j recurrence relation based on finite continued fractions;
  // it is stable forward recurrence algorithm, see William J. Lentz
  // "Generating Bessel functions in Mie scattering calculations
  // using continued fractions" for details
  template <class T>
  inline T hypergeometric_bessel_j_recurrence_next(const T& jvm1, const T& v, const T& z);

  // swap values of Bessel functions:
  template <class T>
  inline void hypergeometric_bessel_j_recurrence_iterate(T& jvm1, T& jv, const T& v, const T& z);

  // next coefficient for 13_3_7
  template <class T>
  inline T hypergeometric_13_3_7_coefficient_next(const T& anm3, const T& anm2, const T& a, const T& b, const unsigned n);

  // next coefficient for 13_3_8
  template <class T>
  inline T hypergeometric_13_3_8_coefficient_next(const T& cnm3, const T& cnm2, const T& cnm1, const T& a, const T& b, const unsigned n);

  // swap values of coefficients:
  template <class T>
  inline void hypergeometric_coefficient_13_3_7_iterate(T& anm3, T& anm2, T& anm1, T& an, const T& a, const T& b, const unsigned n);

  template <class T>
  inline void hypergeometric_coefficient_13_3_8_iterate(T& cnm3, T& cnm2, T& cnm1, T& cn, const T& a, const T& b, const unsigned n);

  // term class of Abramowitz & Stegun 13_3_7 formula
  template <class T>
  struct hypergeometric_1f1_13_3_7_series_term
  {
    typedef T result_type;

    hypergeometric_1f1_13_3_7_series_term(const T& a, const T& b, const T& z):
      a(a), b(b), z(z), n(0u)
    {
      BOOST_MATH_STD_USING

      sqrt_z_pow_n = sqrt_z = sqrt(z);
      sqrt_2b_minus_4a_pow_n = sqrt_2b_minus_4a = sqrt(2 * (b - (2 * a)));
      sqrt_2zb_minus_4za = sqrt_z * sqrt_2b_minus_4a;

      v_current = b - 1;

      anm3 = 1;
      anm2 = 0;
      anm1 = b / 2;
      an = detail::hypergeometric_13_3_7_coefficient_next(anm3, anm2, a, b, 3u);

      jvm1 = boost::math::cyl_bessel_j(v_current, sqrt_2zb_minus_4za);
      v_current++;
      jv = detail::hypergeometric_bessel_j_recurrence_next(jvm1, v_current, sqrt_2zb_minus_4za);

      term = jvm1;
      ++n;
    }

    T operator()()
    {
      const T result = term;

      iterate();
      term = ((anm2 * sqrt_z_pow_n) / sqrt_2b_minus_4a_pow_n) * jv;

      return result;
    }

  private:
    void iterate()
    {
      ++n; ++v_current;
      sqrt_z_pow_n *= sqrt_z;
      sqrt_2b_minus_4a_pow_n *= sqrt_2b_minus_4a;

      detail::hypergeometric_coefficient_13_3_7_iterate(anm3, anm2, anm1, an, a, b, n + 2);
      detail::hypergeometric_bessel_j_recurrence_iterate(jvm1, jv, v_current, sqrt_2zb_minus_4za);
    }

    const T a, b, z;
    unsigned n;
    T sqrt_z, sqrt_2b_minus_4a, sqrt_2zb_minus_4za;
    T sqrt_z_pow_n, sqrt_2b_minus_4a_pow_n;
    T v_current;
    T anm3, anm2, anm1, an;
    T jvm1, jv;
    T term;
  };

  // function for 13_3_7 evaluation
  template <class T, class Policy>
  inline T hypergeometric_1f1_13_3_7_series(const T& a, const T& b, const T& z, const Policy& pol)
  {
    BOOST_MATH_STD_USING

    const T sqrt_bz_div_2_minus_az = sqrt(((b * z) / 2) - (a * z));
    const T prefix = ((boost::math::tgamma(b, pol) * sqrt_bz_div_2_minus_az) /
        pow(sqrt_bz_div_2_minus_az, b)) * exp(z / 2);

    detail::hypergeometric_1f1_13_3_7_series_term<T> s(a, b, z);
    boost::uintmax_t max_iter = boost::math::policies::get_max_series_iterations<Policy>();
#if BOOST_WORKAROUND(__BORLANDC__, BOOST_TESTED_AT(0x582))
    T zero = 0;
    T result = boost::math::tools::sum_series(s, boost::math::policies::get_epsilon<T, Policy>(), max_iter, zero);
#else
    T result = boost::math::tools::sum_series(s, boost::math::policies::get_epsilon<T, Policy>(), max_iter);
#endif
    boost::math::policies::check_series_iterations<T>("boost::math::hypergeometric_1f1_13_3_7_series<%1%>(%1%,%1%,%1%)", max_iter, pol);
    return prefix * result;
  }

  // term class of Abramowitz & Stegun 13_3_8 formula
  template <class T>
  struct hypergeometric_1f1_13_3_8_series_term
  {
    typedef T result_type;

    static const T h;

    hypergeometric_1f1_13_3_8_series_term(const T& a, const T& b, const T& z):
      a(a), b(b), z(z), n(0u)
    {
      BOOST_MATH_STD_USING

      sqrt_minus_az = sqrt(-a * z);
      const T double_sqrt_minus_az = 2 * sqrt_minus_az;

      v_current = b - 1;
      z_pow_n = z; sqrt_minus_az_pow_n = sqrt_minus_az;

      cnm3 = 1;
      cnm2 = -b * h;
      cnm1 = (((1 - (2 * h)) * a) + ((b * (b + 1)) * (h * h))) / 2;
      cn = detail::hypergeometric_13_3_8_coefficient_next(cnm3, cnm2, cnm1, a, b, 3u);

      jvm1 = boost::math::cyl_bessel_j(v_current, double_sqrt_minus_az);
      ++v_current;
      jv = detail::hypergeometric_bessel_j_recurrence_next(jvm1, b, double_sqrt_minus_az);
    }

    T operator()()
    {
      ++n;
      switch (n - 1)
      {
        case 0u:
          return jvm1;
        case 1u:
          detail::hypergeometric_bessel_j_recurrence_iterate(jvm1, jv, v_current, T(2 * sqrt_minus_az));
          return ((cnm2 * z) / sqrt_minus_az) * jvm1;
      }

      ++v_current;
      z_pow_n *= z;
      sqrt_minus_az_pow_n *= sqrt_minus_az;
      const T result = ((cnm1 * z_pow_n) / sqrt_minus_az_pow_n) * jv;

      detail::hypergeometric_coefficient_13_3_8_iterate(cnm3, cnm2, cnm1, cn, a, b, n + 1);
      detail::hypergeometric_bessel_j_recurrence_iterate(jvm1, jv, v_current, T(2 * sqrt_minus_az));
      return result;
    }

  private:
    const T a, b, z;
    unsigned n;
    T sqrt_minus_az;
    T z_pow_n, sqrt_minus_az_pow_n;
    T v_current;
    T cnm3, cnm2, cnm1, cn;
    T jvm1, jv;
  };

  template <class T>
  const T hypergeometric_1f1_13_3_8_series_term<T>::h = -boost::math::constants::pi<T>() / T(10.);

  // function for 13_3_8 evaluation
  template <class T, class Policy>
  inline T hypergeometric_1f1_13_3_8_series(const T& a, const T& b, const T& z, const Policy& pol)
  {
    BOOST_MATH_STD_USING

    const T sqrt_minus_az_pow_b_minus_one = pow(sqrt(-a * z), b - 1);
    const T prefix = (boost::math::tgamma(b, pol) / sqrt_minus_az_pow_b_minus_one) *
                      exp(detail::hypergeometric_1f1_13_3_8_series_term<T>::h * z);

    detail::hypergeometric_1f1_13_3_8_series_term<T> s(a, b, z);
    boost::uintmax_t max_iter = boost::math::policies::get_max_series_iterations<Policy>();
#if BOOST_WORKAROUND(__BORLANDC__, BOOST_TESTED_AT(0x582))
    T zero = 0;
    T result = boost::math::tools::sum_series(s, boost::math::policies::get_epsilon<T, Policy>(), max_iter, zero);
#else
    T result = boost::math::tools::sum_series(s, boost::math::policies::get_epsilon<T, Policy>(), max_iter);
#endif
    boost::math::policies::check_series_iterations<T>("boost::math::hypergeometric_1f1_13_3_8_series<%1%>(%1%,%1%,%1%)", max_iter, pol);
    return prefix * result;
  }

  // definitions of helpers:

  template <class T>
  inline T hypergeometric_bessel_j_recurrence_next(const T& jvm1, const T& v, const T& z)
  {
    BOOST_MATH_STD_USING

    // TODO: we can estimate a possible capacity for these vectors
    std::vector<T> finite_cf_nums;
    std::vector<T> finite_cf_denoms;

    const T a_coef_when_i_is_one = (2 * v) / z;
    const T a_coef_when_i_is_two = (-2 * (v + 1)) / z;

    finite_cf_nums.push_back(a_coef_when_i_is_one);
    finite_cf_nums.push_back(a_coef_when_i_is_two + (1. / a_coef_when_i_is_one));

    finite_cf_denoms.push_back(a_coef_when_i_is_two);

    unsigned i = 3u;

    do {
      const T temp = (2 * (v + (i - 1))) / z;
      const T next_coef = i & 1 ? temp : T(-temp);

      finite_cf_nums.push_back(next_coef + (1. / finite_cf_nums.back()));
      finite_cf_denoms.push_back(next_coef + (1. / finite_cf_denoms.back()));

      ++i;
    } while (fabs(finite_cf_nums.back() - finite_cf_denoms.back()) > boost::math::tools::epsilon<T>());

    T ratio = finite_cf_nums.front();

    for (i = 0u; i < finite_cf_denoms.size(); ++i)
      ratio *= finite_cf_nums[i+1] / finite_cf_denoms[i];

    return jvm1 / ratio;
  }

  template <class T>
  inline void hypergeometric_bessel_j_recurrence_iterate(T& jvm1, T& jv, const T& v, const T& z)
  {
    using std::swap;

    swap(jvm1, jv);
    jv = detail::hypergeometric_bessel_j_recurrence_next(jvm1, v, z);
  }

  template <class T>
  inline T hypergeometric_13_3_7_coefficient_next(const T& anm3, const T& anm2, const T& a, const T& b, const unsigned n)
  {
    const T term_anm2 = ((n + b) - 2) * anm2;
    const T term_anm3 = ((2 * a) - b) * anm3;

    return (term_anm2 + term_anm3) / n;

  }

  template <class T>
  inline T hypergeometric_13_3_8_coefficient_next(const T& cnm3, const T& cnm2, const T& cnm1, const T& a, const T& b, const unsigned n)
  {
    static const T& h = detail::hypergeometric_1f1_13_3_8_series_term<T>::h;
    const T one_minus_two_h = 1 - (2 * h);
    const T h_minus_one = h - 1;

    const T term_cnm1  = ((one_minus_two_h * n) - (b * h)) * cnm1;
    const T term_cnm2 = ((one_minus_two_h * a) - ((h * h_minus_one) * (b + (n - 1)))) * cnm2;
    const T term_cnm3 = ((-h * h_minus_one) * a) * cnm3;

    return ((term_cnm1 + term_cnm2) + term_cnm3) / (n + 1);
  }

  template <class T>
  inline void hypergeometric_coefficient_13_3_7_iterate(T& anm3, T& anm2, T& anm1, T& an, const T& a, const T& b, const unsigned n)
  {
    using std::swap;

    swap(anm3, anm2);
    swap(anm2, anm1);
    swap(anm1, an);
    an = detail::hypergeometric_13_3_7_coefficient_next(anm3, anm2, a, b, n);
  }

  template <class T>
  inline void hypergeometric_coefficient_13_3_8_iterate(T& cnm3, T& cnm2, T& cnm1, T& cn, const T& a, const T& b, const unsigned n)
  {
    using std::swap;

    swap(cnm3, cnm2);
    swap(cnm2, cnm1);
    swap(cnm1, cn);
    cn = detail::hypergeometric_13_3_8_coefficient_next(cnm3, cnm2, cnm1, a, b, n);
  }

  } } } // namespaces

#endif // BOOST_MATH_HYPERGEOMETRIC_1F1_BESSEL_HPP
