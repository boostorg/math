//  Copyright (c) 2013 Christopher Kormanyos
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// This work is based on an earlier work:
// "Algorithm 910: A Portable C++ Multiple-Precision System for Special-Function Calculations",
// in ACM TOMS, {VOL 37, ISSUE 4, (February 2011)} (C) ACM, 2011. http://doi.acm.org/10.1145/1916461.1916469
//
// This header contains implementation details for estimating the zeros
// of the Airy functions airy_ai and airy_bi on the negative real axis.
//
#ifndef _AIRY_AI_BI_ZERO_2013_01_20_HPP_
  #define _AIRY_AI_BI_ZERO_2013_01_20_HPP_

  #include <boost/math/constants/constants.hpp>

  namespace boost { namespace math { namespace detail {
  namespace airy_zero
  {
    template<class T>
    T equation_as_10_4_105(const T& z)
    {
      BOOST_MATH_STD_USING // ADL of std names, needed for pow.

      const T one_over_z         = T(1) / z;
      const T one_over_z_squared = one_over_z * one_over_z;

      const T z_pow_two_thirds(pow(z, T(2) / 3U));

      // Implement the top line of Eq. 10.4.105.
      const T fz(z_pow_two_thirds * (((((                     + (T(162375596875ULL) / 334430208UL)
                                         * one_over_z_squared - (   T(108056875ULL) /   6967296UL))
                                         * one_over_z_squared + (       T(77125UL)  /     82944UL))
                                         * one_over_z_squared - (           T(5U)   /        36U))
                                         * one_over_z_squared + (           T(5U)   /        48U))
                                         * one_over_z_squared + (1)));

      return fz;
    }

    namespace airy_ai_zero_detail
    {
      template<class T>
      T initial_guess(const unsigned m)
      {
        switch(m)
        {
          case 1U:  return T(-2.33810741045976703849);
          case 2U:  return T(-4.08794944413097061664);
          case 3U:  return T(-5.52055982809555105913);
          case 4U:  return T(-6.78670809007175899878);
          case 5U:  return T(-7.94413358712085312314);
          case 6U:  return T(-9.02265085334098038016);
          case 7U:  return T(-10.0401743415580859306);
          case 8U:  return T(-11.0085243037332628932);
          case 9U:  return T(-11.9360155632362625170);
          case 10U: return T(-12.8287767528657572004);
          default:
          {
            const T t = ((boost::math::constants::pi<T>() * 3U) * T((T(m) * 4U) - 1)) / 8U;

            return -boost::math::detail::airy_zero::equation_as_10_4_105(t);
          }
        }
      }
    } // namespace airy_ai_zero_detail

    namespace airy_bi_zero_detail
    {
      template<class T>
      T initial_guess(const unsigned m)
      {
        switch(m)
        {
          case 1U:  return T(-1.17371322270912792492);
          case 2U:  return T(-3.27109330283635271568);
          case 3U:  return T(-4.83073784166201593267);
          case 4U:  return T(-6.16985212831025125983);
          case 5U:  return T(-7.37676207936776371360);
          case 6U:  return T(-8.49194884650938801345);
          case 7U:  return T(-9.53819437934623888663);
          case 8U:  return T(-10.5299135067053579244);
          case 9U:  return T(-11.4769535512787794379);
          case 10U: return T(-12.3864171385827387456);
          default:
          {
            const T t = ((boost::math::constants::pi<T>() * 3U) * T((T(m) * 4U) - 3)) / 8U;

            return -boost::math::detail::airy_zero::equation_as_10_4_105(t);
          }
        }
      }
    } // namespace airy_bi_zero_detail
  } // namespace airy_zero
  } // namespace detail
  } // namespace math
  } // namespaces boost

#endif // _AIRY_AI_BI_ZERO_2013_01_20_HPP_
