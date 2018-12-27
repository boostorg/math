
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
#include <boost/math/special_functions/detail/hypergeometric_pFq_checked_series.hpp>


  namespace boost { namespace math { namespace detail {

#if 1

     template <class T, class Policy>
     T hypergeometric_1F1_divergent_fallback(const T& a, const T& b, const T& z, const Policy& pol, int& log_scaling);

     template <class T>
     bool hypergeometric_1F1_is_tricomi_viable_positive_b(const T& a, const T& b, const T& z)
     {
        if (b <= 100)
           return true;
        // Even though we're in a reasonable domain for Tricomi's approximation, 
        // the arguments to the Bessel functions may be so large that we can't
        // actually evaluate them:
        T x = sqrt(fabs(2 * z * b - 4 * a * z));
        T v = b - 1;
        return log(boost::math::constants::e<T>() * x / (2 * v)) * v > tools::log_min_value<T>();
     }

     template <class T, class Policy>
     T tricomi_Ev(const T& v, const T& z, const Policy& pol)
     {
        T result;
        if (z == 0)
           result = 1 / boost::math::tgamma(v + 1);
        else if (z > 0)
           result = boost::math::cyl_bessel_j(v, 2 * sqrt(z), pol);
        else
           result = boost::math::cyl_bessel_i(v, 2 * sqrt(-z), pol);
        return result;
     }

     template <class T, class Policy>
     struct hypergeometric_1F1_AS_13_3_7_tricomi_series
     {
        //
        // TODO: store and cache Bessel function evaluations via backwards recurrence.
        //
        hypergeometric_1F1_AS_13_3_7_tricomi_series(const T& a, const T& b, const T& z, const Policy& pol_)
           : A_minus_2(1), A_minus_1(0), A(b / 2), mult(z / 2), term(1), b_minus_1_plus_n(b - 1),
            bessel_arg((b / 2 - a) * z),
           two_a_minus_b(2 * a - b), pol(pol_), n(2)
        {
           BOOST_MATH_STD_USING
           term /= pow(fabs(bessel_arg), b_minus_1_plus_n / 2);
           mult /= sqrt(fabs(bessel_arg));
           if ((term == 0) || !(boost::math::isfinite)(term))
           {
              term = -log(fabs(bessel_arg)) * b_minus_1_plus_n / 2;
              log_scale = boost::math::itrunc(term);
              term -= log_scale;
              term = exp(term);
           }
           else
              log_scale = 0;
        }
        T operator()()
        {
           //
           // We return the n-2 term, and do 2 terms at once as every other term can be
           // very small (or zero) when b == 2a:
           //
           BOOST_MATH_STD_USING
           T result = A_minus_2 * term * tricomi_Ev(b_minus_1_plus_n, bessel_arg, pol);
           term *= mult;
           ++n;
           T A_next = ((b_minus_1_plus_n + 2) * A_minus_1 + two_a_minus_b * A_minus_2) / n;
           ++b_minus_1_plus_n;
           A_minus_2 = A_minus_1;
           A_minus_1 = A;
           A = A_next;

           if (A_minus_2 != 0)
           {
              result += A_minus_2 * term * tricomi_Ev(b_minus_1_plus_n, bessel_arg, pol);
           }
           term *= mult;
           ++n;
           A_next = ((b_minus_1_plus_n + 2) * A_minus_1 + two_a_minus_b * A_minus_2) / n;
           ++b_minus_1_plus_n;
           A_minus_2 = A_minus_1;
           A_minus_1 = A;
           A = A_next;

           return result;
        }
        T A_minus_2, A_minus_1, A, mult, term, b_minus_1_plus_n, bessel_arg, two_a_minus_b;
        const Policy& pol;
        int n, log_scale;

        typedef T result_type;
     };

     template <class T, class Policy>
     T hypergeometric_1F1_AS_13_3_7_tricomi(const T& a, const T& b, const T& z, const Policy& pol, int& log_scale)
     {
        BOOST_MATH_STD_USING
        //
        // Works for a < 0, b < 0, z > 0.
        //
        // For convergence we require A * term to be converging otherwise we get
        // a divergent alternating series.  It's actually really hard to analyse this
        // and the best purely heuristic policy we've found is
        // z < fabs((2 * a - b) / (sqrt(fabs(a)))) ; b > 0  or:
        // z < fabs((2 * a - b) / (sqrt(fabs(ab)))) ; b < 0
        //
        T prefix(0);
        int prefix_sgn(0);
        bool use_logs = false;
        int scale = 0;
        try
        {
           prefix = boost::math::tgamma(b, pol);
           prefix *= exp(z / 2);
        }
        catch (const std::runtime_error&)
        {
           use_logs = true;
        }
        if (use_logs || (prefix == 0) || !boost::math::isfinite(prefix))
        {
           use_logs = true;
           prefix = boost::math::lgamma(b, &prefix_sgn, pol) + z / 2;
           scale = boost::math::itrunc(prefix);
           log_scale += scale;
           prefix -= scale;
        }
        hypergeometric_1F1_AS_13_3_7_tricomi_series<T, Policy> s(a, b, z, pol);
        log_scale += s.log_scale;
        boost::uintmax_t max_iter = boost::math::policies::get_max_series_iterations<Policy>();
        bool retry = false;
        T result(0);
        try
        {
           result = boost::math::tools::sum_series(s, boost::math::policies::get_epsilon<T, Policy>(), max_iter);
           if (!(boost::math::isfinite)(result) || (result == 0))
              retry = true;
        }
        catch (const std::overflow_error&)
        {
           retry = true;
        }
        catch (const boost::math::evaluation_error&)
        {
           retry = true;
        }
        if (retry)
        {
           log_scale -= scale;
           log_scale -= s.log_scale;
           return hypergeometric_1F1_divergent_fallback(a, b, z, pol, log_scale);
        }
        boost::math::policies::check_series_iterations<T>("boost::math::hypergeometric_1F1_AS_13_3_7<%1%>(%1%,%1%,%1%)", max_iter, pol);
        if (use_logs)
        {
           int sgn = boost::math::sign(result);
           prefix += log(fabs(result));
           prefix += log_scale;
           log_scale = 0;
           result = sgn * prefix_sgn * exp(prefix);
        }
        else
        {
           result *= prefix;
        }
        return result;
     }


     /****************************************************************************************************************/
     //
     // The following are not used at present:
     //
     /****************************************************************************************************************/

     template <class T, class Policy>
     struct hypergeometric_1F1_AS_13_3_6_series
     {
        //
        // TODO: store and cache Bessel function evaluations via backwards recurrence.
        //
        // Only convergent for a ~ b ?
        //
        hypergeometric_1F1_AS_13_3_6_series(const T& a, const T& b, const T& z, const T& b_minus_a, const Policy& pol_)
           : b_minus_a_minus_half(b_minus_a - 0.5f), half_z(z / 2), poch_1(2 * b_minus_a - 1), poch_2(b_minus_a - a), b_poch(b), term(1), sign(1), n(0), pol(pol_)
        {}
        T operator()()
        {
           T result = term * (b_minus_a_minus_half + n) * sign * boost::math::cyl_bessel_i(b_minus_a_minus_half + n, half_z, pol);
           ++n;
           term *= poch_1++;
           term *= poch_2++;
           term /= n;
           term /= b_poch++;
           sign = -sign;
           return result;
        }
        T b_minus_a_minus_half, half_z, poch_1, poch_2, b_poch, term;
        int sign;
        int n;
        const Policy& pol;

        typedef T result_type;
     };

     template <class T, class Policy>
     T hypergeometric_1F1_AS_13_3_6(const T& a, const T& b, const T& z, const T& b_minus_a, const Policy& pol, int& log_scaling)
     {
        BOOST_MATH_STD_USING
        hypergeometric_1F1_AS_13_3_6_series<T, Policy> s(a, b, z, b_minus_a, pol);
        boost::uintmax_t max_iter = boost::math::policies::get_max_series_iterations<Policy>();
        T result = boost::math::tools::sum_series(s, boost::math::policies::get_epsilon<T, Policy>(), max_iter);
        boost::math::policies::check_series_iterations<T>("boost::math::hypergeometric_1F1_AS_13_3_6<%1%>(%1%,%1%,%1%)", max_iter, pol);
        result *= boost::math::tgamma(b_minus_a - 0.5f) * pow(z / 4, -b_minus_a + 0.5f);
        int scale = itrunc(z / 2);
        log_scaling += scale;
        result *= exp(z / 2 - scale);
        return result;
     }

     template <class T, class Policy>
     struct hypergeometric_1F1_AS_13_3_8_series
     {
        //
        // TODO: store and cache Bessel function evaluations via backwards recurrence.
        //
        // The C term grows by at least an order of magnitude with each iteration, and
        // rate of growth is largely independent of the arguments.  Free parameter h
        // seems to give accurate results for small values (almost zero) or h=1.
        // Convergence and accuracy, only when -a/z > 100, this appears to have no
        // or little benefit over 13.3.7 as it generally requires more iterations?
        //
        hypergeometric_1F1_AS_13_3_8_series(const T& a, const T& b, const T& z, const T& h, const Policy& pol_)
           : C_minus_2(1), C_minus_1(-b * h), C(b * (b + 1) * h * h / 2 - (2 * h - 1) * a / 2),
           bessel_arg(2 * sqrt(-a * z)), bessel_order(b - 1), power_term(pow(-a * z, (1 - b) / 2)), mult(z / sqrt(-a * z)),
           a_(a), b_(b), z_(z), h_(h), n(2), pol(pol_)
        {
        }
        T operator()()
        {
           // we actually return the n-2 term:
           T result = C_minus_2 * power_term * boost::math::cyl_bessel_j(bessel_order, bessel_arg, pol);
           ++bessel_order;
           power_term *= mult;
           ++n;
           T C_next = ((1 - 2 * h_) * (n - 1) - b_ * h_) * C
              + ((1 - 2 * h_) * a_ - h_ * (h_ - 1) *(b_ + n - 2)) * C_minus_1
              - h_ * (h_ - 1) * a_ * C_minus_2;
           C_next /= n;
           C_minus_2 = C_minus_1;
           C_minus_1 = C;
           C = C_next;
           return result;
        }
        T C, C_minus_1, C_minus_2, bessel_arg, bessel_order, power_term, mult, a_, b_, z_, h_;
        const Policy& pol;
        int n;

        typedef T result_type;
     };

     template <class T, class Policy>
     T hypergeometric_1F1_AS_13_3_8(const T& a, const T& b, const T& z, const T& h, const Policy& pol)
     {
        T prefix = exp(h * z) * boost::math::tgamma(b);
        hypergeometric_1F1_AS_13_3_8_series<T, Policy> s(a, b, z, h, pol);
        boost::uintmax_t max_iter = boost::math::policies::get_max_series_iterations<Policy>();
        T result = boost::math::tools::sum_series(s, boost::math::policies::get_epsilon<T, Policy>(), max_iter);
        boost::math::policies::check_series_iterations<T>("boost::math::hypergeometric_1F1_AS_13_3_8<%1%>(%1%,%1%,%1%)", max_iter, pol);
        result *= prefix;
        return result;
     }

     //
     // This is the series from https://dlmf.nist.gov/13.11
     // It appears to be unusable for a,z < 0, and for
     // b < 0 appears to never be better than the defining series
     // for 1F1.
     //
     template <class T, class Policy>
     struct hypergeometric_1f1_13_11_1_series
     {
        typedef T result_type;

        hypergeometric_1f1_13_11_1_series(const T& a, const T& b, const T& z, const Policy& pol_)
           : term(1), two_a_minus_1_plus_s(2 * a - 1), two_a_minus_b_plus_s(2 * a - b), b_plus_s(b), a_minus_half_plus_s(a - 0.5f), half_z(z / 2), s(0), pol(pol_)
        {
        }
        T operator()()
        {
           T result = term * a_minus_half_plus_s * boost::math::cyl_bessel_i(a_minus_half_plus_s, half_z, pol);

           term *= two_a_minus_1_plus_s * two_a_minus_b_plus_s / (b_plus_s * ++s);
           ++two_a_minus_1_plus_s;
           ++a_minus_half_plus_s;
           ++two_a_minus_b_plus_s;
           ++b_plus_s;

           return result;
        }
        T term, two_a_minus_1_plus_s, two_a_minus_b_plus_s, b_plus_s, a_minus_half_plus_s, half_z;
        int s;
        const Policy& pol;
     };

     template <class T, class Policy>
     T hypergeometric_1f1_13_11_1(const T& a, const T& b, const T& z, const Policy& pol, int& log_scaling)
     {
        BOOST_MATH_STD_USING
           bool use_logs = false;
        T prefix;
        int prefix_sgn = 1;
        if (true/*(a < boost::math::max_factorial<T>::value) && (a > 0)*/)
           prefix = boost::math::tgamma(a - 0.5f, pol);
        else
        {
           prefix = boost::math::lgamma(a - 0.5f, &prefix_sgn, pol);
           use_logs = true;
        }
        if (use_logs)
        {
           prefix += z / 2;
           prefix += log(z / 4) * (0.5f - a);
        }
        else if (z > 0)
        {
           prefix *= pow(z / 4, 0.5f - a);
           prefix *= exp(z / 2);
        }
        else
        {
           prefix *= exp(z / 2);
           prefix *= pow(z / 4, 0.5f - a);
        }

        hypergeometric_1f1_13_11_1_series<T, Policy> s(a, b, z, pol);
        boost::uintmax_t max_iter = boost::math::policies::get_max_series_iterations<Policy>();
        T result = boost::math::tools::sum_series(s, boost::math::policies::get_epsilon<T, Policy>(), max_iter);
        boost::math::policies::check_series_iterations<T>("boost::math::hypergeometric_1f1_13_11_1<%1%>(%1%,%1%,%1%)", max_iter, pol);
        if (use_logs)
        {
           int scaling = itrunc(prefix);
           log_scaling += scaling;
           prefix -= scaling;
           result *= exp(prefix) * prefix_sgn;
        }
        else
           result *= prefix;

        return result;
     }

#else

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
  struct hypergeometric_1F1_13_3_7_series_term
  {
    typedef T result_type;

    hypergeometric_1F1_13_3_7_series_term(const T& a, const T& b, const T& z):
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
      ++v_current;
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
  inline T hypergeometric_1F1_13_3_7_series(const T& a, const T& b, const T& z, const Policy& pol, const char* function, int& log_scaling)
  {
    BOOST_MATH_STD_USING

    const T sqrt_bz_div_2_minus_az = sqrt(((b * z) / 2) - (a * z));

    bool use_logs = false;
    if ((b > max_factorial<T>::value) && (b * b > tools::log_max_value<T>()))
       use_logs = true;
    if(z / 2 > tools::log_max_value<T>())
       use_logs = true;

    T prefix(0);
    int sign = 1;
    if (!use_logs)
    {
       try
       {
          prefix = boost::math::tgamma(b, pol);
          if (!(boost::math::isfinite)(prefix))
             use_logs = true;
       }
       catch (std::overflow_error const&)
       {
          use_logs = true;
       }
    }
    if (!use_logs)
    {
       prefix *= (sqrt_bz_div_2_minus_az / pow(sqrt_bz_div_2_minus_az, b)) * exp(z / 2);
    }
    else
    {
       prefix = boost::math::lgamma(b, &sign, pol) + (b - 1) * log(sqrt_bz_div_2_minus_az);
       prefix += z / 2;
    }
    T result(0);
    boost::uintmax_t max_iter = boost::math::policies::get_max_series_iterations<Policy>();
    try {
       detail::hypergeometric_1F1_13_3_7_series_term<T> s(a, b, z);
#if BOOST_WORKAROUND(__BORLANDC__, BOOST_TESTED_AT(0x582))
       T zero = 0;
       result = boost::math::tools::sum_series(s, boost::math::policies::get_epsilon<T, Policy>(), max_iter, zero);
#else
       result = boost::math::tools::sum_series(s, boost::math::policies::get_epsilon<T, Policy>(), max_iter);
#endif
    }
    catch (const std::exception&)
    {
       // The bessel functions can overflow, use a checked series in that case:
       return hypergeometric_1F1_checked_series_impl(a, b, z, pol, log_scaling);
    }
    if (!(boost::math::isfinite)(result))
    {
       // as above:
       return hypergeometric_1F1_checked_series_impl(a, b, z, pol, log_scaling);
    }
    boost::math::policies::check_series_iterations<T>("boost::math::hypergeometric_1F1_13_3_7_series<%1%>(%1%,%1%,%1%)", max_iter, pol);

    if (use_logs)
    {
       if (result < 0)
       {
          result = -result;
          sign = -sign;
       }
       result = log(result) + prefix;
       if (result > tools::log_max_value<T>())
          return sign * policies::raise_overflow_error<T>(function, 0, pol);
       result = sign * exp(result);
    }
    else
       result *= prefix;

    return result;
  }

  // term class of Abramowitz & Stegun 13_3_8 formula
  template <class T>
  struct hypergeometric_1F1_13_3_8_series_term
  {
    typedef T result_type;

    static const T h;

    hypergeometric_1F1_13_3_8_series_term(const T& a, const T& b, const T& z):
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
  const T hypergeometric_1F1_13_3_8_series_term<T>::h = -boost::math::constants::pi<T>() / T(10.);

  // function for 13_3_8 evaluation
  template <class T, class Policy>
  inline T hypergeometric_1F1_13_3_8_series(const T& a, const T& b, const T& z, const Policy& pol)
  {
    BOOST_MATH_STD_USING

    const T sqrt_minus_az_pow_b_minus_one = pow(sqrt(-a * z), b - 1);
    const T prefix = (boost::math::tgamma(b, pol) / sqrt_minus_az_pow_b_minus_one) *
                      exp(detail::hypergeometric_1F1_13_3_8_series_term<T>::h * z);

    detail::hypergeometric_1F1_13_3_8_series_term<T> s(a, b, z);
    boost::uintmax_t max_iter = boost::math::policies::get_max_series_iterations<Policy>();
#if BOOST_WORKAROUND(__BORLANDC__, BOOST_TESTED_AT(0x582))
    T zero = 0;
    T result = boost::math::tools::sum_series(s, boost::math::policies::get_epsilon<T, Policy>(), max_iter, zero);
#else
    T result = boost::math::tools::sum_series(s, boost::math::policies::get_epsilon<T, Policy>(), max_iter);
#endif
    boost::math::policies::check_series_iterations<T>("boost::math::hypergeometric_1F1_13_3_8_series<%1%>(%1%,%1%,%1%)", max_iter, pol);
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
    static const T& h = detail::hypergeometric_1F1_13_3_8_series_term<T>::h;
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

#endif

  } } } // namespaces

#endif // BOOST_MATH_HYPERGEOMETRIC_1F1_BESSEL_HPP
