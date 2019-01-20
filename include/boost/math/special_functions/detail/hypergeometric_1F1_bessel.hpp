
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

     template <class T, class Policy>
     T hypergeometric_1F1_divergent_fallback(const T& a, const T& b, const T& z, const Policy& pol, int& log_scaling);

     template <class T>
     bool hypergeometric_1F1_is_tricomi_viable_positive_b(const T& a, const T& b, const T& z)
     {
        BOOST_MATH_STD_USING
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
        BOOST_MATH_STD_USING
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
              log_scale = itrunc(term);
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
           b_minus_1_plus_n += 1;
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
           b_minus_1_plus_n += 1;
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
        if (use_logs || (prefix == 0) || !(boost::math::isfinite)(prefix))
        {
           use_logs = true;
           prefix = boost::math::lgamma(b, &prefix_sgn, pol) + z / 2;
           scale = itrunc(prefix);
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


     template <class T, class Policy>
     struct hypergeometric_1F1_AS_13_3_6_series
     {
        //
        // TODO: store and cache Bessel function evaluations via backwards recurrence.
        //
        // Only convergent for a ~ b ?
        //
        hypergeometric_1F1_AS_13_3_6_series(const T& a, const T& b, const T& z, const T& b_minus_a, const Policy& pol_)
           : b_minus_a_minus_half(b_minus_a - 0.5f), half_z(z / 2), poch_1(2 * b_minus_a - 1), poch_2(b_minus_a - a), b_poch(b), term(1), last_result(1), sign(1), n(0), pol(pol_)
        {}
        T operator()()
        {
           BOOST_MATH_STD_USING
           T result = term * (b_minus_a_minus_half + n) * sign * boost::math::cyl_bessel_i(b_minus_a_minus_half + n, half_z, pol);
           ++n;
           term *= poch_1;
           poch_1 += 1;
           term *= poch_2;
           poch_2 += 1;
           term /= n;
           term /= b_poch;
           b_poch += 1;
           sign = -sign;
           //if (n % 50)
           //   std::cout << std::setw(5) << std::right << n << std::setw(30) << std::right << result << std::endl;
           if ((fabs(result) > fabs(last_result)) && (n > 100))
              return 0;  // series has gone divergent!
           last_result = result;
           return result;
        }
        T b_minus_a_minus_half, half_z, poch_1, poch_2, b_poch, term, last_result;
        int sign;
        int n;
        const Policy& pol;

        typedef T result_type;
     };

     template <class T, class Policy>
     T hypergeometric_1F1_AS_13_3_6(const T& a, const T& b, const T& z, const T& b_minus_a, const Policy& pol, int& log_scaling)
     {
        BOOST_MATH_STD_USING
        if(b_minus_a == 0)
        {
           // special case: M(a,a,z) == exp(z)
           int scale = itrunc(z, pol);
           log_scaling += scale;
           return exp(z - scale);
        }
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

     /****************************************************************************************************************/
     //
     // The following are not used at present:
     //
     /****************************************************************************************************************/

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
           bessel_arg(2 * sqrt(-a * z)), bessel_order(b - 1), power_term(std::pow(-a * z, (1 - b) / 2)), mult(z / std::sqrt(-a * z)),
           a_(a), b_(b), z_(z), h_(h), n(2), pol(pol_)
        {
        }
        T operator()()
        {
           // we actually return the n-2 term:
           T result = C_minus_2 * power_term * boost::math::cyl_bessel_j(bessel_order, bessel_arg, pol);
           bessel_order += 1;
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
        BOOST_MATH_STD_USING
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
           two_a_minus_1_plus_s += 1;
           a_minus_half_plus_s += 1;
           two_a_minus_b_plus_s += 1;
           b_plus_s += 1;

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

  } } } // namespaces

#endif // BOOST_MATH_HYPERGEOMETRIC_1F1_BESSEL_HPP
