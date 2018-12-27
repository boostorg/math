
///////////////////////////////////////////////////////////////////////////////
//  Copyright 2018 John Maddock
//  Distributed under the Boost
//  Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_HYPERGEOMETRIC_1F1_BY_RATIOS_HPP_
#define BOOST_HYPERGEOMETRIC_1F1_BY_RATIOS_HPP_

#include <boost/math/tools/recurrence.hpp>
#include <boost/math/policies/error_handling.hpp>

  namespace boost { namespace math { namespace detail {

     template <class T, class Policy>
     T hypergeometric_1F1_imp(const T& a, const T& b, const T& z, const Policy& pol, int& log_scaling);

     /*
      Evaluation by method of ratios for domain b < 0 < a,z

      We first convert the recurrence relation into a ratio
      of M(a+1, b+1, z) / M(a, b, z) using Shintan's equivalence
      between solving a recurrence relation using Miller's method
      and continued fractions.  The continued fraction is VERY rapid
      to converge (typically < 10 terms), but may converge to entirely
      the wrong value if we're in a bad part of the domain.  Strangely
      it seems to matter not whether we use recurrence on a, b or a and b
      they all work or not work about the same, so we might as well make
      life easy for ourselves and use the a and b recurrence to avoid
      having to apply one extra recurrence to convert from an a or b
      recurrence to an a+b one.

      See: H. SHINTAN, Note on Miller’s recurrence algorithm, J. Sci. Hiroshima Univ. Ser. A-I
      Math., 29 (1965), pp. 121-133.
      Also: Computational Aspects of Three Term Recurrence Relations, SIAM Review, January 1967.

      The following table lists by experiment, how large z can be in order to
      ensure the continued fraction converges to the correct value:

          a         b    max  z
         13,      -130,      22
         13,     -1300,     335
         13,    -13000,    3585
        130,      -130,       8
        130,     -1300,     263
        130,   - 13000,    3420
       1300,      -130,       1
       1300,     -1300,      90
       1300,    -13000,    2650
      13000,       -13,       -
      13000,      -130,       -
      13000,     -1300,      13
      13000,    -13000,     750

      So try z_limit = -b / (4 - 5 * sqrt(log(a)) * a / b);
      or     z_limit = -b / (4 - 5 * (log(a)) * a / b)  for a < 100
      
      This still isn't quite right for both a and b small, but we'll be using a Bessel approximation
      in that region anyway.

      Normalization using wronksian {1,2} from A&S 13.1.20 (also 13.1.12, 13.1.13):

      W{ M(a,b,z), z^(1-b)M(1+a-b, 2-b, z) } = (1-b)z^-b e^z

       = M(a,b,z) M2'(a,b,z) - M'(a,b,z) M2(a,b,z)
       = M(a,b,z) [(a-b+1)z^(1-b)/(2-b) M2(a+1,b+1,z) + (1-b)z^-b M2(a,b,z)] - a/b M(a+1,b+1,z) z^(1-b)M2(a,b,z)
       = M(a,b,z) [(a-b+1)z^(1-b)/(2-b) M2(a+1,b+1,z) + (1-b)z^-b M2(a,b,z)] - a/b R(a,b,z) M(a,b,z) z^(1-b)M2(a,b,z)
       = M(a,b,z) [(a-b+1)z^(1-b)/(2-b) M2(a+1,b+1,z) + (1-b)z^-b M2(a,b,z) - a/b R(a,b,z) z^(1-b)M2(a,b,z) ]
       so:
       (1-b)e^z = M(a,b,z) [(a-b+1)z/(2-b) M2(a+1,b+1,z) + (1-b) M2(a,b,z) - a/b z R(a,b,z) M2(a,b,z) ]

      */

     template <class T>
     inline bool is_in_hypergeometric_1F1_from_function_ratio_negative_b_region(const T& a, const T& b, const T& z)
     {
        BOOST_MATH_STD_USING
        if (a < 100)
           return z < -b / (4 - 5 * (log(a)) * a / b);
        else
           return z < -b / (4 - 5 * sqrt(log(a)) * a / b);
     }

     template <class T, class Policy>
     T hypergeometric_1F1_from_function_ratio_negative_b(const T& a, const T& b, const T& z, const Policy& pol, int& log_scaling)
     {
        //
        // Get the function ratio, M(a+1, b+1, z)/M(a,b,z):
        //
        boost::uintmax_t max_iter = boost::math::policies::get_max_series_iterations<Policy>();
        boost::math::detail::hypergeometric_1F1_recurrence_a_and_b_coefficients<T> coef(a + 1, b + 1, z);
        T ratio = boost::math::tools::function_ratio_from_backwards_recurrence(coef, boost::math::policies::get_epsilon<T, Policy>(), max_iter);
        boost::math::policies::check_series_iterations<T>("boost::math::hypergeometric_1F1_from_function_ratio_negative_b_positive_a<%1%>(%1%,%1%,%1%)", max_iter, pol);
        //
        // Let M2 = M(1+a-b, 2-b, z)
        // This is going to be a mighty big number:
        //
        int local_scaling = 0;
        T M2 = boost::math::detail::hypergeometric_1F1_imp(1 + a - b, 2 - b, z, pol, local_scaling);
        log_scaling -= local_scaling; // all the M2 terms are in the denominator
        //
        // Since a, b and z are all likely to be large we need the Wronksian
        // calculation below to not overflow, so scale everything right down:
        //
        if (fabs(M2) > 1)
        {
           int s = boost::math::itrunc(log(fabs(M2)));
           log_scaling -= s;  // M2 will be in the denominator, so subtract the scaling!
           M2 *= exp(T(-s));
        }
        //
        // Let M3 = M(1+a-b + 1, 2-b + 1, z)
        // we can get to this from the ratio which is cheaper to calculate:
        //
        max_iter = boost::math::policies::get_max_series_iterations<Policy>();
        boost::math::detail::hypergeometric_1F1_recurrence_a_and_b_coefficients<T> coef2(1 + a - b + 1, 2 - b + 1, z);
        T M3 = boost::math::tools::function_ratio_from_backwards_recurrence(coef2, boost::math::policies::get_epsilon<T, Policy>(), max_iter) * M2;
        boost::math::policies::check_series_iterations<T>("boost::math::hypergeometric_1F1_from_function_ratio_negative_b_positive_a<%1%>(%1%,%1%,%1%)", max_iter, pol);
        //
        // Get the RHS of the Wronksian:
        //
        int fz = itrunc(z);
        log_scaling += fz;
        T rhs = (1 - b) * exp(z - fz);
        //
        // We need to divide that by:
        // [(a-b+1)z/(2-b) M2(a+1,b+1,z) + (1-b) M2(a,b,z) - a/b z^b R(a,b,z) M2(a,b,z) ]
        // Note that at this stage, both M3 and M2 are scaled by exp(local_scaling).
        //
        T lhs = (a - b + 1) * z * M3 / (2 - b) + (1 - b) * M2 - a * z * ratio * M2 / b;

        return rhs / lhs;
     }
     //
     // This next version is largely the same as above, but calculates the ratio for the b recurrence relation
     // which has a larger area of stability than the ab recurrence when a,b < 0.  We can then use a single
     // recurrence step to convert this to the ratio for the ab recursion and proceed largly as before.
     // The routine is quite insensitive to the size of z, but requires |a| < |5b| for accuracy.
     // Fortunately the accuracy outside this domain falls off steadily rather than suddenly switching
     // to a different behaviour.
     //
     template <class T, class Policy>
     T hypergeometric_1F1_from_function_ratio_negative_ab(const T& a, const T& b, const T& z, const Policy& pol, int& log_scaling)
     {
        //
        // Get the function ratio, M(a+1, b+1, z)/M(a,b,z):
        //
        boost::uintmax_t max_iter = boost::math::policies::get_max_series_iterations<Policy>();
        boost::math::detail::hypergeometric_1F1_recurrence_b_coefficients<T> coef(a, b + 1, z);
        T ratio = boost::math::tools::function_ratio_from_backwards_recurrence(coef, boost::math::policies::get_epsilon<T, Policy>(), max_iter);
        boost::math::policies::check_series_iterations<T>("boost::math::hypergeometric_1F1_from_function_ratio_negative_b_positive_a<%1%>(%1%,%1%,%1%)", max_iter, pol);
        //
        // We need to use A&S 13.4.3 to convert a ratio for M(a, b + 1, z) / M(a, b, z)
        // to M(a+1, b+1, z) / M(a, b, z)
        //
        // We have:        (a-b)M(a, b+1, z) - aM(a+1, b+1, z) + bM(a, b, z) = 0
        // and therefore:  M(a + 1, b + 1, z) / M(a, b, z) = ((a - b)M(a, b + 1, z) / M(a, b, z) + b) / a
        //
        ratio = ((a - b) * ratio + b) / a;
        //
        // Let M2 = M(1+a-b, 2-b, z)
        // This is going to be a mighty big number:
        //
        int local_scaling = 0;
        T M2 = boost::math::detail::hypergeometric_1F1_imp(1 + a - b, 2 - b, z, pol, local_scaling);
        log_scaling -= local_scaling; // all the M2 terms are in the denominator
        //
        // Let M3 = M(1+a-b + 1, 2-b + 1, z)
        // We don't use the ratio to get this as it's not clear that it's reliable:
        //
        int local_scaling2 = 0;
        T M3 = boost::math::detail::hypergeometric_1F1_imp(2 + a - b, 3 - b, z, pol, local_scaling2);
        //
        // M2 and M3 must be identically scaled:
        //
        if (local_scaling != local_scaling2)
        {
           M3 *= exp(local_scaling2 - local_scaling);
        }
        //
        // Get the RHS of the Wronksian:
        //
        int fz = itrunc(z);
        log_scaling += fz;
        T rhs = (1 - b) * exp(z - fz);
        //
        // We need to divide that by:
        // [(a-b+1)z/(2-b) M2(a+1,b+1,z) + (1-b) M2(a,b,z) - a/b z^b R(a,b,z) M2(a,b,z) ]
        // Note that at this stage, both M3 and M2 are scaled by exp(local_scaling).
        //
        T lhs = (a - b + 1) * z * M3 / (2 - b) + (1 - b) * M2 - a * z * ratio * M2 / b;

        return rhs / lhs;
     }


  } } } // namespaces

#endif // BOOST_HYPERGEOMETRIC_1F1_BY_RATIOS_HPP_
