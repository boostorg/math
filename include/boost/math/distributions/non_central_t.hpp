// boost\math\distributions\non_central_t.hpp

// Copyright John Maddock 2008.

// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_SPECIAL_NON_CENTRAL_T_HPP
#define BOOST_MATH_SPECIAL_NON_CENTRAL_T_HPP

#include <boost/math/distributions/fwd.hpp>
#include <boost/math/distributions/non_central_beta.hpp> // for nc beta
#include <boost/math/distributions/normal.hpp> // for normal CDF and quantile
#include <boost/math/distributions/detail/generic_quantile.hpp> // quantile

namespace boost
{
   namespace math
   {

      template <class RealType, class Policy>
      class non_central_t_distribution;

      namespace detail{

         template <class T, class Policy>
         T non_central_t2_p(T n, T delta, T x, T y, const Policy& pol, T init_val)
         {
            BOOST_MATH_STD_USING
            //
            // Variables come first:
            //
            boost::uintmax_t max_iter = policies::get_max_series_iterations<Policy>();
            T errtol = ldexp(1.0f, -boost::math::policies::digits<T, Policy>());
            T d2 = delta * delta / 2;
            //
            // k is the starting point for iteration, and is the
            // maximum of the poisson weighting term:
            //
            int k = boost::math::itrunc(d2);
            // Starting Poisson weight:
            T pois = gamma_p_derivative(T(k+1), d2, pol) 
               * tgamma_delta_ratio(T(k + 1), T(0.5f))
               * delta / constants::root_two<T>();
            // Starting beta term:
            T beta = x < y
               ? ibeta(T(k + 1), n / 2, x, pol)
               : ibetac(n / 2, T(k + 1), y, pol);
            // Recurance term:
            T xterm = x < y
               ? ibeta_derivative(T(k + 1), n / 2, x, pol)
               : ibeta_derivative(n / 2, T(k + 1), y, pol);
            xterm *= y / (n / 2 + k);
            T poisf(pois), betaf(beta), xtermf(xterm);
            T sum = init_val;

            //
            // Backwards recursion first, this is the stable
            // direction for recursion:
            //
            int count = 0;
            for(int i = k; i >= 0; --i)
            {
               T term = beta * pois;
               sum += term;
               if(fabs(term/sum) < errtol)
                  break;
               pois *= (i + 0.5f) / d2;
               beta += xterm;
               xterm *= (i) / (x * (n / 2 + i - 1));
               ++count;
            }
            for(int i = k + 1; i - k < max_iter; ++i)
            {
               poisf *= d2 / (i + 0.5f);
               xtermf *= (x * (n / 2 + i - 1)) / (i);
               betaf -= xtermf;
               T term = poisf * betaf;
               sum += term;
               if(fabs(term/sum) < errtol)
                  break;
               ++count;
            }
            return sum;
         }

         template <class T, class Policy>
         T non_central_t2_q(T n, T delta, T x, T y, const Policy& pol, T init_val)
         {
            BOOST_MATH_STD_USING
            //
            // Variables come first:
            //
            boost::uintmax_t max_iter = policies::get_max_series_iterations<Policy>();
            T errtol = ldexp(1.0f, -boost::math::policies::digits<T, Policy>());
            T d2 = delta * delta / 2;
            //
            // k is the starting point for iteration, and is the
            // maximum of the poisson weighting term:
            //
            int k = boost::math::itrunc(d2);
            // Starting Poisson weight:
            T pois = gamma_p_derivative(T(k+1), d2, pol) 
               * tgamma_delta_ratio(T(k + 1), T(0.5f))
               * delta / constants::root_two<T>();
            // Starting beta term:
            T beta = x < y 
               ? ibetac(T(k + 1), n / 2, x, pol)
               : ibeta(n / 2, T(k + 1), y, pol);
            // Recurance term:
            T xterm = x < y
               ? ibeta_derivative(T(k + 1), n / 2, x, pol)
               : ibeta_derivative(n / 2, T(k + 1), y, pol);
            xterm *= y / (n / 2 + k);
            T poisf(pois), betaf(beta), xtermf(xterm);
            T sum = init_val;

            //
            // Forward recursion first, this is the stable direction:
            //
            int count = 0;
            for(int i = k + 1; i - k < max_iter; ++i)
            {
               poisf *= d2 / (i + 0.5f);
               xtermf *= (x * (n / 2 + i - 1)) / (i);
               betaf += xtermf;

               T term = poisf * betaf;
               sum += term;
               if(fabs(term/sum) < errtol)
                  break;
               ++count;
            }
            //
            // Backwards recursion:
            //
            for(int i = k; i >= 0; --i)
            {
               T term = beta * pois;
               sum += term;
               if(fabs(term/sum) < errtol)
                  break;
               pois *= (i + 0.5f) / d2;
               beta -= xterm;
               xterm *= (i) / (x * (n / 2 + i - 1));
               ++count;
            }
            return sum;
         }

         template <class T, class Policy>
         T non_central_t_cdf(T n, T delta, T t, bool invert, const Policy& pol)
         {
            //
            // For t < 0 we have to use reflect:
            //
            if(t < 0)
            {
               t = -t;
               delta = -delta;
               invert = !invert;
            }
            //
            // x and y are the corresponding random
            // variables for the noncentral beta distribution,
            // with y = 1 - x:
            //
            T x = t * t / (n + t * t);
            T y = n / (n + t * t);
            T d2 = delta * delta;
            T a = 0.5f;
            T b = n / 2;
            T c = a + b + d2 / 2;
            //
            // Crossover point for calculating p or q is the same
            // as for the noncentral beta:
            //
            T cross = 1 - (b / c) * (1 + d2 / (2 * c * c));
            T result;
            if(x < cross)
            {
               //
               // Calculate p:
               //
               result = non_central_beta_p(a, b, d2, x, y, pol);
               result = non_central_t2_p(n, delta, x, y, pol, result);
               result /= 2;
               result += cdf(boost::math::normal_distribution<T, Policy>(), -delta);
            }
            else
            {
               //
               // Calculate q:
               //
               invert = !invert;
               result = non_central_beta_q(a, b, d2, x, y, pol);
               result = non_central_t2_q(n, delta, x, y, pol, result);
               result /= 2;
            }
            if(invert)
               result = 1 - result;
            return result;
         }

         template <class T, class Policy>
         T non_central_t_quantile(T v, T delta, T p, T q, const Policy&)
         {
            static const char* function = "quantile(non_central_t_distribution<%1%>, %1%)";
            typedef typename policies::evaluation<T, Policy>::type value_type;
            typedef typename policies::normalise<
               Policy, 
               policies::promote_float<false>, 
               policies::promote_double<false>, 
               policies::discrete_quantile<>,
               policies::assert_undefined<> >::type forwarding_policy;

               T r;
               if(!detail::check_df(
                  function,
                  v, &r, Policy())
                  ||
               !detail::check_finite(
                  function,
                  delta,
                  &r,
                  Policy())
                  ||
               !detail::check_probability(
                  function,
                  p,
                  &r,
                  Policy()))
                     return r;

            value_type mean = delta * sqrt(v / 2) * tgamma_delta_ratio((v - 1) * 0.5f, T(0.5f));
            value_type var = ((delta * delta + 1) * v) / (v - 2) - mean * mean;
            value_type guess;
            if(p < q)
               guess = quantile(normal_distribution<value_type, forwarding_policy>(mean, var), p);
            else
               guess = quantile(complement(normal_distribution<value_type, forwarding_policy>(mean, var), q));
            value_type result = detail::generic_quantile(
               non_central_t_distribution<value_type, forwarding_policy>(v, delta), 
               (p < q ? p : q), 
               guess, 
               (p >= q), 
               function);
            return policies::checked_narrowing_cast<T, forwarding_policy>(
               result, 
               function);
         }

         template <class T, class Policy>
         T non_central_t2_pdf(T n, T delta, T x, T y, const Policy& pol, T init_val)
         {
            BOOST_MATH_STD_USING
            //
            // Variables come first:
            //
            boost::uintmax_t max_iter = policies::get_max_series_iterations<Policy>();
            T errtol = ldexp(1.0f, -boost::math::policies::digits<T, Policy>());
            T d2 = delta * delta / 2;
            //
            // k is the starting point for iteration, and is the
            // maximum of the poisson weighting term:
            //
            int k = boost::math::itrunc(d2);
            // Starting Poisson weight:
            T pois = gamma_p_derivative(T(k+1), d2, pol) 
               * tgamma_delta_ratio(T(k + 1), T(0.5f))
               * delta / constants::root_two<T>();
            // Starting beta term:
            T xterm = x < y
               ? ibeta_derivative(T(k + 1), n / 2, x, pol)
               : ibeta_derivative(n / 2, T(k + 1), y, pol);
            T poisf(pois), xtermf(xterm);
            T sum = init_val;

            //
            // Backwards recursion first, this is the stable
            // direction for recursion:
            //
            int count = 0;
            for(int i = k; i >= 0; --i)
            {
               T term = xterm * pois;
               sum += term;
               if(fabs(term/sum) < errtol)
                  break;
               pois *= (i + 0.5f) / d2;
               xterm *= (i) / (x * (n / 2 + i));
               ++count;
            }
            for(int i = k + 1; i - k < max_iter; ++i)
            {
               poisf *= d2 / (i + 0.5f);
               xtermf *= (x * (n / 2 + i)) / (i);
               T term = poisf * xtermf;
               sum += term;
               if(fabs(term/sum) < errtol)
                  break;
               ++count;
            }
            return sum;
         }

         template <class T, class Policy>
         T non_central_t_pdf(T n, T delta, T t, const Policy& pol)
         {
            //
            // For t < 0 we have to use reflect:
            //
            if(t < 0)
            {
               t = -t;
               delta = -delta;
            }
            //
            // x and y are the corresponding random
            // variables for the noncentral beta distribution,
            // with y = 1 - x:
            //
            T x = t * t / (n + t * t);
            T y = n / (n + t * t);
            T a = 0.5f;
            T b = n / 2;
            T d2 = delta * delta;
            //
            // Calculate pdf:
            //
            T dt = 2 * n * t / (n * n + 2 * n * t * t + t * t * t * t);
            T result = non_central_beta_pdf(a, b, d2, x, pol);
            result = non_central_t2_pdf(n, delta, x, y, pol, result);
            result *= dt / 2;
            return result;
         }

         template <class T, class Policy>
         T mean(T v, T delta, const Policy& pol)
         {
            return delta * sqrt(v / 2) * tgamma_delta_ratio((v - 1) * 0.5f, T(0.5f), pol);
         }

         template <class T, class Policy>
         T variance(T v, T delta, const Policy& pol)
         {
            T result = ((delta * delta + 1) * v) / (v - 2);
            T m = mean(v, delta, pol);
            result -= m * m;
            return result;
         }

         template <class T, class Policy>
         T skewness(T v, T delta, const Policy& pol)
         {
            T mean = boost::math::detail::mean(v, delta, pol);
            T l2 = delta * delta;
            T var = ((l2 + 1) * v) / (v - 2) - mean * mean;
            T result = -2 * var;
            result += v * (l2 + 2 * v - 3) / ((v - 3) * (v - 2));
            result *= mean;
            result /= pow(var, T(1.5f));
            return result;
         }

         template <class T, class Policy>
         T kurtosis_excess(T v, T delta, const Policy& pol)
         {
            T mean = boost::math::detail::mean(v, delta, pol);
            T l2 = delta * delta;
            T var = ((l2 + 1) * v) / (v - 2) - mean * mean;
            T result = -3 * var;
            result += v * (l2 * (v + 1) + 3 * (3 * v - 5)) / ((v - 3) * (v - 2));
            result *= -mean * mean;
            result += v * v * (l2 * l2 + 6 * l2 + 3) / ((v - 4) * (v - 2));
            result /= var * var;
            return result;
         }

      } // namespace detail

      template <class RealType = double, class Policy = policies::policy<> >
      class non_central_t_distribution
      {
      public:
         typedef RealType value_type;
         typedef Policy policy_type;

         non_central_t_distribution(RealType v_, RealType lambda) : v(v_), ncp(lambda)
         { 
            const char* function = "boost::math::non_central_t_distribution<%1%>::non_central_t_distribution(%1%,%1%)";
            RealType r;
            detail::check_df(
               function,
               v, &r, Policy());
            detail::check_finite(
               function,
               lambda,
               &r,
               Policy());
         } // non_central_t_distribution constructor.

         RealType degrees_of_freedom() const
         { // Private data getter function.
            return v;
         }
         RealType non_centrality() const
         { // Private data getter function.
            return ncp;
         }
      private:
         // Data member, initialized by constructor.
         RealType v;   // degrees of freedom
         RealType ncp; // non-centrality parameter
      }; // template <class RealType, class Policy> class non_central_t_distribution

      typedef non_central_t_distribution<double> non_central_t; // Reserved name of type double.

      // Non-member functions to give properties of the distribution.

      template <class RealType, class Policy>
      inline const std::pair<RealType, RealType> range(const non_central_t_distribution<RealType, Policy>& /* dist */)
      { // Range of permissible values for random variable k.
         using boost::math::tools::max_value;
         return std::pair<RealType, RealType>(-max_value<RealType>(), max_value<RealType>());
      }

      template <class RealType, class Policy>
      inline const std::pair<RealType, RealType> support(const non_central_t_distribution<RealType, Policy>& /* dist */)
      { // Range of supported values for random variable k.
         // This is range where cdf rises from 0 to 1, and outside it, the pdf is zero.
         using boost::math::tools::max_value;
         return std::pair<RealType, RealType>(-max_value<RealType>(), max_value<RealType>());
      }

      template <class RealType, class Policy>
      inline RealType mode(const non_central_t_distribution<RealType, Policy>& dist)
      { // mode.
         static const char* function = "mode(non_central_t_distribution<%1%> const&)";
         RealType v = dist.degrees_of_freedom();
         RealType l = dist.non_centrality();
         RealType r;
         if(!detail::check_df(
            function,
            v, &r, Policy())
            ||
         !detail::check_finite(
            function,
            l,
            &r,
            Policy()))
               return (RealType)r;
         return detail::generic_find_mode(
            dist, 
            detail::mean(v, l, Policy()), 
            function);
      }

      template <class RealType, class Policy>
      inline RealType mean(const non_central_t_distribution<RealType, Policy>& dist)
      { 
         BOOST_MATH_STD_USING
         const char* function = "mean(const non_central_t_distribution<%1%>&)";
         typedef typename policies::evaluation<RealType, Policy>::type value_type;
         typedef typename policies::normalise<
            Policy, 
            policies::promote_float<false>, 
            policies::promote_double<false>, 
            policies::discrete_quantile<>,
            policies::assert_undefined<> >::type forwarding_policy;
         RealType v = dist.degrees_of_freedom();
         RealType l = dist.non_centrality();
         RealType r;
         if(!detail::check_df(
            function,
            v, &r, Policy())
            ||
         !detail::check_finite(
            function,
            l,
            &r,
            Policy()))
               return (RealType)r;
         if(v <= 1)
            return policies::raise_domain_error<RealType>(
               function, 
               "The non central t distribution has no defined mean for degrees of freedom <= 1: got v=%1%.", v, Policy());
         // return l * sqrt(v / 2) * tgamma_delta_ratio((v - 1) * 0.5f, RealType(0.5f));
         return policies::checked_narrowing_cast<RealType, forwarding_policy>(
            detail::mean(static_cast<value_type>(v), static_cast<value_type>(l), forwarding_policy()), function);

      } // mean

      template <class RealType, class Policy>
      inline RealType variance(const non_central_t_distribution<RealType, Policy>& dist)
      { // variance.
         const char* function = "variance(const non_central_t_distribution<%1%>&)";
         typedef typename policies::evaluation<RealType, Policy>::type value_type;
         typedef typename policies::normalise<
            Policy, 
            policies::promote_float<false>, 
            policies::promote_double<false>, 
            policies::discrete_quantile<>,
            policies::assert_undefined<> >::type forwarding_policy;
         BOOST_MATH_STD_USING
         RealType v = dist.degrees_of_freedom();
         RealType l = dist.non_centrality();
         RealType r;
         if(!detail::check_df(
            function,
            v, &r, Policy())
            ||
         !detail::check_finite(
            function,
            l,
            &r,
            Policy()))
               return (RealType)r;
         if(v <= 2)
            return policies::raise_domain_error<RealType>(
               function, 
               "The non central t distribution has no defined variance for degrees of freedom <= 2: got v=%1%.", v, Policy());
         return policies::checked_narrowing_cast<RealType, forwarding_policy>(
            detail::variance(static_cast<value_type>(v), static_cast<value_type>(l), forwarding_policy()), function);
      }

      // RealType standard_deviation(const non_central_t_distribution<RealType, Policy>& dist)
      // standard_deviation provided by derived accessors.

      template <class RealType, class Policy>
      inline RealType skewness(const non_central_t_distribution<RealType, Policy>& dist)
      { // skewness = sqrt(l).
         const char* function = "skewness(const non_central_t_distribution<%1%>&)";
         typedef typename policies::evaluation<RealType, Policy>::type value_type;
         typedef typename policies::normalise<
            Policy, 
            policies::promote_float<false>, 
            policies::promote_double<false>, 
            policies::discrete_quantile<>,
            policies::assert_undefined<> >::type forwarding_policy;
         RealType v = dist.degrees_of_freedom();
         RealType l = dist.non_centrality();
         RealType r;
         if(!detail::check_df(
            function,
            v, &r, Policy())
            ||
         !detail::check_finite(
            function,
            l,
            &r,
            Policy()))
               return (RealType)r;
         if(v <= 3)
            return policies::raise_domain_error<RealType>(
               function, 
               "The non central t distribution has no defined skewness for degrees of freedom <= 3: got v=%1%.", v, Policy());;
         return policies::checked_narrowing_cast<RealType, forwarding_policy>(
            detail::skewness(static_cast<value_type>(v), static_cast<value_type>(l), forwarding_policy()), function);
      }

      template <class RealType, class Policy>
      inline RealType kurtosis_excess(const non_central_t_distribution<RealType, Policy>& dist)
      { 
         const char* function = "kurtosis_excess(const non_central_t_distribution<%1%>&)";
         typedef typename policies::evaluation<RealType, Policy>::type value_type;
         typedef typename policies::normalise<
            Policy, 
            policies::promote_float<false>, 
            policies::promote_double<false>, 
            policies::discrete_quantile<>,
            policies::assert_undefined<> >::type forwarding_policy;
         RealType v = dist.degrees_of_freedom();
         RealType l = dist.non_centrality();
         RealType r;
         if(!detail::check_df(
            function,
            v, &r, Policy())
            ||
         !detail::check_finite(
            function,
            l,
            &r,
            Policy()))
               return (RealType)r;
         if(v <= 4)
            return policies::raise_domain_error<RealType>(
               function, 
               "The non central t distribution has no defined kurtosis for degrees of freedom <= 4: got v=%1%.", v, Policy());;
         return policies::checked_narrowing_cast<RealType, forwarding_policy>(
            detail::kurtosis_excess(static_cast<value_type>(v), static_cast<value_type>(l), forwarding_policy()), function);
      } // kurtosis_excess

      template <class RealType, class Policy>
      inline RealType kurtosis(const non_central_t_distribution<RealType, Policy>& dist)
      {
         return kurtosis_excess(dist) + 3;
      }

      template <class RealType, class Policy>
      inline RealType pdf(const non_central_t_distribution<RealType, Policy>& dist, const RealType& t)
      { // Probability Density/Mass Function.
         const char* function = "cdf(non_central_t_distribution<%1%>, %1%)";
         typedef typename policies::evaluation<RealType, Policy>::type value_type;
         typedef typename policies::normalise<
            Policy, 
            policies::promote_float<false>, 
            policies::promote_double<false>, 
            policies::discrete_quantile<>,
            policies::assert_undefined<> >::type forwarding_policy;

         RealType v = dist.degrees_of_freedom();
         RealType l = dist.non_centrality();
         RealType r;
         if(!detail::check_df(
            function,
            v, &r, Policy())
            ||
         !detail::check_finite(
            function,
            l,
            &r,
            Policy())
            ||
         !detail::check_x(
            function,
            t,
            &r,
            Policy()))
               return (RealType)r;
         return policies::checked_narrowing_cast<RealType, forwarding_policy>(
            detail::non_central_t_pdf(static_cast<value_type>(v), 
               static_cast<value_type>(l), 
               static_cast<value_type>(t), 
               Policy()),
            function);
      } // pdf

      template <class RealType, class Policy>
      RealType cdf(const non_central_t_distribution<RealType, Policy>& dist, const RealType& x)
      { 
         const char* function = "boost::math::non_central_t_distribution<%1%>::cdf(%1%)";
         typedef typename policies::evaluation<RealType, Policy>::type value_type;
         typedef typename policies::normalise<
            Policy, 
            policies::promote_float<false>, 
            policies::promote_double<false>, 
            policies::discrete_quantile<>,
            policies::assert_undefined<> >::type forwarding_policy;

         RealType v = dist.degrees_of_freedom();
         RealType l = dist.non_centrality();
         RealType r;
         if(!detail::check_df(
            function,
            v, &r, Policy())
            ||
         !detail::check_finite(
            function,
            l,
            &r,
            Policy())
            ||
         !detail::check_x(
            function,
            x,
            &r,
            Policy()))
               return (RealType)r;

         return policies::checked_narrowing_cast<RealType, forwarding_policy>(
            detail::non_central_t_cdf(
               static_cast<value_type>(v), 
               static_cast<value_type>(l), 
               static_cast<value_type>(x), 
               false, Policy()),
            function);
      } // cdf

      template <class RealType, class Policy>
      RealType cdf(const complemented2_type<non_central_t_distribution<RealType, Policy>, RealType>& c)
      { // Complemented Cumulative Distribution Function
         const char* function = "boost::math::non_central_t_distribution<%1%>::cdf(%1%)";
         typedef typename policies::evaluation<RealType, Policy>::type value_type;
         typedef typename policies::normalise<
            Policy, 
            policies::promote_float<false>, 
            policies::promote_double<false>, 
            policies::discrete_quantile<>,
            policies::assert_undefined<> >::type forwarding_policy;

         non_central_t_distribution<RealType, Policy> const& dist = c.dist;
         RealType x = c.param;
         RealType v = dist.degrees_of_freedom();
         RealType l = dist.non_centrality();
         RealType r;
         if(!detail::check_df(
            function,
            v, &r, Policy())
            ||
         !detail::check_finite(
            function,
            l,
            &r,
            Policy())
            ||
         !detail::check_x(
            function,
            x,
            &r,
            Policy()))
               return (RealType)r;

         return policies::checked_narrowing_cast<RealType, forwarding_policy>(
            detail::non_central_t_cdf(
               static_cast<value_type>(v), 
               static_cast<value_type>(l), 
               static_cast<value_type>(x), 
               true, Policy()),
            function);
      } // ccdf

      template <class RealType, class Policy>
      inline RealType quantile(const non_central_t_distribution<RealType, Policy>& dist, const RealType& p)
      { // Quantile (or Percent Point) function.
         RealType v = dist.degrees_of_freedom();
         RealType l = dist.non_centrality();
         return detail::non_central_t_quantile(v, l, p, 1-p, Policy());
      } // quantile

      template <class RealType, class Policy>
      inline RealType quantile(const complemented2_type<non_central_t_distribution<RealType, Policy>, RealType>& c)
      { // Quantile (or Percent Point) function.
         non_central_t_distribution<RealType, Policy> const& dist = c.dist;
         RealType q = c.param;
         RealType v = dist.degrees_of_freedom();
         RealType l = dist.non_centrality();
         return detail::non_central_t_quantile(v, l, 1-q, q, Policy());
      } // quantile complement.

   } // namespace math
} // namespace boost

// This include must be at the end, *after* the accessors
// for this distribution have been defined, in order to
// keep compilers that support two-phase lookup happy.
#include <boost/math/distributions/detail/derived_accessors.hpp>

#endif // BOOST_MATH_SPECIAL_NON_CENTRAL_T_HPP

