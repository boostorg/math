// boost\math\distributions\non_central_chi_squared.hpp

// Copyright John Maddock 2008.

// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_SPECIAL_NON_CENTRAL_CHI_SQUARE_HPP
#define BOOST_MATH_SPECIAL_NON_CENTRAL_CHI_SQUARE_HPP

#include <boost/math/distributions/fwd.hpp>
#include <boost/math/special_functions/gamma.hpp> // for incomplete gamma. gamma_q
#include <boost/math/special_functions/bessel.hpp> // for cyl_bessel_i
#include <boost/math/special_functions/round.hpp> // for iround
#include <boost/math/distributions/complement.hpp> // complements
#include <boost/math/distributions/chi_squared.hpp> // central distribution
#include <boost/math/distributions/detail/common_error_handling.hpp> // error checks
#include <boost/math/special_functions/fpclassify.hpp> // isnan.
#include <boost/math/tools/roots.hpp> // for root finding.

namespace boost
{
   namespace math
   {

      template <class RealType, class Policy>
      class non_central_chi_squared_distribution;

      namespace detail{

         template <class T, class Policy>
         T non_central_chi_square_q(T x, T f, T theta, const Policy& pol)
         {
            //
            // Computes the complement of the Non-Central Chi-Square
            // Distribution CDF by summing a weighted sum of complements
            // of the central-distributions.  The weighting factor is
            // a Poisson Distribution.
            //
            // This is an application of the technique described in:
            //
            // Computing discrete mixtures of continuous
            // distributions: noncentral chisquare, noncentral t
            // and the distribution of the square of the sample
            // multiple correlation coeficient.
            // D. Benton, K. Krishnamoorthy.
            // Computational Statistics & Data Analysis 43 (2003) 249 – 267
            //
            BOOST_MATH_STD_USING
            using namespace boost::math;

            //
            // Initialize the variables we'll be using:
            //
            T lambda = theta / 2;
            T del = f / 2;
            T y = x / 2;
            boost::uintmax_t max_iter = policies::get_max_series_iterations<Policy>();
            T errtol = ldexp(1.0, -boost::math::policies::digits<T, Policy>());
            T sum = 0;
            //
            // k is the starting location for iteration, we'll
            // move both forwards and backwards from this point.
            // k is chosen as the peek of the Poisson weights, which
            // will occur *before* the largest term.
            //
            int k = iround(lambda);
            // Forwards and backwards Poisson weights:
            T poisf = boost::math::gamma_p_derivative(1 + k, lambda, pol);
            T poisb = poisf * k / lambda;
            // Initial forwards central chi squared term:
            T gamf = boost::math::gamma_q(del + k, y, pol);
            // Forwards and backwards recursion terms on the central chi squared:
            T xtermf = boost::math::gamma_p_derivative(del + 1 + k, y, pol);
            T xtermb = xtermf * (del + k) / y;
            // Initial backwards central chi squared term:
            T gamb = gamf - xtermb;

            //
            // Forwards iteration first, this is the
            // stable direction for the gamma function
            // recurrences:
            //
            for(int i = k; i < max_iter; ++i)
            {
               T term = poisf * gamf;
               sum += term;
               poisf *= lambda / (i + 1);
               gamf += xtermf;
               xtermf *= y / (del + i + 1);
               if((sum == 0) || (term / sum < errtol))
                  break;
            }
            //
            // Now backwards iteration: the gamma
            // function recurrences are unstable in this
            // direction, we rely on the terms deminishing in size
            // faster than we introduce cancellation errors.  
            // For this reason it's very important that we start
            // *before* the largest term so that backwards iteration
            // is strictly converging.
            //
            for(int i = k - 1; i >= 0; --i)
            {
               T term = poisb * gamb;
               sum += term;
               poisb *= i / lambda;
               xtermb *= (del + i) / y;
               gamb -= xtermb;
               if((sum == 0) || (term / sum < errtol))
                  break;
            }

            return sum;
         }

         template <class T, class Policy>
         T non_central_chi_square_p_ding(T x, T f, T theta, const Policy& pol)
         {
            //
            // This is an implementation of:
            //
            // Algorithm AS 275: 
            // Computing the Non-Central #2 Distribution Function
            // Cherng G. Ding
            // Applied Statistics, Vol. 41, No. 2. (1992), pp. 478-482.
            //
            // This uses a stable forward iteration to sum the
            // CDF, unfortunately this can not be used for large
            // values of the non-centrality parameter because:
            // * The first term may underfow to zero.
            // * We may need an extra-ordinary number of terms
            //   before we reach the first *significant* term.
            //
            BOOST_MATH_STD_USING
            using namespace boost::math;
            T tk = boost::math::gamma_p_derivative(f/2 + 1, x/2, pol);
            T lambda = theta / 2;
            T vk = exp(-lambda);
            T uk = vk;
            T sum = tk * vk;
            if(sum == 0)
               return sum;

            boost::uintmax_t max_iter = policies::get_max_series_iterations<Policy>();
            T errtol = ldexp(1.0, -boost::math::policies::digits<T, Policy>());

            for(int i = 1; i < max_iter; ++i)
            {
               tk = tk * x / (f + 2 * i);
               uk = uk * lambda / i;
               vk = vk + uk;
               T term = vk * tk;
               sum += term;
               if(term / sum < errtol)
                  break;
            }
            return sum;
         }


         template <class T, class Policy>
         T non_central_chi_square_p(T y, T n, T lambda, const Policy& pol)
         {
            //
            // This is taken more or less directly from:
            //
            // Computing discrete mixtures of continuous
            // distributions: noncentral chisquare, noncentral t
            // and the distribution of the square of the sample
            // multiple correlation coeficient.
            // D. Benton, K. Krishnamoorthy.
            // Computational Statistics & Data Analysis 43 (2003) 249 – 267
            //
            // We're summing a Poisson weighting term multiplied by
            // a central chi squared distribution.
            //
            using namespace boost::math;
            BOOST_MATH_STD_USING
            boost::uintmax_t max_iter = policies::get_max_series_iterations<Policy>();
            T errtol = ldexp(1.0, -boost::math::policies::digits<T, Policy>());
            T errorf, errorb;

            T x = y / 2;
            T del = lambda / 2;
            //
            // Starting location for the iteration, we'll iterate
            // both forwards and backwards from this point.  The
            // location chosen is the maximum of the Poisson weight
            // function, which ocurrs *after* the largest term in the
            // sum.
            //
            int k = iround(del);
            T a = n / 2 + k;
            // Central chi squared term for forward iteration:
            T gamkf = boost::math::gamma_p(a, x, pol);

            if(lambda == 0)
               return gamkf;
            // Central chi squared term for backward iteration:
            T gamkb = gamkf;
            // Forwards Poisson weight:
            T poiskf = gamma_p_derivative(k+1, del, pol);
            // Backwards Poisson weight:
            T poiskb = poiskf;
            // Forwards gamma function recursion term:
            T xtermf = boost::math::gamma_p_derivative(a, x, pol);
            // Backwards gamma function recursion term:
            T xtermb = xtermf * x / a;
            T sum = poiskf * gamkf;
            if(sum == 0)
               return sum;
            int i = 1;
            //
            // Backwards recursion first, this is the stable
            // direction for gamma function recurrences:
            //
            while(i <= k) 
            {
               xtermb *= (a - i + 1) / x;
               gamkb += xtermb;
               poiskb = poiskb * (k - i + 1) / del;
               errorb = gamkb * poiskb;
               sum += errorb;
               if(errorb / sum < errtol)
                  break;
               ++i;
            }
            i = 1;
            //
            // Now forwards recursion, the gamma function
            // recurrence relation is unstable in this direction,
            // so we rely on the magnitude of successive terms
            // decreasing faster than we introduce cancellation error.
            // For this reason it's vital that k is chosen to be *after*
            // the largest term, so that successive forward iterations
            // are strictly (and rapidly) converging.
            //
            do
            {
               xtermf = xtermf * x / (a + i - 1);
               gamkf = gamkf - xtermf;
               poiskf = poiskf * del / (k + i);
               errorf = poiskf * gamkf;
               sum += errorf;
               ++i;
            }while((errorf / sum > errtol) && (i < max_iter));

            return sum;
         }

         template <class RealType, class Policy>
         inline RealType non_central_chi_squared_cdf(RealType x, RealType k, RealType l, bool invert, const Policy&)
         {
            typedef typename policies::evaluation<RealType, Policy>::type value_type;
            typedef typename policies::normalise<
               Policy, 
               policies::promote_float<false>, 
               policies::promote_double<false>, 
               policies::discrete_quantile<>,
               policies::assert_undefined<> >::type forwarding_policy;

            value_type result;
            if(l == 0)
               result = cdf(boost::math::chi_squared_distribution<RealType, Policy>(k), x);
            else if(x > k + l)
            {
               // Complement is the smaller of the two:
               result = detail::non_central_chi_square_q(
                  static_cast<value_type>(x), 
                  static_cast<value_type>(k), 
                  static_cast<value_type>(l), forwarding_policy());
               invert = !invert;
            }
            else if(l < 200)
            {
               // For small values of the non-centrality parameter
               // we can use Ding's method:
               result = detail::non_central_chi_square_p_ding(
                  static_cast<value_type>(x), 
                  static_cast<value_type>(k), 
                  static_cast<value_type>(l), forwarding_policy());
            }
            else
            {
               // For largers values of the non-centrality
               // parameter Ding's method will consume an
               // extra-ordinary number of terms, and worse
               // may return zero when the result is in fact
               // finite, use Krishnamoorthy's method instead:
               result = detail::non_central_chi_square_p(
                  static_cast<value_type>(x), 
                  static_cast<value_type>(k), 
                  static_cast<value_type>(l), forwarding_policy());
            }
            if(invert)
               result = 1 - result;
            return policies::checked_narrowing_cast<RealType, forwarding_policy>(
               result, 
               "boost::math::non_central_chi_squared_cdf<%1%>(%1%, %1%, %1%)");
         }

         template <class T, class Policy>
         struct nccs_quantile_functor
         {
            nccs_quantile_functor(const non_central_chi_squared_distribution<T,Policy>& d, T t, bool c)
               : dist(d), target(t), comp(c) {}

            T operator()(const T& x)
            {
               return comp ?
                  target - cdf(complement(dist, x))
                  : cdf(dist, x) - target;
            }

         private:
            non_central_chi_squared_distribution<T,Policy> dist;
            T target;
            bool comp;
         };

         template <class RealType, class Policy>
         RealType nccs_quantile(const non_central_chi_squared_distribution<RealType, Policy>& dist, const RealType& p, bool comp)
         {
            static const char* function = "quantile(non_central_chi_squared_distribution<%1%>, %1%)";
            typedef typename policies::evaluation<RealType, Policy>::type value_type;
            typedef typename policies::normalise<
               Policy, 
               policies::promote_float<false>, 
               policies::promote_double<false>, 
               policies::discrete_quantile<>,
               policies::assert_undefined<> >::type forwarding_policy;

            value_type k = dist.degrees_of_freedom();
            value_type l = dist.non_centrality();
            value_type r;
            if(!detail::check_df(
               function,
               k, &r, Policy())
               ||
            !detail::check_non_centrality(
               function,
               l,
               &r,
               Policy())
               ||
            !detail::check_probability(
               function,
               static_cast<value_type>(p),
               &r,
               Policy()))
                  return (RealType)r;
            //
            // Special cases first:
            //
            if(p == 0)
               return comp
               ? policies::raise_overflow_error<RealType>(function, 0, Policy())
               : 0;
            if(p == 1)
               return !comp
               ? policies::raise_overflow_error<RealType>(function, 0, Policy())
               : 0;

            value_type b = (l * l) / (k + 3 * l);
            value_type c = (k + 3 * l) / (k + 2 * l);
            value_type ff = (k + 2 * l) / (c * c);
            value_type guess;
            if(comp)
               guess = b + c * quantile(complement(chi_squared_distribution<value_type, forwarding_policy>(ff), p));
            else
               guess = b + c * quantile(chi_squared_distribution<value_type, forwarding_policy>(ff), p);

            if(guess < 0)
               guess = tools::min_value<value_type>();

            detail::nccs_quantile_functor<value_type, Policy>
               f(non_central_chi_squared_distribution<value_type, Policy>(k, l), p, comp);
            tools::eps_tolerance<value_type> tol(policies::digits<RealType, Policy>());
            boost::uintmax_t max_iter = policies::get_max_root_iterations<Policy>();
            std::pair<value_type, value_type> ir = tools::bracket_and_solve_root(
               f, guess, value_type(2), true, tol, max_iter, Policy());
            value_type result = ir.first + (ir.second - ir.first) / 2;
            if(max_iter == policies::get_max_root_iterations<Policy>())
            {
               policies::raise_evaluation_error<value_type>(function, "Unable to locate solution in a reasonable time:"
                  " either there is no answer to quantile of the non central chi squared distribution"
                  " or the answer is infinite.  Current best guess is %1%", result, Policy());
            }
            return policies::checked_narrowing_cast<RealType, forwarding_policy>(
               result, 
               function);
         }

      }

      template <class RealType = double, class Policy = policies::policy<> >
      class non_central_chi_squared_distribution
      {
      public:
         typedef RealType value_type;
         typedef Policy policy_type;

         non_central_chi_squared_distribution(RealType df_, RealType lambda) : df(df_), ncp(lambda)
         { 
            const char* function = "boost::math::non_central_chi_squared_distribution<%1%>::non_central_chi_squared_distribution(%1%,%1%)";
            RealType r;
            detail::check_df(
               function,
               df, &r, Policy());
            detail::check_non_centrality(
               function,
               ncp,
               &r,
               Policy());
         } // non_central_chi_squared_distribution constructor.

         RealType degrees_of_freedom() const
         { // Private data getter function.
            return df;
         }
         RealType non_centrality() const
         { // Private data getter function.
            return ncp;
         }
      private:
         // Data member, initialized by constructor.
         RealType df; // degrees of freedom.
         RealType ncp; // non-centrality parameter
      }; // template <class RealType, class Policy> class non_central_chi_squared_distribution

      typedef non_central_chi_squared_distribution<double> non_central_chi_squared; // Reserved name of type double.

      // Non-member functions to give properties of the distribution.

      template <class RealType, class Policy>
      inline const std::pair<RealType, RealType> range(const non_central_chi_squared_distribution<RealType, Policy>& /* dist */)
      { // Range of permissible values for random variable k.
         using boost::math::tools::max_value;
         return std::pair<RealType, RealType>(0, max_value<RealType>()); // Max integer?
      }

      template <class RealType, class Policy>
      inline const std::pair<RealType, RealType> support(const non_central_chi_squared_distribution<RealType, Policy>& /* dist */)
      { // Range of supported values for random variable k.
         // This is range where cdf rises from 0 to 1, and outside it, the pdf is zero.
         using boost::math::tools::max_value;
         return std::pair<RealType, RealType>(0,  max_value<RealType>());
      }

      template <class RealType, class Policy>
      inline RealType mean(const non_central_chi_squared_distribution<RealType, Policy>& dist)
      { // Mean of poisson distribution = lambda.
         const char* function = "boost::math::non_central_chi_squared_distribution<%1%>::mean()";
         RealType k = dist.degrees_of_freedom();
         RealType l = dist.non_centrality();
         RealType r;
         if(!detail::check_df(
            function,
            k, &r, Policy())
            ||
         !detail::check_non_centrality(
            function,
            l,
            &r,
            Policy()))
               return r;
         return k + l;
      } // mean

      template <class RealType, class Policy>
      inline RealType mode(const non_central_chi_squared_distribution<RealType, Policy>& dist)
      { // mode.
         // TODO!
         return 0;
      }

      template <class RealType, class Policy>
      inline RealType variance(const non_central_chi_squared_distribution<RealType, Policy>& dist)
      { // variance.
         const char* function = "boost::math::non_central_chi_squared_distribution<%1%>::variance()";
         RealType k = dist.degrees_of_freedom();
         RealType l = dist.non_centrality();
         RealType r;
         if(!detail::check_df(
            function,
            k, &r, Policy())
            ||
         !detail::check_non_centrality(
            function,
            l,
            &r,
            Policy()))
               return r;
         return 2 * (2 * l + k);
      }

      // RealType standard_deviation(const non_central_chi_squared_distribution<RealType, Policy>& dist)
      // standard_deviation provided by derived accessors.

      template <class RealType, class Policy>
      inline RealType skewness(const non_central_chi_squared_distribution<RealType, Policy>& dist)
      { // skewness = sqrt(l).
         const char* function = "boost::math::non_central_chi_squared_distribution<%1%>::skewness()";
         RealType k = dist.degrees_of_freedom();
         RealType l = dist.non_centrality();
         RealType r;
         if(!detail::check_df(
            function,
            k, &r, Policy())
            ||
         !detail::check_non_centrality(
            function,
            l,
            &r,
            Policy()))
               return r;
         BOOST_MATH_STD_USING
            return pow(2 / (k + 2 * l), RealType(3)/2) * (k + 3 * l);
      }

      template <class RealType, class Policy>
      inline RealType kurtosis_excess(const non_central_chi_squared_distribution<RealType, Policy>& dist)
      { 
         const char* function = "boost::math::non_central_chi_squared_distribution<%1%>::kurtosis_excess()";
         RealType k = dist.degrees_of_freedom();
         RealType l = dist.non_centrality();
         RealType r;
         if(!detail::check_df(
            function,
            k, &r, Policy())
            ||
         !detail::check_non_centrality(
            function,
            l,
            &r,
            Policy()))
               return r;
         return 12 * (k + 4 * l) / ((k + 2 * l) * (k + 2 * l));
      } // kurtosis_excess

      template <class RealType, class Policy>
      inline RealType kurtosis(const non_central_chi_squared_distribution<RealType, Policy>& dist)
      {
         return kurtosis_excess(dist) + 3;
      }

      template <class RealType, class Policy>
      RealType pdf(const non_central_chi_squared_distribution<RealType, Policy>& dist, const RealType& x)
      { // Probability Density/Mass Function.
         const char* function = "boost::math::non_central_chi_squared_distribution<%1%>::pdf(%1%)";
         RealType k = dist.degrees_of_freedom();
         RealType l = dist.non_centrality();
         RealType r;
         if(!detail::check_df(
            function,
            k, &r, Policy())
            ||
         !detail::check_non_centrality(
            function,
            l,
            &r,
            Policy()))
               return r;

         if(l == 0)
            return pdf(boost::math::chi_squared_distribution<RealType, Policy>(k), x);

         r = -(x + l) / 2 + log(x / l) * (k / 4 - 0.5f);

         return 0.5f * exp(r)
            * boost::math::cyl_bessel_i(k/2 - 1, sqrt(l * x), Policy());
      } // pdf

      template <class RealType, class Policy>
      RealType cdf(const non_central_chi_squared_distribution<RealType, Policy>& dist, const RealType& x)
      { 
         const char* function = "boost::math::non_central_chi_squared_distribution<%1%>::cdf(%1%)";
         RealType k = dist.degrees_of_freedom();
         RealType l = dist.non_centrality();
         RealType r;
         if(!detail::check_df(
            function,
            k, &r, Policy())
            ||
         !detail::check_non_centrality(
            function,
            l,
            &r,
            Policy()))
               return r;

         return detail::non_central_chi_squared_cdf(x, k, l, false, Policy());
      } // cdf

      template <class RealType, class Policy>
      RealType cdf(const complemented2_type<non_central_chi_squared_distribution<RealType, Policy>, RealType>& c)
      { // Complemented Cumulative Distribution Function
         const char* function = "boost::math::non_central_chi_squared_distribution<%1%>::cdf(%1%)";
         non_central_chi_squared_distribution<RealType, Policy> const& dist = c.dist;
         RealType x = c.param;
         RealType k = dist.degrees_of_freedom();
         RealType l = dist.non_centrality();
         RealType r;
         if(!detail::check_df(
            function,
            k, &r, Policy())
            ||
         !detail::check_non_centrality(
            function,
            l,
            &r,
            Policy()))
               return r;

         return detail::non_central_chi_squared_cdf(x, k, l, true, Policy());
      } // ccdf

      template <class RealType, class Policy>
      inline RealType quantile(const non_central_chi_squared_distribution<RealType, Policy>& dist, const RealType& p)
      { // Quantile (or Percent Point) function.
         return detail::nccs_quantile(dist, p, false);
      } // quantile

      template <class RealType, class Policy>
      inline RealType quantile(const complemented2_type<non_central_chi_squared_distribution<RealType, Policy>, RealType>& c)
      { // Quantile (or Percent Point) function.
         return detail::nccs_quantile(c.dist, c.param, true);
      } // quantile complement.

   } // namespace math
} // namespace boost

// This include must be at the end, *after* the accessors
// for this distribution have been defined, in order to
// keep compilers that support two-phase lookup happy.
#include <boost/math/distributions/detail/derived_accessors.hpp>

#endif // BOOST_MATH_SPECIAL_NON_CENTRAL_CHI_SQUARE_HPP



