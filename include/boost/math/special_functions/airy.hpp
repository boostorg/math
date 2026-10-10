// Copyright John Maddock 2012.
// Copyright Matt Borland 2024.
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_AIRY_HPP
#define BOOST_MATH_AIRY_HPP

#include <boost/math/tools/config.hpp>
#include <boost/math/tools/numeric_limits.hpp>
#include <boost/math/tools/precision.hpp>
#include <boost/math/tools/cstdint.hpp>
#include <boost/math/special_functions/math_fwd.hpp>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/cbrt.hpp>
#include <boost/math/special_functions/detail/airy_ai_bi_zero.hpp>
#include <boost/math/tools/roots.hpp>
#include <boost/math/tools/series.hpp>
#include <boost/math/policies/error_handling.hpp>
#include <boost/math/constants/constants.hpp>

BOOST_MATH_NAMESPACE_BEGIN

namespace detail{

// Whether Ai, Bi and their derivatives at zero are available as correctly rounded constants:
template <class T>
BOOST_MATH_GPU_ENABLED constexpr bool airy_has_exact_zero_values()
{
   return BOOST_MATH_NAMESPACE::numeric_limits<T>::radix != 10 && BOOST_MATH_NAMESPACE::tools::digits<T>() <= 390;
}

// Ai(0), Bi(0), Ai'(0) and Bi'(0).  Types of up to 390 bits get constants to 120 digits (from Arb),
// evaluating them via tgamma loses a few bits as the arguments 1/3 and 2/3 are not exact:
template <class T, class Policy>
BOOST_MATH_GPU_ENABLED T airy_ai_zero_value(const Policy& pol)
{
   BOOST_MATH_STD_USING
   BOOST_MATH_IF_CONSTEXPR(airy_has_exact_zero_values<T>())
      return BOOST_MATH_BIG_CONSTANT(T, 390, 0.355028053887817239260063186004183176397979174199177240583326510300810042450126712957174246054040271688420448730349495840);
   else
      return 1 / (pow(T(3), constants::twothirds<T>()) * BOOST_MATH_NAMESPACE::tgamma(constants::twothirds<T>(), pol));
}

template <class T, class Policy>
BOOST_MATH_GPU_ENABLED T airy_bi_zero_value(const Policy& pol)
{
   BOOST_MATH_STD_USING
   BOOST_MATH_IF_CONSTEXPR(airy_has_exact_zero_values<T>())
      return BOOST_MATH_BIG_CONSTANT(T, 390, 0.614926627446000735150922369093613553594728188648596505040878753014296519305520640529387343345267569240728438782242516725);
   else
      return 1 / (sqrt(BOOST_MATH_NAMESPACE::cbrt(T(3), pol)) * BOOST_MATH_NAMESPACE::tgamma(constants::twothirds<T>(), pol));
}

// Returns -Ai'(0):
template <class T, class Policy>
BOOST_MATH_GPU_ENABLED T airy_ai_prime_zero_value(const Policy& pol)
{
   BOOST_MATH_STD_USING
   BOOST_MATH_IF_CONSTEXPR(airy_has_exact_zero_values<T>())
      return BOOST_MATH_BIG_CONSTANT(T, 390, 0.258819403792806798405183560189203963479091138354934582210001813856102772676790280654196405827275384313371193211789133381);
   else
      return 1 / (BOOST_MATH_NAMESPACE::cbrt(T(3), pol) * BOOST_MATH_NAMESPACE::tgamma(constants::third<T>(), pol));
}

template <class T, class Policy>
BOOST_MATH_GPU_ENABLED T airy_bi_prime_zero_value(const Policy& pol)
{
   BOOST_MATH_STD_USING
   BOOST_MATH_IF_CONSTEXPR(airy_has_exact_zero_values<T>())
      return BOOST_MATH_BIG_CONSTANT(T, 390, 0.448288357353826357914823710398828390866226799212262061082808778372330755009780647185046574400736362878496329249031699774);
   else
      return sqrt(BOOST_MATH_NAMESPACE::cbrt(T(3), pol)) / BOOST_MATH_NAMESPACE::tgamma(constants::third<T>(), pol);
}

// Terms of the Maclaurin series of Ai and Bi, which are combinations of
//   f(x) = 1 + x^3 / 6 + ...,  g(x) = x + x^4 / 12 + ...  and their derivatives
//   f'(x) = x^2 / 2 + ..., g'(x) = 1 + x^3 / 3 + ...
// with each term being the previous one times x^3 / ((3k + a)(3k + b)):
//
template <class T>
struct airy_maclaurin_series
{
   typedef T result_type;
   BOOST_MATH_GPU_ENABLED airy_maclaurin_series(T first, T x, int a_, int b_)
      : term(first), x3(x * x * x), a(a_), b(b_), k(0) {}
   BOOST_MATH_GPU_ENABLED T operator()()
   {
      T result = term;
      ++k;
      term *= x3 / (T(3 * k + a) * T(3 * k + b));
      return result;
   }
private:
   T term, x3;
   int a, b, k;
};

template <class T, class Policy>
BOOST_MATH_GPU_ENABLED T airy_maclaurin_sum(T first, T x, int a, int b, const Policy& pol)
{
   airy_maclaurin_series<T> series(first, x, a, b);
   BOOST_MATH_NAMESPACE::uintmax_t max_iter = policies::get_max_series_iterations<Policy>();
   T result = tools::sum_series(series, policies::get_epsilon<T, Policy>(), max_iter);
   policies::check_series_iterations<T>("boost::math::airy<%1%>(%1%)", max_iter, pol);
   return result;
}

template <class T, class Policy>
BOOST_MATH_GPU_ENABLED T airy_f(T x, const Policy& pol) { return airy_maclaurin_sum(T(1), x, -1, 0, pol); }
template <class T, class Policy>
BOOST_MATH_GPU_ENABLED T airy_g(T x, const Policy& pol) { return airy_maclaurin_sum(x, x, 0, 1, pol); }
template <class T, class Policy>
BOOST_MATH_GPU_ENABLED T airy_f_prime(T x, const Policy& pol) { return airy_maclaurin_sum(T(x * x / 2), x, 0, 2, pol); }
template <class T, class Policy>
BOOST_MATH_GPU_ENABLED T airy_g_prime(T x, const Policy& pol) { return airy_maclaurin_sum(T(1), x, -2, 0, pol); }

template <class T, class Policy>
BOOST_MATH_GPU_ENABLED T airy_ai_imp(T x, const Policy& pol)
{
   BOOST_MATH_STD_USING

   if(x > -3 && x < 0.5f)
   {
      return airy_ai_zero_value<T>(pol) * airy_f(x, pol) - airy_ai_prime_zero_value<T>(pol) * airy_g(x, pol);
   }

   if(x < 0)
   {
      T p = (-x * sqrt(-x) * 2) / 3;
      T v = T(1) / 3;
      T j1 = BOOST_MATH_NAMESPACE::cyl_bessel_j(v, p, pol);
      T j2 = BOOST_MATH_NAMESPACE::cyl_bessel_j(-v, p, pol);
      T ai = sqrt(-x) * (j1 + j2) / 3;
      //T bi = sqrt(-x / 3) * (j2 - j1);
      return ai;
   }
   else
   {
      T p = 2 * x * sqrt(x) / 3;
      T v = T(1) / 3;
      //T j1 = boost::math::cyl_bessel_i(-v, p, pol);
      //T j2 = boost::math::cyl_bessel_i(v, p, pol);
      //
      // Note that although we can calculate ai from j1 and j2, the accuracy is horrible
      // as we're subtracting two very large values, so use the Bessel K relation instead:
      //
      T ai = cyl_bessel_k(v, p, pol) * sqrt(x / 3) / BOOST_MATH_NAMESPACE::constants::pi<T>();  //sqrt(x) * (j1 - j2) / 3;
      //T bi = sqrt(x / 3) * (j1 + j2);
      return ai;
   }
}

template <class T, class Policy>
BOOST_MATH_GPU_ENABLED T airy_bi_imp(T x, const Policy& pol)
{
   BOOST_MATH_STD_USING

   if(x > -3 && x < 8)
   {
      return airy_bi_zero_value<T>(pol) * airy_f(x, pol) + airy_bi_prime_zero_value<T>(pol) * airy_g(x, pol);
   }

   if(x < 0)
   {
      T p = (-x * sqrt(-x) * 2) / 3;
      T v = T(1) / 3;
      T j1 = BOOST_MATH_NAMESPACE::cyl_bessel_j(v, p, pol);
      T j2 = BOOST_MATH_NAMESPACE::cyl_bessel_j(-v, p, pol);
      //T ai = sqrt(-x) * (j1 + j2) / 3;
      T bi = sqrt(-x / 3) * (j2 - j1);
      return bi;
   }
   else
   {
      T p = 2 * x * sqrt(x) / 3;
      T v = T(1) / 3;
      T j1 = BOOST_MATH_NAMESPACE::cyl_bessel_i(-v, p, pol);
      T j2 = BOOST_MATH_NAMESPACE::cyl_bessel_i(v, p, pol);
      T bi = sqrt(x / 3) * (j1 + j2);
      return bi;
   }
}

template <class T, class Policy>
BOOST_MATH_GPU_ENABLED T airy_ai_prime_imp(T x, const Policy& pol)
{
   BOOST_MATH_STD_USING

   if(x > -3 && x < 0.5f)
   {
      return airy_ai_zero_value<T>(pol) * airy_f_prime(x, pol) - airy_ai_prime_zero_value<T>(pol) * airy_g_prime(x, pol);
   }

   if(x < 0)
   {
      T p = (-x * sqrt(-x) * 2) / 3;
      T v = T(2) / 3;
      T j1 = BOOST_MATH_NAMESPACE::cyl_bessel_j(v, p, pol);
      T j2 = BOOST_MATH_NAMESPACE::cyl_bessel_j(-v, p, pol);
      T aip = -x * (j1 - j2) / 3;
      return aip;
   }
   else
   {
      T p = 2 * x * sqrt(x) / 3;
      T v = T(2) / 3;
      //T j1 = boost::math::cyl_bessel_i(-v, p, pol);
      //T j2 = boost::math::cyl_bessel_i(v, p, pol);
      //
      // Note that although we can calculate ai from j1 and j2, the accuracy is horrible
      // as we're subtracting two very large values, so use the Bessel K relation instead:
      //
      T aip = -cyl_bessel_k(v, p, pol) * x / (BOOST_MATH_NAMESPACE::constants::root_three<T>() * BOOST_MATH_NAMESPACE::constants::pi<T>());
      return aip;
   }
}

template <class T, class Policy>
BOOST_MATH_GPU_ENABLED T airy_bi_prime_imp(T x, const Policy& pol)
{
   BOOST_MATH_STD_USING

   if(x > -3 && x < 8)
   {
      return airy_bi_zero_value<T>(pol) * airy_f_prime(x, pol) + airy_bi_prime_zero_value<T>(pol) * airy_g_prime(x, pol);
   }

   if(x < 0)
   {
      T p = (-x * sqrt(-x) * 2) / 3;
      T v = T(2) / 3;
      T j1 = BOOST_MATH_NAMESPACE::cyl_bessel_j(v, p, pol);
      T j2 = BOOST_MATH_NAMESPACE::cyl_bessel_j(-v, p, pol);
      T aip = -x * (j1 + j2) / constants::root_three<T>();
      return aip;
   }
   else
   {
      T p = 2 * x * sqrt(x) / 3;
      T v = T(2) / 3;
      T j1 = BOOST_MATH_NAMESPACE::cyl_bessel_i(-v, p, pol);
      T j2 = BOOST_MATH_NAMESPACE::cyl_bessel_i(v, p, pol);
      T aip = x * (j1 + j2) / BOOST_MATH_NAMESPACE::constants::root_three<T>();
      return aip;
   }
}

template <class T, class Policy>
BOOST_MATH_GPU_ENABLED T airy_ai_zero_imp(int m, const Policy& pol)
{
   BOOST_MATH_STD_USING // ADL of std names, needed for log, sqrt.

   // Handle cases when a negative zero (negative rank) is requested.
   if(m < 0)
   {
      return policies::raise_domain_error<T>("boost::math::airy_ai_zero<%1%>(%1%, int)",
         "Requested the %1%'th zero, but the rank must be 1 or more !", static_cast<T>(m), pol);
   }

   // Handle case when the zero'th zero is requested.
   if(m == 0U)
   {
      return policies::raise_domain_error<T>("boost::math::airy_ai_zero<%1%>(%1%,%1%)",
        "The requested rank of the zero is %1%, but must be 1 or more !", static_cast<T>(m), pol);
   }

   // Set up the initial guess for the upcoming root-finding.
   const T guess_root = BOOST_MATH_NAMESPACE::detail::airy_zero::airy_ai_zero_detail::initial_guess<T>(m, pol);

   // Select the maximum allowed iterations based on the number
   // of decimal digits in the numeric type T, being at least 12.
   const int my_digits10 = static_cast<int>(static_cast<float>(policies::digits<T, Policy>() * 0.301F));

   const std::uintmax_t iterations_allowed = static_cast<std::uintmax_t>(BOOST_MATH_GPU_SAFE_MAX(12, my_digits10 * 2));

   std::uintmax_t iterations_used = iterations_allowed;

   // Use a dynamic tolerance because the roots get closer the higher m gets.
   T tolerance;  // LCOV_EXCL_LINE

   if     (m <=   10) { tolerance = T(0.3F); }
   else if(m <=  100) { tolerance = T(0.1F); }
   else if(m <= 1000) { tolerance = T(0.05F); }
   else               { tolerance = T(1) / sqrt(T(m)); }

   // Perform the root-finding using Newton-Raphson iteration from Boost.Math.
   const T am =
      BOOST_MATH_NAMESPACE::tools::newton_raphson_iterate(
         BOOST_MATH_NAMESPACE::detail::airy_zero::airy_ai_zero_detail::function_object_ai_and_ai_prime<T, Policy>(pol),
         guess_root,
         T(guess_root - tolerance),
         T(guess_root + tolerance),
         policies::digits<T, Policy>(),
         iterations_used);

   static_cast<void>(iterations_used);

   return am;
}

template <class T, class Policy>
BOOST_MATH_GPU_ENABLED T airy_bi_zero_imp(int m, const Policy& pol)
{
   BOOST_MATH_STD_USING // ADL of std names, needed for log, sqrt.

   // Handle cases when a negative zero (negative rank) is requested.
   if(m < 0)
   {
      return policies::raise_domain_error<T>("boost::math::airy_bi_zero<%1%>(%1%, int)",
         "Requested the %1%'th zero, but the rank must 1 or more !", static_cast<T>(m), pol);
   }

   // Handle case when the zero'th zero is requested.
   if(m == 0U)
   {
      return policies::raise_domain_error<T>("boost::math::airy_bi_zero<%1%>(%1%,%1%)",
        "The requested rank of the zero is %1%, but must be 1 or more !", static_cast<T>(m), pol);
   }
   // Set up the initial guess for the upcoming root-finding.
   const T guess_root = BOOST_MATH_NAMESPACE::detail::airy_zero::airy_bi_zero_detail::initial_guess<T>(m, pol);

   // Select the maximum allowed iterations based on the number
   // of decimal digits in the numeric type T, being at least 12.
   const int my_digits10 = static_cast<int>(static_cast<float>(policies::digits<T, Policy>() * 0.301F));

   const std::uintmax_t iterations_allowed = static_cast<std::uintmax_t>(BOOST_MATH_GPU_SAFE_MAX(12, my_digits10 * 2));

   std::uintmax_t iterations_used = iterations_allowed;

   // Use a dynamic tolerance because the roots get closer the higher m gets.
   T tolerance; // LCOV_EXCL_LINE

   if     (m <=   10) { tolerance = T(0.3F); }
   else if(m <=  100) { tolerance = T(0.1F); }
   else if(m <= 1000) { tolerance = T(0.05F); }
   else               { tolerance = T(1) / sqrt(T(m)); }

   // Perform the root-finding using Newton-Raphson iteration from Boost.Math.
   const T bm =
      BOOST_MATH_NAMESPACE::tools::newton_raphson_iterate(
         BOOST_MATH_NAMESPACE::detail::airy_zero::airy_bi_zero_detail::function_object_bi_and_bi_prime<T, Policy>(pol),
         guess_root,
         T(guess_root - tolerance),
         T(guess_root + tolerance),
         policies::digits<T, Policy>(),
         iterations_used);

   static_cast<void>(iterations_used);

   return bm;
}

} // namespace detail

BOOST_MATH_EXPORT template <class T, class Policy>
BOOST_MATH_GPU_ENABLED inline typename tools::promote_args<T>::type airy_ai(T x, const Policy&)
{
   BOOST_FPU_EXCEPTION_GUARD
   typedef typename tools::promote_args<T>::type result_type;
   typedef typename policies::evaluation<result_type, Policy>::type value_type;
   typedef typename policies::normalise<
      Policy, 
      policies::promote_float<false>, 
      policies::promote_double<false>, 
      policies::discrete_quantile<>,
      policies::assert_undefined<> >::type forwarding_policy;

   return policies::checked_narrowing_cast<result_type, Policy>(detail::airy_ai_imp<value_type>(static_cast<value_type>(x), forwarding_policy()), "boost::math::airy<%1%>(%1%)");
}

BOOST_MATH_EXPORT template <class T>
BOOST_MATH_GPU_ENABLED inline typename tools::promote_args<T>::type airy_ai(T x)
{
   return airy_ai(x, policies::policy<>());
}

BOOST_MATH_EXPORT template <class T, class Policy>
BOOST_MATH_GPU_ENABLED inline typename tools::promote_args<T>::type airy_bi(T x, const Policy&)
{
   BOOST_FPU_EXCEPTION_GUARD
   typedef typename tools::promote_args<T>::type result_type;
   typedef typename policies::evaluation<result_type, Policy>::type value_type;
   typedef typename policies::normalise<
      Policy, 
      policies::promote_float<false>, 
      policies::promote_double<false>, 
      policies::discrete_quantile<>,
      policies::assert_undefined<> >::type forwarding_policy;

   return policies::checked_narrowing_cast<result_type, Policy>(detail::airy_bi_imp<value_type>(static_cast<value_type>(x), forwarding_policy()), "boost::math::airy<%1%>(%1%)");
}

BOOST_MATH_EXPORT template <class T>
BOOST_MATH_GPU_ENABLED inline typename tools::promote_args<T>::type airy_bi(T x)
{
   return airy_bi(x, policies::policy<>());
}

BOOST_MATH_EXPORT template <class T, class Policy>
BOOST_MATH_GPU_ENABLED inline typename tools::promote_args<T>::type airy_ai_prime(T x, const Policy&)
{
   BOOST_FPU_EXCEPTION_GUARD
   typedef typename tools::promote_args<T>::type result_type;
   typedef typename policies::evaluation<result_type, Policy>::type value_type;
   typedef typename policies::normalise<
      Policy, 
      policies::promote_float<false>, 
      policies::promote_double<false>, 
      policies::discrete_quantile<>,
      policies::assert_undefined<> >::type forwarding_policy;

   return policies::checked_narrowing_cast<result_type, Policy>(detail::airy_ai_prime_imp<value_type>(static_cast<value_type>(x), forwarding_policy()), "boost::math::airy<%1%>(%1%)");
}

BOOST_MATH_EXPORT template <class T>
BOOST_MATH_GPU_ENABLED inline typename tools::promote_args<T>::type airy_ai_prime(T x)
{
   return airy_ai_prime(x, policies::policy<>());
}

BOOST_MATH_EXPORT template <class T, class Policy>
BOOST_MATH_GPU_ENABLED inline typename tools::promote_args<T>::type airy_bi_prime(T x, const Policy&)
{
   BOOST_FPU_EXCEPTION_GUARD
   typedef typename tools::promote_args<T>::type result_type;
   typedef typename policies::evaluation<result_type, Policy>::type value_type;
   typedef typename policies::normalise<
      Policy, 
      policies::promote_float<false>, 
      policies::promote_double<false>, 
      policies::discrete_quantile<>,
      policies::assert_undefined<> >::type forwarding_policy;

   return policies::checked_narrowing_cast<result_type, Policy>(detail::airy_bi_prime_imp<value_type>(static_cast<value_type>(x), forwarding_policy()), "boost::math::airy<%1%>(%1%)");
}

BOOST_MATH_EXPORT template <class T>
BOOST_MATH_GPU_ENABLED inline typename tools::promote_args<T>::type airy_bi_prime(T x)
{
   return airy_bi_prime(x, policies::policy<>());
}

BOOST_MATH_EXPORT template <class T, class Policy>
BOOST_MATH_GPU_ENABLED inline T airy_ai_zero(int m, const Policy& /*pol*/)
{
   BOOST_FPU_EXCEPTION_GUARD
   typedef typename policies::evaluation<T, Policy>::type value_type;
   typedef typename policies::normalise<
      Policy, 
      policies::promote_float<false>, 
      policies::promote_double<false>, 
      policies::discrete_quantile<>,
      policies::assert_undefined<> >::type forwarding_policy;

   static_assert(    false == std::numeric_limits<T>::is_specialized
                           || (   true  == std::numeric_limits<T>::is_specialized
                               && false == std::numeric_limits<T>::is_integer),
                           "Airy value type must be a floating-point type.");

   return policies::checked_narrowing_cast<T, Policy>(detail::airy_ai_zero_imp<value_type>(m, forwarding_policy()), "boost::math::airy_ai_zero<%1%>(unsigned)");
}

BOOST_MATH_EXPORT template <class T>
BOOST_MATH_GPU_ENABLED inline T airy_ai_zero(int m)
{
   return airy_ai_zero<T>(m, policies::policy<>());
}

BOOST_MATH_EXPORT template <class T, class OutputIterator, class Policy>
BOOST_MATH_GPU_ENABLED inline OutputIterator airy_ai_zero(
                         int start_index,
                         unsigned number_of_zeros,
                         OutputIterator out_it,
                         const Policy& pol)
{
   typedef T result_type;

   static_assert(    false == std::numeric_limits<T>::is_specialized
                           || (   true  == std::numeric_limits<T>::is_specialized
                               && false == std::numeric_limits<T>::is_integer),
                           "Airy value type must be a floating-point type.");

   for(unsigned i = 0; i < number_of_zeros; ++i)
   {
      *out_it = BOOST_MATH_NAMESPACE::airy_ai_zero<result_type>(start_index + i, pol);
      ++out_it;
   }
   return out_it;
}

BOOST_MATH_EXPORT template <class T, class OutputIterator>
BOOST_MATH_GPU_ENABLED inline OutputIterator airy_ai_zero(
                         int start_index,
                         unsigned number_of_zeros,
                         OutputIterator out_it)
{
   return airy_ai_zero<T>(start_index, number_of_zeros, out_it, policies::policy<>());
}

BOOST_MATH_EXPORT template <class T, class Policy>
BOOST_MATH_GPU_ENABLED inline T airy_bi_zero(int m, const Policy& /*pol*/)
{
   BOOST_FPU_EXCEPTION_GUARD
   typedef typename policies::evaluation<T, Policy>::type value_type;
   typedef typename policies::normalise<
      Policy, 
      policies::promote_float<false>, 
      policies::promote_double<false>, 
      policies::discrete_quantile<>,
      policies::assert_undefined<> >::type forwarding_policy;

   static_assert(    false == std::numeric_limits<T>::is_specialized
                           || (   true  == std::numeric_limits<T>::is_specialized
                               && false == std::numeric_limits<T>::is_integer),
                           "Airy value type must be a floating-point type.");

   return policies::checked_narrowing_cast<T, Policy>(detail::airy_bi_zero_imp<value_type>(m, forwarding_policy()), "boost::math::airy_bi_zero<%1%>(unsigned)");
}

BOOST_MATH_EXPORT template <typename T>
BOOST_MATH_GPU_ENABLED inline T airy_bi_zero(int m)
{
   return airy_bi_zero<T>(m, policies::policy<>());
}

BOOST_MATH_EXPORT template <class T, class OutputIterator, class Policy>
BOOST_MATH_GPU_ENABLED inline OutputIterator airy_bi_zero(
                         int start_index,
                         unsigned number_of_zeros,
                         OutputIterator out_it,
                         const Policy& pol)
{
   typedef T result_type;

   static_assert(    false == std::numeric_limits<T>::is_specialized
                           || (   true  == std::numeric_limits<T>::is_specialized
                               && false == std::numeric_limits<T>::is_integer),
                           "Airy value type must be a floating-point type.");

   for(unsigned i = 0; i < number_of_zeros; ++i)
   {
      *out_it = BOOST_MATH_NAMESPACE::airy_bi_zero<result_type>(start_index + i, pol);
      ++out_it;
   }
   return out_it;
}

BOOST_MATH_EXPORT template <class T, class OutputIterator>
BOOST_MATH_GPU_ENABLED inline OutputIterator airy_bi_zero(
                         int start_index,
                         unsigned number_of_zeros,
                         OutputIterator out_it)
{
   return airy_bi_zero<T>(start_index, number_of_zeros, out_it, policies::policy<>());
}

BOOST_MATH_NAMESPACE_END

#endif // BOOST_MATH_AIRY_HPP
