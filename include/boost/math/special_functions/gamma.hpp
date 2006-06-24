//  (C) Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_SF_GAMMA_HPP
#define BOOST_MATH_SF_GAMMA_HPP

#include <boost/lexical_cast.hpp>
#include <boost/config.hpp>
#include <boost/math/tools/series.hpp>
#include <boost/math/tools/fraction.hpp>
#include <boost/math/tools/precision.hpp>
#include <boost/math/tools/real_cast.hpp>
#include <boost/math/tools/error_handling.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/log1p.hpp>
#include <boost/math/special_functions/powm1.hpp>
#include <boost/math/special_functions/sqrtp1m1.hpp>
#include <boost/math/special_functions/lanczos.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/type_traits/is_convertible.hpp>
#include <boost/assert.hpp>
#include <cmath>
#include <algorithm>

#ifdef BOOST_MATH_INSTRUMENT
#include <iostream>
#include <iomanip>
#include <typeinfo>
#endif


namespace boost{ namespace math{ 
   
//
// Forward declarations:
//
template <class T>
T tgamma(T z);
template <class T>
T lgamma(T z, int* sign);
template <class T>
T lgamma(T x);
template <class T>
T tgammap1m1(T z);template <class T>
T tgamma(T a, T z);
template <class T>
T tgamma_lower(T a, T z);
template <class T>
T gamma_Q(T a, T z);
template <class T>
T gamma_P(T a, T z);
   
namespace detail{

template <class T>
inline bool is_odd(T v, const boost::true_type&)
{
   int i = static_cast<int>(v);
   return i&1;
}
template <class T>
inline bool is_odd(T v, const boost::false_type&)
{
   // Oh dear can't cast T to int!
   using namespace std;
   T modulus = v - 2 * floor(v/2);
   return static_cast<bool>(modulus != 0);
}
template <class T>
inline bool is_odd(T v)
{
   return is_odd(v, ::boost::is_convertible<T, int>());
}

template <class T>
T sinpx(T z)
{
   // add hoc function calculates x*sin(pi*x),
   // taking extra care near when x is near a whole number.
   using namespace std;
   int sign = 1;
   if(z < 0)
   {
      z = -z;
   }
   else
   {
      sign = -sign;
   }
   T fl = floor(z);
   T dist;
   if(is_odd(fl))
   {
      fl += 1;
      dist = fl - z;
      sign = -sign;
   }
   else
   {
      dist = z - fl;
   }
   BOOST_ASSERT(fl >= 0);
   if(dist > 0.5)
      dist = 1 - dist;
   T result = sin(dist*boost::math::constants::pi<T>());
   return sign*z*result;
}
//
// tgamma(z), with Lanczos support:
//
template <class T, class L>
T gamma_imp(T z, const L& l)
{
   using namespace std;

   T result = 1;

#ifdef BOOST_MATH_INSTRUMENT
   static bool b = false;
   if(!b)
   {
      std::cout << "tgamma_imp called with " << typeid(z).name() << " " << typeid(l).name() << std::endl;
      b = true;
   }
#endif

   if((z <= 0) && (floor(z) == z))
      return tools::pole_error<T>(BOOST_CURRENT_FUNCTION, "Evaluation of tgamma at a negative integer.");
   if(z <= -20)
   {
      result = gamma_imp(-z, l) * sinpx(z);
      if((fabs(result) < 1) && (tools::max_value<T>() * fabs(result) < boost::math::constants::pi<T>()))
         return tools::overflow_error<T>(BOOST_CURRENT_FUNCTION, "Result of tgamma is too large to represent.");
      result = -boost::math::constants::pi<T>() / result;
      if(result == 0)
         return tools::underflow_error<T>(BOOST_CURRENT_FUNCTION, "Result of tgamma is too small to represent.");
      if(boost::math::fpclassify(result) == FP_SUBNORMAL)
         return tools::denorm_error<T>(result, BOOST_CURRENT_FUNCTION, "Result of tgamma is denormalized.");
      return result;
   }

   // shift z to > 1:
   while(z < 1)
   {
      result /= z;
      z += 1;
   }
   result *= L::lanczos_sum(z);
   if(z * log(z) > tools::log_max_value<T>())
   {
      // we're going to overflow unless this is done with care:
      T zgh = (z + L::g() - boost::math::constants::half<T>());
      T hp = pow(zgh, (z / 2) - T(0.25));
      result *= hp / exp(zgh);
      if(tools::max_value<T>() / hp < result)
         return tools::overflow_error<T>(BOOST_CURRENT_FUNCTION, "Result of tgamma is too large to represent.");
      result *= hp;
   }
   else
   {
      T zgh = (z + L::g() - boost::math::constants::half<T>());
      result *= pow(zgh, z - boost::math::constants::half<T>()) / exp(zgh);
   }
   return result;
}
//
// lgamma(z) with Lanczos support:
//
template <class T, class L>
T lgamma_imp(T z, const L& l, int* sign = 0)
{
#ifdef BOOST_MATH_INSTRUMENT
   static bool b = false;
   if(!b)
   {
      std::cout << "lgamma_imp called with " << typeid(z).name() << " " << typeid(l).name() << std::endl;
      b = true;
   }
#endif

   using namespace std;

   T result = 0;
   int sresult = 1;
   if(z <= 0)
   {
      // reflection formula:
      if(floor(z) == z)
         return tools::pole_error<T>(BOOST_CURRENT_FUNCTION, "Evaluation of lgamma at a negative integer.");

      T t = sinpx(z);
      z = -z;
      if(t < 0)
      {
         t = -t;
      }
      else
      {
         sresult = -sresult;
      }
      result = log(boost::math::constants::pi<T>()) - lgamma_imp(z, l) - log(t);
   }
   else if(z < 0.5)
   {
      // taking the log of tgamma reduces the error, no danger of overflow here:
      result = log(gamma_imp(z, l));
   }
   else if((z >= 3) && (z < 100))
   {
      // taking the log of tgamma reduces the error, no danger of overflow here:
      result = log(gamma_imp(z, l));
   }
   else if(z > 100)
   {
      // regular evaluation:
      T zgh = (z + L::g() - boost::math::constants::half<T>()) / boost::math::constants::e<T>();
      T l = L::lanczos_sum_expG_scaled(z);
      T p = z - boost::math::constants::half<T>();
      l *= zgh;
      p -= 1;
      result = log(l) + p * log(zgh);
   }
   else if(z >= 1.5)
   {
      // special case near 2:
      T dz = z - 2;
      result = dz * log((z + L::g() - T(0.5)) / boost::math::constants::e<T>());
      result += boost::math::log1p(dz / (L::g() + T(1.5))) * T(1.5);
      result += boost::math::log1p(L::lanczos_sum_near_2(z));
   }
   else
   {
      // special case near 1:
      T dz = z - 1;
      result = dz * log((z + L::g() - T(0.5)) / boost::math::constants::e<T>());
      result += boost::math::log1p(dz / (L::g() + T(0.5))) / 2;
      result += boost::math::log1p(L::lanczos_sum_near_1(dz));
   }
   if(sign)
      *sign = sresult;
   return result;
}

//
// Incomplete gamma functions follow:
//
template <class T>
struct upper_incomplete_gamma_fract
{
private:
   T z, a;
   int k;
public:
   typedef std::pair<T,T> result_type;

   upper_incomplete_gamma_fract(T a1, T z1)
      : z(z1-a1+1), a(a1), k(0)
   {
   }

   result_type operator()()
   {
      ++k;
      z += 2;
      return result_type(k * (a - k), z);
   }
};

template <class T>
T upper_gamma_fraction(T a, T z, int bits)
{
   // Multiply result by z^a * e^-z to get the full
   // upper incomplete integral.  Divide by tgamma(z)
   // to normalise.
   upper_incomplete_gamma_fract<T> f(a, z);
   return 1 / (z - a + 1 + boost::math::tools::continued_fraction_a(f, bits));
}

template <class T>
struct lower_incomplete_gamma_series
{
private:
   T a, z, result;
public:
   typedef T result_type;
   lower_incomplete_gamma_series(T a1, T z1) : a(a1), z(z1), result(1){}

   T operator()()
   {
      T r = result;
      a += 1;
      result *= z/a;
      return r;
   }
};

template <class T>
T lower_gamma_series(T a, T z, int bits)
{
   // Multiply result by ((z^a) * (e^-z) / a) to get the full
   // lower incomplete integral. Then divide by tgamma(a)
   // to get the normalised value.
   lower_incomplete_gamma_series<T> s(a, z);
   return boost::math::tools::sum_series(s, bits);
}

//
// Fully generic tgamma and lgamma use the incomplete partial
// sums added together:
//
template <class T>
T gamma_imp(T z, const lanczos::undefined_lanczos& l)
{
   using namespace std;
   if((z <= 0) && (floor(z) == z))
      return tools::pole_error<T>(BOOST_CURRENT_FUNCTION, "Evaluation of tgamma at a negative integer.");
   if(z <= -20)
   {
      T result = gamma_imp(-z, l) * sinpx(z);
      if((fabs(result) < 1) && (tools::max_value<T>() * fabs(result) < boost::math::constants::pi<T>()))
         return tools::overflow_error<T>(BOOST_CURRENT_FUNCTION, "Result of tgamma is too large to represent.");
      result = -boost::math::constants::pi<T>() / result;
      if(result == 0)
         return tools::underflow_error<T>(BOOST_CURRENT_FUNCTION, "Result of tgamma is too small to represent.");
      if(boost::math::fpclassify(result) == FP_SUBNORMAL)
         return tools::denorm_error<T>(result, BOOST_CURRENT_FUNCTION, "Result of tgamma is denormalized.");
      return result;
   }
   //
   // The upper gamma fraction is *very* slow for z < 6, actually it's very
   // slow to converge everywhere but recursing until z > 6 gets rid of the 
   // worst of it's behaviour.
   //
   T prefix = 1;
   while(z < 6)
   {
      prefix /= z;
      z += 1;
   }
   prefix = prefix * pow(z / boost::math::constants::e<T>(), z);
   T sum = detail::lower_gamma_series(z, z, ::boost::math::tools::digits<T>()) / z;
   sum += detail::upper_gamma_fraction(z, z, ::boost::math::tools::digits<T>());
   if(fabs(tools::max_value<T>() / prefix) < fabs(sum))
      return tools::overflow_error<T>(BOOST_CURRENT_FUNCTION, "Result of tgamma is too large to represent.");
   return sum * prefix;
}

template <class T>
T lgamma_imp(T z, const lanczos::undefined_lanczos&, int*sign)
{
   using namespace std;

   T result = 0;
   int sresult = 1;
   if(z <= 0)
   {
      if(floor(z) == z)
         return tools::pole_error<T>(BOOST_CURRENT_FUNCTION, "Evaluation of tgamma at a negative integer.");
      T t = detail::sinpx(z);
      z = -z;
      if(t < 0)
      {
         t = -t;
      }
      else
      {
         sresult = -sresult;
      }
      result = log(boost::math::constants::pi<T>()) - lgamma(z) - log(t);
   }
   else if((z != 1) && (z != 2))
   {
      T limit = (std::max)(z+1, T(10));
      T prefix = z * log(limit) - limit;
      T sum = detail::lower_gamma_series(z, limit, ::boost::math::tools::digits<T>()) / z;
      sum += detail::upper_gamma_fraction(z, limit, ::boost::math::tools::digits<T>());
      result = log(sum) + prefix;
   }
   if(sign)
      *sign = sresult;
   return result;
}
//
// This helper calculates tgamma(dz+1)-1 without cancellation errors,
// used by the upper incomplete gamma with z < 1:
//
template <class T, class L>
T tgammap1m1_imp(T dz, const L&)
{
   // 
   using namespace std;

   T zgh = (L::g() + T(0.5) + dz) / boost::math::constants::e<T>();
   T A = boost::math::powm1(zgh, dz);
   T B = boost::math::sqrtp1m1(dz / (L::g() + T(0.5)));
   T C = L::lanczos_sum_near_1(dz);
   T Ap1 = pow(zgh, dz);
   T Bp1 = sqrt(1 + (dz / (L::g() + T(0.5))) );

   return Bp1 * (A + C * Ap1) + B;
}

template <class T>
T tgammap1m1_imp(T dz, 
                 const ::boost::math::lanczos::undefined_lanczos& l)
{
   //
   // TODO FIXME!!!!
   // There should be a better solution than this, but the
   // algebra isn't easy for the general case....
   //
   T result = gamma_imp(1 + dz, l) - 1;
   if(std::pow(2.0, boost::math::tools::digits<T>()) * result 
      < std::pow(2.0, std::numeric_limits<long double>::digits))
   {
      // Cancellation errors mean that we have 
      // fewer good digits left in the result than
      // there are digits in a long double.  To limit
      // the damage, call the long double version:
      typedef boost::math::lanczos::lanczos_traits<long double>::evaluation_type eval_type;
      result = tgammap1m1_imp(boost::math::tools::real_cast<long double>(dz), eval_type());
   }
   return result;
}

//
// Series representation for upper fraction when z is small:
//
template <class T>
struct small_gamma2_series
{
   typedef T result_type;

   small_gamma2_series(T a_, T x_) : result(-x_), x(-x_), apn(a_+1), n(1){}

   T operator()()
   {
      T r = result / (apn);
      result *= x;
      result /= ++n;
      apn += 1;
      return r;
   }

private:
   T result, x, apn;
   int n;
};
//
// calculate power term prefix (z^a)(e^-z) used in the non-normalised
// incomplete gammas:
//
template <class T>
T full_igamma_prefix(T a, T z)
{
   using namespace std;

   T prefix;
   T alz = a * log(z);

   if(z >= 1)
   {
      if((alz < tools::log_max_value<T>()) && (-z > tools::log_min_value<T>()))
      {
         prefix = pow(z, a) * exp(-z);
      }
      else if(a >= 1)
      {
         prefix = pow(z / exp(z/a), a);
      }
      else
      {
         prefix = exp(alz - z);
      }
   }
   else
   {
      if(alz > tools::log_min_value<T>())
      {
         prefix = pow(z, a) * exp(-z);
      }
      else if(z/a < tools::log_max_value<T>())
      {
         prefix = pow(z / exp(z/a), a);
      }
      else
      {
         prefix = exp(alz - z);
      }
   }
   //
   // This error handling isn't very good: it happens after the fact
   // rather than before it...
   //
   if(boost::math::fpclassify(prefix) == FP_INFINITE)
      tools::overflow_error<T>(BOOST_CURRENT_FUNCTION, "Result of incomplete gamma function is too large to represent.");

   return prefix;
}
//
// Full lower integral:
//
template <class T>
T tgamma_lower_part(T a, T z)
{
   using namespace std;
   T prefix = full_igamma_prefix(a, z);
   T r2 = detail::lower_gamma_series(a, z, boost::math::tools::digits<T>()) / a;
   if(tools::max_value<T>() / prefix < r2)
      tools::overflow_error<T>(BOOST_CURRENT_FUNCTION, "Result of incomplete gamma function is too large to represent.");
   prefix *= r2;
   return prefix;
}
//
// Full upper integral:
//
template <class T>
T tgamma_upper_part(T a, T z)
{
   using namespace std;
   if((z < 1) && ((a <= 0.75) || (z < a)))
   {
      T result = tgammap1m1(a) - powm1(z, a);
      result /= a;
      detail::small_gamma2_series<T> s(a, z);
      result -= pow(z, a) * tools::sum_series(s, boost::math::tools::digits<T>());
      return result;
   }
   else
   {
      T prefix = full_igamma_prefix(a, z);
      return (prefix == 0) ? 0 : prefix * detail::upper_gamma_fraction(a, z, boost::math::tools::digits<T>());
   }
}
//
// Full lower integral, computing via lower integral as required:
//
template <class T>
T tgamma_imp(T a, T z)
{
   if((a <= 0) || (z < 0))
      tools::domain_error<T>(BOOST_CURRENT_FUNCTION, "Negative argument to the incomplete gamma function.");

   using namespace std;

   // TODO: Find a better heuristic here:
   if((a < 0.75) && (z < 1) && ((a*a < z) || (a < 0.125)))
   {
      return tgamma_upper_part(a, z);
   }
   else if((z < a+1) && (a > 1))
   {
      T result = ::boost::math::tgamma(a);
      if(std::numeric_limits<T>::is_specialized && (result >= (std::numeric_limits<T>::max)()))
         return result;
      return result - tgamma_lower_part(a, z);
   }
   else
   {
      return tgamma_upper_part(a, z);
   }
}
//
// Full lower integral, computing via upper integral as required:
//
template <class T>
T tgamma_lower_imp(T a, T z)
{
   if((a <= 0) || (z < 0))
      tools::domain_error<T>(BOOST_CURRENT_FUNCTION, "Negative argument to the incomplete gamma function.");

   if(z < (std::max)(a+1, T(10)))
   {
      return tgamma_lower_part(a, z);
   }
   else
   {
      T result = ::boost::math::tgamma(a);
      if(std::numeric_limits<T>::is_specialized && (result >= (std::numeric_limits<T>::max)()))
         return result;
      T r2 = tgamma_upper_part(a, z);
      return result - r2;
   }
}
//
// Helper to compute log(1+x)-x:
//
template <class T>
T log1pmx(T x)
{
   boost::math::detail::log1p_series<T> s(x);
   s();
   return boost::math::tools::sum_series(s, tools::digits<T>() + 2);
}
//
// Compute (z^a)(e^-z)/tgamma(a)
// most if the error occurs in this function:
//
template <class T, class L>
T regularised_gamma_prefix(T a, T z, const L& l)
{
   using namespace std;
   T agh = a + L::g() - T(0.5);
   T prefix;
   T d = ((z - a) - L::g() + T(0.5)) / agh;

   if(a < 1)
   {
      //
      // We have to treat a < 1 as a special case because our Lanczos
      // approximations are optimised against the factorials with a > 1,
      // and for high precision types especially (128-bit reals for example)
      // very small values of a can give rather eroneous results for gamma
      // unless we do this:
      //
      // TODO: is this still required?  Lanczos approx should be better now?
      //
      if(z <= tools::log_min_value<T>())
      {
         // Oh dear, have to use logs, should be free of cancellation errors though:
         return exp(a * log(z) - z - lgamma_imp(a, l));
      }
      else
      {
         // direct calculation, no danger of overflow as gamma(a) < 1/a 
         // for small a.
         return pow(z, a) * exp(-z) / gamma_imp(a, l);
      }
   }
   else if((fabs(d*d*a) <= 100) && (a > 150))
   {
      // special case for large a and a ~ z.
      prefix = a * log1pmx(d) + z * (0.5 - L::g()) / agh;
      prefix = exp(prefix);
   }
   else
   {
      //
      // general case.
      // direct computation is most accurate, but use various fallbacks
      // for different parts of the problem domain:
      //
      T alz = a * log(z / agh);
      T amz = a - z;
      if(((std::min)(alz, amz) <= tools::log_min_value<T>()) || ((std::max)(alz, amz) >= tools::log_max_value<T>()))
      {
         T amza = amz / a;
         if(((std::min)(alz, amz)/2 > tools::log_min_value<T>()) && ((std::max)(alz, amz)/2 < tools::log_max_value<T>()))
         {
            // compute square root of the result and then square it:
            T sq = pow(z / agh, a / 2) * exp(amz / 2);
            prefix = sq * sq;
         }
         else if(((std::min)(alz, amz)/4 > tools::log_min_value<T>()) && ((std::max)(alz, amz)/4 < tools::log_max_value<T>()) && (z > a))
         {
            // compute the 4th root of the result then square it twice:
            T sq = pow(z / agh, a / 4) * exp(amz / 4);
            prefix = sq * sq;
            prefix *= prefix;
         }
         else if((amza > tools::log_min_value<T>()) && (amza < tools::log_max_value<T>()))
         {
            prefix = pow((z * exp(amza)) / agh, a);
         }
         else
         {
            prefix = exp(alz + amz);
         }
      }
      else
      {
         prefix = pow(z / agh, a) * exp(amz);
      }
   }
   prefix *= sqrt(agh / boost::math::constants::e<T>()) / L::lanczos_sum_expG_scaled(a);
   return prefix;
}
//
// And again, without Lanczos support:
//
template <class T>
T regularised_gamma_prefix(T a, T z, const lanczos::undefined_lanczos&)
{
   using namespace std;

   T limit = (std::max)(T(10), a);
   T sum = detail::lower_gamma_series(a, limit, ::boost::math::tools::digits<T>()) / a;
   sum += detail::upper_gamma_fraction(a, limit, ::boost::math::tools::digits<T>());

   if(a < 10)
   {
      // special case for small a:
      T prefix = pow(z / 10, a);
      prefix *= exp(10-z);
      if(0 == prefix)
      {
         prefix = pow((z * exp((10-z)/a)) / 10, a);
      }
      prefix /= sum;
      return prefix;
   }

   T zoa = z / a;
   T amz = a - z;
   T alzoa = a * log(zoa);
   T prefix;
   if(((std::min)(alzoa, amz) <= tools::log_min_value<T>()) || ((std::max)(alzoa, amz) >= tools::log_max_value<T>()))
   {
      T amza = amz / a;
      if((amza <= tools::log_min_value<T>()) || (amza >= tools::log_max_value<T>()))
      {
         prefix = exp(alzoa + amz);
      }
      else
      {
         prefix = pow(zoa * exp(amza), a);
      }
   }
   else
   {
      prefix = pow(zoa, a) * exp(amz);
   }
   prefix /= sum;
   return prefix;
}

//
// regularised incomplete gamma, divide through by Lanczos 
// approximation and combine terms:
//
template <class T, class L>
T gamma_Q_imp(T a, T z, const L& l)
{
   if((a <= 0) || (z < 0))
      tools::domain_error<T>(BOOST_CURRENT_FUNCTION, "Negative argument to the incomplete gamma function.");

   using namespace std;

   if((a < 0.75) && (z < 1))
   {
      T tga = gamma_imp(a, l);
      T result = tgamma_upper_part(a, z);
      return result / tga;
   }
   else 
   {
      T prefix = regularised_gamma_prefix(a, z, l);
      if(z < a+1)
      {
         T result = 1;
         T r2 = prefix * detail::lower_gamma_series(a, z, boost::math::tools::digits<T>()) / a;
         return result - r2;
      }
      else
         return (prefix == 0) ? 0 : prefix * detail::upper_gamma_fraction(a, z, boost::math::tools::digits<T>());
   }
}
template <class T, class L>
T gamma_P_imp(T a, T z, const L& l)
{
   if((a <= 0) || (z < 0))
      tools::domain_error<T>(BOOST_CURRENT_FUNCTION, "Negative argument to the incomplete gamma function.");

   using namespace std;

   T prefix = regularised_gamma_prefix(a, z, l);
   if(z < a+1)
   {
      return (prefix == 0) ? 0 : prefix * detail::lower_gamma_series(a, z, boost::math::tools::digits<T>()) / a;
   }
   else
   {
      T result = 1;
      T r2 = prefix * detail::upper_gamma_fraction(a, z, boost::math::tools::digits<T>());
      return result - r2;
   }
}

template <class T, class L>
T tgamma_delta_ratio_imp(T z, T delta, const L&)
{
   using namespace std;

   if((z <= 0) || (z + delta <= 0))
      tools::domain_error<T>(BOOST_CURRENT_FUNCTION, "Gamma function ratios only implemented for positive arguments.");

   T zgh = z + L::g() - constants::half<T>();
   T result;
   if(fabs(delta) < 10)
   {
      result = exp((constants::half<T>() - z) * boost::math::log1p(delta / zgh));
   }
   else
   {
      result = pow(zgh / (zgh + delta), z - constants::half<T>());
   }
   result *= pow(constants::e<T>() / (zgh + delta), delta);
   result *= L::lanczos_sum(z) / L::lanczos_sum(z + delta);
   return result;
}

template <class T>
T tgamma_delta_ratio_imp(T z, T delta, const lanczos::undefined_lanczos&)
{
   using namespace std;

   if((z <= 0) || (z + delta <= 0))
      tools::domain_error<T>(BOOST_CURRENT_FUNCTION, "Gamma function ratios only implemented for positive arguments.");

   //
   // The upper gamma fraction is *very* slow for z < 6, actually it's very
   // slow to converge everywhere but recursing until z > 6 gets rid of the 
   // worst of it's behaviour.
   //
   T prefix = 1;
   T zd = z + delta;
   while((zd < 6) && (z < 6))
   {
      prefix /= z;
      prefix *= zd;
      z += 1;
      zd += 1;
   }
   if(delta < 10)
   {
      prefix *= exp(-z * boost::math::log1p(delta / z));
   }
   else
   {
      prefix *= pow(z / zd, z);
   }
   prefix *= pow(constants::e<T>() / zd, delta);
   T sum = detail::lower_gamma_series(z, z, ::boost::math::tools::digits<T>()) / z;
   sum += detail::upper_gamma_fraction(z, z, ::boost::math::tools::digits<T>());
   T sumd = detail::lower_gamma_series(zd, zd, ::boost::math::tools::digits<T>()) / zd;
   sumd += detail::upper_gamma_fraction(zd, zd, ::boost::math::tools::digits<T>());
   sum /= sumd;
   if(fabs(tools::max_value<T>() / prefix) < fabs(sum))
      return tools::overflow_error<T>(BOOST_CURRENT_FUNCTION, "Result of tgamma is too large to represent.");
   return sum * prefix;
}

} // namespace detail

template <class T>
inline T tgamma(T z)
{
   typedef typename lanczos::lanczos_traits<T>::value_type value_type;
   typedef typename lanczos::lanczos_traits<T>::evaluation_type evaluation_type;
   return tools::checked_narrowing_cast<T>(detail::gamma_imp(static_cast<value_type>(z), evaluation_type()), BOOST_CURRENT_FUNCTION);
}

template <class T>
inline T lgamma(T z, int* sign)
{
   typedef typename lanczos::lanczos_traits<T>::value_type value_type;
   typedef typename lanczos::lanczos_traits<T>::evaluation_type evaluation_type;
   return tools::checked_narrowing_cast<T>(detail::lgamma_imp(static_cast<value_type>(z), evaluation_type(), sign), BOOST_CURRENT_FUNCTION);
}

template <class T>
inline T lgamma(T x)
{
   return ::boost::math::lgamma(x, 0);
}

template <class T>
inline T tgammap1m1(T z)
{
   typedef typename lanczos::lanczos_traits<T>::value_type value_type;
   typedef typename lanczos::lanczos_traits<T>::evaluation_type evaluation_type;
   return tools::checked_narrowing_cast<T>(detail::tgammap1m1_imp(static_cast<value_type>(z), evaluation_type()), BOOST_CURRENT_FUNCTION);
}

//
// Full upper incomplete gamma:
//
template <class T>
inline T tgamma(T a, T z)
{
   typedef typename boost::math::tools::evaluation<T>::type eval_type;
   return tools::checked_narrowing_cast<T>(detail::tgamma_imp(static_cast<eval_type>(a), static_cast<eval_type>(z)), BOOST_CURRENT_FUNCTION);
}

template <class T>
inline T tgamma_lower(T a, T z)
{
   typedef typename boost::math::tools::evaluation<T>::type eval_type;
   return tools::checked_narrowing_cast<T>(detail::tgamma_lower_imp(static_cast<eval_type>(a), static_cast<eval_type>(z)), BOOST_CURRENT_FUNCTION);
}

template <class T>
inline T gamma_Q(T a, T z)
{
   typedef typename lanczos::lanczos_traits<T>::value_type value_type;
   typedef typename lanczos::lanczos_traits<T>::evaluation_type evaluation_type;
   return tools::checked_narrowing_cast<T>(
      detail::gamma_Q_imp(static_cast<value_type>(a), 
      static_cast<value_type>(z), 
      evaluation_type()), BOOST_CURRENT_FUNCTION);
}

template <class T>
inline T gamma_P(T a, T z)
{
   typedef typename lanczos::lanczos_traits<T>::value_type value_type;
   typedef typename lanczos::lanczos_traits<T>::evaluation_type evaluation_type;
   return tools::checked_narrowing_cast<T>(
      detail::gamma_P_imp(static_cast<value_type>(a), 
      static_cast<value_type>(z), 
      evaluation_type()), BOOST_CURRENT_FUNCTION);
}

// ratios of gamma functions:
template <class T>
inline T tgamma_delta_ratio(T z, T delta)
{
   typedef typename lanczos::lanczos_traits<T>::value_type value_type;
   typedef typename lanczos::lanczos_traits<T>::evaluation_type evaluation_type;
   return tools::checked_narrowing_cast<T>(detail::tgamma_delta_ratio_imp(static_cast<value_type>(z), static_cast<value_type>(delta), evaluation_type()), BOOST_CURRENT_FUNCTION);
}
template <class T>
inline T tgamma_ratio(T a, T b)
{
   typedef typename lanczos::lanczos_traits<T>::value_type value_type;
   typedef typename lanczos::lanczos_traits<T>::evaluation_type evaluation_type;
   return tools::checked_narrowing_cast<T>(detail::tgamma_delta_ratio_imp(static_cast<value_type>(a), static_cast<value_type>(b - a), evaluation_type()), BOOST_CURRENT_FUNCTION);
}

} // namespace math
} // namespace boost

#include <boost/math/special_functions/igamma_inverse.hpp>

#endif // BOOST_MATH_SF_GAMMA_HPP


