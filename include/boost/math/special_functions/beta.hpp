//  (C) Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_SPECIAL_BETA_HPP
#define BOOST_MATH_SPECIAL_BETA_HPP

#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/special_functions/erf.hpp>
#include <boost/math/special_functions/log1p.hpp>
#include <boost/math/special_functions/expm1.hpp>
#include <boost/math/tools/roots.hpp>

namespace boost{ namespace math{

template <class T>
T beta(T a, T b);
template <class T>
T ibeta(T a, T b, T x);
template <class T>
T ibetac(T a, T b, T x);

namespace detail{

//
// Implementation of Beta(a,b) using the Lanczos approximation:
//
template <class T, class L>
T beta_imp(T a, T b, const L& l)
{
   using namespace std;

   if((a <= 0) || (b <= 0))
      tools::domain_error<T>(BOOST_CURRENT_FUNCTION, "Negative argument to the beta function.");
   
   T result;

   T prefix = 1;
   T c = a + b;

   // special cases:
   if(c == a)
      return boost::math::tgamma(b);
   else if(c == b)
      return boost::math::tgamma(a);
   if(b == 1)
      return 1/a;
   else if(a == 1)
      return 1/b;

   // if a or b are less than 1, shift to greater than 1:
   if(a < 1)
   {
      prefix *= c / a;
      c += 1;
      a += 1;
   }
   if(b < 1)
   {
      prefix *= c / b;
      c += 1;
      b += 1;
   }
   if(a < b)
      std::swap(a, b);

   // Lanczos calculation:
   T agh = a + L::g() - T(0.5);
   T bgh = b + L::g() - T(0.5);
   T cgh = c + L::g() - T(0.5);
   result = L::lanczos_sum_expG_scaled(a) * L::lanczos_sum_expG_scaled(b) / L::lanczos_sum_expG_scaled(c);
   T ambh = a - T(0.5) - b;
   if((fabs(b * ambh) < (cgh * 100)) && (a > 100))
   {
      // special case where the base of the power term is close to 1
      // compute (1+x)^y instead:
      result *= exp(ambh * boost::math::log1p(-b / cgh));
   }
   else
   {
      result *= pow(agh / cgh, a - T(0.5) - b);
   }
   result *= pow((agh * bgh) / (cgh * cgh), b);
   result *= sqrt(boost::math::constants::e<T>() / bgh);

   // if a and b were originally less than 1 we need to scale the result:
   result *= prefix;

   return result;
}

//
// Generic implementation of Beta(a,b) without Lanczos approximation support
// (Caution this is slow!!!):
//
template <class T>
T beta_imp(T a, T b, const lanczos::undefined_lanczos& l)
{
   using namespace std;

   if((a <= 0) || (b <= 0))
      tools::domain_error<T>(BOOST_CURRENT_FUNCTION, "Negative argument to the beta function.");

   T result;

   T prefix = 1;
   T c = a + b;

   // special cases:
   if(c == a)
      return boost::math::tgamma(b);
   else if(c == b)
      return boost::math::tgamma(a);
   if(b == 1)
      return 1/a;
   else if(a == 1)
      return 1/b;

   // shift to a and b > 1 if required:
   if(a < 1)
   {
      prefix *= c / a;
      c += 1;
      a += 1;
   }
   if(b < 1)
   {
      prefix *= c / b;
      c += 1;
      b += 1;
   }
   if(a < b)
      std::swap(a, b);

   // set integration limits:
   T la = (std::max)(T(10), a);
   T lb = (std::max)(T(10), b);
   T lc = (std::max)(T(10), a+b);

   // calculate the fraction parts:
   T sa = detail::lower_gamma_series(a, la, ::boost::math::tools::digits(a)) / a;
   sa += detail::upper_gamma_fraction(a, la, ::boost::math::tools::digits(a));
   T sb = detail::lower_gamma_series(b, lb, ::boost::math::tools::digits(b)) / b;
   sb += detail::upper_gamma_fraction(b, lb, ::boost::math::tools::digits(b));
   T sc = detail::lower_gamma_series(c, lc, ::boost::math::tools::digits(c)) / c;
   sc += detail::upper_gamma_fraction(c, lc, ::boost::math::tools::digits(c));

   // and the exponent part:
   result = exp(lc - la - lb) * pow(la/lc, a) * pow(lb/lc, b);
   
   // and combine:
   result *= sa * sb / sc;

   // if a and b were originally less than 1 we need to scale the result:
   result *= prefix;

   return result;
}

//
// Series approximation to the incomplete beta:
//
template <class T>
struct ibeta_series_t
{
   typedef T result_type;
   ibeta_series_t(T a_, T b_, T x_, T mult) : result(mult), x(x_), apn(a_), poch(1-b_), n(1) {}
   T operator()()
   {
      T r = result / apn;
      apn += 1;
      result *= poch * x / n;
      ++n;
      poch += 1;
      return r;
   }
private:
   T result, x, apn, poch;
   int n;
};

template <class T, class L>
T ibeta_series(T a, T b, T x, T s0, const L&, bool normalised)
{
   using namespace std;

   T result;

   if(normalised)
   {
      T prefix = 1;
      T c = a + b;

      // incomplete beta power term, combined with the Lanczos approximation:
      T agh = a + L::g() - T(0.5);
      T bgh = b + L::g() - T(0.5);
      T cgh = c + L::g() - T(0.5);
      result = L::lanczos_sum_expG_scaled(c) / (L::lanczos_sum_expG_scaled(a) * L::lanczos_sum_expG_scaled(b));
      result *= pow(cgh / bgh, b - 0.5);
      result *= pow(x * cgh / agh, a);
      result *= sqrt(agh / boost::math::constants::e<T>());
      result *= prefix;
   }
   else
   {
      // Non-normalised, just compute the power:
      result = pow(x, a);
   }
   ibeta_series_t<T> s(a, b, x, result);
   return boost::math::tools::sum_series(s, boost::math::tools::digits(x), s0);
}
//
// Incomplete Beta series again, this time without Lanczos support:
//
template <class T>
T ibeta_series(T a, T b, T x, T s0, const boost::math::lanczos::undefined_lanczos&, bool normalised)
{
   using namespace std;

   T result;

   if(normalised)
   {
      T prefix = 1;
      T c = a + b;

      // figure out integration limits for the gamma function:
      T la = (std::max)(T(10), a);
      T lb = (std::max)(T(10), b);
      T lc = (std::max)(T(10), a+b);

      // calculate the gamma parts:
      T sa = detail::lower_gamma_series(a, la, ::boost::math::tools::digits(a)) / a;
      sa += detail::upper_gamma_fraction(a, la, ::boost::math::tools::digits(a));
      T sb = detail::lower_gamma_series(b, lb, ::boost::math::tools::digits(b)) / b;
      sb += detail::upper_gamma_fraction(b, lb, ::boost::math::tools::digits(b));
      T sc = detail::lower_gamma_series(c, lc, ::boost::math::tools::digits(c)) / c;
      sc += detail::upper_gamma_fraction(c, lc, ::boost::math::tools::digits(c));

      // and their combined power-terms:
      T b1 = (x * lc) / la;
      T b2 = lc/lb;
      T e1 = lc - la - lb;
      T lb1 = a * log(b1);
      T lb2 = b * log(b2);

      if((lb1 >= tools::log_max_value(b1)) 
         || (lb1 <= tools::log_min_value(b1))
         || (lb2 >= tools::log_max_value(b1)) 
         || (lb2 <= tools::log_min_value(b1))
         || (e1 >= tools::log_max_value(b1)) 
         || (e1 <= tools::log_min_value(b1)) )
      {
         T p = lb1 + lb2 - e1;
         result = exp(p);
      }
      else
      {
         result = pow(b1, a) * pow(b2, b) / exp(e1);
      }
      // and combine the results:
      result /= sa * sb / sc;
   }
   else
   {
      // non-normalised, just compute the power:
      result = pow(x, a);
   }
   ibeta_series_t<T> s(a, b, x, result);
   return boost::math::tools::sum_series(s, boost::math::tools::digits(x), s0);
}
//
// Compute the leading power terms in the incomplete Beta:
//
// (x^a)(y^b)/Beta(a,b) when normalised, and
// (x^a)(y^b) otherwise.
//
// Almost all of the error in the incomplete beta comes from this
// function: particularly when a and b are large. Computing large
// powers are *hard* though, and using logarithms just leads to
// horrendous cancellation errors.
//
template <class T, class L>
T ibeta_power_terms(T a, 
                        T b, 
                        T x, 
                        T y, 
                        const L&,
                        bool normalised)
{
   using namespace std;

   if(!normalised)
   {
      // can we do better here?
      return pow(x, a) * pow(y, b);
   }

   T result;

   T prefix = 1;
   T c = a + b;

   // combine power terms with Lanczos approximation:
   T agh = a + L::g() - T(0.5);
   T bgh = b + L::g() - T(0.5);
   T cgh = c + L::g() - T(0.5);
   result = L::lanczos_sum_expG_scaled(c) / (L::lanczos_sum_expG_scaled(a) * L::lanczos_sum_expG_scaled(b));

   // l1 and l2 are the base of the exponents minus one:
   T l1 = (x * b - y * agh) / agh;
   T l2 = (y * a - x * bgh) / bgh;
   if(((std::min)(a, b) > 5) 
      && ((std::min)(fabs(l1), fabs(l2)) < 0.2))
   {
      // when the base of the exponent is very near 1 we get really
      // gross errors unless extra care is taken:
      if((std::max)(fabs(l1), fabs(l2)) < 0.5)
      {
         // both bases near 1:
         if(a < b)
         {
            T l = boost::math::expm1((b/a) * boost::math::log1p(l2));
            l = l1 + l + l * l1;
            l = a * boost::math::log1p(l);
            result *= exp(l);
         }
         else
         {
            T l = boost::math::expm1((a/b) * boost::math::log1p(l1));
            l = l2 + l + l * l2;
            l = b * boost::math::log1p(l);
            result *= exp(l);
         }
      }
      else if(fabs(l1) < fabs(l2))
      {
         // first base near 1:
         T l = a * boost::math::log1p(l1) 
            + b * log((y * cgh) / bgh);
         result *= exp(l);
      }
      else
      {
         // second base near 1:
         T l = b * boost::math::log1p(l2) 
            + a * log((x * cgh) / agh);
         result *= exp(l);
      }
   }
   else
   {
      // general case:
      T b1 = (x * cgh) / agh;
      T b2 = (y * cgh) / bgh;
      T l1 = a * log(b1);
      T l2 = b * log(b2);
      if((l1 >= tools::log_max_value(l1))
         || (l1 <= tools::log_min_value(l1)) 
         || (l2 >= tools::log_max_value(l1))
         || (l2 <= tools::log_min_value(l1)) 
         )
      {
         // Oops, overflow, sidestep:
         if(a < b)
            result *= pow(pow(b2, b/a) * b1, a);
         else
            result *= pow(pow(b1, a/b) * b2, b);
      }
      else
      {
         // finally the normal case:
         result *= pow(b1, a) * pow(b2, b);
      }
   }
   // combine with the leftover terms from the Lanczos approximation:
   result *= sqrt(bgh / boost::math::constants::e<T>());
   result *= sqrt(agh / cgh);
   result *= prefix;

   return result;
}
//
// Compute the leading power terms in the incomplete Beta:
//
// (x^a)(y^b)/Beta(a,b) when normalised, and
// (x^a)(y^b) otherwise.
//
// Almost all of the error in the incomplete beta comes from this
// function: particularly when a and b are large. Computing large
// powers are *hard* though, and using logarithms just leads to
// horrendous cancellation errors.
//
// This version is generic, slow, and does not use the Lanczos approximation.
//
template <class T>
T ibeta_power_terms(T a, 
                        T b, 
                        T x, 
                        T y, 
                        const boost::math::lanczos::undefined_lanczos&,
                        bool normalised)
{
   using namespace std;

   if(!normalised)
   {
      return pow(x, a) * pow(y, b);
   }

   T result;

   T prefix = 1;
   T c = a + b;

   // integration limits for the gamma functions:
   T la = (std::max)(T(10), a);
   T lb = (std::max)(T(10), b);
   T lc = (std::max)(T(10), a+b);
   // gamma function partials:
   T sa = detail::lower_gamma_series(a, la, ::boost::math::tools::digits(a)) / a;
   sa += detail::upper_gamma_fraction(a, la, ::boost::math::tools::digits(a));
   T sb = detail::lower_gamma_series(b, lb, ::boost::math::tools::digits(b)) / b;
   sb += detail::upper_gamma_fraction(b, lb, ::boost::math::tools::digits(b));
   T sc = detail::lower_gamma_series(c, lc, ::boost::math::tools::digits(c)) / c;
   sc += detail::upper_gamma_fraction(c, lc, ::boost::math::tools::digits(c));
   // gamma function powers combined with incomplete beta powers:

   T b1 = (x * lc) / la;
   T b2 = (y * lc) / lb;
   T e1 = lc - la - lb;
   T lb1 = a * log(b1);
   T lb2 = b * log(b2);

   if((lb1 >= tools::log_max_value(lb1))
      || (lb1 <= tools::log_min_value(lb1))
      || (lb2 >= tools::log_max_value(lb1))
      || (lb2 <= tools::log_min_value(lb1))
      || (e1 >= tools::log_max_value(lb1))
      || (e1 <= tools::log_min_value(lb1))
      )
   {
      result = exp(lb1 + lb2 - e1);
   }
   else
   {
      T p1 = pow(b1, a);
      T p2 = pow(b2, b);
      T p3 = exp(e1);
      result = p1 * p2 / p3;
   }
   // and combine with the remaining gamma function components:
   result /= sa * sb / sc;

   return result;
}
//
// Continued fraction for the incomplete beta:
//
template <class T>
struct ibeta_fraction2_t
{
   typedef std::pair<T, T> result_type;

   ibeta_fraction2_t(T a_, T b_, T x_) : a(a_), b(b_), x(x_), m(0) {}

   result_type operator()()
   {
      T aN = (a + m - 1) * (a + b + m - 1) * m * (b - m) * x * x;
      T denom = (a + 2 * m - 1);
      aN /= denom * denom;

      T bN = m;
      bN += (m * (b - m) * x) / (a + 2*m - 1);
      bN += ((a + m) * (a - (a + b) * x + 1 + m *(2 - x))) / (a + 2*m + 1);

      ++m;
      
      return std::make_pair(aN, bN);
   }

private:
   T a, b, x;
   int m;
};
//
// evaluate the incomplete beta via the continued fraction representation:
//
template <class T, class L>
T ibeta_fraction2(T a, T b, T x, T y, const L& l, bool normalised)
{
   using namespace std;
   T result = ibeta_power_terms(a, b, x, y, l, normalised);
   if(result == 0)
      return result;

   ibeta_fraction2_t<T> f(a, b, x);
   T fract = boost::math::tools::continued_fraction_b(f, boost::math::tools::digits(x));
   return result / fract;
}
//
// Computes the difference between ibeta(a,b,x) and ibeta(a+k,b,x):
//
template <class T, class L>
T ibeta_a_step(T a, T b, T x, T y, int k, const L& l, bool normalised)
{
   T prefix = ibeta_power_terms(a, b, x, y, l, normalised)/a;
   if(prefix == 0) 
      return prefix;
   T sum = 1;
   T term = 1;
   // series summation from 0 to k-1:
   for(int i = 0; i < k-1; ++i)
   {
      term *= (a+b+i) * x / (a+i+1);
      sum += term;
   }
   prefix *= sum;

   return prefix;
}
//
// This function is only needed for the non-regular incomplete beta,
// it computes the delta in:
// beta(a,b,x) = prefix + delta * beta(a+k,b,x)
// it is currently only called for small k.
//
template <class T>
T rising_factorial_ratio(T a, T b, int k)
{
   // calculate:
   // (a)(a+1)(a+2)...(a+k-1)
   // _______________________
   // (b)(b+1)(b+2)...(b+k-1)

   // This is only called with small k, for large k
   // it is grossly inefficient, do not use outside it's
   // intended purpose!!!
   if(k == 0)
      return 1;
   T result = 1;
   for(int i = 0; i < k; ++i)
      result *= (a+i) / (b+i);
   return result;
}
//
// Routine for a > 15, b < 1
//
template <class T, class L>
T beta_small_b_large_a_series(T a, T b, T x, T y, T s0, T mult, const L& l, bool normalised)
{
   using namespace std;
   //
   // This is DiDonato and Morris's BGRAT routine, see Eq's 9 through 9.6.
   //
   // Some values we'll need later, these are Eq 9.1:
   //
   T bm1 = b - 1;
   T t = a + bm1 / 2;
   T lx, u;
   if(y < 0.35)
      lx = boost::math::log1p(-y);
   else
      lx = log(x);
   u = -t * lx;
   // and from from 9.2:
   T prefix;
   T h = regularised_gamma_prefix(b, u, l);
   if(h == 0)
      return s0;
   if(normalised)
   {
      prefix = h / tgamma_delta_ratio_imp(a, b, l);
      prefix /= pow(t, b);
   }
   else
   {
      prefix = full_igamma_prefix(b, u) / pow(t, b);
   }
   prefix *= mult;
   //
   // now we need the quantity Pn, unfortunatately this is computed
   // recursively, and requires a full history of all the previous values
   // so no choice but to declare a big table and hope it's big enough...
   //
   T p[30] = { 1 };  // see 9.3.
   //
   // Now an initial value for J, see 9.6:
   //
   T j = gamma_Q(b, u) / h;
   //
   // Now we can start to pull things together and evaluate the sum in Eq 9:
   //
   T sum = s0 + prefix * j;  // Value at N = 0
   // some variables we'll need:
   unsigned tnp1 = 1; // 2*N+1 
   T lx2 = lx / 2;
   lx2 *= lx2;
   T lxp = 1;
   T t4 = 4 * t * t;
   T b2n = b; 

   for(unsigned n = 1; n < 30; ++n)
   {
      //
      // begin by evaluating the next Pn from Eq 9.4:
      //
      tnp1 += 2;
      p[n] = 0;
      T mbn = b - n;
      unsigned tmp1 = 3;
      for(unsigned m = 1; m < n; ++m)
      {
         mbn = m * b - n;
         p[n] += mbn * p[n-m] / unchecked_factorial<T>(tmp1);
         tmp1 += 2;
      }
      p[n] /= n;
      p[n] += bm1 / unchecked_factorial<T>(tnp1);
      //
      // Now we want Jn from Jn-1 using Eq 9.6:
      //
      j = (b2n * (b2n + 1) * j + (u + b2n + 1) * lxp) / t4;
      lxp *= lx2;
      b2n += 2;
      //
      // pull it together with Eq 9:
      //
      T r = prefix * p[n] * j;
      sum += r;
      if(fabs(r) < tools::epsilon(r) * sum)
         break;
   }
   return sum;
}

//
// The incomplete beta function implementation:
// This is just a big bunch of spagetti code to divide up the
// input range and select the right implementation method for
// each domain:
//
template <class T, class L>
T ibeta_imp(T a, T b, T x, const L& l, bool inv, bool normalised)
{
   bool invert = inv;
   T fract;
   T y = 1 - x;

   if((a <= 0) || (b <= 0))
      tools::domain_error<T>(BOOST_CURRENT_FUNCTION, "Negative argument to the incomplete beta function.");
   if((x < 0) || (x > 1))
      tools::domain_error<T>(BOOST_CURRENT_FUNCTION, "Parameter x outside the range [0,1] in the incomplete beta function.");

   if(x == 0)
      return (invert ? (normalised ? 1 : beta_imp(a, b, l)) : 0);
   if(x == 1)
      return (invert == 0 ? (normalised ? 1 : beta_imp(a, b, l)) : 0);

   if((std::min)(a, b) <= 1)
   {
      if(x > 0.5)
      {
         std::swap(a, b);
         std::swap(x, y);
         invert = !invert;
      }
      if((std::max)(a, b) <= 1)
      {
         // Both a,b < 1:
         if((a >= (std::min)(T(0.2), b)) || (pow(x, a) <= 0.9))
         {
            if(!invert)
               fract = ibeta_series(a, b, x, T(0), l, normalised);
            else
            {
               fract = -(normalised ? 1 : beta_imp(a, b, l));
               invert = false;
               fract = -ibeta_series(a, b, x, fract, l, normalised);
            }
         }
         else 
         {
            std::swap(a, b);
            std::swap(x, y);
            invert = !invert;
            if(y >= 0.3)
            {
               if(!invert)
                  fract = ibeta_series(a, b, x, T(0), l, normalised);
               else
               {
                  fract = -(normalised ? 1 : beta_imp(a, b, l));
                  invert = false;
                  fract = -ibeta_series(a, b, x, fract, l, normalised);
               }
            }
            else
            {
               // sidestep on a, and then use the series representation:
               T prefix;
               if(!normalised)
               {
                  prefix = rising_factorial_ratio(a+b, a, 20);
               }
               else
               {
                  prefix = 1;
               }
               fract = ibeta_a_step(a, b, x, y, 20, l, normalised);
               if(!invert)
                  fract = beta_small_b_large_a_series(a + 20, b, x, y, fract, prefix, l, normalised);
               else
               {
                  fract -= (normalised ? 1 : beta_imp(a, b, l));
                  invert = false;
                  fract = -beta_small_b_large_a_series(a + 20, b, x, y, fract, prefix, l, normalised);
               }
            }
         }
      }
      else
      {
         // One of a, b < 1 only:
         if((b <= 1) || ((x < 0.1) && (pow(b * x, a) <= 0.7)))
         {
            if(!invert)
               fract = ibeta_series(a, b, x, T(0), l, normalised);
            else
            {
               fract = -(normalised ? 1 : beta_imp(a, b, l));
               invert = false;
               fract = -ibeta_series(a, b, x, fract, l, normalised);
            }
         }
         else 
         {
            std::swap(a, b);
            std::swap(x, y);
            invert = !invert;

            if(y >= 0.3)
            {
               if(!invert)
                  fract = ibeta_series(a, b, x, T(0), l, normalised);
               else
               {
                  fract = -(normalised ? 1 : beta_imp(a, b, l));
                  invert = false;
                  fract = -ibeta_series(a, b, x, fract, l, normalised);
               }
            }
            else if(a >= 15)
            {
               if(!invert)
                  fract = beta_small_b_large_a_series(a, b, x, y, T(0), T(1), l, normalised);
               else
               {
                  fract = -(normalised ? 1 : beta_imp(a, b, l));
                  invert = false;
                  fract = -beta_small_b_large_a_series(a, b, x, y, fract, T(1), l, normalised);
               }
            }
            else
            {
               // sidestep to improve errors:
               T prefix;
               if(!normalised)
               {
                  prefix = rising_factorial_ratio(a+b, a, 20);
               }
               else
               {
                  prefix = 1;
               }
               fract = ibeta_a_step(a, b, x, y, 20, l, normalised);
               if(!invert)
                  fract = beta_small_b_large_a_series(a + 20, b, x, y, fract, prefix, l, normalised);
               else
               {
                  fract -= (normalised ? 1 : beta_imp(a, b, l));
                  invert = false;
                  fract = -beta_small_b_large_a_series(a + 20, b, x, y, fract, prefix, l, normalised);
               }
            }
         }
      }
   }
   else
   {
      // Both a,b >= 1:
      T lambda;
      if(a < b)
      {
         lambda = a - (a + b) * x;
      }
      else
      {
         lambda = (a + b) * y - b;
      }
      if(lambda < 0)
      {
         std::swap(a, b);
         std::swap(x, y);
         invert = !invert;
      }

      if(b < 40)
      {
         if(b * x <= 0.7)
         {
            if(!invert)
               fract = ibeta_series(a, b, x, T(0), l, normalised);
            else
            {
               fract = -(normalised ? 1 : beta_imp(a, b, l));
               invert = false;
               fract = -ibeta_series(a, b, x, fract, l, normalised);
            }
         }
         else if(a > 15)
         {
            // sidestep so we can use the series representation:
            int n = static_cast<int>(boost::math::tools::real_cast<long double>(floor(b)));
            if(n == b)
               --n;
            T bbar = b - n;
            T prefix;
            if(!normalised)
            {
               prefix = rising_factorial_ratio(a+bbar, bbar, n);
            }
            else
            {
               prefix = 1;
            }
            fract = ibeta_a_step(bbar, a, y, x, n, l, normalised);
            fract += ibeta_series(a, bbar, x, T(0), l, normalised);
            fract /= prefix;
         }
         else if(normalised)
         {
            // the formula here for the non-normalised case is tricky to figure
            // out (for me!!), and requires two pochhammer calculations rather 
            // than one, so leave it for now....
            int n = static_cast<int>(boost::math::tools::real_cast<long double>(floor(b)));
            if(n == b)
               --n;
            T bbar = b - n;
            fract = ibeta_a_step(bbar, a, y, x, n, l, normalised);
            fract += ibeta_a_step(a, bbar, x, y, 20, l, normalised);
            if(invert)
               fract -= (normalised ? 1 : beta_imp(a, b, l));
            fract = ibeta_series(a+20, bbar, x, fract, l, normalised);
            if(invert)
            {
               fract = -fract;
               invert = false;
            }
         }
         else
            fract = ibeta_fraction2(a, b, x, y, l, normalised);
      }
      else
         fract = ibeta_fraction2(a, b, x, y, l, normalised);
   }
   return invert ? (normalised ? 1 : beta_imp(a, b, l)) - fract : fract;
}

} // namespace detail

//
// The actual function entrypoints now follow, these just figure out
// which Lanczos approximation to use and forward to the implementation
// functions:
//
template <class T>
T beta(T a, T b)
{
   typedef typename lanczos::lanczos_traits<T>::value_type value_type;
   typedef typename lanczos::lanczos_traits<T>::evaluation_type evaluation_type;
   return tools::checked_narrowing_cast<T>(detail::beta_imp(static_cast<value_type>(a), static_cast<value_type>(b), evaluation_type()), BOOST_CURRENT_FUNCTION);
}

template <class T>
T beta(T a, T b, T x)
{
   typedef typename lanczos::lanczos_traits<T>::value_type value_type;
   typedef typename lanczos::lanczos_traits<T>::evaluation_type evaluation_type;
   return tools::checked_narrowing_cast<T>(detail::ibeta_imp(static_cast<value_type>(a), static_cast<value_type>(b), static_cast<value_type>(x), evaluation_type(), false, false), BOOST_CURRENT_FUNCTION);
}

template <class T>
T betac(T a, T b, T x)
{
   typedef typename lanczos::lanczos_traits<T>::value_type value_type;
   typedef typename lanczos::lanczos_traits<T>::evaluation_type evaluation_type;
   return tools::checked_narrowing_cast<T>(detail::ibeta_imp(static_cast<value_type>(a), static_cast<value_type>(b), static_cast<value_type>(x), evaluation_type(), true, false), BOOST_CURRENT_FUNCTION);
}

template <class T>
T ibeta(T a, T b, T x)
{
   typedef typename lanczos::lanczos_traits<T>::value_type value_type;
   typedef typename lanczos::lanczos_traits<T>::evaluation_type evaluation_type;
   return tools::checked_narrowing_cast<T>(detail::ibeta_imp(static_cast<value_type>(a), static_cast<value_type>(b), static_cast<value_type>(x), evaluation_type(), false, true), BOOST_CURRENT_FUNCTION);
}

template <class T>
T ibetac(T a, T b, T x)
{
   typedef typename lanczos::lanczos_traits<T>::value_type value_type;
   typedef typename lanczos::lanczos_traits<T>::evaluation_type evaluation_type;
   return tools::checked_narrowing_cast<T>(detail::ibeta_imp(static_cast<value_type>(a), static_cast<value_type>(b), static_cast<value_type>(x), evaluation_type(), true, true), BOOST_CURRENT_FUNCTION);
}

}} // namespaces

#include <boost/math/special_functions/ibeta_inverse.hpp>

#endif // BOOST_MATH_SPECIAL_BETA_HPP



